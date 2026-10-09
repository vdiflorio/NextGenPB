/*
 *  Copyright (C) 2019-2025 Carlo de Falco
 *  Copyright (C) 2020-2021 Martina Politi
 *  Copyright (C) 2021-2025 Vincenzo Di Florio
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef HAVE_PB_CLASS_H
#define HAVE_PB_CLASS_H

#include <bim_timing.h>
#include <tmesh_3d.h>

#define TIC() MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
#define TOC(S) MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc (S); }
#include <cmath>

const double p4esttol = 1 / std::pow (2, P8EST_QMAXLEVEL);

#include <algorithm>
#include <array>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <memory>

#include "raytracer.h"
#include <nanoshaper.h>


#include "wrapper_search.h"


// Problem parameters
constexpr double e_0 = 8.85418781762e-12; //Dielectric void const [F/m]
constexpr double kb = 1.380649e-23; //Boltzmann constant [J/K]
// constexpr double T = 273.15 + 25; //Temperature [K]
// constexpr double T = 273.15 + 25.00005445539216;
constexpr double e = 1.602176634e-19; //Charge of an electron [C]
constexpr double N_av = 6.022e23; //Avogadro Number [mol^-1]
constexpr double Angs = 1e-10; //Angstrom [m]
constexpr double pi = 3.14159265358979323846;


// Mobile-ion model for the nonlinear PBE (1:1 salt, potential u in kT/e).
// The ionic charge density in the solvent is  rho_ion = -eps_out k^2 g(u),
// i.e. the Newton residual/Jacobian use C*g(u) and C*g'(u), C = reaction_nodes.
//
//   nu == 0 : ideal Boltzmann ions        g = sinh u
//   nu  > 0 : Size modified lattice gas (finite ion size a):
//               D(u) = 1 + nu (cosh u - 1),   nu = 2 a^3 n_b  (bulk packing fraction)
//               g = sinh u / D,   g' = ((1 - nu) cosh u + nu) / D^2
//             for nu < 1, g is bounded (|g| -> 1/nu) and 0 < g' <= cosh u with
//             g' -> 0 at large |u|: counterions saturate at the close-packing
//             density 1/a^3 instead of growing as exp(|u|).
//
// Free energy (per 2 n_b kT, i.e. per C/(4 pi) in code units):
//   osm(u)   = (P - P0) / (2 n_b kT)      osmotic pressure excess
//            = cosh u - 1                 (ideal)
//            = log(D(u)) / nu             (steric; -> cosh u - 1 as nu -> 0)
//   f_exc(u) = -1/2 rho_s phi / (2 n_b kT) - osm(u) = u/2 g(u) - osm(u)
//            excess ionic free energy density; it cancels to O(u^4) in the
//            linear limit:
//              f_exc = (1 - 3 nu) u^4/24 + (1 - 15 nu + 30 nu^2) u^6/360
//                    + (1 - 63 nu + 420 nu^2 - 630 nu^3) u^8/13440 + O(u^10)
//            The closed form loses ~12 eps/u^2 in relative precision (it is the
//            difference of two O(u^2) terms, each accurate to eps once cosh u - 1
//            is evaluated as 2 sinh^2(u/2)); the series truncates at ~u^10/1e6.
//            They break even at |u| ~ 0.05, where the series takes over.
//            f_exc >= 0 only for nu = 0; for nu > 0 it goes like -u/(2 nu) at large |u|.
struct ion_model_t {
  double nu = 0.0;

  // cosh u - 1 without cancellation (cosh(u) - 1.0 has absolute error ~eps).
  static double
  coshm1 (double u)
  {
    const double s = std::sinh (0.5 * u);
    return 2.0 * s * s;
  }

  double
  D (double u) const
  { return 1.0 + nu * coshm1 (u); }

  double
  g (double u) const
  { return std::sinh (u) / D (u); }

  double
  dg (double u) const
  {
    const double d = D (u);
    return ((1.0 - nu) * std::cosh (u) + nu) / (d * d);
  }

  double
  osm (double u) const
  {
    const double c = coshm1 (u);
    return nu > 0.0 ? std::log1p (nu * c) / nu : c;
  }

  double
  f_exc (double u) const
  {
    if (std::fabs (u) < 0.05) {
      const double u2 = u * u, u4 = u2 * u2, nu2 = nu * nu;
      return u4 * ((1.0 - 3.0 * nu) / 24.0
                   + u2 * ((1.0 - 15.0 * nu + 30.0 * nu2) / 360.0
                           + u2 * (1.0 - 63.0 * nu + 420.0 * nu2 - 630.0 * nu2 * nu) / 13440.0));
    }
    return 0.5 * u * g (u) - osm (u);
  }

  const char *
  name () const
  { return nu > 0.0 ? "steric" : "sinh"; }
};

// ------------------------------------------------------------
//  Size-modified model with non-uniform sizes, 1:1 salt: cation volume v_p,
//  anion volume v_m and solvent volume v_w all different (Chu 2007,
//  Li 2009b, Zhou 2011; manuscript eqs. smpb_density, s_equation).
//  Same conventions as ion_model_t: u = e phi / kT, rho_ion = -2 e n_b g(u).
//
//  Unknown at each node: s = log(theta_w / theta_w^b). The ion densities are
//      n_p = n_b e^{x_p},  x_p = -u + r_p s,      n_m = n_b e^{x_m},  x_m = u + r_m s,
//  with r_p = v_p / v_w, r_m = v_m / v_w, and eq. s_equation reads
//      F(s) = thb e^s + vf_p e^{x_p} + vf_m e^{x_m} = 1,
//  thb = theta_w^b = 1 - vf_p - vf_m,  vf_p = v_p n_b,  vf_m = v_m n_b.
//  Newton is applied to G(s) = log F(s), which has the same root: a log of a
//  sum of exponentials of s is convex and increasing, so Newton converges
//  from any start in a few steps, with no bisection (Newton on theta_w can
//  step to theta_w < 0, and Newton on F is slow at large |u|).
//  Start: s0 = min(sigma u, asymptotes): sigma u is the first-order root
//  (manuscript eq. ds_lin), (u - log vf_p) / r_p and (-u - log vf_m) / r_m are
//  the points where the cation or the anion alone fills the volume; the root
//  lies left of both.
//  Stop when |ds| <= 1e-8 |s| (|u| <= 1) or 1e-8 max(1,|s|) (|u| > 1):
//  convergence is quadratic, so the remaining error is ~ds^2, i.e. round-off.
//
//  Quantities:
//    g     = (n_m - n_p) / (2 n_b) = (expm1(x_m) - expm1(x_p)) / 2
//    g'    = [v_w theta (n_p + n_m) + n_p n_m (v_p + v_m)^2]
//            / [2 n_b (v_w theta + v_p^2 n_p + v_m^2 n_m)]
//            (eq. c_models rewritten with the Lagrange identity: all terms
//            >= 0, no cancellation when the counterions saturate, g' -> 0)
//    osm   = (P - P0) / (2 n_b kT)
//          = [osm_term(x_p) + osm_term(x_m)] / 2 + thb / (2 n_b v_w) osm_term(s)
//    f_exc = u g / 2 - osm
//          = [exc_term(x_p) + exc_term(x_m)] / 2 + thb / (2 n_b v_w) exc_term(s)
//  osm_term(x) = e^x - 1 - x and exc_term(x) = (x/2 - 1) e^x + 1 + x/2 come
//  from eqs. smpb_pressure and G_exc_smpb rewritten with the definition of
//  theta_w and electroneutrality: the solvent enters as one more species of
//  bulk density thb / v_w with x = s. Each term is O(u^2) (osm) or O(u^3)
//  (f_exc), so the large cancellations of the direct formulas are gone; the
//  one left inside osm_term and exc_term is handled by their Taylor series.
//  For v_p = v_m = v_w this is the model of ion_model_t (which is what the
//  solver uses in that case); for v_p = v_m = 0 it is the ideal (sinh) model.
// ------------------------------------------------------------
struct ion_model_nonuniform_t {
  static constexpr int max_iter = 50;

  // input
  double n_b = 0.0;      // bulk density of each ion [1/A^3]
  double v_p = 0.0;      // cation volume [A^3], 0 = point-like
  double v_m = 0.0;      // anion volume [A^3], 0 = point-like
  double v_w = 0.0;      // solvent volume [A^3]

  // derived in init ()
  double vf_p = 0.0, vf_m = 0.0;          // bulk volume fractions v_p n_b, v_m n_b
  double log_vf_p = 0.0, log_vf_m = 0.0;  // their logs (only used if v > 0)
  double r_p = 0.0, r_m = 0.0;            // v_p / v_w, v_m / v_w
  double thb = 1.0, log_thb = 0.0;        // bulk free volume fraction theta_w^b
  double w_s = 0.0;                       // thb / (2 n_b v_w): weight of the solvent terms
  double sigma = 0.0;                     // ds/du at u = 0

  // e^x - 1 - x: the term of each species in the osmotic pressure
  static double
  osm_term (double x)
  {
    if (std::fabs (x) < 0.5) {
      double t = 0.0;                     // Horner: t = sum_{k>=1} x^k / (k+1)!
      for (int m = 20; m >= 2; --m)
        t = (t + 1.0) * x / m;
      return t * x;
    }
    return std::expm1 (x) - x;
  }

  // (x/2 - 1) e^x + 1 + x/2 = sum_{m>=3} (m - 2) x^m / (2 m!):
  // the term of each species in the excess free energy
  static double
  exc_term (double x)
  {
    if (std::fabs (x) < 0.5) {
      double t = 0.0, f = 1.0, xm = 1.0;
      for (int m = 1; m <= 20; ++m) {
        f *= m;
        xm *= x;
        if (m >= 3)
          t += (m - 2) * xm / (2.0 * f);
      }
      return t;
    }
    return (0.5 * x - 1.0) * std::expm1 (x) + x;
  }

  // Bulk density n_b [1/A^3] and volumes [A^3]; false if the bulk has no
  // free volume left (theta_w^b <= 0) or the input is not valid.
  bool
  init (double nb, double vp, double vm, double vw)
  {
    n_b = nb;  v_p = vp;  v_m = vm;  v_w = vw;
    if (!(n_b > 0.0) || !(v_w > 0.0) || v_p < 0.0 || v_m < 0.0)
      return false;
    vf_p = v_p * n_b;
    vf_m = v_m * n_b;
    thb = 1.0 - vf_p - vf_m;
    if (!(thb > 0.0))
      return false;
    log_vf_p = v_p > 0.0 ? std::log (vf_p) : 0.0;
    log_vf_m = v_m > 0.0 ? std::log (vf_m) : 0.0;
    r_p = v_p / v_w;
    r_m = v_m / v_w;
    log_thb = std::log (thb);
    w_s = thb / (2.0 * n_b * v_w);
    sigma = v_w * n_b * (v_p - v_m) / (v_w * thb + n_b * (v_p * v_p + v_m * v_m));
    return true;
  }

  // g'(0) = kappa_eff^2 / kappa^2 (manuscript eq. keff_sym), diagnostics only.
  double
  dg0 () const
  {
    const double vs = v_p + v_m;
    return (v_w * thb * 2.0 * n_b + n_b * n_b * vs * vs)
           / ((v_w * thb + n_b * (v_p * v_p + v_m * v_m)) * 2.0 * n_b);
  }

  // Root s(u) of G; returns the number of iterations, or -1 if Newton did
  // not converge in max_iter (s is then the last iterate).
  //   |u| <= 1: G = log1p(thb expm1(s) + vf_p expm1(x_p) + vf_m expm1(x_m)),
  //             i.e. log F with thb + vf_p + vf_m = 1 used exactly: accurate
  //             to relative round-off in s as u -> 0 (s = O(u), O(u^2) for
  //             v_p = v_m). All exponents are O(1) there: no overflow.
  //   |u| > 1 : G = log-sum-exp of log thb + s, log vf_p + x_p, log vf_m + x_m,
  //             safe from overflow, accurate to absolute round-off in s.
  int
  solve_s (double u, double &s) const
  {
    s = sigma * u;
    if (v_p > 0.0)
      s = std::min (s, (u - log_vf_p) / r_p);
    if (v_m > 0.0)
      s = std::min (s, (-u - log_vf_m) / r_m);

    const bool small = std::fabs (u) <= 1.0;
    for (int it = 1; it <= max_iter; ++it) {
      const double x_p = -u + r_p * s, x_m = u + r_m * s;
      double G, dG;
      if (small) {
        const double f = thb * std::expm1 (s) + vf_p * std::expm1 (x_p)
                         + vf_m * std::expm1 (x_m);
        dG = (thb * std::exp (s) + r_p * vf_p * std::exp (x_p)
              + r_m * vf_m * std::exp (x_m)) / (1.0 + f);
        G  = std::log1p (f);
      }
      else {
        // log-sum-exp: subtract the largest exponent before exp. A point-like
        // ion (v = 0) has no term; no infinities are used (-Ofast assumes none).
        const double L_w = log_thb + s;
        const double L_p = log_vf_p + x_p, L_m = log_vf_m + x_m;
        double m = L_w;
        if (v_p > 0.0) m = std::max (m, L_p);
        if (v_m > 0.0) m = std::max (m, L_m);
        const double e_w = std::exp (L_w - m);
        const double e_p = v_p > 0.0 ? std::exp (L_p - m) : 0.0;
        const double e_m = v_m > 0.0 ? std::exp (L_m - m) : 0.0;
        const double W = e_w + e_p + e_m;
        G  = m + std::log (W);
        dG = (e_w + r_p * e_p + r_m * e_m) / W;
      }
      const double ds = -G / dG;
      s += ds;
      // Quadratic convergence: after a step ds the error is ~ds^2.
      // Relative to |s| near u = 0 (floor 1e-16 |u| where s crosses zero),
      // relative to max(1, |s|) elsewhere.
      const double scale = small ? std::max (std::fabs (s), 1.0e-8 * std::fabs (u))
                                 : std::max (1.0, std::fabs (s));
      if (std::fabs (ds) <= 1.0e-8 * scale)
        return it;
    }
    return -1;
  }

  // g and g' for the Newton system; false if the root was not found.
  bool
  newton (double u, double &g, double &dg) const
  {
    double s;
    const bool ok = solve_s (u, s) > 0;
    const double x_p = -u + r_p * s, x_m = u + r_m * s;
    const double n_p = n_b * std::exp (x_p), n_m = n_b * std::exp (x_m);
    const double theta = thb * std::exp (s);
    const double vs = v_p + v_m;
    g  = 0.5 * (std::expm1 (x_m) - std::expm1 (x_p));
    dg = (v_w * theta * (n_p + n_m) + n_p * n_m * vs * vs)
         / (2.0 * n_b * (v_w * theta + v_p * v_p * n_p + v_m * v_m * n_m));
    return ok;
  }

  // Everything the energy and Gamma+- need, from one root.
  struct obs_t {
    double g, osm, f_exc;
    double dn_p, dn_m;           // n_p / n_b - 1, n_m / n_b - 1
  };

  bool
  observables (double u, obs_t &o) const
  {
    double s;
    const bool ok = solve_s (u, s) > 0;
    const double x_p = -u + r_p * s, x_m = u + r_m * s;
    o.dn_p  = std::expm1 (x_p);
    o.dn_m  = std::expm1 (x_m);
    o.g     = 0.5 * (o.dn_m - o.dn_p);
    o.osm   = 0.5 * (osm_term (x_p) + osm_term (x_m)) + w_s * osm_term (s);
    o.f_exc = 0.5 * (exc_term (x_p) + exc_term (x_m)) + w_s * exc_term (s);
    return ok;
  }
};

// ------------------------------------------------------------
//  Stern layer (ion-free shell): union of the spheres |x - r_i| < R_i + s,
//  the APBS (ion radius) / Amber pbsa (iprob) convention.
//
//  Uniform cell list with cell side L = max(R_i) + s: a sphere that contains
//  x has its centre in the cell of x or in one of the 26 neighbours, so a
//  query visits a few tens of atoms whatever the size of the molecule, and
//  returns at the first sphere that contains x.
// ------------------------------------------------------------
struct stern_grid_t {
  double L = 1.0;
  std::array<double,3> lo {}, hi {};       // atom box enlarged by L
  std::array<int,3> n {};                  // cells per direction
  std::vector<int> start;                  // atoms of cell c: [start[c], start[c+1])
  std::vector<std::array<double,4>> at;    // x, y, z, (R_i + s)^2, sorted by cell

  void
  build (const std::vector<std::array<double,3>> &pos,
         const std::vector<double> &rad, double s);

  // Does the box [a, b] overlap the region where inside() can be true?
  bool
  box_near (const std::array<double,3> &a, const std::array<double,3> &b) const
  {
    for (int d = 0; d < 3; ++d)
      if (b[d] < lo[d] || a[d] > hi[d])
        return false;

    return true;
  }

  bool
  inside (double x, double y, double z) const;
};


struct
  poisson_boltzmann {

  p4est_topidx_t simple_conn_num_vertices;
  p4est_topidx_t simple_conn_num_trees;
  std::unique_ptr<double[]> simple_conn_p;
  std::unique_ptr<p4est_topidx_t[]> simple_conn_t;
  // double *simple_conn_p;
  // p4est_topidx_t *simple_conn_t;
  std::vector<std::pair<p4est_topidx_t, p4est_topidx_t>> bcells;

  std::vector<NS::Atom> atoms;
  std::vector<std::array<double,3>> pos_atoms;
  std::vector<double> charge_atoms;
  std::vector<double> r_atoms;
  std::vector<int> index_atoms;


  // input parameters
  std::string filetype;
  std::string radiusfilename;
  std::string chargefilename;
  std::string name_pqr;
  int write_pqr;
  // Center of the system
  double cc[3];

  //Cubic mesh:
  double ll[3]; //min value between all the coordinates
  double rr[3]; //max value between all the coordinates

  //Stretched mesh:
  double l_c[3]; //min x, y, z value
  double r_c[3]; //max x, y, z value
  double l_cr[3]; //refined box min x, y, z value
  double r_cr[3]; //refined box max x, y, z value
  double l_box[3]; //refined box min x, y, z value focusing
  double r_box[3]; //refined box max x, y, z value focusing
  double pot_bc = 0.0;
  //number of trees
  p4est_topidx_t num_trees[3];
  double len;

  //Focusing mesh:
  double cc_focusing[3];
  int n_grid;

  //mesh:
  int maxlevel;
  int minlevel;
  int unilevel;
  int outlevel;
  int loc_refinement;
  int mesh_shape;
  int refine_box;
  int rand_center;
  int rand_seed;         // rand_center: 0 = random seed, n > 0 = fixed seed
  int scale_level;
  int scale_level_min_box;
  double scale, scale_min, scale_max;
  double perfil1, perfil2;
  int loc_ref = 0;
  int aligned = 0;
  //model:
  int linearized;
  int bc;
  double e_in, e_out, ionic_strength; //[M]
  double ion_size;        // [Angstrom] finite ion size a; 0 -> ideal (sinh) ions
  ion_model_t ion_model;  // nu = 2 a^3 n_b, set right after reading the options
  double T;
  int calc_energy;
  double energy_pol = 0.0;
  double energy_react = 0.0;
  double coul_energy = 0.0;
  double energy_exc = 0.0;

  int calc_coulombic;
  int calc_potential_term;
  int calc_field_term;



  //surface:
  NS::surface_type surf_type;
  int surf_type_num = 0;
  double surf_param;
  double prb_radius = 1.4; //typical prob radious;
  int stern_layer_surf;
  double stern_layer;
  unsigned num_threads;

  //algorithm:
  std::string linear_solver_name;
  std::string linear_solver_options;
  int newton_compress;   // 1 (default): after Newton it 0 (== linear solve) map solvent phi -> 2 asinh(phi/2)

  MPI_Comm mpicomm;
  tmesh_3d tmsh;

  std::string optionsfilename;
  std::string pqrfilename;
  std::string p4estfilename;
  std::string surffilename;
  std::string markerfilename;

  std::string pqr_atoms;

  //post_processing
  int atoms_write;
  int surf_write;
  std::string map_type;
  int potential_map;
  int eps_map;
  int dataset_write;    // 1 = export AI/ML vertex dataset to CSV
  std::map<int, tmesh_3d::quadrant_t> lookup_table;


  std::vector<double> marker;
  std::vector<double> epsilon;
  std::vector<double> epsilon_in;
  std::vector<double> epsilon_out;

  std::vector<int> border_quad;

  std::vector<double> const_ones;

  std::unique_ptr<distributed_vector> markn;
  std::unique_ptr<distributed_vector> epsilon_nodes;
  std::unique_ptr<distributed_vector> reaction_nodes;
  double net_charge;

  std::unique_ptr<distributed_vector> phi;
  std::unique_ptr<distributed_vector> rho_fixed;
  std::unique_ptr<distributed_vector> ones;
  std::unique_ptr<distributed_vector> rhs;
  std::unique_ptr<distributed_sparse_matrix> A;

  static constexpr
  std::array<int, 12> edge_axis = {0,1,0,1,0,1,0,1,2,2,2,2};

  static constexpr
  std::array<int, 24> edge2nodes = {
    0, 1,
    1, 3,
    2, 3,
    0, 2,
    4, 5,
    5, 7,
    6, 7,
    4, 6,
    0, 4,
    1, 5,
    3, 7,
    2, 6
  };



  int edgeTable[256]= {
    0x0, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x99, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x33, 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0xaa, 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x66, 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff, 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55, 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0xcc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x55, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0xff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x66, 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa, 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33, 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99, 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0
  };
  int triTable[256][16] = {
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
    {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
    {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
    {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
    {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
    {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
    {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
    {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
    {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
    {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
    {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
    {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
    {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
    {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
    {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
    {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
    {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
    {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
    {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
    {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
    {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
    {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
    {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
    {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
    {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
    {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
    {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
    {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
    {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
    {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
    {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
    {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
    {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
    {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
    {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
    {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
    {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
    {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
    {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
    {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
    {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
    {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
    {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
    {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
    {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
    {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
    {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
    {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
    {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
    {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
    {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
    {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
    {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
    {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
    {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
    {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
    {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
    {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
    {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
    {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
    {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
    {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
    {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
    {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
    {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
    {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
    {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
    {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
    {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
    {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
    {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
    {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
    {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
    {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
    {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
    {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
    {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
    {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
    {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
    {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
    {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
    {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
    {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
    {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
    {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
    {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
    {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
    {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
    {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
    {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
    {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
    {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
    {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
    {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
    {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
    {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
    {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
  };

  std::array<std::array<int,3>,5> triangles;

  poisson_boltzmann (int maxlevel_ = 9, int minlevel_ = 3, int unilevel_ = 5, int mesh_shape_ = 1,
                     int bc_ = 1, int linearized_ = 1,
                     double e_in_ = 2.0, double e_out_ = 80.0, double ionic_strength_ = 0.145,
                     std::string linear_solver_name_ = "mumps", std::string linear_solver_options_ = "",
                     MPI_Comm mpicomm_ = MPI_COMM_WORLD)
    : maxlevel (maxlevel_),
      minlevel (minlevel_),
      unilevel (unilevel_),
      mesh_shape (mesh_shape_),
      bc (bc_),
      linearized (linearized_),
      e_in (e_in_),
      e_out (e_out_),
      ionic_strength (ionic_strength_),
      linear_solver_name (linear_solver_name_),
      linear_solver_options (linear_solver_options_),
      mpicomm (mpicomm_),
      tmsh (mpicomm)
  { };

  double
  levelsetfun (double x, double y, double z);

  /**
   * @brief Determines whether a point is inside a molecular surface along a specified direction.
   *
   * This function evaluates whether a point defined by coordinates `(x, y, z)` lies
   * inside, outside, or on the boundary of a molecular surface based on ray intersections
   * retrieved from a ray cache.
   *
   * @param ray_cache A reference to the ray tracing cache that stores intersection
   *                  data and manages ray operations.
   * @param x The x-coordinate of the point.
   * @param y The y-coordinate of the point.
   * @param z The z-coordinate of the point.
   * @param dir The direction of the evaluation:
   *            - `0`: Evaluate in the yz-plane.
   *            - `1`: Evaluate in the xz-plane.
   *            - `2`: Evaluate in the xy-plane.
   *
   * @return A value indicating the position of the point relative to the molecular surface:
   *         - `0.0`: The point is outside the surface.
   *         - `1.0`: The point is inside the surface.
   *         - `-1.0`: The point requires additional ray tracing or data is unavailable.
   *
   * ### Algorithm Details
   * - Coordinates `(x, y, z)` are reordered based on the specified evaluation direction (`dir`).
   * - Intersections along the specified direction are retrieved from the ray cache.
   * - If no intersections are found, or the point lies before the first intersection, it is marked as outside.
   * - Iteratively evaluates whether the point alternates between inside and outside based on intersection crossings.
   * - Returns `1.0` if the number of intersections passed is odd (inside) and `0.0` if even (outside).
   *
   * ### Notes
   * - Requires the `ray_cache` to be properly initialized and populated with
   *   intersection data.
   * - Assumes that intersections are sorted and stored in the ray cache.
   * - This function is intended to be used within an MPI-based parallel environment.
   */
  double
  is_in_ns_surf (ray_cache_t & ray_cache, double x, double y, double z, int dir);

  static int
  uniform_refinement (tmesh_3d::quadrant_iterator quadrant)
  {
    return 1;
  }


  void
  create_mesh ();

  void
  create_mesh_ns ();

  void
  create_mesh_scale ();

  int
  parse_options (int argc, char **argv);

  void
  print_options ();

  void
  read_atoms_from_pqr (std::basic_istream<char> &inputfile);

  void
  read_atoms_from_pdb (std::basic_istream<char> &inputfile);

  void
  read_atoms_from_class ();

  void
  broadcast_vectors ();

  void
  write_atoms_to_pqr (std::basic_ostream<char> &outputfile);

  friend std::basic_istream<char>&
  operator>> (std::basic_istream<char>& inputfile, NS::Atom &a);

  friend std::basic_istream<char>&
  operator>> (std::basic_istream<char>& inputfile, std::array<float,5> &a);

  void
  write_potential_on_atoms ();

  void
  write_potential_on_atoms_fast ();

  void
  init_tmesh ();

  void
  init_tmesh_with_refine_scale ();

  void
  init_tmesh_with_refine_box_scale ();

  bool
  is_in (const NS::Atom& i, tmesh_3d::quadrant_iterator q);

  bool
  is_in_ref (const NS::Atom& i, tmesh_3d::quadrant_iterator q);

  void
  refine_surface (ray_cache_t & ray_cache);

  void
  refine_only_surface (ray_cache_t & ray_cache);

  /**
   * @brief Initializes and updates markers for quadrants in a forest mesh.
   *
   * This function sets up markers for quadrants within a forest mesh to distinguish
   * whether they are inside, outside, or on the boundary of a molecule. It also
   * computes dielectric properties and reaction rates for the nodes based on their
   * location relative to the molecule and optional Stern layer.
   *
   * @param ray_cache A reference to a ray tracing cache structure that holds
   *                  information on required rays and their intersections.
   *
   * This function performs the following:
   * - Initializes the dielectric constants (`epsilon`) for the inside and outside
   *   regions.
   * - Computes reaction rates based on ionic strength and other physical parameters.
   * - Loops over quadrants in the mesh, evaluating the position of nodes relative to
   *   the molecule and Stern layer (if present).
   * - Updates markers for quadrants and nodes based on their location:
   *     - `0.0`: Inside the molecule.
   *     - `0.5`: On the boundary of the molecule.
   *     - `1.0`: Outside the molecule.
   * - Updates reaction and dielectric properties for nodes inside the molecule,
   *   and the reaction for nodes inside the Stern layer.
   * - Handles MPI-based parallelism, including barrier synchronization and data
   *   exchange.
   * - Ensures rays are calculated and cached for points near molecular boundaries.
   *
   * ### Stern Layer Handling
   * If `stern_layer_surf` is set to `1`, `reaction_nodes` is also set to zero on
   * the nodes inside the union of the spheres R_i + `stern_layer` (stern_grid_t),
   * so every solver and post-processing step that uses the nodal reaction
   * coefficient sees the ion-free shell.
   *
   * ### Constants
   * - Dielectric constants for inside (`eps_in`) and outside (`eps_out`) regions.
   * - Reaction rate constant based on ionic strength (`k2`).
   *
   * ### Notes
   * - The function assumes that the mesh and ray tracing cache are correctly
   *   initialized before calling.
   * - It performs two cycles of refinement for multi-process configurations and
   *   one cycle for single-process configurations.
   */
  void
  create_markers (ray_cache_t & ray_cache);

  void
  export_tmesh (ray_cache_t & ray_cache);

  void
  export_potential_map (ray_cache_t & ray_cache);

  void
  export_marked_tmesh ();

  void
  export_p4est ();

  void
  assemple_system_matrix (ray_cache_t & ray_cache);

  void 
  newton_solve (ray_cache_t & ray_cache);

  void 
  assemble_newton_system (ray_cache_t & ray_cache, distributed_vector & phi_cur);

  void
  create_density_map (ray_cache_t & ray_cache);

  /**
   * @brief New distributed vector on the mesh nodes, ghosts already set up.
   *
   * With MPI it is a copy of epsilon_nodes, whose ghost entries were built
   * once by bim3a_solution_with_ghosts in create_markers: overwrite the owned
   * values, then v->assemble (op) updates the ghosts with one exchange,
   * without the sweep over the mesh and the remap. The ghost values of the
   * copy are those of epsilon_nodes until then. Valid after create_markers
   * (the mesh does not change afterwards).
   */
  std::unique_ptr<distributed_vector>
  new_node_vector ();

  void
  mumps_compute_electric_potential (ray_cache_t & ray_cache);

  void
  lis_compute_electric_potential (ray_cache_t & ray_cache);

  int
  classifyCube (tmesh_3d::quadrant_iterator& quadrant,double isolevel);

  int
  classifyCube_fast (tmesh_3d::quadrant_iterator& quadrant,double isolevel);

  std::tuple<std::array<double,8>, std::array<double,8>, std::vector<int>,std::vector<int>>
  classifyCube_flux (tmesh_3d::quadrant_iterator& quadrant,
                     std::array<double,8>& tmp_phi,
                     std::array<double,8>& tmp_eps);

  std::tuple<std::array<double,8>, std::array<double,8>, std::vector<int>,std::vector<int> >
  classifyCube_flux_fast (tmesh_3d::quadrant_iterator& quadrant,
                          std::array<double,8>& tmp_phi,
                          std::array<double,8>& tmp_eps);
  void
  energy_pot_field (ray_cache_t & ray_cache);

  void
  energy_excess_nonlinear (ray_cache_t & ray_cache);

  void
  write_potential_on_surface (ray_cache_t & ray_cache);

  double
  coulomb_boundary_conditions (double x, double y, double z);

  double
  analytic_solution (double x, double y, double z);

  void
  analitic_potential ();

  void
  abs_value_field (distributed_vector &phi);

  std::array<double,12>
  cube_fraction_intersection (tmesh_3d::quadrant_iterator& quadrant,
                              const ray_cache_t & ray_cache);

  void
  normal_intersection (tmesh_3d::quadrant_iterator& quadrant,
                       const ray_cache_t & ray_cache,
                       int edge, std::array<double,3> &norm,double &frac);

  int
  getTriangles (int cubeindex,
                std::array<std::array<int,3>,5> &triangles);

  bool
  controlla_coordinate (int i, const p8est_quadrant_t *quadrant);

  int
  cerca_atomo (p8est_t * p4est,
               p4est_topidx_t which_tree,
               p8est_quadrant_t * quadrant,
               p4est_locidx_t local_num,
               void *point);

  void
  search_points ();

  void
  write_dataset (ray_cache_t & ray_cache);
};

std::basic_istream<char>&
operator>> (std::basic_istream<char>& inputfile, NS::Atom &a);
std::basic_istream<char>&
operator>> (std::basic_istream<char>& inputfile, std::array<float,5> &a);


struct EdgeKey {
  int a, b;
  bool operator==(const EdgeKey& other) const {
    return a == other.a && b == other.b;
  }
};

struct EdgeHash {
  std::size_t operator()(const EdgeKey& k) const {
    return std::hash<int>()(k.a) ^ (std::hash<int>()(k.b) << 1);
  }
};

static constexpr
std::array<int, 24> edge2nodes_nn1 = {
  2, 4,  3, 5,  0, 6,  4, 1,  6, 0,  1, 4,
  4, 2,  0, 5,  1, 2,  0, 3,  2, 1,  3, 0
};
static constexpr
std::array<int, 24> edge2index_nn1 = {
  0, 2,  0, 3,  1, 2,  0, 2,  0, 3,  1, 3,
  1, 3,  1, 2,  0, 2,  1, 2,  1, 3,  0, 3
};
static constexpr
std::array<int, 24> edge2nodes_nn2 = {
  3, 5,  7, 2,  1, 7,  6, 3,  7, 1,  3, 6,
  5, 3,  2, 7,  5, 6,  4, 7,  6, 5,  7, 4
};
static constexpr
std::array<int, 24> edge2index_nn2 = {
  0, 2,  0, 3,  1, 2,  0, 2,  0, 3,  1, 3,
  1, 3,  1, 2,  0, 2,  1, 2,  1, 3,  0, 3
};

struct VertexData {
  int axis = -1;
  double phi0 = 0.0;
  double alpha = 0.0;
  double phi1 = 0.0;
  double phi2 = 0.0;
  double eps1 = 0.0;
  double eps2 = 0.0;
  std::array<double, 3> pos1 = {0.0, 0.0, 0.0};
  std::array<double, 3> pos0 = {0.0, 0.0, 0.0};
  std::array<double, 3> N    = {0.0, 0.0, 0.0};
  std::array<double, 4> phi1_nn = {0.0, 0.0, 0.0, 0.0};
  std::array<double, 4> eps1_nn = {0.0, 0.0, 0.0, 0.0};
  std::array<double, 4> phi2_nn = {0.0, 0.0, 0.0, 0.0};
  std::array<double, 4> eps2_nn = {0.0, 0.0, 0.0, 0.0};
};

#endif
