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

#include "pb_class.h"
#include "GetPot"

#include <bim_distributed_vector.h>
#include <quad_operators_3d.h>
#include <mumps_class.h>
#include <lis_class.h>
// #include <lis_class_distributed.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdio.h>

#include <p8est.h>
#include <random>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <regex>
#include <cctype>
#include <unordered_map>
#include <cstdint>
#include <cstring>

// NaN/Inf test that survives -Ofast. The optimized builds (local_setting/*_fast.mk,
// recipe.def) use -Ofast, whose -ffinite-math-only lets the compiler assume
// that no NaN or Inf exists and fold std::isfinite (x) to true (checked with
// GCC 11.5). A double is NaN or Inf iff its 11 exponent bits are all set;
// testing them is integer arithmetic, which that flag does not touch.
static bool
is_finite_bits (double x)
{
  std::uint64_t b;
  std::memcpy (&b, &x, sizeof b);
  return (b & 0x7ff0000000000000ULL) != 0x7ff0000000000000ULL;
}


void
poisson_boltzmann::create_mesh ()
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);

  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/ (e*e); //adim e_out
  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/ (e_0*e_out*kb*T);
  double k = std::sqrt (k2);

  int nx, ny, nz;
  double scale_tmp, scale_x, scale_y, scale_z;
  double lmax = 0;
  double max_len = 0;
  double l[3];
  scale_level = unilevel;
  double maxradius = *std::max_element (r_atoms.begin (), r_atoms.end ());

  auto comp_pos_x = [] (const std::array<double, 3>& a1, const std::array<double, 3>& a2) -> bool {
    return a1[0] < a2[0];
  };
  auto comp_pos_y = [] (const std::array<double, 3>& a1, const std::array<double, 3>& a2) -> bool {
    return a1[1] < a2[1];
  };
  auto comp_pos_z = [] (const std::array<double, 3>& a1, const std::array<double, 3>& a2) -> bool {
    return a1[2] < a2[2];
  };


  auto minmax_x = std::minmax_element (pos_atoms.begin (), pos_atoms.end (), comp_pos_x);
  auto minmax_y = std::minmax_element (pos_atoms.begin (), pos_atoms.end (), comp_pos_y);
  auto minmax_z = std::minmax_element (pos_atoms.begin (), pos_atoms.end (), comp_pos_z);



  net_charge = std::accumulate (charge_atoms.begin(), charge_atoms.end(), 0.0);
  int num_atoms = charge_atoms.size ();

  if (rank == 0) {
    std::cout << "\n========== [ System Information ] ==========\n";
    std::cout << "  Number of atoms    : " << num_atoms << '\n';
    std::cout << "  Size protein [Å]   : ";
    std::cout << "[" << (*minmax_x.second)[0] - (*minmax_x.first)[0] + 2*maxradius << ", "
              << (*minmax_y.second)[1] - (*minmax_y.first)[1] + 2*maxradius << ", "
              << (*minmax_z.second)[2] - (*minmax_z.first)[2] + 2*maxradius << "]\n";

    if (std::fabs (net_charge - std::round (net_charge)) > 1.e-5)
      std::cerr << "  [WARNING] Net charge is not an integer: " << net_charge << '\n';

    std::cout << "  Solute epsilon     : " << e_in << '\n';
    std::cout << "  Solvent epsilon    : " << e_out << '\n';
    std::cout << "  Temperature        : " << T << " [K] \n";
    std::cout << "  Ionic strength     : " << ionic_strength << " [mol/L] \n";
    if (linearized == 0) {
      std::cout << "  Ion model          : "
                << (nonuniform ? "steric non-uniform (free volume solved at each node)"
                               : ion_model.name ()) << '\n';
      if (per_species_sizes) {
        std::cout << "  Ion sizes [A]      : a+ = " << ion_size_pos << ", a- = " << ion_size_neg
                  << ", solvent = " << solvent_size;
        if (ionic_strength == 0.0)
          std::cout << " (no ions: sizes have no effect)";
        else if (nonuniform && ion_size_pos == ion_size_neg)
          std::cout << " (equal ions, different solvent: kappa_eff = kappa)";
        else if (! nonuniform && ion_size > 0.0)
          std::cout << " (all equal: uniform model, same as ion_size = " << ion_size << ")";
        else if (! nonuniform)
          std::cout << " (point-like ions: solvent_size has no effect)";
        std::cout << '\n';
      }
      if (nonuniform) {
        // Linear response at u = 0: kappa_eff^2 = kappa^2 g'(0) (manuscript
        // eq. keff_sym), I_eq = I g'(0). Printed only, not used by the solver.
        const ion_model_nonuniform_t & m = ion_model_nu;
        const double dg0 = m.dg0 ();
        std::cout << "  Ion volumes [A^3]  : v+ = " << m.v_p << ", v- = " << m.v_m
                  << ", v_w = " << m.v_w << " (v+/v_w = " << m.r_p << ", v-/v_w = "
                  << m.r_m << ")\n";
        std::cout << "  Packing fraction   : " << 1.0 - m.thb
                  << " (1 - theta_w^b, bulk ion volume fraction)\n";
        std::cout << "  Linear response    : kappa_eff^2/kappa^2 = " << dg0
                  << ", I_eq = " << ionic_strength * dg0 << " [mol/L] (diagnostic only)\n";
        std::cout << "  Debye length       : 1/kappa_eff = " << 1.0 / (k * std::sqrt (dg0))
                  << " [Å], 1/kappa = " << 1.0 / k << " [Å]\n";
      } else if (ion_model.nu > 0.0) {
        std::cout << "  Ion size           : " << ion_size << " [Å] \n";
        std::cout << "  Packing fraction   : " << ion_model.nu << " (nu = 2 a^3 n_b)\n";
      }
    }
    std::cout << "============================================\n\n";
  }

  if (mesh_shape !=2) {
    l_c[0] = (*minmax_x.first)[0] - maxradius - 2*prb_radius;
    l_c[1] = (*minmax_y.first)[1] - maxradius - 2*prb_radius;
    l_c[2] = (*minmax_z.first)[2] - maxradius - 2*prb_radius;
    r_c[0] = (*minmax_x.second)[0] + maxradius + 2*prb_radius;
    r_c[1] = (*minmax_y.second)[1] + maxradius + 2*prb_radius;
    r_c[2] = (*minmax_z.second)[2] + maxradius + 2*prb_radius;

    for (int kk = 0; kk < 3; ++kk) {
      l[kk] = (r_c[kk] - l_c[kk]);
      cc[kk] = (r_c[kk] + l_c[kk])*0.5;
      lmax = l[kk] > lmax ? l[kk] : lmax;
    }

    // For random displacement of the grid. The shift is drawn on rank 0 and
    // broadcast: every rank must build the mesh around the same centre.
    if (rand_center == 1) {
      double shift[3] = {0.0, 0.0, 0.0};

      if (rank == 0) {
        // rand_seed > 0: fixed seed, reproducible shift; 0: random seed.
        std::random_device rd;
        std::mt19937 gen (rand_seed > 0 ? static_cast<unsigned> (rand_seed) : rd ());
        std::uniform_real_distribution<> dis (-1./scale*0.5, 1./scale*0.5);

        for (int n = 0; n < 3; ++n)
          shift[n] = dis (gen);
      }

      MPI_Bcast (shift, 3, MPI_DOUBLE, 0, mpicomm);

      for (int n = 0; n < 3; ++n)
        cc[n] += shift[n];

      if (rank == 0)
        std::cout << "  Random centre shift [Å]: [" << shift[0] << ", " << shift[1]
                  << ", " << shift[2] << "]"
                  << (rand_seed > 0 ? "  (seed " + std::to_string (rand_seed) + ")" : "")
                  << "\n";
    }

    for (int kk = 0; kk < 3; ++kk) {
      ll[kk] = cc[kk] - lmax*0.5;
      rr[kk] = cc[kk] + lmax*0.5;
    }

    //stretched box with perfil1
    l_c[0] -= l[0]*0.5* (1.0/perfil1 - 1);
    l_c[1] -= l[1]*0.5* (1.0/perfil1 - 1);
    l_c[2] -= l[2]*0.5* (1.0/perfil1 - 1);
    r_c[0] += l[0]*0.5* (1.0/perfil1 - 1);
    r_c[1] += l[1]*0.5* (1.0/perfil1 - 1);
    r_c[2] += l[2]*0.5* (1.0/perfil1 - 1);

    //cubic box with perfil1
    ll[0] -= lmax*0.5* (1.0/perfil1 - 1);
    ll[1] -= lmax*0.5* (1.0/perfil1 - 1);
    ll[2] -= lmax*0.5* (1.0/perfil1 - 1);
    rr[0] += lmax*0.5* (1.0/perfil1 - 1);
    rr[1] += lmax*0.5* (1.0/perfil1 - 1);
    rr[2] += lmax*0.5* (1.0/perfil1 - 1);

    l_cr[0] = l_c[0];
    r_cr[0] = r_c[0];
    l_cr[1] = l_c[1];
    r_cr[1] = r_c[1];
    l_cr[2] = l_c[2];
    r_cr[2] = r_c[2];
  }

  if (mesh_shape == 0) {

    //cubic box with max perfil2
    double size = 1.0/scale;
    scale_level = 0;

    for (int ii = 0; ii<29; ++ii) {
      size *= 2;
      scale_level ++;

      if (lmax/size < perfil2)
        break;
    }

    ll[0] = cc[0] - size/2;
    ll[1] = cc[1] - size/2;
    ll[2] = cc[2] - size/2;
    rr[0] = cc[0] + size/2;
    rr[1] = cc[1] + size/2;
    rr[2] = cc[2] + size/2;


    nx = (int) ( (r_cr[0] - cc[0])*scale + 0.5);
    ny = (int) ( (r_cr[1] - cc[1])*scale + 0.5);
    nz = (int) ( (r_cr[2] - cc[2])*scale + 0.5);

    //refined box
    l_cr[0] = cc[0] - nx*1.0/scale;
    l_cr[1] = cc[1] - ny*1.0/scale;
    l_cr[2] = cc[2] - nz*1.0/scale;
    r_cr[0] = cc[0] + nx*1.0/scale;
    r_cr[1] = cc[1] + ny*1.0/scale;
    r_cr[2] = cc[2] + nz*1.0/scale;

    double dist = size/2 - lmax*0.5;
    pot_bc = std::exp (-k*dist)/ (dist*eps_out);

    if (rank == 0) {
      std::cout << "========== [ Domain Information ] ==========\n";
      std::cout << "  Scale:  " << scale << "\n";

      std::cout << "  Center of the System [Å]:";
      std::cout << "  [" << cc[0] << ", " << cc[1] << ", " << cc[2] << "]\n";

      std::cout << "  Perfil outer box:  " << perfil2 << "\n";
      std::cout << "  Complete Domain Box Size [Å]:\n";
      std::cout << "      x = [" << ll[0] << ", " << rr[0] << "]\n";
      std::cout << "      y = [" << ll[1] << ", " << rr[1] << "]\n";
      std::cout << "      z = [" << ll[2] << ", " << rr[2] << "]\n";

      std::cout << "  Perfil uniform grid:  " << perfil1 << "\n";
      std::cout << "  Uniform grid Size [Å]:\n";
      std::cout << "      x = [" << l_cr[0] << ", " << r_cr[0] << "]\n";
      std::cout << "      y = [" << l_cr[1] << ", " << r_cr[1] << "]\n";
      std::cout << "      z = [" << l_cr[2] << ", " << r_cr[2] << "]\n";

      std::cout << "  Number of Subdivisions in the Uniform grid:";
      std::cout << "  nx = " << nx * 2 << "  ny = " << ny * 2 << "  nz = " << nz * 2 << '\n';

      std::cout << "============================================\n";
    }

    simple_conn_num_vertices = 8;
    simple_conn_num_trees = 1;

    simple_conn_p = std::make_unique<double[]> (simple_conn_num_vertices*3);
    simple_conn_t = std::make_unique<p4est_topidx_t[]> (simple_conn_num_vertices);

    auto tmp_p = {ll[0], ll[1], ll[2],
                  rr[0], ll[1], ll[2],
                  ll[0], rr[1], ll[2],
                  rr[0], rr[1], ll[2],
                  ll[0], ll[1], rr[2],
                  rr[0], ll[1], rr[2],
                  ll[0], rr[1], rr[2],
                  rr[0], rr[1], rr[2]
                 };
    auto tmp_t = {1, 2, 3, 4, 5, 6, 7, 8};

    std::copy (tmp_p.begin (), tmp_p.end (), simple_conn_p.get ());
    std::copy (tmp_t.begin (), tmp_t.end (), simple_conn_t.get ());

    for (int i = 0; i<6; i++)
      bcells.push_back (std::make_pair (0, i));
  } else if (mesh_shape == 1) {
    l_cr[0] = ll[0];
    l_cr[1] = ll[1];
    l_cr[2] = ll[2];
    r_cr[0] = rr[0];
    r_cr[1] = rr[1];
    r_cr[2] = rr[2];
    scale = (1<<unilevel)/ (rr[0]-ll[0]);

    if (refine_box == 1) {
      //cubic box with perfil2
      ll[0] = cc[0] - lmax*0.5*1.0/perfil2;
      ll[1] = cc[1] - lmax*0.5*1.0/perfil2;
      ll[2] = cc[2] - lmax*0.5*1.0/perfil2;
      rr[0] = cc[0] + lmax*0.5*1.0/perfil2;
      rr[1] = cc[1] + lmax*0.5*1.0/perfil2;
      rr[2] = cc[2] + lmax*0.5*1.0/perfil2;
      scale_tmp = (1<<unilevel)/ (rr[0]-ll[0]);

      nx = (int) ( (r_cr[0] - cc[0])*scale_tmp + 0.5);
      ny = (int) ( (r_cr[1] - cc[1])*scale_tmp + 0.5);
      nz = (int) ( (r_cr[2] - cc[2])*scale_tmp + 0.5);

      //refined stretched box
      l_cr[0] = cc[0] - nx*1.0/scale_tmp;
      l_cr[1] = cc[1] - ny*1.0/scale_tmp;
      l_cr[2] = cc[2] - nz*1.0/scale_tmp;
      r_cr[0] = cc[0] + nx*1.0/scale_tmp;
      r_cr[1] = cc[1] + ny*1.0/scale_tmp;
      r_cr[2] = cc[2] + nz*1.0/scale_tmp;
    }

    double dist = ((rr[0]-ll[0]) - lmax)*0.5;
    pot_bc = std::exp (-k*dist)/ (dist*eps_out);

    if (rank == 0) {
      std::cout << "========== [ Domain Information ] ==========\n";
      std::cout << "  Scale:  " << scale << "\n";

      std::cout << "  Center of the System [Å]:";
      std::cout << "  [" << cc[0] << ", " << cc[1] << ", " << cc[2] << "]\n";

      std::cout << "  Perfil box:  " << perfil1 << "\n";
      std::cout << "  Complete Domain Box Size [Å]:\n";
      std::cout << "      x = [" << ll[0] << ", " << rr[0] << "]\n";
      std::cout << "      y = [" << ll[1] << ", " << rr[1] << "]\n";
      std::cout << "      z = [" << ll[2] << ", " << rr[2] << "]\n";

      std::cout << "  Number of Subdivisions:";
      std::cout << "  nx = " << nx * 2 << "  ny = " << ny * 2 << "  nz = " << nz * 2 << '\n';

      std::cout << "============================================\n";
    }

    simple_conn_num_vertices = 8;
    simple_conn_num_trees = 1;

    simple_conn_p = std::make_unique<double[]> (simple_conn_num_vertices*3);
    simple_conn_t = std::make_unique<p4est_topidx_t[]> (simple_conn_num_vertices);

    auto tmp_p = {ll[0], ll[1], ll[2],
                  rr[0], ll[1], ll[2],
                  ll[0], rr[1], ll[2],
                  rr[0], rr[1], ll[2],
                  ll[0], ll[1], rr[2],
                  rr[0], ll[1], rr[2],
                  ll[0], rr[1], rr[2],
                  rr[0], rr[1], rr[2]
                 };
    auto tmp_t = {1, 2, 3, 4, 5, 6, 7, 8};

    std::copy (tmp_p.begin (), tmp_p.end (), simple_conn_p.get ());
    std::copy (tmp_t.begin (), tmp_t.end (), simple_conn_t.get ());

    for (int i = 0; i<6; i++)
      bcells.push_back (std::make_pair (0, i));
  } else if (mesh_shape == 2) {

    l_cr[0] = l_c[0];
    l_cr[1] = l_c[1];
    l_cr[2] = l_c[2];
    r_cr[0] = r_c[0];
    r_cr[1] = r_c[1];
    r_cr[2] = r_c[2];
    scale = (1<<unilevel)/ (r_c[0]-l_c[0]);

    double dist = ((r_cr[0]-l_cr[0]) - lmax)*0.5;
    pot_bc = std::exp (-k*dist)/ (dist*eps_out);

    if (rank == 0) {
      std::cout << "========== [ Domain Information ] ==========\n";

      std::cout << "  Scale:  " << scale << "\n";

      std::cout << "  Center of the System [Å]:";
      std::cout << "  [" << (r_c[0] + l_c[0])*0.5 << ", "
                << (r_c[1] + l_c[1])*0.5 << ", "
                << (r_c[2] + l_c[2])*0.5 << "]\n";

      std::cout << "  Complete Domain Box Size [Å]:\n";
      std::cout << "      x = [" << l_cr[0] << ", " << r_cr[0] << "]\n";
      std::cout << "      y = [" << l_cr[1] << ", " << r_cr[1] << "]\n";
      std::cout << "      z = [" << l_cr[2] << ", " << r_cr[2] << "]\n";

      std::cout << "  Number of Subdivisions:";
      std::cout << "  nx = " << (1<<unilevel) << "  ny = " << (1<<unilevel) << "  nz = " << (1<<unilevel) << '\n';

      std::cout << "============================================\n";
    }

    if (refine_box == 1) {
      double nx_tmp, ny_tmp, nz_tmp;
      double cc_tmp[3];

      scale_x = (1<<unilevel)/ (r_c[0]-l_c[0]);
      scale_y = (1<<unilevel)/ (r_c[1]-l_c[1]);
      scale_z = (1<<unilevel)/ (r_c[2]-l_c[2]);

      cc[0] = (r_c[0]+l_c[0])*0.5;
      cc[1] = (r_c[1]+l_c[1])*0.5;
      cc[2] = (r_c[2]+l_c[2])*0.5;

      cc_tmp[0] = (r_cr[0]+l_cr[0])*0.5;
      cc_tmp[1] = (r_cr[1]+l_cr[1])*0.5;
      cc_tmp[2] = (r_cr[2]+l_cr[2])*0.5;

      nx_tmp = (int) ( (cc[0] - cc_tmp[0])*scale_x + 0.5);
      ny_tmp = (int) ( (cc[1] - cc_tmp[1])*scale_y + 0.5);
      nz_tmp = (int) ( (cc[2] - cc_tmp[2])*scale_z + 0.5);

      cc_tmp[0] = cc[0] + nx_tmp*1.0/scale_x;
      cc_tmp[1] = cc[1] + ny_tmp*1.0/scale_y;
      cc_tmp[2] = cc[2] + nz_tmp*1.0/scale_z;

      nx = (int) ( (r_cr[0] - cc_tmp[0])*scale_x + 0.5);
      ny = (int) ( (r_cr[1] - cc_tmp[1])*scale_y + 0.5);
      nz = (int) ( (r_cr[2] - cc_tmp[2])*scale_z + 0.5);
      //refined stretched box
      l_cr[0] = cc_tmp[0] - nx*1.0/scale_x;
      l_cr[1] = cc_tmp[1] - ny*1.0/scale_y;
      l_cr[2] = cc_tmp[2] - nz*1.0/scale_z;
      r_cr[0] = cc_tmp[0] + nx*1.0/scale_x;
      r_cr[1] = cc_tmp[1] + ny*1.0/scale_y;
      r_cr[2] = cc_tmp[2] + nz*1.0/scale_z;

      std::cout << "xb: " << l_cr[0] << ", " << r_cr[0] << std::endl;
      std::cout << "yb: " << l_cr[1] << ", " << r_cr[1] << std::endl;
      std::cout << "zb: " << l_cr[2] << ", " << r_cr[2] << "\n" << std::endl;
    }

    simple_conn_num_vertices = 8;
    simple_conn_num_trees = 1;

    simple_conn_p = std::make_unique<double[]> (simple_conn_num_vertices*3);
    simple_conn_t = std::make_unique<p4est_topidx_t[]> (simple_conn_num_vertices);

    auto tmp_p = {l_c[0], l_c[1], l_c[2],
                  r_c[0], l_c[1], l_c[2],
                  l_c[0], r_c[1], l_c[2],
                  r_c[0], r_c[1], l_c[2],
                  l_c[0], l_c[1], r_c[2],
                  r_c[0], l_c[1], r_c[2],
                  l_c[0], r_c[1], r_c[2],
                  r_c[0], r_c[1], r_c[2]
                 };
    auto tmp_t = {1, 2, 3, 4, 5, 6, 7, 8};

    std::copy (tmp_p.begin (), tmp_p.end (), simple_conn_p.get ());
    std::copy (tmp_t.begin (), tmp_t.end (), simple_conn_t.get ());

    for (int i = 0; i<6; i++)
      bcells.push_back (std::make_pair (0, i));

  } else if (mesh_shape == 3) {
    //cubic box with max perfil2
    double size = 1.0/scale;
    double scale_min_box = 0.5;
    scale_level = 0;
    scale_level_min_box = 0;


    for (int ii = 0; ii<28; ++ii) {
      size *= 2;
      scale_level ++;

      if (lmax/size < perfil2)
        break;
    }

    for (int ii = 0; ii<scale_level; ++ii) {
      scale_level_min_box ++;

      if ( (std::pow (2,scale_level_min_box)+1)/size >= scale_min_box)
        break;
    }

    ll[0] = cc[0] - size/2;
    ll[1] = cc[1] - size/2;
    ll[2] = cc[2] - size/2;
    rr[0] = cc[0] + size/2;
    rr[1] = cc[1] + size/2;
    rr[2] = cc[2] + size/2;

    //refined box NS
    nx = (int) ( (r_cr[0] - cc[0])*scale + 0.5);
    ny = (int) ( (r_cr[1] - cc[1])*scale + 0.5);
    nz = (int) ( (r_cr[2] - cc[2])*scale + 0.5);

    l_cr[0] = cc[0] - nx*1.0/scale;
    l_cr[1] = cc[1] - ny*1.0/scale;
    l_cr[2] = cc[2] - nz*1.0/scale;
    r_cr[0] = cc[0] + nx*1.0/scale;
    r_cr[1] = cc[1] + ny*1.0/scale;
    r_cr[2] = cc[2] + nz*1.0/scale;

    //refined box FOCUS

    double err_l;
    l_box[0] = cc_focusing[0] - std::round (n_grid/2.)/scale;
    l_box[1] = cc_focusing[1] - std::round (n_grid/2.)/scale;
    l_box[2] = cc_focusing[2] - std::round (n_grid/2.)/scale;
    r_box[0] = cc_focusing[0] + std::round (n_grid/2.)/scale;
    r_box[1] = cc_focusing[1] + std::round (n_grid/2.)/scale;
    r_box[2] = cc_focusing[2] + std::round (n_grid/2.)/scale;


    if ( (r_box[0]- l_box[0])>= (r_cr[0] - l_cr[0])) {
      r_box[0] = r_cr[0];
      l_box[0] = l_cr[0];
    } else {
      if (l_box[0] <= l_cr[0]) {
        err_l = l_cr[0]-l_box[0];
        l_box[0] = l_box[0] + err_l;
        r_box[0] = r_box[0] + err_l;
      }

      if (r_box[0] >= r_cr[0]) {
        err_l = r_box[0]-r_cr[0];
        l_box[0] = l_box[0] - err_l;
        r_box[0] = r_box[0] - err_l;
      }
    }

    if ( (r_box[1]- l_box[1])>= (r_cr[1] - l_cr[1])) {
      r_box[1] = r_cr[1];
      l_box[1] = l_cr[1];
    } else {
      if (l_box[1] <= l_cr[1]) {
        err_l = l_cr[1]-l_box[1];
        l_box[1] = l_box[1] + err_l;
        r_box[1] = r_box[1] + err_l;
      }

      if (r_box[1] >= r_cr[1]) {
        err_l = r_box[1]-r_cr[1];
        l_box[1] = l_box[1] - err_l;
        r_box[1] = r_box[1] - err_l;
      }
    }

    if ( (r_box[2]- l_box[2])>= (r_cr[2] - l_cr[2])) {
      r_box[2] = r_cr[2];
      l_box[2] = l_cr[2];
    } else {
      if (l_box[2] <= l_cr[2]) {
        err_l = l_cr[2]-l_box[2];
        l_box[2] = l_box[2] + err_l;
        r_box[2] = r_box[2] + err_l;
      }

      if (r_box[2] >= r_cr[2]) {
        err_l = r_box[2]-r_cr[2];
        r_box[2] = r_box[2] - err_l;
      }
    }

    //calcolo outlevel
    unsigned int ratio_l_b;
    double len_box_foc = n_grid/scale;
    double len_box = rr[0]-ll[0];

    for (int ii = 1; ii<=5; ++ii) {
      len_box /= 2.0;

      if (len_box <= len_box_foc) {
        ratio_l_b =ii;
        break;
      }
    }

    outlevel = ratio_l_b;
    // outlevel = 1;
    ////////////////////////////////////////


    if (rank == 0) {
      std::cout << "========== [ Domain Information ] ==========\n";


      std::cout << "  Center of the System [Å]:";
      std::cout << "  [" << cc[0] << ", " << cc[1] << ", " << cc[2] << "]\n";

      std::cout << "  Center of the focusing [Å]:";
      std::cout << "  [" << cc_focusing[0] << ", " << cc_focusing[1] << ", " << cc_focusing[2] << "]\n";

      std::cout << "  Scale in the box focusing:  " << scale << "\n";

      std::cout << "  Complete Domain Box Size [Å]:\n";
      std::cout << "      x = [" << ll[0] << ", " << rr[0] << "]\n";
      std::cout << "      y = [" << ll[1] << ", " << rr[1] << "]\n";
      std::cout << "      z = [" << ll[2] << ", " << rr[2] << "]\n";

      std::cout << "  Focusing Box Size [Å]:\n";
      std::cout << "      x = [" << l_box[0] << ", " << r_box[0] << "]\n";
      std::cout << "      y = [" << l_box[1] << ", " << r_box[1] << "]\n";
      std::cout << "      z = [" << l_box[2] << ", " << r_box[2] << "]\n";


      std::cout << "  Number of Subdivisions in the focusing Box:";
      std::cout << "  nx = " << n_grid << "  ny = " << n_grid << "  nz = " << n_grid << '\n';

      std::cout << "============================================\n";
    }


    simple_conn_num_vertices = 8;
    simple_conn_num_trees = 1;

    simple_conn_p = std::make_unique<double[]> (simple_conn_num_vertices*3);
    simple_conn_t = std::make_unique<p4est_topidx_t[]> (simple_conn_num_vertices);

    auto tmp_p = {ll[0], ll[1], ll[2],
                  rr[0], ll[1], ll[2],
                  ll[0], rr[1], ll[2],
                  rr[0], rr[1], ll[2],
                  ll[0], ll[1], rr[2],
                  rr[0], ll[1], rr[2],
                  ll[0], rr[1], rr[2],
                  rr[0], rr[1], rr[2]
                 };
    auto tmp_t = {1, 2, 3, 4, 5, 6, 7, 8};

    std::copy (tmp_p.begin (), tmp_p.end (), simple_conn_p.get ());
    std::copy (tmp_t.begin (), tmp_t.end (), simple_conn_t.get ());

    for (int i = 0; i<6; i++)
      bcells.push_back (std::make_pair (0, i));
  } else if (mesh_shape == 4) {

    //cubic box with max perfil2

    double size = 1.0/scale_max;

    scale = scale_max;
    maxlevel = 0;

    for (int ii = 0; ii<28; ++ii) {
      std::cout << size << "  " << maxlevel << std::endl;
      size *= 2;
      maxlevel ++;

      if (lmax/size < perfil2)
        break;
    }

    minlevel = (int) (maxlevel - std::sqrt (scale_max/scale_min));

    unilevel = (int) ( (maxlevel + minlevel)*0.5);
    unilevel = unilevel + 1;
    outlevel = minlevel;
    scale_level = maxlevel;

    ll[0] = cc[0] - size/2;
    ll[1] = cc[1] - size/2;
    ll[2] = cc[2] - size/2;
    rr[0] = cc[0] + size/2;
    rr[1] = cc[1] + size/2;
    rr[2] = cc[2] + size/2;

    //refined box NS
    nx = (int) ( (r_cr[0] - cc[0])*scale + 0.5);
    ny = (int) ( (r_cr[1] - cc[1])*scale + 0.5);
    nz = (int) ( (r_cr[2] - cc[2])*scale + 0.5);

    l_cr[0] = cc[0] - nx*1.0/scale;
    l_cr[1] = cc[1] - ny*1.0/scale;
    l_cr[2] = cc[2] - nz*1.0/scale;
    r_cr[0] = cc[0] + nx*1.0/scale;
    r_cr[1] = cc[1] + ny*1.0/scale;
    r_cr[2] = cc[2] + nz*1.0/scale;


    if (rank == 0) {
      std::cout << "cx: " << cc[0]
                << " , cy: " << cc[1]
                << " , cz: " << cc[2] << std::endl;

      std::cout << "x: " << ll[0] << ", " << rr[0] << std::endl;
      std::cout << "y: " << ll[1] << ", " << rr[1] << std::endl;
      std::cout << "z: " << ll[2] << ", " << rr[2] << std::endl;

      std::cout << minlevel <<" "<<maxlevel << "  " << unilevel <<std::endl;
    }

    simple_conn_num_vertices = 8;
    simple_conn_num_trees = 1;

    simple_conn_p = std::make_unique<double[]> (simple_conn_num_vertices*3);
    simple_conn_t = std::make_unique<p4est_topidx_t[]> (simple_conn_num_vertices);

    auto tmp_p = {ll[0], ll[1], ll[2],
                  rr[0], ll[1], ll[2],
                  ll[0], rr[1], ll[2],
                  rr[0], rr[1], ll[2],
                  ll[0], ll[1], rr[2],
                  rr[0], ll[1], rr[2],
                  ll[0], rr[1], rr[2],
                  rr[0], rr[1], rr[2]
                 };
    auto tmp_t = {1, 2, 3, 4, 5, 6, 7, 8};

    std::copy (tmp_p.begin (), tmp_p.end (), simple_conn_p.get ());
    std::copy (tmp_t.begin (), tmp_t.end (), simple_conn_t.get ());

    for (int i = 0; i<6; i++)
      bcells.push_back (std::make_pair (0, i));
  } else if (mesh_shape == 5) {
    l_cr[0] = cc[0] - len/2;
    l_cr[1] = cc[1] - len/2;
    l_cr[2] = cc[2] - len/2;
    r_cr[0] = cc[0] + len/2;
    r_cr[1] = cc[1] + len/2;
    r_cr[2] = cc[2] + len/2;


    if (unilevel == 0)
      scale = (num_trees[0])/ (len);
    else
      scale = (num_trees[0]* (1<<unilevel))/ (len);

    //////////////////////////////
    num_trees[1] = num_trees[0];
    num_trees[2] = num_trees[0];

    /////////////////////////////
    if (rank == 0) {
      std::cout << "x: " << l_cr[0] << ", " << r_cr[0] << std::endl;
      std::cout << "y: " << l_cr[1] << ", " << r_cr[1] << std::endl;
      std::cout << "z: " << l_cr[2] << ", " << r_cr[2] << "\n" << std::endl;
      std::cout << "Number of trees: " << num_trees[0] << std::endl;
    }

    double bound_x = std::abs (r_cr[0]-l_cr[0]);
    double bound_y = std::abs (r_cr[1]-l_cr[1]);
    double bound_z = std::abs (r_cr[2]-l_cr[2]);
    double step[3] = { bound_x/num_trees[0],
                       bound_y/num_trees[1],
                       bound_z/num_trees[2]
                     };
    double bound[3] = {bound_x, bound_y, bound_z};
    make_connectivity_3d (num_trees, step, simple_conn_p,
                          simple_conn_num_vertices, simple_conn_t,
                          simple_conn_num_trees, bcells);

    for (p4est_topidx_t i =0; i < simple_conn_num_vertices; ++i) {
      p4est_topidx_t j = 0;
      simple_conn_p[3*i + j++] += l_cr[0];
      simple_conn_p[3*i + j++] += l_cr[1];
      simple_conn_p[3*i + j] += l_cr[2];
    }
  } else {
    if (rank == 0) {
      std::cout << "x: " << ll[0] << ", " << rr[0] << std::endl;
      std::cout << "y: " << ll[1] << ", " << rr[1] << std::endl;
      std::cout << "z: " << ll[2] << ", " << rr[2] << "\n" << std::endl;
    }

    simple_conn_num_vertices = 8;
    simple_conn_num_trees = 1;

    simple_conn_p = std::make_unique<double[]> (simple_conn_num_vertices*3);
    simple_conn_t = std::make_unique<p4est_topidx_t[]> (simple_conn_num_vertices);

    auto tmp_p = {ll[0], ll[1], ll[2],
                  rr[0], ll[1], ll[2],
                  ll[0], rr[1], ll[2],
                  rr[0], rr[1], ll[2],
                  ll[0], ll[1], rr[2],
                  rr[0], ll[1], rr[2],
                  ll[0], rr[1], rr[2],
                  rr[0], rr[1], rr[2]
                 };
    auto tmp_t = {1, 2, 3, 4, 5, 6, 7, 8};

    std::copy (tmp_p.begin (), tmp_p.end (), simple_conn_p.get ());
    std::copy (tmp_t.begin (), tmp_t.end (), simple_conn_t.get ());

    for (int i = 0; i<6; i++)
      bcells.push_back (std::make_pair (0, i));
  }

  tmsh.read_connectivity (simple_conn_p.get (), simple_conn_num_vertices,
                          simple_conn_t.get (), simple_conn_num_trees);

}



double
poisson_boltzmann::levelsetfun (double x, double y, double z)
{
  double dist = 0.0;

  for (const NS::Atom& i : atoms) {

    dist += std::exp (surf_param * ( (std::pow (x - i.pos[0], 2) +
                                      std::pow (y - i.pos[1], 2) +
                                      std::pow (z - i.pos[2], 2)) /
                                     (i.radius*i.radius) - 1.0));

    if (dist > 1.5)
      break;

  }

  return dist;
}


double
poisson_boltzmann::is_in_ns_surf (ray_cache_t & ray_cache, double x, double y, double z, int dir)
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);
  double x1 = x;
  double x2 = y;
  double x3 = z;

  if (dir == 0) {
    x1 = y;
    x2 = z;
    x3 = x;
  } else if (dir == 1) {
    x1 = x;
    x2 = z;
    x3 = y;
  }

  crossings_t & ct = ray_cache (x1, x2, dir);

  if (!ct.init && rank != 0) {
    std::array<double, 2> ray = {x1, x2};
    ray_cache.rays[dir].erase (ray);
    return -1.;
  }


  int i = 0;

  if (ct.inters.size () == 0 || x3 < ct.inters[i])
    return 0; //if there are no inters or y_the coord is before the first intersection, the point is outside.

  while (i < ct.inters.size () && x3 >= ct.inters[i]) //go on until the inters is passed
    i++;

  return (i % 2);
}


int
poisson_boltzmann::parse_options (int argc, char **argv)
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);

  GetPot g (argc, argv);

  // =============================
  // 1. Leggi file dei parametri
  // =============================
  if (!g.search ("--prmfile") && !g.search ("--potfile")) {
    if (rank == 0)
      std::cout << "Warning: No parameters file selected, using the default one."
                << "\nTo select one use --prmfile or --potfile option followed by the desired file.\n";
  }

  // Cerca il file, dando precedenza a --prmfile se entrambi sono presenti
  if (g.search ("--prmfile")) {
    optionsfilename = g.next ("../../data/options.prm");
  } else if (g.search ("--potfile")) {
    optionsfilename = g.next ("../../data/options.prm");
  }

  if (rank == 0)
    std::cout << "Selected parameters file: " << optionsfilename << std::endl;

  std::ifstream optionsfile (optionsfilename);

  if (!optionsfile) {
    if (rank == 0)
      std::cerr << "Cannot find the options file" << std::endl;

    return 1;
  }

  GetPot g2 (optionsfilename.c_str());

  // =============================
  // 2. Leggi parametri input/
  // =============================
  const std::string input_section = "input/";
  std::string filename_from_file;

  filetype = g2 ((input_section + "filetype").c_str(), "pqr");
  pqrfilename = g2 ((input_section + "filename").c_str(), "../../data/1CCM.pqr");
  radiusfilename = g2 ((input_section + "radius_file").c_str(), "../../data/radius.siz");
  chargefilename = g2 ((input_section + "charge_file").c_str(), "../../data/charge.crg");
  write_pqr = g2 ((input_section + "write_pqr").c_str(), 0);
  name_pqr = g2 ((input_section + "name_pqr").c_str(), "output.pqr");



  // =============================
  // 3. Override da riga di comando: --pqrfile
  // =============================
  if (g.search ("--pqrfile")) {
    pqrfilename = g.next ("");
    filetype = "pqr"; // Forza il filetype a pqr se viene specificato un file
  }

  if (rank == 0)
    std::cout << "Selected molecule file:   " << pqrfilename << std::endl;

  std::ifstream pqrfile (pqrfilename);

  if (!pqrfile) {
    if (rank == 0) {
      std::cerr << "Cannot find the pqr file" << std::endl;
      return 1;
    }
  }

  const std::string mesh_options = "mesh/";
  maxlevel = g2 ( (mesh_options + "maxlevel").c_str (), 9);
  minlevel = g2 ( (mesh_options + "minlevel").c_str (), 3);
  unilevel = g2 ( (mesh_options + "unilevel").c_str (), 5);
  outlevel = g2 ( (mesh_options + "outlevel").c_str (), 1);
  loc_refinement = g2 ( (mesh_options + "loc_refinement").c_str (), 0);
  mesh_shape = g2 ( (mesh_options + "mesh_shape").c_str (), 1);
  refine_box = g2 ( (mesh_options + "refine_box").c_str (), 0);
  rand_center = g2 ( (mesh_options + "rand_center").c_str (), 0);
  rand_seed = g2 ( (mesh_options + "rand_seed").c_str (), 0);
  aligned = g2 ( (mesh_options + "aligned").c_str (), 0);

  if (mesh_shape < 2) {
    perfil1 = g2 ( (mesh_options + "perfil1").c_str (), 0.8);
    perfil2 = g2 ( (mesh_options + "perfil2").c_str (), 0.2);
    scale = g2 ( (mesh_options + "scale").c_str (), 2.0);
  }

  if (mesh_shape == 2) {
    l_c[0] = g2 ( (mesh_options + "x1").c_str (), -128.0);
    r_c[0] = g2 ( (mesh_options + "x2").c_str (), 128.0);
    l_c[1] = g2 ( (mesh_options + "y1").c_str (), -128.0);
    r_c[1] = g2 ( (mesh_options + "y2").c_str (), 128.0);
    l_c[2] = g2 ( (mesh_options + "z1").c_str (), -128.0);
    r_c[2] = g2 ( (mesh_options + "z2").c_str (), 128.0);
    l_cr[0] = g2 ( (mesh_options + "refine_x1").c_str (), -64.0);
    r_cr[0] = g2 ( (mesh_options + "refine_x2").c_str (), 64.0);
    l_cr[1] = g2 ( (mesh_options + "refine_y1").c_str (), -64.0);
    r_cr[1] = g2 ( (mesh_options + "refine_y2").c_str (), 64.0);
    l_cr[2] = g2 ( (mesh_options + "refine_z1").c_str (), -64.0);
    r_cr[2] = g2 ( (mesh_options + "refine_z2").c_str (), 64.0);
  }

  if (mesh_shape == 3) {
    perfil1 = g2 ( (mesh_options + "perfil1").c_str (), 0.8);
    perfil2 = g2 ( (mesh_options + "perfil2").c_str (), 0.2);
    scale = g2 ( (mesh_options + "scale").c_str (), 2.0);
    cc_focusing[0] = g2 ( (mesh_options + "cx_foc").c_str (), 0.0);
    cc_focusing[1] = g2 ( (mesh_options + "cy_foc").c_str (), 0.0);
    cc_focusing[2] = g2 ( (mesh_options + "cz_foc").c_str (), 0.0);
    n_grid = g2 ( (mesh_options + "n_grid").c_str (), 8);
  }

  if (mesh_shape == 4) {
    perfil1 = g2 ( (mesh_options + "perfil1").c_str (), 0.8);
    perfil2 = g2 ( (mesh_options + "perfil2").c_str (), 0.2);
    scale_min = g2 ( (mesh_options + "scale_min").c_str (), 0.5);
    scale_max = g2 ( (mesh_options + "scale_max").c_str (), 2.0);
  }

  if (mesh_shape == 5) {
    num_trees [0] = g2 ( (mesh_options + "num_trees_x").c_str (), 10);
    num_trees [1] = g2 ( (mesh_options + "num_trees_y").c_str (), 10);
    num_trees [2] = g2 ( (mesh_options + "num_trees_z").c_str (), 10);
    len = g2 ( (mesh_options + "lato").c_str (), 50.0);
    perfil1 = g2 ( (mesh_options + "perfil1").c_str (), 0.8);
  }

  const std::string model_options = "model/";
  linearized = g2 ( (model_options + "linearized").c_str (), 1);
  bc = g2 ( (model_options + "bc_type").c_str (), 1);
  ionic_strength = g2 ( (model_options + "ionic_strength").c_str (), 0.145);
  ion_size = g2 ( (model_options + "ion_size").c_str (), 0.0);

  // Negative values passed silently before: ion_size < 0 gives nu < 0 and
  // D(u) = 1 + nu (cosh u - 1) can vanish in the Newton solver.
  if (ionic_strength < 0.0 || ion_size < 0.0) {
    if (rank == 0)
      std::cerr << "ERROR: ionic_strength = " << ionic_strength << " M, ion_size = "
                << ion_size << " A: both must be >= 0.\n";
    return 1;
  }

  // Second input form, for different sizes: checked below, after the
  // unknown-key test.
  ion_size_pos = g2 ( (model_options + "ion_size_pos").c_str (), 0.0);
  ion_size_neg = g2 ( (model_options + "ion_size_neg").c_str (), 0.0);
  solvent_size = g2 ( (model_options + "solvent_size").c_str (), 0.0);

  e_in = g2 ( (model_options + "molecular_dielectric_constant").c_str (), 2.);
  e_out = g2 ( (model_options + "solvent_dielectric_constant").c_str (), 80.);
  T = g2 ( (model_options + "T").c_str (), 298.15);
  calc_energy = g2 ( (model_options + "calc_energy").c_str (), 2);
  calc_coulombic = g2 ( (model_options + "calc_coulombic").c_str (), 0);
  calc_potential_term = g2 ( (model_options + "calc_potential_terms").c_str (), 0);
  calc_field_term = g2 ( (model_options + "calc_field_terms").c_str (), 0);
  atoms_write = g2 ( (model_options + "atoms_write").c_str (), 0);
  surf_write = g2 ( (model_options + "surf_write").c_str (), 0);
  surf_write = g2 ( (model_options + "surf_write").c_str (), 0);
  map_type = g2 ( (model_options + "map_type").c_str (), "vtu");
  potential_map = g2 ( (model_options + "potential_map").c_str (), 0);
  eps_map = g2 ( (model_options + "eps_map").c_str (), 0);
  dataset_write = g2 ( (model_options + "dataset_write").c_str (), 0);

  // Every key of [model] has been read above. GetPot returns the default for
  // a key it cannot find, so a misspelt key (e.g. ion_sise) would be ignored
  // and the run would go on with the default model: list the keys of [model]
  // that were never read and stop. The other sections are not checked: some
  // of their keys are read only for some mesh_shape values.
  {
    std::vector<std::string> unknown;
    for (const auto & name : g2.unidentified_variables ())
      if (name.compare (0, model_options.size (), model_options) == 0)
        unknown.push_back (name);
    if (!unknown.empty ()) {
      if (rank == 0) {
        std::cerr << "ERROR: unknown key(s) in [model] of " << optionsfilename << ":";
        for (const auto & name : unknown)
          std::cerr << " " << name.substr (model_options.size ());
        std::cerr << "\n       Check the spelling (the valid keys are listed in data/options.prm).\n";
      }
      return 1;
    }
  }

  // Second input form, for different sizes (1:1 salt): cation, anion and
  // solvent, all cube sides (diameters) in A, volumes v = a^3. The three keys
  // go together and exclude ion_size. Checked only now, after the unknown-key
  // test, so that a misspelt size key is reported as unknown, not as missing.
  // vector_variable_size () is 0 for a key that is not in the file.
  {
    const auto given = [&] (const char * key)
      { return g2.vector_variable_size ( (model_options + key).c_str ()) > 0; };
    const bool given_pos = given ("ion_size_pos");
    const bool given_neg = given ("ion_size_neg");
    const bool given_w   = given ("solvent_size");
    per_species_sizes = given_pos || given_neg || given_w;

    if (per_species_sizes) {
      const char * hint = "       Sizes are cube sides (diameters) in A. Give either ion_size alone"
                          " (one size for cations, anions and solvent)\n"
                          "       or all of ion_size_pos, ion_size_neg (0 = point-like ion) and"
                          " solvent_size (3.1 = water;\n"
                          "       equal to the ion sizes for the uniform model).\n";
      if (given ("ion_size")) {
        if (rank == 0)
          std::cerr << "ERROR: ion_size cannot be used together with"
                    << " ion_size_pos, ion_size_neg, solvent_size.\n" << hint;
        return 1;
      }
      if (! (given_pos && given_neg && given_w)) {
        if (rank == 0)
          std::cerr << "ERROR: missing key(s) in [model]:"
                    << (given_pos ? "" : " ion_size_pos") << (given_neg ? "" : " ion_size_neg")
                    << (given_w ? "" : " solvent_size") << "\n" << hint;
        return 1;
      }
      if (ion_size_pos < 0.0 || ion_size_neg < 0.0 || ! (solvent_size > 0.0)) {
        if (rank == 0)
          std::cerr << "ERROR: ion_size_pos = " << ion_size_pos << " A, ion_size_neg = "
                    << ion_size_neg << " A, solvent_size = " << solvent_size
                    << " A: ion sizes must be >= 0 and solvent_size > 0.\n";
        return 1;
      }

      // Bulk free volume theta_w^b = 1 - n_b (v+ + v-) (manuscript eq. s_equation
      // at u = 0) must be positive, as nu < 1 for the uniform model below.
      const double n_b = 1.0e3 * N_av * ionic_strength * Angs * Angs * Angs;
      const double v_p = ion_size_pos * ion_size_pos * ion_size_pos;
      const double v_m = ion_size_neg * ion_size_neg * ion_size_neg;
      const double v_w = solvent_size * solvent_size * solvent_size;
      const double packing = n_b * (v_p + v_m);
      if (! (packing < 1.0)) {
        if (rank == 0)
          std::cerr << "ERROR: ion_size_pos = " << ion_size_pos << " A, ion_size_neg = "
                    << ion_size_neg << " A at ionic_strength = " << ionic_strength
                    << " M: bulk volume fractions v+ n_b = " << v_p * n_b << ", v- n_b = "
                    << v_m * n_b << ",\n       theta_w^b = 1 - v+ n_b - v- n_b = "
                    << 1.0 - packing << " <= 0: no bulk solvent left, reduce the ion sizes.\n";
        return 1;
      }

      // Warnings, only for this form (the runs with ion_size are unchanged).
      if (rank == 0) {
        const std::pair<const char *, double> sizes[3] =
          {{"ion_size_pos", ion_size_pos}, {"ion_size_neg", ion_size_neg},
           {"solvent_size", solvent_size}};
        for (const auto & s : sizes) {
          if (s.second > 15.0)
            std::cerr << "WARNING: " << s.first << " = " << s.second
                      << " A is above 15 A: sizes are cube sides (diameters) in A,"
                      << " not volumes in A^3 or lengths in pm.\n";
          else if (s.second > 0.0 && s.second < 1.0)
            std::cerr << "WARNING: " << s.first << " = " << s.second
                      << " A is below 1 A: sizes are cube sides (diameters) in A, not nm.\n";
        }
        if (packing > 0.5)
          std::cerr << "WARNING: bulk packing fraction 1 - theta_w^b = " << packing
                    << " > 0.5: the ions fill more than half of the bulk,"
                    << " where the lattice-gas model is crude.\n";
        if (ionic_strength == 0.0)
          std::cerr << "WARNING: ionic_strength = 0 (no ions): the ion sizes have no effect.\n";
        else if (ion_size_pos == 0.0 && ion_size_neg == 0.0)
          std::cerr << "WARNING: point-like ions (ion_size_pos = ion_size_neg = 0):"
                    << " ideal (sinh) model, solvent_size has no effect.\n";
      }

      // Model choice, by exact equality of the input values.
      if (ionic_strength == 0.0 || (ion_size_pos == 0.0 && ion_size_neg == 0.0))
        ion_size = 0.0;                    // ideal model
      else if (ion_size_pos == ion_size_neg && ion_size_neg == solvent_size)
        ion_size = ion_size_pos;           // uniform model, same as ion_size = a
      else {
        nonuniform = true;
        if (! ion_model_nu.init (n_b, v_p, v_m, v_w)) {
          if (rank == 0)
            std::cerr << "ERROR: invalid non-uniform ion model (n_b = " << n_b << " A^-3, v+ = "
                      << v_p << ", v- = " << v_m << ", v_w = " << v_w << " A^3).\n";
          return 1;
        }
      }
    }
  }

  // Bulk packing fraction nu = 2 a^3 n_b (1:1 salt, n_b in 1/Angs^3).
  // nu >= 1 leaves no solvent in the bulk (theta_FS = 1 - nu <= 0) and makes
  // g'(u) change sign at large |u|: the lattice-gas model is undefined there.
  {
    const double n_b = 1.0e3 * N_av * ionic_strength * Angs * Angs * Angs;
    ion_model.nu = 2.0 * ion_size * ion_size * ion_size * n_b;
    if (ion_model.nu >= 1.0) {
      if (rank == 0)
        std::cerr << "ERROR: ion_size = " << ion_size << " A at ionic_strength = "
                  << ionic_strength << " M gives packing fraction nu = 2 a^3 n_b = "
                  << ion_model.nu << " >= 1: no bulk solvent left, reduce ion_size.\n";
      return 1;
    }
  }

  const std::string surf_options = "surface/";
  surf_type_num = g2 ( (surf_options + "surface_type").c_str (), 0);

  if (surf_type_num == 1) surf_type = NS::skin;
  else if (surf_type_num == 0) surf_type = NS::ses;
  // else if (surf_type_num == 2) surf_type = NS::blobby;
  else surf_type = NS::ses;

  surf_param = g2 ( (surf_options + "surface_parameter").c_str (), 0.45);
  stern_layer_surf = g2 ( (surf_options + "stern_layer_surf").c_str (), 0);
  stern_layer = g2 ( (surf_options + "stern_layer_thickness").c_str (), 2.);
  num_threads = g2 ( (surf_options + "number_of_threads").c_str (), 1);

  const std::string alg_options = "algorithm/";
  linear_solver_name = g2 ( (alg_options + "linear_solver").c_str (), "lis");
  linear_solver_options = g2 ( (alg_options + "solver_options").c_str (), "-p ssor -ssor_omega 0.51 -i cgs -tol 1.e-6 -print 2 -conv_cond 2 -tol_w 0");
  newton_compress = g2 ( (alg_options + "newton_compress").c_str (), 1);

  const std::string out_options = "output/";
  p4estfilename = g2 ( (out_options + "p4estfilename").c_str (), "poisson_boltzmann_p4est");
  markerfilename = g2 ( (out_options + "markerfilename").c_str (), "poisson_boltzmann_marker_0");
  surffilename = g2 ( (out_options + "surffilename").c_str (), "poisson_boltzmann_surface_0");

  return 0;
}

// ====================================
// Autovettore dominante per matrice 3x3 simmetrica
// Usa il metodo di Jacobi (iterativo, robusto)
// ====================================
void compute_dominant_eigenvector (double cov[3][3], double axis[3])
{
  // Inizializza axis = (1,0,0)
  axis[0] = 1.0;
  axis[1] = 0.0;
  axis[2] = 0.0;

  // Potenza iterativa per il vettore principale
  for (int iter = 0; iter < 20; ++iter) {
    double x = cov[0][0]*axis[0] + cov[0][1]*axis[1] + cov[0][2]*axis[2];
    double y = cov[1][0]*axis[0] + cov[1][1]*axis[1] + cov[1][2]*axis[2];
    double z = cov[2][0]*axis[0] + cov[2][1]*axis[1] + cov[2][2]*axis[2];

    double norm = std::sqrt (x*x + y*y + z*z);

    if (norm < 1e-12) break;

    axis[0] = x / norm;
    axis[1] = y / norm;
    axis[2] = z / norm;
  }
}

void align_atoms_to_Z (std::vector<NS::Atom> &atoms)
{
  if (atoms.empty()) return;

  // 1. Centro geometrico (NON pesato)
  double center[3] = {0.0, 0.0, 0.0};

  for (const auto &a : atoms) {
    center[0] += a.pos[0];
    center[1] += a.pos[1];
    center[2] += a.pos[2];
  }

  center[0] /= atoms.size();
  center[1] /= atoms.size();
  center[2] /= atoms.size();

  // Traslazione
  for (auto &a : atoms) {
    a.pos[0] -= center[0];
    a.pos[1] -= center[1];
    a.pos[2] -= center[2];
  }

  // 2. Matrice di covarianza simmetrica
  double cov[3][3] = {{0.0}};

  for (const auto &a : atoms) {
    cov[0][0] += a.pos[0] * a.pos[0];
    cov[0][1] += a.pos[0] * a.pos[1];
    cov[0][2] += a.pos[0] * a.pos[2];
    cov[1][1] += a.pos[1] * a.pos[1];
    cov[1][2] += a.pos[1] * a.pos[2];
    cov[2][2] += a.pos[2] * a.pos[2];
  }

  cov[1][0] = cov[0][1];
  cov[2][0] = cov[0][2];
  cov[2][1] = cov[1][2];

  // 3. Autovettore principale (metodo potenza)
  double axis[3];
  compute_dominant_eigenvector (cov, axis);

  // Normalizza
  double norm = std::sqrt (axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);

  if (norm > 1e-12) {
    axis[0] /= norm;
    axis[1] /= norm;
    axis[2] /= norm;
  }

  // 4. Rotazione: porta axis su (0,0,1) usando rotazione di Rodrigues
  double z_axis[3] = {0.0, 0.0, 1.0};
  double v[3] = {
    axis[1]*z_axis[2] - axis[2]*z_axis[1],
    axis[2]*z_axis[0] - axis[0]*z_axis[2],
    axis[0]*z_axis[1] - axis[1]*z_axis[0]
  };
  double s = std::sqrt (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  double c = axis[0]*z_axis[0] + axis[1]*z_axis[1] + axis[2]*z_axis[2];

  double R[3][3];

  if (s < 1e-8) {
    R[0][0] = R[1][1] = R[2][2] = 1.0;
    R[0][1] = R[0][2] = R[1][0] = R[1][2] = R[2][0] = R[2][1] = 0.0;
  } else {
    double vx = v[0]/s, vy = v[1]/s, vz = v[2]/s;
    double k = 1.0 - c;
    R[0][0] = c + vx*vx*k;
    R[0][1] = vx*vy*k - vz*s;
    R[0][2] = vx*vz*k + vy*s;
    R[1][0] = vy*vx*k + vz*s;
    R[1][1] = c + vy*vy*k;
    R[1][2] = vy*vz*k - vx*s;
    R[2][0] = vz*vx*k - vy*s;
    R[2][1] = vz*vy*k + vx*s;
    R[2][2] = c + vz*vz*k;
  }

  // 5. Applica rotazione
  for (auto &a : atoms) {
    double x = a.pos[0], y = a.pos[1], z = a.pos[2];
    a.pos[0] = R[0][0]*x + R[0][1]*y + R[0][2]*z;
    a.pos[1] = R[1][0]*x + R[1][1]*y + R[1][2]*z;
    a.pos[2] = R[2][0]*x + R[2][1]*y + R[2][2]*z;
  }
}

// Parse a single PQR atomic line ("ATOM"/"HETATM") robustly, tolerating the
// quirks of fixed-width PQR that break a naive ">>"-based reader:
//   - record name merged with the serial when serial >= 10000 (e.g. "HETATM10936");
//   - coordinates/charges merged together when a negative value overflows its
//     field (e.g. "-0.412-100.217");
//   - chain merged with the residue sequence number (e.g. "E1589").
// Returns true if the line was a well-formed atomic record and 'a' was filled.
static bool
parse_pqr_atom_line (const std::string &line, NS::Atom &a)
{
  // Identify and strip the record name; the serial may be glued to it.
  std::string rest;

  if (line.compare (0, 6, "HETATM") == 0)
    rest = line.substr (6);
  else if (line.compare (0, 4, "ATOM") == 0)
    rest = line.substr (4);
  else
    return false;   // not an atomic record (TER, END, REMARK, MODEL, #, ...)

  // The five trailing floats (x, y, z, charge, radius) are the only tokens with
  // a decimal point. A regex that also matches a leading sign splits values that
  // got glued together on their sign (e.g. "-0.412-100.217").
  static const std::regex float_re ("[-+]?[0-9]+\\.[0-9]+");

  std::vector<std::smatch> matches;
  for (std::sregex_iterator it (rest.begin (), rest.end (), float_re), end;
       it != end; ++it)
    matches.push_back (*it);

  if (matches.size () < 5)
    return false;   // malformed: not enough numeric fields

  // Keep the last five floats; everything before the first of them is the head.
  const std::smatch &first_float = matches[matches.size () - 5];
  std::string head = rest.substr (0, first_float.position ());

  a.pos[0]   = std::stod (matches[matches.size () - 5].str ());
  a.pos[1]   = std::stod (matches[matches.size () - 4].str ());
  a.pos[2]   = std::stod (matches[matches.size () - 3].str ());
  a.charge   = std::stod (matches[matches.size () - 2].str ());
  a.radius   = std::stod (matches[matches.size () - 1].str ());

  if (a.radius < 1.e-5)
    a.radius = 1.0;

  // Head layout: serial name resName [chain] resNum   (chain optional / glued).
  std::istringstream hs (head);
  std::string serial;
  std::vector<std::string> tok;
  std::string t;

  while (hs >> t)
    tok.push_back (t);

  if (tok.size () < 4)
    return false;   // need at least serial, name, resName, resNum

  // serial = tok[0] (discarded, ngpb renumbers atoms itself)
  a.ai.name    = tok[1];
  a.ai.resName = tok[2];

  if (tok.size () >= 5) {
    // chain and resNum given as separate tokens (e.g. "B 4")
    a.ai.chain  = tok[3];
    a.ai.resNum = std::atoi (tok[4].c_str ());
  } else {
    // single trailing token: either a bare resNum ("22194") or chain glued to
    // resNum ("E1589").
    const std::string &last = tok[3];
    bool is_number = !last.empty () &&
                     (std::isdigit (static_cast<unsigned char> (last[0])) ||
                      last[0] == '-' || last[0] == '+');

    if (is_number) {
      a.ai.chain.clear ();
      a.ai.resNum = std::atoi (last.c_str ());
    } else {
      // split leading non-digits (chain) from trailing digits (resNum)
      std::size_t p = 0;
      while (p < last.size () &&
             !std::isdigit (static_cast<unsigned char> (last[p])))
        ++p;

      a.ai.chain  = last.substr (0, p);
      a.ai.resNum = (p < last.size ()) ? std::atoi (last.c_str () + p) : -1;
    }
  }

  return true;
}

void
poisson_boltzmann::read_atoms_from_pqr (std::basic_istream<char> &inputfile)
{
  NS::Atom a;
  atoms.clear ();

  // Line-based, fault-tolerant parsing: a malformed line is skipped (with a
  // warning) instead of silently aborting the whole read, as the previous
  // ">>"-based loop did on the first unexpected token.
  std::string line;
  std::size_t n_read = 0, n_skipped = 0;

  while (std::getline (inputfile, line)) {
    // strip a trailing CR (Windows line endings)
    if (!line.empty () && line.back () == '\r')
      line.pop_back ();

    // only ATOM/HETATM records carry atoms; skip everything else silently
    if (line.compare (0, 6, "HETATM") != 0 &&
        line.compare (0, 4, "ATOM") != 0)
      continue;

    if (parse_pqr_atom_line (line, a)) {
      atoms.push_back (a);
      ++n_read;
    } else {
      ++n_skipped;
      std::cerr << "  [WARNING] skipping malformed PQR atom line: "
                << line << '\n';
    }
  }

  std::cout << "  PQR reader: " << n_read << " atomic records read, "
            << n_skipped << " skipped\n";

  if (aligned == 1) {
    align_atoms_to_Z (atoms);
  }

  if (atoms.size() < 4) {
    auto comp = [] (const NS::Atom &a1, const NS::Atom &a2) -> bool {
      return a1.radius < a2.radius;
    };

    auto max_iter = std::max_element (atoms.begin(), atoms.end(), comp);
    std::array<double, 3> max_pos = {
      max_iter->pos[0],
      max_iter->pos[1],
      max_iter->pos[2]
    };

    const double epsilon = 0.001;

    std::array<std::array<int, 3>, 6> directions = {{
        {{-1, 0, 0}},
        {{+1, 0, 0}},
        {{0, -1, 0}},
        {{0, +1, 0}},
        {{0, 0, -1}},
        {{0, 0, +1}}
      }
    };

    for (int ii = 0; ii < 6; ++ii) {
      std::array<double, 3> new_pos = {
        max_pos[0] + epsilon * directions[ii][0],
        max_pos[1] + epsilon * directions[ii][1],
        max_pos[2] + epsilon * directions[ii][2]
      };
      NS::Atom dummy;
      dummy.pos[0] = new_pos[0];
      dummy.pos[1] = new_pos[1];
      dummy.pos[2] = new_pos[2];
      dummy.charge = 0.0;
      dummy.radius = 0.0;
      atoms.push_back (dummy);
    }
  }
}

void
poisson_boltzmann::read_atoms_from_class ()
{
  static std::array<double,3> pos;

  if (atoms_write == 1) {
    int atom_number = 1;

    for (const NS::Atom& i : atoms) {
      pos[0] = i.pos[0];
      pos[1] = i.pos[1];
      pos[2] = i.pos[2];
      pos_atoms.push_back (pos);
      index_atoms.push_back (atom_number);
      charge_atoms.push_back (i.charge);
      r_atoms.push_back (i.radius);
      atom_number++;
    }
  } else {
    for (const NS::Atom& i : atoms) {
      pos[0] = i.pos[0];
      pos[1] = i.pos[1];
      pos[2] = i.pos[2];
      pos_atoms.push_back (pos);
      charge_atoms.push_back (i.charge);
      r_atoms.push_back (i.radius);
    }
  }
}

void
poisson_boltzmann::broadcast_vectors ()
{
  int size_vec = charge_atoms.size ();
  MPI_Bcast (&size_vec, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Ogni processo alloca il vettore
  pos_atoms.resize (size_vec);
  charge_atoms.resize (size_vec);
  r_atoms.resize (size_vec);
  index_atoms.resize (size_vec);

  // Effettuare il broadcast del vettore
  MPI_Bcast (pos_atoms.data (), size_vec * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (charge_atoms.data (), size_vec, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (r_atoms.data (), size_vec, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (index_atoms.data (), size_vec, MPI_INT, 0, MPI_COMM_WORLD);
}

void
poisson_boltzmann::write_atoms_to_pqr (std::basic_ostream<char> &outputfile)
{
  int Atom_number = 1;

  outputfile << std::setw (10) << std::left << "fieldname" << std::setw (12)
             << std::left <<"Atom_number" << std::setw (12) << std::left << "Atom_name" << std::setw (16) << std::left
             << "Residue_name" << std::setw (16) << std::left << "Residue_number" << std::setw (10) << std::left << "X"
             << std::setw (10) << std::left << "Y" << std::setw (10) << std::left
             << "Z" << std::setw (10) << std::left << "Charge" << std::setw (10) << std::left << "Radius" << std::endl;

  for (auto & ii : atoms)
    outputfile << std::setw (10) << std::left << "ATOM" << std::setw (12)
               << std::left << Atom_number++ << std::setw (12) << std::left << ii.ai.name << std::setw (16) << std::left
               << ii.ai.resName << std::setw (16) << std::left << ii.ai.resNum << std::setw (10) << std::left << ii.pos[0]
               << std::setw (10) << std::left << ii.pos[1] << std::setw (10) << std::left
               << ii.pos[2] << std::setw (10) << std::left << ii.charge << std::setw (10) << std::left << ii.radius << std::endl;

}



std::basic_istream<char>&
operator>> (std::basic_istream<char>& inputfile, NS::Atom &a)
{
  int Atom_number;
  std::string Field_name;

  inputfile >> Field_name
            >> Atom_number
            >> a.ai.name
            >> a.ai.resName;

  std::string token;
  inputfile >> token;

  // Verifica se è un numero (resNum) oppure una stringa (chain)
  bool is_number = !token.empty() &&
                   (std::isdigit (token[0]) || token[0] == '-' || token[0] == '+');

  if (is_number) {
    // Era resNum - metti solo il token indietro nello stream
    for (auto it = token.rbegin(); it != token.rend(); ++it) {
      inputfile.putback (*it);
    }

    a.ai.chain.clear(); // Nessuna catena specificata
  } else {
    // Era chain - leggi il resNum successivo
    a.ai.chain = token;
  }

  // Ora leggi i dati numerici
  inputfile >> a.ai.resNum >> a.pos[0] >> a.pos[1] >> a.pos[2]
            >> a.charge >> a.radius;

  if (a.radius < 1.e-5)
    a.radius = 1.0;

  return inputfile;
}

std::basic_istream<char>&
operator>> (std::basic_istream<char>& inputfile, std::array<float,5> &a)
{
  int Atom_number;
  std::string Field_name;
  std::string name;
  std::string resName;
  int resNum;

  inputfile >> Field_name
            >> Atom_number
            >> name
            >> resName
            >> resNum
            >> a[0] // x_pos
            >> a[1] // y_pos
            >> a[2] // z_pos
            >> a[3] // charge
            >> a[4]; // radius

  if (a[4] < 1.e-5)
    a[4] = 1.0;

  return inputfile;
}

void
poisson_boltzmann::init_tmesh ()
{
  for (auto i = 0; i < unilevel; ++i) {
    tmsh.set_refine_marker (uniform_refinement);
    tmsh.refine (0, 1);
  }
}


void
poisson_boltzmann::init_tmesh_with_refine_scale ()
{

  for (auto i = 0; i < outlevel; ++i) {
    tmsh.set_refine_marker (uniform_refinement);
    tmsh.refine (0, 1);
  }

  auto refinement = [this]
  (tmesh_3d::quadrant_iterator q) -> int {
    int currentlevel = static_cast<int> (q->the_quadrant->level);
    int retval = 0;
    double x1, y1, z1;

    if (currentlevel >= this->scale_level)
      retval = 0;
    else {
      for (int ii = 0; ii < 8; ++ii) {
        if (! q->is_hanging (ii)) {
          x1 = q -> p (0, ii);
          y1 = q -> p (1, ii);
          z1 = q -> p (2, ii);

          if ( (x1 >= this->l_c[0]) && (x1 <= this->r_c[0])
               && (y1 >= this->l_c[1]) && (y1 <= this->r_c[1])
               && (z1 >= this->l_c[2]) && (z1 <= this->r_c[2])) {
            retval = 1;
            break;
          }
        }
      }
    }
    return (retval);
  };
  for (auto i = 0; i < scale_level; ++i) {
    tmsh.set_refine_marker (refinement);
    tmsh.refine (0, 1);
  }
  
}


void
poisson_boltzmann::init_tmesh_with_refine_box_scale ()
{

  for (auto i = 0; i < outlevel; ++i) {
    tmsh.set_refine_marker (uniform_refinement);
    tmsh.refine (0, 1);
  }

  auto refinement = [this]
  (tmesh_3d::quadrant_iterator q) -> int {
    int currentlevel = static_cast<int> (q->the_quadrant->level);
    int retval = 0;
    double x1, y1, z1;

    if (currentlevel >= this->scale_level)
      retval = 0;
    else {
      for (int ii = 0; ii < 8; ++ii) {
        if (! q->is_hanging (ii)) {
          x1 = q -> p (0, ii);
          y1 = q -> p (1, ii);
          z1 = q -> p (2, ii);

          if ( (x1 >= this->l_box[0]) && (x1 <= this->r_box[0])
               && (y1 >= this->l_box[1]) && (y1 <= this->r_box[1])
               && (z1 >= this->l_box[2]) && (z1 <= this->r_box[2])) {
            retval = 1;
            break;
          }

          // else if ((x1 >= this->l_c[0]) && (x1 <= this->r_c[0])
          // && (y1 >= this->l_c[1]) && (y1 <= this->r_c[1])
          // && (z1 >= this->l_c[2]) && (z1 <= this->r_c[2])
          // &&currentlevel < this->scale_level_min_box) {
          // retval = 1;
          // break;
          // }
        }
      }
    }

    return (retval);
  };

  for (auto i = 0; i < scale_level; ++i) {
    tmsh.set_refine_marker (refinement);
    tmsh.refine (0, 1);
  }
}

bool
poisson_boltzmann::is_in (const NS::Atom& i,
                          tmesh_3d::quadrant_iterator q)
{
  double tol = p4esttol * (rr[0]-ll[0]);

  if (mesh_shape == 2)
    tol = p4esttol * (r_c[0]-l_c[0]);

  bool retval = false;
  double l, r, t, b, f, bk;

  l = q->p (0, 0);
  r = q->p (0, 7);

  f = q->p (1, 0);
  bk = q->p (1, 7);

  b = q->p (2, 0);
  t = q->p (2, 7);


  // for (int ii = 1; ii < 8; ++ii) {
  // l = q->p (0, ii) < l ? q->p (0, ii) : l;
  // r = q->p (0, ii) > r ? q->p (0, ii) : r;
  // f = q->p (1, ii) < f ? q->p (1, ii) : f;
  // bk = q->p (1, ii) > bk ? q->p (1, ii) : bk;
  // b = q->p (2, ii) < b ? q->p (2, ii) : b;
  // t = q->p (2, ii) > t ? q->p (2, ii) : t;

  // }

  retval = (i.pos[0] > l - tol) && (i.pos[0] <= r - tol); //make sure that the charge is assigned only once
  retval = retval && (i.pos[1] > f - tol) && (i.pos[1] <= bk - tol);
  retval = retval && (i.pos[2] > b - tol) && (i.pos[2] <= t - tol);

  return retval;
}

bool
poisson_boltzmann::is_in_ref (const NS::Atom& i,
                              tmesh_3d::quadrant_iterator q)
{

  bool retval = false;
  double l, r, t, b, f, bk;

  l = q->p (0, 0);
  r = q->p (0, 7);

  f = q->p (1, 0);
  bk = q->p (1, 7);

  b = q->p (2, 0);
  t = q->p (2, 7);


  // for (int ii = 1; ii < 8; ++ii) {
  // l = q->p (0, ii) < l ? q->p (0, ii) : l;
  // r = q->p (0, ii) > r ? q->p (0, ii) : r;
  // f = q->p (1, ii) < f ? q->p (1, ii) : f;
  // bk = q->p (1, ii) > bk ? q->p (1, ii) : bk;
  // b = q->p (2, ii) < b ? q->p (2, ii) : b;
  // t = q->p (2, ii) > t ? q->p (2, ii) : t;

  // }

  retval = (i.pos[0] >= l) && (i.pos[0] <= r);
  retval = retval && (i.pos[1] >= f) && (i.pos[1] <= bk);
  retval = retval && (i.pos[2] >= b) && (i.pos[2] <= t);

  return retval;
}

void
poisson_boltzmann::refine_surface (ray_cache_t & ray_cache)
{
  int rank, size;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);

  int num_cycles = 2;

  if (size == 1 || surf_type_num == 2)
    num_cycles = 1;

  int coars_ref_cycles = (maxlevel - unilevel) > (unilevel - minlevel) ? (maxlevel - unilevel) : (unilevel - minlevel);

  for (int kk = 0; kk < coars_ref_cycles; ++kk) {
    // REFINEMENT
    {
      distributed_vector rcoeff (tmsh.num_owned_nodes ());

      for (int jj = 0; jj < num_cycles; jj++) {
        ray_cache.num_req_rays[0] = 0; //zero at each ref/coarsen cycle
        ray_cache.num_req_rays[1] = 0; //zero at each ref/coarsen cycle
        ray_cache.num_req_rays[2] = 0; //zero at each ref/coarsen cycle
        ray_cache.rays_list[0].clear ();
        ray_cache.rays_list[1].clear ();
        ray_cache.rays_list[2].clear ();

        if (rank == 0 && jj == 0)
          std::cout << "Refinement: " << kk << std::endl;

        for (auto quadrant = tmsh.begin_quadrant_sweep ();
             quadrant != tmsh.end_quadrant_sweep ();
             ++quadrant) {

          for (int ii = 0; ii < 8; ++ii) {

            if (! quadrant->is_hanging (ii)) {
              if (surf_type_num == 2)
                rcoeff[quadrant->gt (ii)] = levelsetfun (quadrant->p (0, ii),
                                            quadrant->p (1, ii),
                                            quadrant->p (2, ii));
              else {
                for (int idir = 0; idir < 3; ++idir) {
                  rcoeff[quadrant->gt (ii)] = is_in_ns_surf (ray_cache,
                                              quadrant->p (0, ii),
                                              quadrant->p (1, ii),
                                              quadrant->p (2, ii), idir);

                  if (rcoeff[quadrant->gt (ii)] < -0.5) {
                    ray_cache.num_req_rays[idir]++;

                    std::array<double, 2> ray;

                    std::vector<int> direzioni = {0,1,2};
                    direzioni.erase (direzioni.begin ()+idir);

                    for (unsigned i = 0; i < direzioni.size (); ++i) {
                      ray[i] = quadrant->p (direzioni[i], ii);
                    }

                    ray_cache.rays_list[idir].insert (ray);
                  }
                }
              }
            } else {
              for (int idir = 0; idir < 3; ++idir) {
                double pippo = is_in_ns_surf (ray_cache,
                                              quadrant->p (0, ii),
                                              quadrant->p (1, ii),
                                              quadrant->p (2, ii), idir);

                if (pippo < -0.5) {
                  ray_cache.num_req_rays[idir]++;

                  std::array<double, 2> ray;

                  std::vector<int> direzioni = {0,1,2};
                  direzioni.erase (direzioni.begin ()+idir);

                  for (unsigned i = 0; i < direzioni.size (); ++i) {
                    ray[i] = quadrant->p (direzioni[i], ii);
                  }

                  ray_cache.rays_list[idir].insert (ray);
                }
              }

              for (int jj = 0; jj < quadrant->num_parents (ii); ++jj) {
                rcoeff[quadrant->gparent (jj, ii)] += 0.;
              }
            }
          }
        }

        MPI_Barrier (mpicomm);
        ray_cache.fill_cache ();
      }

      auto refinement = [&rcoeff,this]
      (tmesh_3d::quadrant_iterator q) -> int {
        int currentlevel = static_cast<int> (q->the_quadrant->level);
        int retval = 0;
        double min = 2.0;
        double max = 0.0;
        double tmp = 0.0;

        if (currentlevel >= this->maxlevel)
          retval = 0;
        else {
          for (int ii = 0; ii < 8; ++ii) {

            if (! q->is_hanging (ii)) {
              tmp = rcoeff[q->gt (ii)];

              if (tmp > max) max = tmp;

              if (tmp < min) min = tmp;
            }

          }

          if (this->surf_type_num == 2) {
            if (max > 1 && min < 1)
              retval = this->maxlevel - currentlevel;
            else
              for (const NS::Atom& i : atoms)
                if (is_in_ref (i, q)) {
                  retval = this->maxlevel - currentlevel;
                  break;
                }
          } else {
            if (max > 0.5 && min < 0.5)
              retval = this->maxlevel - currentlevel;
            else
              for (const NS::Atom& i : atoms)
                if (is_in_ref (i, q)) {
                  retval = this->maxlevel - currentlevel;
                  break;
                }
          }
        }

        return (retval);
      };

      tmsh.set_refine_marker (refinement);
      tmsh.refine (0, 1);
    }

    // COARSENING
    {
      distributed_vector rcoeff (tmsh.num_owned_nodes ());

      for (int jj = 0; jj < num_cycles; jj++) {
        ray_cache.num_req_rays[0] = 0; //zero at each ref/coarsen cycle
        ray_cache.num_req_rays[1] = 0; //zero at each ref/coarsen cycle
        ray_cache.num_req_rays[2] = 0; //zero at each ref/coarsen cycle
        ray_cache.rays_list[0].clear ();
        ray_cache.rays_list[1].clear ();
        ray_cache.rays_list[2].clear ();

        if (rank == 0 && jj == 0)
          std::cout << "Coarsening: " << kk << std::endl;

        for (auto quadrant = tmsh.begin_quadrant_sweep ();
             quadrant != tmsh.end_quadrant_sweep ();
             ++quadrant) {

          for (int ii = 0; ii < 8; ++ii) {

            if (! quadrant->is_hanging (ii)) {
              if (surf_type_num == 2)
                rcoeff[quadrant->gt (ii)] = levelsetfun (quadrant->p (0, ii),
                                            quadrant->p (1, ii),
                                            quadrant->p (2, ii));
              else {
                for (int idir = 0; idir < 3; ++idir) {
                  rcoeff[quadrant->gt (ii)] = is_in_ns_surf (ray_cache,
                                              quadrant->p (0, ii),
                                              quadrant->p (1, ii),
                                              quadrant->p (2, ii), idir);

                  if (rcoeff[quadrant->gt (ii)] < -0.5) {
                    ray_cache.num_req_rays[idir]++;

                    std::array<double, 2> ray;

                    std::vector<int> direzioni = {0,1,2};
                    direzioni.erase (direzioni.begin ()+idir);

                    for (unsigned i = 0; i < direzioni.size (); ++i) {
                      ray[i] = quadrant->p (direzioni[i], ii);
                    }

                    // ray_cache.rays[idir].erase(ray);
                    ray_cache.rays_list[idir].insert (ray);
                  }
                }
              }
            } else {
              for (int idir = 0; idir < 3; ++idir) {
                double pippo = is_in_ns_surf (ray_cache,
                                              quadrant->p (0, ii),
                                              quadrant->p (1, ii),
                                              quadrant->p (2, ii), idir);

                if (pippo < -0.5) {
                  ray_cache.num_req_rays[idir]++;

                  std::array<double, 2> ray;

                  std::vector<int> direzioni = {0,1,2};
                  direzioni.erase (direzioni.begin ()+idir);

                  for (unsigned i = 0; i < direzioni.size (); ++i) {
                    ray[i] = quadrant->p (direzioni[i], ii);
                  }

                  ray_cache.rays_list[idir].insert (ray);
                }
              }

              for (int jj = 0; jj < quadrant->num_parents (ii); ++jj) {
                rcoeff[quadrant->gparent (jj, ii)] += 0.;
              }
            }
          }
        }

        MPI_Barrier (mpicomm);
        ray_cache.fill_cache ();
      }

      auto coarsening = [&rcoeff,this]
      (tmesh_3d::quadrant_iterator q) -> int {
        int currentlevel = static_cast<int> (q->the_quadrant->level);
        int retval = 0;
        double min = 2.0;
        double max = 0.0;
        double tmp = 0.0;

        if (currentlevel <= this->minlevel)
          retval = 0;
        else {
          for (int ii = 0; ii < 8; ++ii) {

            if (! q->is_hanging (ii)) {
              tmp = rcoeff[q->gt (ii)];

              if (tmp > max) max = tmp;

              if (tmp < min) min = tmp;
            }

          }

          if (this->surf_type_num == 2) {
            if (min > 1 || max < 1)
              retval = currentlevel - this->minlevel;
          } else {
            if (min > 0.5 || max < 0.5)
              retval = currentlevel - this->minlevel;
          }

          for (const NS::Atom& i : atoms)
            if (is_in_ref (i, q)) {
              retval = 0;
              break;
            }
        }

        return (retval);
      };

      tmsh.set_coarsen_marker (coarsening);
      tmsh.coarsen (0, 1);
    }
  }
}



void
stern_grid_t::build (const std::vector<std::array<double,3>> &pos,
                     const std::vector<double> &rad, double s)
{
  const std::size_t na = pos.size ();
  at.clear ();
  start.assign (1, 0);
  n = {0, 0, 0};

  if (na == 0)
    return;

  L = *std::max_element (rad.begin (), rad.end ()) + s;

  for (int d = 0; d < 3; ++d) {
    lo[d] = hi[d] = pos[0][d];

    for (const auto &p : pos) {
      lo[d] = std::min (lo[d], p[d]);
      hi[d] = std::max (hi[d], p[d]);
    }

    lo[d] -= L;
    hi[d] += L;
    n[d] = std::max (1, static_cast<int> (std::ceil ((hi[d] - lo[d]) / L)));
  }

  auto cell_of = [&] (const std::array<double,3> &p) {
    int c[3];

    for (int d = 0; d < 3; ++d)
      c[d] = std::min (n[d] - 1, static_cast<int> ((p[d] - lo[d]) / L));

    return (c[0] * n[1] + c[1]) * n[2] + c[2];
  };

  // Counting sort of the atoms by cell (CSR layout).
  start.assign (static_cast<std::size_t> (n[0]) * n[1] * n[2] + 1, 0);
  std::vector<int> cell (na);

  for (std::size_t i = 0; i < na; ++i) {
    cell[i] = cell_of (pos[i]);
    ++start[cell[i] + 1];
  }

  for (std::size_t c = 1; c < start.size (); ++c)
    start[c] += start[c - 1];

  at.resize (na);
  std::vector<int> fill (start.begin (), start.end () - 1);

  for (std::size_t i = 0; i < na; ++i) {
    const double r = rad[i] + s;
    at[fill[cell[i]]++] = {pos[i][0], pos[i][1], pos[i][2], r * r};
  }
}

bool
stern_grid_t::inside (double x, double y, double z) const
{
  if (at.empty () || x < lo[0] || x > hi[0] || y < lo[1] || y > hi[1]
      || z < lo[2] || z > hi[2])
    return false;

  const int c0 = std::min (n[0] - 1, static_cast<int> ((x - lo[0]) / L));
  const int c1 = std::min (n[1] - 1, static_cast<int> ((y - lo[1]) / L));
  const int c2 = std::min (n[2] - 1, static_cast<int> ((z - lo[2]) / L));

  for (int i = std::max (0, c0 - 1); i <= std::min (n[0] - 1, c0 + 1); ++i)
    for (int j = std::max (0, c1 - 1); j <= std::min (n[1] - 1, c1 + 1); ++j)
      for (int k = std::max (0, c2 - 1); k <= std::min (n[2] - 1, c2 + 1); ++k) {
        const int c = (i * n[1] + j) * n[2] + k;

        for (int a = start[c]; a < start[c + 1]; ++a) {
          const double dx = x - at[a][0];
          const double dy = y - at[a][1];
          const double dz = z - at[a][2];

          if (dx * dx + dy * dy + dz * dz < at[a][3])
            return true;
        }
      }

  return false;
}


void
poisson_boltzmann::create_markers (ray_cache_t & ray_cache)
{

  int size, rank;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);

  double eps_in = 4.0*pi*e_0*e_in*kb*T*Angs/ (e*e); //adim e_in
  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/ (e*e); //adim e_out

  /////////////////////////////////////////////////////////
  //reactions
  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/ (e_0*e_out*kb*T);

  this->marker.assign (this->tmsh.num_local_quadrants (), 0.0); //marker = 0 -> in

  // Stern layer: reaction_nodes = 0 also on the nodes inside the spheres
  // R_i + stern_layer (see stern_grid_t). It is nodal like the molecule
  // interior, so the linear and Newton solvers, the excess energy and
  // Gamma+- all see it. Only tested in the last cycle, when every node
  // inside the molecule is already known (and has reaction_nodes = 0).
  stern_grid_t stern;

  if (stern_layer_surf == 1) {
    stern.build (pos_atoms, r_atoms, stern_layer);
    std::vector<double> ().swap (r_atoms);
  }

  epsilon_nodes = std::make_unique<distributed_vector> (tmsh.num_owned_nodes (),mpicomm);
  epsilon_nodes->get_owned_data ().assign (tmsh.num_owned_nodes (), eps_out);

  reaction_nodes = std::make_unique<distributed_vector> (tmsh.num_owned_nodes (),mpicomm);
  reaction_nodes->get_owned_data ().assign (tmsh.num_owned_nodes (), eps_out*k2);

  int num_cycles = 2;

  if (size == 1) {
    num_cycles = 1;
  }

  int local_num;

  for (int jj = 0; jj < num_cycles; ++jj) {
    ray_cache.num_req_rays[0] = 0; //zero at each ref/coarsen cycle
    ray_cache.num_req_rays[1] = 0; //zero at each ref/coarsen cycle
    ray_cache.num_req_rays[2] = 0; //zero at each ref/coarsen cycle
    ray_cache.rays_list[0].clear ();
    ray_cache.rays_list[1].clear ();
    ray_cache.rays_list[2].clear ();


    for (auto quadrant = this->tmsh.begin_quadrant_sweep ();
         quadrant != this->tmsh.end_quadrant_sweep ();
         ++quadrant) {
      int num_int_nodes = 0;
      int num_hanging[3] = {0, 0, 0};
      double x,y,z;

      const bool test_stern = stern_layer_surf == 1 && (jj != 0 || num_cycles == 1)
                              && stern.box_near ({quadrant->p (0, 0), quadrant->p (1, 0), quadrant->p (2, 0)},
                                                 {quadrant->p (0, 7), quadrant->p (1, 7), quadrant->p (2, 7)});

      for (int ii = 0; ii < 8; ++ii) {
        local_num =quadrant->gt (ii);
        x = quadrant->p (0, ii);
        y = quadrant->p (1, ii);
        z = quadrant->p (2, ii);

        if (! quadrant->is_hanging (ii)) {
          const bool in_mol = this->is_in_ns_surf (ray_cache, x, y, z, 2) > 0.5;

          if (in_mol) { //inside the molecule
            ++num_int_nodes;
            (*epsilon_nodes)[local_num] = eps_in;
            (*reaction_nodes)[local_num] = 0.0;
          } else if (this->is_in_ns_surf (ray_cache, x, y, z, 2) < -0.5) {
            ray_cache.num_req_rays[2]++;

            // std::array<double, 2> ray {x,y};

            ray_cache.rays_list[2].insert ({x,y});
          }

          if (this->is_in_ns_surf (ray_cache, x, y, z, 0) < -0.5) {
            ray_cache.num_req_rays[0]++;

            // std::array<double, 2> ray {y,z};

            ray_cache.rays_list[0].insert ({y,z});
          }

          if (this->is_in_ns_surf (ray_cache, x, y, z, 1) < -0.5) {
            ray_cache.num_req_rays[1]++;

            std::array<double, 2> ray {x,z};

            ray_cache.rays_list[1].insert ({x,z});
          }

          if (test_stern && !in_mol && stern.inside (x, y, z)) //inside the stern layer
            (*reaction_nodes)[local_num] = 0.0;

        } else
          for (int idir = 0; idir < 3; ++idir) {
            ++num_hanging[idir];

            if (this->is_in_ns_surf (ray_cache, x, y, z, idir) < -0.5) {
              ray_cache.num_req_rays[idir]++;
              std::array<double, 2> ray;

              std::vector<int> direzioni {0,1,2};
              direzioni.erase (direzioni.begin ()+idir);

              for (unsigned i = 0; i < direzioni.size (); ++i) {
                ray[i] = quadrant->p (direzioni[i], ii);
              }

              ray_cache.rays_list[idir].insert (ray);

            }
          }
      }

      if (jj != 0 || num_cycles == 1) {
        if (num_int_nodes == 0) { //if there's no node inside the molecule
          this->marker[quadrant->get_forest_quad_idx ()] = 1.0; //quadrant is out
        } else if (num_int_nodes < (8 - num_hanging[2])) { //if the non hanging nodes are not all inside
          this->marker[quadrant->get_forest_quad_idx ()] = 1.0/2.0; //"border"
          border_quad.push_back (quadrant->get_forest_quad_idx ());
        }

        //else: all the nodes are inside: the quadrant is inside and the marker value is 0
      }

    }

    MPI_Barrier (mpicomm);
    ray_cache.fill_cache ();
  }

  if (size >1) {
    bim3a_solution_with_ghosts (tmsh, *epsilon_nodes, replace_op);

    // reaction_nodes gets the ghost entries of epsilon_nodes (new_node_vector)
    // instead of a second bim3a_solution_with_ghosts, whose sweep over the
    // mesh and its neighbours costs seconds on large meshes. The old vector
    // is freed first, so the two never coexist.
    std::vector<double> c;
    c.swap (reaction_nodes->get_owned_data ());
    reaction_nodes.reset ();
    reaction_nodes = new_node_vector ();
    reaction_nodes->get_owned_data ().swap (c);
    reaction_nodes->assemble (replace_op);
  }

  if (stern_layer_surf == 1 && k2 > 0.0) {
    // Stern nodes: solvent (eps_out) without ions.
    long loc = 0, glob = 0;
    const auto & ed = epsilon_nodes->get_owned_data ();
    const auto & cd = reaction_nodes->get_owned_data ();

    for (std::size_t i = 0; i < cd.size (); ++i)
      loc += (cd[i] == 0.0 && ed[i] == eps_out);

    MPI_Reduce (&loc, &glob, 1, MPI_LONG, MPI_SUM, 0, mpicomm);

    if (rank == 0)
      std::cout << "  Stern layer: ion-free spheres R_i + " << stern_layer
                << " A, " << glob << " solvent nodes without ions\n";
  }
}

std::unique_ptr<distributed_vector>
poisson_boltzmann::new_node_vector ()
{
  int size;
  MPI_Comm_size (mpicomm, &size);

  if (size > 1)
    return std::make_unique<distributed_vector> (*epsilon_nodes);

  return std::make_unique<distributed_vector> (tmsh.num_owned_nodes (), mpicomm);
}

void
poisson_boltzmann::create_density_map (ray_cache_t & ray_cache)
{
  int size, rank;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);

  // ------------------------------------------------------------
  // Early exit: if rho_fixed is already initialized, skip the whole procedure.
  // The density map must be created only once.
  // ------------------------------------------------------------
  if (this->rho_fixed) {
    if (rank == 0)
      std::cout << "[INFO] Density map already initialized. Skipping.\n";

    return;
  }

  // ------------------------------------------------------------
  // Allocate and initialize nodal density vector (rho)
  // ------------------------------------------------------------
  // Ghost entries from new_node_vector; the exchange sets them to zero
  // before the charges are accumulated below.
  this->rho_fixed = new_node_vector ();
  this->rho_fixed->get_owned_data().assign (tmsh.num_owned_nodes(), 0.0);

  if (size > 1)
    this->rho_fixed->assemble (replace_op);

  const int rho_non_local = this->rho_fixed->non_local_size ();

  // ------------------------------------------------------------
  // Vector of ones at nodes
  // ------------------------------------------------------------
  this->ones = new_node_vector ();
  this->ones->get_owned_data().assign (tmsh.num_owned_nodes(), 1.0);

  // ------------------------------------------------------------
  // Per-cell constant field used for computing patch volumes
  // ------------------------------------------------------------
  this->const_ones.assign (tmsh.num_local_quadrants(), 1.0);

  // ------------------------------------------------------------
  // vol_patch will store the nodal patch volumes obtained via BIM
  // ------------------------------------------------------------
  std::unique_ptr<distributed_vector> vol_patch =
    std::make_unique<distributed_vector> (tmsh.num_owned_nodes(), mpicomm);

  // ------------------------------------------------------------
  // Update ghost values if running in parallel
  // ------------------------------------------------------------
  if (size > 1)
    ones->assemble (replace_op);

  // ------------------------------------------------------------
  // Compute patch volumes by integrating the constant field = 1
  // ------------------------------------------------------------
  bim3a_rhs (tmsh, const_ones, *ones, *vol_patch);

  if (size > 1)
    vol_patch->assemble();

  // ------------------------------------------------------------
  // Locate the atom positions inside the adaptive mesh
  // ------------------------------------------------------------
  search_points();

  // ------------------------------------------------------------
  // Distribute the atomic charges to the surrounding nodes
  // using a linear volumetric approximation.
  // ------------------------------------------------------------
  for (auto it = lookup_table.begin(); it != lookup_table.end(); ++it) {

    // Compute cell volume (Cartesian and axis-aligned)
    double volume =
      (it->second.p (0, 7) - it->second.p (0, 0)) *
      (it->second.p (1, 7) - it->second.p (1, 0)) *
      (it->second.p (2, 7) - it->second.p (2, 0));

    for (int ii = 0; ii < 8; ++ii) {

      // Linear interpolation weight based on opposite corner distances
      double weight = std::abs (
                        (pos_atoms[it->first][0] - it->second.p (0, 7 - ii)) *
                        (pos_atoms[it->first][1] - it->second.p (1, 7 - ii)) *
                        (pos_atoms[it->first][2] - it->second.p (2, 7 - ii))) / volume;

      // Regular node (not hanging)
      if (!it->second.is_hanging (ii)) {
        (*rho_fixed)[it->second.gt (ii)] +=
          charge_atoms[it->first] * 4.0 * pi * weight /
          (*vol_patch)[it->second.gt (ii)];
      }
      // Hanging node → distribute to parents
      else {
        for (int jj = 0; jj < it->second.num_parents (ii); ++jj) {
          double denom =
            it->second.num_parents (ii) *
            (*vol_patch)[it->second.gparent (jj, ii)];
          (*rho_fixed)[it->second.gparent (jj, ii)] +=
            charge_atoms[it->first] * 4.0 * pi * weight / denom;
        }
      }
    }
  }

  // Release temporary vector
  vol_patch.reset();

  // ------------------------------------------------------------
  // Sync ghost nodes of rho_fixed in parallel runs: the ghost contributions
  // are summed into their owners, then sent back to the ghosts. If a charge
  // touched a node outside the known ghost entries (the map grew), rebuild
  // them with the full bim3a_solution_with_ghosts so nothing is lost.
  // ------------------------------------------------------------
  if (size > 1) {
    if (rho_fixed->non_local_size () == rho_non_local)
      rho_fixed->assemble ();
    else
      bim3a_solution_with_ghosts (tmsh, *rho_fixed);
  }
}

void
poisson_boltzmann::assemple_system_matrix (ray_cache_t & ray_cache)
{
  int size, rank;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);

  // ------------------------------------------------------------
  // 1) CHECK / BUILD REQUIRED DATA STRUCTURES
  // ------------------------------------------------------------

  // If the density map has not been created yet,
  // the system matrix cannot be assembled. Build it first.
  if (!rho_fixed) {
    create_density_map (ray_cache);
  }

  // Function computing fractional volume intersections (for cut-cells)
  auto func_frac = [&] (tmesh_3d::quadrant_iterator& quadrant) {
    return cube_fraction_intersection (quadrant, ray_cache);
  };

  // Allocate sparse matrix A and RHS vector
  A = std::make_unique<distributed_sparse_matrix> (mpicomm);
  A->set_ranges (tmsh.num_owned_nodes());

  rhs = std::make_unique<distributed_vector> (tmsh.num_owned_nodes(), mpicomm);

  // ------------------------------------------------------------
  // 2) ASSEMBLE THE LINEAR SYSTEM (RHS + STIFFNESS MATRIX)
  // ------------------------------------------------------------

  // Assemble RHS using previously computed fixed charge density
  bim3a_rhs (tmsh, const_ones, *rho_fixed, *rhs);

  // rho_fixed and const_ones are no longer needed after building RHS
  rho_fixed.reset();
  std::vector<double>().swap (const_ones);

  // Assemble Laplace operator (with fractional cell treatment)
  bim3a_laplacian_frac (tmsh, *epsilon_nodes, *A, func_frac);

  // Reaction term, fractional on the cut cells (zero inside the molecule
  // and in the Stern layer)
  bim3a_reaction_frac (tmsh, (*reaction_nodes), *ones, *A, func_frac);

  // Reaction-related vectors are no longer required
  reaction_nodes.reset();
  ones.reset();
  std::vector<double>().swap (marker);

  // ------------------------------------------------------------
  // 3) APPLY BOUNDARY CONDITIONS
  // ------------------------------------------------------------
  // Build Dirichlet boundary list
  dirichlet_bcs3 bcs;

  if (bc == 1) { // Homogeneous Dirichlet BC
    if (std::fabs (pot_bc) > 1.e-5 && rank == 0)
      std::cerr << "[WARNING] Boundary conditions may be inaccurate!!\n";

    for (auto const& ibc : bcells) {
      bcs.emplace_back (ibc.first, ibc.second,
      [] (double, double, double) {
        return 0.0;
      });
    }

    bim3a_dirichlet_bc (tmsh, bcs, *A, *rhs);
  }

  if (bc == 2) { // Coulombic Dirichlet BC
    for (auto const& ibc : bcells) {
      bcs.emplace_back (ibc.first, ibc.second,
      [&] (double x, double y, double z) {
        return coulomb_boundary_conditions (x, y, z);
      });
    }

    bim3a_dirichlet_bc (tmsh, bcs, *A, *rhs);
  }

  if (bc == 3) { // Analytic Dirichlet BC (sphere test case)
    for (auto const& ibc : bcells) {
      bcs.emplace_back (ibc.first, ibc.second,
      [&] (double x, double y, double z) {
        return analytic_solution (x, y, z);
      });
    }

    bim3a_dirichlet_bc (tmsh, bcs, *A, *rhs);
  }

  // ------------------------------------------------------------
  // 4) PARALLEL ASSEMBLY (MPI)
  // ------------------------------------------------------------

  if (size > 1) {
    A->assemble();
    rhs->assemble();
  }
}

// ============================================================
//  Newton assembly: build J(phi) and RHS for ONE Newton step
//  Nonlinear PBE:  -div(eps grad phi) + C * g(phi) = rho_fixed
//    where C = reaction_nodes (== 0 inside the molecule) and g is the
//    ion model (ion_model_t): g = sinh for ideal ions, g = sinh/D for the
//    steric model.
//  "Full" Newton form (BCs identical to the linear case):
//    [ A_stiff + M[C*g'(phi)] ] phi_new
//        = rho_load + M[ C*(phi*g'(phi) - g(phi)) ]
//  NOTE: does NOT free rho_fixed / reaction_nodes / ones / const_ones,
//        because the Newton loop reuses them every iteration.
// ============================================================
void
poisson_boltzmann::assemble_newton_system (ray_cache_t & ray_cache,
                                           distributed_vector & phi_cur)
{
  int size, rank;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);

  // The operators read all 8 local nodes (incl. ghosts): phi_cur's ghost
  // values must be up to date. newton_solve syncs phi after every change.

  // Same cut-cell helper as the linear assembly.
  auto func_frac = [&] (tmesh_3d::quadrant_iterator& quadrant) {
    return cube_fraction_intersection (quadrant, ray_cache);
  };

  // Fresh matrix + rhs for this Newton iteration.
  A = std::make_unique<distributed_sparse_matrix> (mpicomm);
  A->set_ranges (tmsh.num_owned_nodes());
  rhs = std::make_unique<distributed_vector> (tmsh.num_owned_nodes(), mpicomm);



  // --- build nodal Newton coefficients from the current iterate ---
  // C = reaction_nodes (nodal reaction coefficient, 0 inside molecule)
  const std::size_t n = tmsh.num_owned_nodes();

  // Ghost entries from new_node_vector (no sweep over the mesh).
  auto cosh_coeff = new_node_vector (); // C*g'(phi)
  auto rhs_extra  = new_node_vector (); // C*(phi*g'(phi)-g(phi))

  auto & C_data     = reaction_nodes->get_owned_data ();
  auto & phi_data   = phi_cur.get_owned_data ();
  auto & cosh_data  = cosh_coeff->get_owned_data ();
  auto & extra_data = rhs_extra->get_owned_data ();

  for (std::size_t i = 0; i < n; ++i) {
    const double Ci = C_data[i];
    const double ui = phi_data[i];

    // The ionic term lives only where C != 0 (the solvent). Inside the
    // molecule C == 0, and the potential there can be huge near the point
    // charge -> cosh/sinh would overflow and 0*inf = NaN. So skip them.
    if (Ci == 0.0) {
      cosh_data[i]  = 0.0;
      extra_data[i] = 0.0;
      continue;
    }

    const double dg = ion_model.dg (ui);
    const double g  = ion_model.g (ui);
    cosh_data[i]  = Ci * dg;
    extra_data[i] = Ci * (ui * dg - g);
  }

  if (size > 1) {
    cosh_coeff->assemble (replace_op);
    rhs_extra->assemble (replace_op);
  }

   // Add extra RHS term  M_frac * [ C (phi g'(phi) - g(phi)) ].
  // Must use the SAME fractional cut-cell quadrature as the LHS reaction
  // term (bim3a_reaction_frac); the full-cell bim3a_rhs makes the two
  // disagree at cut cells and scales the Newton step by m_full/m_frac.
  bim3a_rhs_frac (tmsh, *ones, *rhs_extra, *rhs, func_frac);

  // Then accumulate the fixed-charge load  b = M * rho_fixed.
  bim3a_rhs (tmsh, const_ones, *rho_fixed, *rhs);

  // --- stiffness  -div(eps grad .), identical to the linear case ---
  bim3a_laplacian_frac (tmsh, *epsilon_nodes, *A, func_frac);

  // --- Jacobian reaction term with coefficient C*g'(phi) ---
  bim3a_reaction_frac (tmsh, *cosh_coeff, *ones, *A, func_frac);

  // --- Dirichlet BCs: identical to assemple_system_matrix ---
  dirichlet_bcs3 bcs;
  if (bc == 1) {
    for (auto const& ibc : bcells)
      bcs.emplace_back (ibc.first, ibc.second,
                        [] (double, double, double) { return 0.0; });
    bim3a_dirichlet_bc (tmsh, bcs, *A, *rhs);
  }
  if (bc == 2) {
    for (auto const& ibc : bcells)
      bcs.emplace_back (ibc.first, ibc.second,
                        [&] (double x, double y, double z) {
                          return coulomb_boundary_conditions (x, y, z); });
    bim3a_dirichlet_bc (tmsh, bcs, *A, *rhs);
  }
  if (bc == 3) {
    for (auto const& ibc : bcells)
      bcs.emplace_back (ibc.first, ibc.second,
                        [&] (double x, double y, double z) {
                          return analytic_solution (x, y, z); });
    bim3a_dirichlet_bc (tmsh, bcs, *A, *rhs);
  }

  if (size > 1) {
    A->assemble();
    rhs->assemble();
  }
}

// ============================================================
//  Newton driver for the nonlinear PBE
//
//    -div(eps grad phi) + C g(phi) = rho_fixed ,   C = reaction_nodes
//
//  g is the ion model (ion_model_t): sinh for ideal ions, sinh/D for the
//  steric model (ion_size > 0). C == 0 inside the molecule, so the
//  equation is LINEAR there: the entire nonlinearity lives in the solvent
//  (C != 0). Convergence of this iteration is therefore governed by the
//  solvent nodes alone.
//
//  Each iteration solves the "full Newton form"
//    [ A_stiff + M[C g'(phi_k)] ] phi_{k+1}
//        = rho_load + M[ C (phi_k g'(phi_k) - g(phi_k)) ]
//  which is algebraically identical to solving for the correction
//  du = phi_{k+1} - phi_k; du is recovered explicitly below so it can be
//  clamped (globalization).
//
//  phi^0 = 0 => g'(0)=1 and the extra RHS term vanishes, so iteration 0
//  reproduces the linear solve exactly (for both ion models). This is
//  deliberate: it is a free regression check against the linear solver.
//
//  Globalization: from iteration 1 onward the step is scaled so that
//  |du|_inf over the ion-accessible nodes (C != 0) is at most maxdu.
//  Iteration 0 is unclamped so the jump to the linear solution is taken whole.
//
//  Stopping: |du|_inf <= newton_tol * max(1, |phi|_inf), both norms over
//  the ion-accessible nodes; or stagnation at the accuracy of the linear
//  solver (du small and no longer decreasing).
//
//  newton_compress = 1 (default): right after iteration 0, i.e. on the linear
//  solution, solvent nodes are mapped phi -> 2 asinh(phi/2). This is the
//  Grahame relation for a 1:1 electrolyte (sigma ~ phi in DH vs
//  sigma ~ 2 sinh(phi/2) nonlinear): same surface charge, nonlinear surface
//  potential. It is ~identity where |phi| << 1 and only tames the thin layer
//  |phi| >~ 1 at the surface, which otherwise costs one Newton iteration per
//  unit of |phi| (the step on sinh is -+1 for large |phi|). Applied ONCE, at
//  iteration 0 only; it is just a better initial guess for iteration 1.
// ============================================================
void
poisson_boltzmann::newton_solve (ray_cache_t & ray_cache)
{
  int size, rank;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);

  const int    newton_max_iter = 100;
  const double newton_tol      = 1.0e-6; // on |du|_inf,solvent, see the test below
  const double newton_stall    = 100.0;  // stagnation test active below newton_stall * tol
  const double maxdu           = 2.0;    // clamp: max |du|_inf per iteration
  const bool   newton_verbose  = true;   // per-iteration diagnostics

  // --- preconditions: fail loudly instead of segfaulting or lying ---
  if (!reaction_nodes) {
    if (rank == 0)
      std::cerr << "  [Newton] ERROR: reaction_nodes not allocated. create_markers() "
                   "must run first, and assemple_system_matrix() must NOT have run "
                   "(it frees reaction_nodes).\n";
    return;
  }
  if (!rho_fixed || !ones) {
    if (rank == 0)
      std::cerr << "  [Newton] ERROR: rho_fixed/ones not allocated. "
                   "create_density_map() must run before newton_solve().\n";
    return;
  }

  const std::size_t n = tmsh.num_owned_nodes ();

  // Initial guess phi^0 = 0.
  phi = new_node_vector ();
  phi->get_owned_data ().assign (n, 0.0);
  if (size > 1)
    phi->assemble (replace_op);

  std::vector<double> phi_old (n, 0.0);
  std::vector<double> du (n, 0.0);

  bool   converged  = false;
  int    it_done    = 0;
  double dsolv_prev = 0.0;   // |du|_inf,solvent of the previous iteration

  for (int it = 0; it < newton_max_iter; ++it) {
    it_done = it;

    // Save the iterate the system is about to be built from.
    {
      auto & phi_cur = phi->get_owned_data ();
      for (std::size_t i = 0; i < n; ++i)
        phi_old[i] = phi_cur[i];
    }

    assemble_newton_system (ray_cache, *phi);

    // Solve J * phi_new = rhs.
    if (linear_solver_name == "mumps")
      mumps_compute_electric_potential (ray_cache);
    else
      lis_compute_electric_potential (ray_cache);

    // IMPORTANT: both solvers REPLACE the phi object (make_unique), so this
    // reference must be taken AFTER the solve, never before.
    auto & phi_new = phi->get_owned_data ();
    auto & Cd      = reaction_nodes->get_owned_data ();

    // --- Newton correction du = phi_new - phi_old ---
    // Norms over all nodes (anorm) and over the ion-accessible nodes, C != 0
    // (dsolv: correction, usolv: iterate).
    double loc[3] = {0.0, 0.0, 0.0};   // anorm, dsolv, usolv
    for (std::size_t i = 0; i < n; ++i) {
      du[i] = phi_new[i] - phi_old[i];
      const double ad = std::fabs (du[i]);
      if (ad > loc[0]) loc[0] = ad;
      if (Cd[i] != 0.0) {
        if (ad > loc[1]) loc[1] = ad;
        const double au = std::fabs (phi_old[i]);
        if (au > loc[2]) loc[2] = au;
      }
    }

    double glob[3] = {loc[0], loc[1], loc[2]};
    if (size > 1)
      MPI_Allreduce (loc, glob, 3, MPI_DOUBLE, MPI_MAX, mpicomm);
    const double anorm = glob[0], dsolv = glob[1], usolv = glob[2];
    const double tol_s = newton_tol * std::max (1.0, usolv);

    // ---------------- diagnostics ----------------
    // Every line tagged "[Newton]" so one grep catches the whole trace.
    // NOTE: rank 0's OWNED nodes only -- indicative, not a global reduction.
    if (newton_verbose && rank == 0) {
      auto report = [&] (const char * tag, const double * v) {
        double m_all = 0.0, m_solv = 0.0, m_mol = 0.0;
        std::size_t i_all = 0, i_solv = 0;
        for (std::size_t i = 0; i < n; ++i) {
          const double a = std::fabs (v[i]);
          if (a > m_all)    { m_all  = a; i_all  = i; }
          if (Cd[i] != 0.0) { if (a > m_solv) { m_solv = a; i_solv = i; } }
          else              { if (a > m_mol)  { m_mol  = a; } }
        }
        std::cout << "  [Newton]   " << tag
                  << " all=" << std::setprecision (17) << m_all
                  << " (node " << i_all << ", C=" << Cd[i_all] << ")"
                  << "  solvent=" << m_solv << " (node " << i_solv << ")"
                  << "  molecule=" << m_mol
                  << std::setprecision (6) << std::endl;
      };

      if (it == 0) {
        std::size_t nsolv = 0;
        for (std::size_t i = 0; i < n; ++i) if (Cd[i] != 0.0) ++nsolv;
        std::cout << "  [Newton]   solvent nodes (rank 0): "
                  << nsolv << " / " << n << std::endl;
      }

      report ("phi_old:", phi_old.data ());   // system was built from this; == 0 at it==0
      report ("phi_new:", phi_new.data ());   // what the solve just returned
      report ("du     :", du.data ());

      double extramax = 0.0;
      for (std::size_t i = 0; i < n; ++i)
        if (Cd[i] != 0.0) {
          const double e = std::fabs (Cd[i] * (phi_old[i] * ion_model.dg (phi_old[i])
                                               - ion_model.g (phi_old[i])));
          if (e > extramax) extramax = e;
        }
      std::cout << "  [Newton]   max|rhs_extra| from phi_old = "
                << std::setprecision (17) << extramax
                << std::setprecision (6) << std::endl;
    }

    if (rank == 0)
      std::cout << "  [Newton] iter " << it
                << "   ||du||_inf = "         << std::setprecision (17) << anorm
                << "   ||du||_inf,solvent = " << dsolv
                << "   ||phi||_inf,solvent = " << usolv
                << "   tol = " << tol_s
                << std::setprecision (6) << std::endl;

    // --- divergence check: scan phi itself, not just the norm ---
    // (is_finite_bits, not std::isfinite: see its definition.) The node scan
    // is the test that matters: a NaN never wins the max in anorm.
    bool bad = !is_finite_bits (anorm);
    for (std::size_t i = 0; i < n && !bad; ++i)
      if (!is_finite_bits (phi_new[i])) bad = true;
    if (size > 1) {
      int b = bad ? 1 : 0, gb = 0;
      MPI_Allreduce (&b, &gb, 1, MPI_INT, MPI_MAX, mpicomm);
      bad = (gb != 0);
    }
    if (bad) {
      if (rank == 0)
        std::cout << "  [Newton] DIVERGED (non-finite solution) at iter "
                  << it << std::endl;
      break;
    }

    // --- convergence test, on the ion-accessible nodes only ---
    // |du|_s <= tol * max(1, |phi|_s): absolute OR relative, whichever is
    // looser. The norms over all nodes would be dominated by the potential
    // at the point charges inside the molecule, which grows as 1/h and has
    // nothing to do with the nonlinearity.
    if (it >= 1 && dsolv <= tol_s) {
      converged = true;
      if (rank == 0)
        std::cout << "  [Newton] converged in " << it + 1
                  << " iterations." << std::endl;
      break;
    }

    // --- stagnation: stop at the accuracy of the linear solver ---
    // Each solve starts from zero, so du cannot go below the error of one
    // linear solve. If du is already small and no longer decreases, stop.
    if (it >= 2 && dsolv < newton_stall * tol_s && dsolv > 0.5 * dsolv_prev) {
      converged = true;
      if (rank == 0)
        std::cout << "  [Newton] WARNING: stopped at the accuracy of the linear "
                     "solver after " << it + 1 << " iterations, ||du||_inf,solvent = "
                  << std::setprecision (17) << dsolv << " (tol = " << tol_s << ")"
                  << std::setprecision (6) << std::endl;
      break;
    }
    dsolv_prev = dsolv;

    // --- clamping globalization (iteration 0 unclamped) ---
    // The scale is set by the ion-accessible nodes only (dsolv). The next
    // Newton system depends on phi only where C != 0, so the step at the
    // other nodes (molecule, Stern) is overwritten by the next solve and
    // needs no damping. Clamping on all nodes let the linear jump inside the
    // molecule after the compression (~25 kT/e) throttle the whole step.
    double scale = 1.0;
    if (it >= 1 && dsolv > maxdu) {
      scale = maxdu / dsolv;
      if (rank == 0)
        std::cout << "  [Newton] CLAMP active: scale = " << scale << std::endl;
    }

    for (std::size_t i = 0; i < n; ++i)
      phi_new[i] = phi_old[i] + scale * du[i];

    // --- asinh compression of the linear solution (iteration 0 ONLY) ---
    // phi_new here is exactly the linear PB solution (scale == 1 at it 0).
    // Solvent nodes only (C != 0); the molecule is linear and re-adjusts
    // exactly at iteration 1.
    if (it == 0 && newton_compress == 1) {
      double loc_before = 0.0, loc_after = 0.0;   // diagnostics only
      for (std::size_t i = 0; i < n; ++i)
        if (Cd[i] != 0.0) {
          loc_before = std::max (loc_before, std::fabs (phi_new[i]));
          phi_new[i] = 2.0 * std::asinh (0.5 * phi_new[i]);
          loc_after  = std::max (loc_after,  std::fabs (phi_new[i]));
        }
      double g_before = loc_before, g_after = loc_after;
      if (size > 1) {
        MPI_Allreduce (&loc_before, &g_before, 1, MPI_DOUBLE, MPI_MAX, mpicomm);
        MPI_Allreduce (&loc_after,  &g_after,  1, MPI_DOUBLE, MPI_MAX, mpicomm);
      }
      if (rank == 0)
        std::cout << "  [Newton] COMPRESS (it 0, solvent): max|phi|_solv "
                  << std::setprecision (17) << g_before << " -> " << g_after
                  << std::setprecision (6) << std::endl;
    }

    // phi comes from new_node_vector (via the solver): only the values move.
    if (size > 1)
      phi->assemble (replace_op);
  }

  if (!converged && rank == 0)
    std::cout << "  [Newton] WARNING: did NOT converge in " << it_done + 1
              << " iterations. The reported potential is the last iterate."
              << std::endl;
}

void
poisson_boltzmann::export_tmesh (ray_cache_t & ray_cache)
{
  int size, rank;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);
  bim3a_solution_with_ghosts (tmsh, (*epsilon_nodes), replace_op);
  tmsh.octbin_export ("eps_map_0", (*epsilon_nodes));
}

void
poisson_boltzmann::export_potential_map (ray_cache_t & ray_cache)
{
  int size, rank;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);
  tmsh.octbin_export ("potential_map_0", (*phi));
}

void
poisson_boltzmann::export_marked_tmesh ()
{
  tmsh.octbin_export_quadrant (markerfilename.c_str (), marker);
}

void
poisson_boltzmann::export_p4est ()
{
  tmsh.save (p4estfilename.c_str ());
}



void
poisson_boltzmann::mumps_compute_electric_potential (ray_cache_t & ray_cache)
{
  int rank, size;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);


  mumps mumps_solver;

  std::vector<double> vals;
  std::vector<int> irow, jcol;

  (*A).aij (vals, irow, jcol, mumps_solver.get_index_base ());

  mumps_solver.set_lhs_distributed ();
  mumps_solver.set_distributed_lhs_structure (tmsh.num_global_nodes (), irow, jcol);
  mumps_solver.set_distributed_lhs_data (vals);
  mumps_solver.set_rhs_distributed (*rhs);
  rhs.reset ();
  A.reset ();
  std::cout << "mumps_solver.analyze () = "
            << mumps_solver.analyze ()
            << std::endl;
  std::cout << "mumps_solver.factorize () = "
            << mumps_solver.factorize ()
            << std::endl;
  std::cout << "mumps_solver.solve () = "
            << mumps_solver.solve ()
            << std::endl;

  // MUMPS gathers the solution on rank 0 (it owns every entry; the other
  // ranks hold their rows as non-local entries): read our rows by global
  // index, so phi keeps the distribution of the mesh.
  {
    const distributed_vector sol = mumps_solver.get_distributed_solution ();
    phi = new_node_vector ();
    auto & od = phi->get_owned_data ();
    const int is = phi->get_range_start ();

    for (std::size_t i = 0; i < od.size (); ++i)
      od[i] = sol[is + i];
  }

  if (size > 1)
    phi->assemble (replace_op);

  ///////

  mumps_solver.cleanup ();
}



void
poisson_boltzmann::lis_compute_electric_potential (ray_cache_t & ray_cache)
{
  int rank, size;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);

  //CSR
  std::vector<double> vals;
  std::vector<int> irow, jcol;


  (*A).csr (vals, jcol, irow);

  // lis RHS
  LIS_INT i, is, ie, n_rhs, ln;
  LIS_VECTOR rhs_lis;
  //n_rhs = tmsh.num_global_nodes();
  ln = rhs->get_owned_data ().size ();

  lis_vector_create (mpicomm, &rhs_lis);
  lis_vector_set_size (rhs_lis, ln, 0);
  lis_vector_get_range (rhs_lis, &is, &ie);

  double rhs_max = 0.0;   // for the zero-solution check after the solve

  for (i=is; i<ie; i++) {
    const double v = rhs->get_owned_data ()[i-is];
    lis_vector_set_value (LIS_INS_VALUE, i, v, rhs_lis);
    rhs_max = std::max (rhs_max, std::fabs (v));
  }

  //cleaning of rhs
  rhs.reset ();

  //lis_vector_print(rhs_lis);
  // lis PHI
  LIS_VECTOR phi_lis;

  lis_vector_create (mpicomm, &phi_lis);
  lis_vector_set_size (phi_lis, ln, 0);
  //lis_vector_set_size(phi_lis, 0, n_rhs);
  lis_vector_get_range (phi_lis, &is, &ie);

  // lis MATRIX
  LIS_INT n, nnz; //n: matrix dim ; nnz: numb of non zero elems
  LIS_INT *index; //array of integer containing the col index of non zero elems
  LIS_INT *ptr; //array of integer with starting points of rows
  LIS_SCALAR *value; //array of double stores non-zero elements of matrix A along the row
  LIS_MATRIX A_lis; //array of integer containing the col index of non zero elems

  nnz = (*A).owned_nnz ();
  n = tmsh.num_owned_nodes ();

  A.reset ();

  lis_matrix_create (mpicomm, &A_lis);

  lis_matrix_set_size (A_lis, n, 0);
  ptr = &irow[0];
  index = &jcol[0];
  value = &vals[0];

  lis_matrix_set_csr (nnz, ptr, index, value, A_lis);


  lis_matrix_assemble (A_lis);

  //Solve linear system
  LIS_SOLVER solver;

  lis_solver_create (&solver);

  std::string opts = linear_solver_options;

  lis_solver_set_option (&opts[0], solver);

  lis_solve (A_lis, rhs_lis, phi_lis, solver);

  lis_solver_destroy (solver);
  lis_vector_destroy (rhs_lis);

  phi = new_node_vector ();

  lis_vector_get_values (phi_lis, is, ln, phi->get_owned_data ().data ());

  lis_vector_destroy (phi_lis);

  // LIS starts from x = 0 and, with an absolute stopping test (conv_cond 2,
  // ||b - Ax||_1 <= tol) and tol >= ||b||_1, accepts it without iterating.
  // A zero solution with a nonzero right-hand side is never correct: stop.
  {
    double loc[2] = {rhs_max, 0.0};
    for (const double v : phi->get_owned_data ())
      loc[1] = std::max (loc[1], std::fabs (v));

    double glob[2] = {loc[0], loc[1]};
    if (size > 1)
      MPI_Allreduce (loc, glob, 2, MPI_DOUBLE, MPI_MAX, mpicomm);

    if (glob[0] > 0.0 && glob[1] == 0.0) {
      if (rank == 0)
        std::cerr << "ERROR: the linear solver returned the zero vector for a "
                     "nonzero right-hand side.\n       Check solver_options "
                     "(with -conv_cond 2 the -tol threshold is absolute).\n";
      MPI_Abort (mpicomm, 1);
    }
  }

  if (size > 1)
    phi->assemble (replace_op);

}


void
poisson_boltzmann::write_potential_on_atoms_fast ()
{
  int rank, size;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  std::vector<std::string> local_lines;

  double phi_on_atom;
  double phi_hang_nodes = 0.0;

  // Costruisci stringhe localmente
  for (auto it = lookup_table.begin(); it != lookup_table.end(); ++it) {
    phi_on_atom = 0.0;
    double volume = (it->second.p (0, 7) - it->second.p (0, 0)) *
                    (it->second.p (1, 7) - it->second.p (1, 0)) *
                    (it->second.p (2, 7) - it->second.p (2, 0));

    for (int ii = 0; ii < 8; ++ii) {
      double weigth = std::abs ((pos_atoms[it->first][0] - it->second.p (0, 7-ii)) *
                                (pos_atoms[it->first][1] - it->second.p (1, 7-ii)) *
                                (pos_atoms[it->first][2] - it->second.p (2, 7-ii))) / volume;

      if (!it->second.is_hanging (ii)) {
        phi_on_atom += (*phi)[it->second.gt (ii)] * weigth;
      } else {
        phi_hang_nodes = 0.0;

        for (int jj = 0; jj < it->second.num_parents (ii); ++jj) {
          phi_hang_nodes += (*phi)[it->second.gparent (jj, ii)] / it->second.num_parents (ii);
        }

        phi_on_atom += phi_hang_nodes * weigth;
      }
    }

    std::ostringstream oss;
    oss << std::setw (8) << index_atoms[it->first]
        << std::fixed << std::setprecision (3)
        << std::setw (8) << pos_atoms[it->first][0]
        << std::setw (8) << pos_atoms[it->first][1]
        << std::setw (8) << pos_atoms[it->first][2]
        << std::fixed << std::setprecision (4)
        << "  " << phi_on_atom << "\n";

    local_lines.push_back (oss.str());
  }

  // Serializzazione delle stringhe
  std::string local_data;

  for (const auto& line : local_lines)
    local_data += line;

  int local_size = local_data.size();
  std::vector<int> all_sizes (size);

  MPI_Gather (&local_size, 1, MPI_INT, all_sizes.data(), 1, MPI_INT, 0, mpicomm);

  std::vector<int> displs (size);
  std::string global_data;

  if (rank == 0) {
    int total_size = 0;

    for (int i = 0; i < size; ++i) {
      displs[i] = total_size;
      total_size += all_sizes[i];
    }

    global_data.resize (total_size);
  }

  MPI_Gatherv (local_data.data(), local_size, MPI_CHAR,
               rank == 0 ? &global_data[0] : nullptr,
               all_sizes.data(), displs.data(), MPI_CHAR,
               0, mpicomm);

  // Solo rank 0 scrive sul file
  if (rank == 0) {
    std::ofstream phi_atoms ("phi_on_atoms.txt");
    phi_atoms << global_data;
    phi_atoms.close();
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////

std::array<double,12>
poisson_boltzmann::cube_fraction_intersection (tmesh_3d::quadrant_iterator& quadrant,
    const ray_cache_t & ray_cache)


// v6_________e7_________v7
// /|                  /|
// e8 / |                 / |
// /  |             e6 /  |
// /   | e12           /   | e11
// v4/____|_____e5_______/v5  |
// |    |              |    |
// |  v2|______e3______|____|v3
// e9  |   /               |   /
// |  /            e10 |  /
// | /  e4             | / e2
// |/                  |/
// v0/_________e1________/v1
{
  std::array<double,12> fraction = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};

  int dir;
  int i1, i2;
  double x1, x2;
  std::array<double,2> ray;

  if (marker[quadrant->get_forest_quad_idx ()] == 0.5) {
    for (int j: {
           0,2,4,6
         })
      //x-axis edge
    {
      dir = 0;
      i1 = edge2nodes[2 * j ];
      i2 = edge2nodes[2 * j + 1];
      x1 = quadrant->p (dir, i1);
      x2 = quadrant->p (dir, i2);
      std::vector<int> direzioni {0,1,2};
      direzioni.erase (direzioni.begin ()+dir);

      for (unsigned i = 0; i < direzioni.size (); ++i) {
        ray[i] = quadrant->p (direzioni[i], i1);
      }

      auto it0 = ray_cache.rays[dir].find (ray);
      auto inters = it0->second.inters;

      for (int ii =0; ii<inters.size (); ii++) {
        if (inters[ii]>= x1 && inters[ii] <=x2) {
          fraction[j] = (inters[ii] - x1)/ (x2 - x1);
        }
      }
    }

    for (int j: {
           1,3,5,7
         })
      //y-axis edge
    {
      dir = 1;
      i1 = edge2nodes[2 * j ];
      i2 = edge2nodes[2 * j + 1];
      x1 = quadrant->p (dir, i1);
      x2 = quadrant->p (dir, i2);
      std::vector<int> direzioni {0,1,2};
      direzioni.erase (direzioni.begin ()+dir);

      for (unsigned i = 0; i < direzioni.size (); ++i) {
        ray[i] = quadrant->p (direzioni[i], i1);
      }

      auto it0 = ray_cache.rays[dir].find (ray);
      auto inters = it0->second.inters;

      for (int ii =0; ii<inters.size (); ii++) {
        if (inters[ii]>= x1 && inters[ii] <=x2) {
          fraction[j] = (inters[ii] - x1)/ (x2 - x1);
        }
      }
    }

    for (int j: {
           8,9,10,11
         })
      //z-axis edge
    {
      dir = 2;
      i1 = edge2nodes[2 * j ];
      i2 = edge2nodes[2 * j + 1];
      x1 = quadrant->p (dir, i1);
      x2 = quadrant->p (dir, i2);
      std::vector<int> direzioni {0,1,2};
      direzioni.erase (direzioni.begin ()+dir);

      for (unsigned i = 0; i < direzioni.size (); ++i) {
        ray[i] = quadrant->p (direzioni[i], i1);
      }

      auto it0 = ray_cache.rays[dir].find (ray);
      auto inters = it0->second.inters;

      for (int ii =0; ii<inters.size (); ii++) {
        if (inters[ii]>= x1 && inters[ii] <=x2) {
          fraction[j] = (inters[ii] - x1)/ (x2 - x1);
        }
      }
    }
  }

  return fraction;
}

void
poisson_boltzmann::normal_intersection (tmesh_3d::quadrant_iterator& quadrant,
                                        const ray_cache_t & ray_cache,
                                        int edge, std::array<double,3> &norm,
                                        double &frac)
{
  int dir = edge_axis[edge];
  int i1 = edge2nodes[2*edge ];
  int i2 = edge2nodes[2*edge + 1];
  double x1 = quadrant->p (dir, i1),
         x2 = quadrant->p (dir, i2);

  std::array<double,2> ray;
  std::vector<int> direzioni = {0,1,2};
  direzioni.erase (direzioni.begin ()+dir);

  for (unsigned i = 0; i < direzioni.size (); ++i) {
    ray[i] = quadrant->p (direzioni[i], i1);
  }

  auto it0 = ray_cache.rays[dir].find (ray);
  auto normali = it0->second.normals;
  auto inters = it0->second.inters;

  frac = 0.5;

  for (int ii =0; ii<inters.size (); ii++) {
    if (inters[ii]>= x1 && inters[ii] <=x2) {
      norm[0] = normali[0 + 3*ii];
      norm[1] = normali[1 + 3*ii];
      norm[2] = normali[2 + 3*ii];
      frac = (inters[ii] - x1)/ (x2 - x1);
    }
  }
}

int
poisson_boltzmann::classifyCube (tmesh_3d::quadrant_iterator& quadrant,
                                 double isolevel)
{
  int cubeindex = 0;
  int index = 1;
  double tmp = 0;
  constexpr double EPSILON = 1e-10;

  for (int ii : {
         0,1,3,2,4,5,7,6
       }) {
    if (! quadrant->is_hanging (ii)) {
      if ( (*epsilon_nodes)[quadrant->gt (ii)] < (isolevel - EPSILON)) cubeindex |= index;
    } else {
      for (int jj = 0; jj < quadrant->num_parents (ii); ++jj) {
        tmp += (*epsilon_nodes)[quadrant->gparent (jj, ii)] / quadrant->num_parents (ii);
      }

      if (tmp < (isolevel - EPSILON)) cubeindex |= index;
    }

    tmp = 0;
    index *= 2;
  }

  // Cube is entirely in/out of the surface
  if (edgeTable[cubeindex] == 0)
    return -1;

  return cubeindex;
}

int
poisson_boltzmann::classifyCube_fast (tmesh_3d::quadrant_iterator& quadrant,
                                      double isolevel)
{
  int cubeindex = 0;
  int index = 1;
  double tmp = 0;
  constexpr double EPSILON = 1e-10;

  for (int ii : {
         0,1,3,2,4,5,7,6
       }) {
    if ( (*epsilon_nodes)[quadrant->gt (ii)] < (isolevel - EPSILON)) cubeindex |= index;

    index *= 2;
  }

  // Cube is entirely in/out of the surface
  if (edgeTable[cubeindex] == 0)
    return -1;

  return cubeindex;
}

std::tuple<std::array<double,8>, std::array<double,8>, std::vector<int>,std::vector<int> >
poisson_boltzmann::classifyCube_flux (tmesh_3d::quadrant_iterator& quadrant,
                                      std::array<double,8>& tmp_phi,
                                      std::array<double,8>& tmp_eps)
{
  std::vector<int> edges {};
  std::vector<int> flux {};

  for (int ii = 0; ii < 8; ++ii) {
    if (! quadrant->is_hanging (ii)) {
      tmp_eps[ii]= (*epsilon_nodes)[quadrant->gt (ii)];
      tmp_phi[ii]= (*phi)[quadrant->gt (ii)];
    } else {
      tmp_eps[ii] = 0.0;
      tmp_phi[ii] = 0.0;

      for (int jj = 0; jj < quadrant->num_parents (ii); ++jj) {
        tmp_eps[ii] += (*epsilon_nodes)[quadrant->gparent (jj, ii)] / quadrant->num_parents (ii);
        tmp_phi[ii] += (*phi)[quadrant->gparent (jj, ii)] / quadrant->num_parents (ii);
      }
    }
  }

  for (int ii = 0; ii < 12; ++ii) {
    if (tmp_eps[edge2nodes[2*ii]] < tmp_eps[edge2nodes[2*ii +1]]) {
      flux.push_back (1);
      edges.push_back (ii);
    } else if (tmp_eps[edge2nodes[2*ii]] > tmp_eps[edge2nodes[2*ii +1]]) {
      flux.push_back (-1);
      edges.push_back (ii);
    }
  }

  return make_tuple (tmp_phi, tmp_eps, edges, flux);
}

std::tuple<std::array<double,8>, std::array<double,8>, std::vector<int>,std::vector<int> >
poisson_boltzmann::classifyCube_flux_fast (tmesh_3d::quadrant_iterator& quadrant,
    std::array<double,8>& tmp_phi,
    std::array<double,8>& tmp_eps)
{
  std::vector<int> edges {};
  std::vector<int> flux {};


  for (int ii = 0; ii < 8; ++ii) {
    tmp_eps[ii]= (*epsilon_nodes)[quadrant->gt (ii)];
    tmp_phi[ii]= (*phi)[quadrant->gt (ii)];
  }

  for (int ii = 0; ii < 12; ++ii) {
    if (tmp_eps[edge2nodes[2*ii]] < tmp_eps[edge2nodes[2*ii +1]]) {
      flux.push_back (1);
      edges.push_back (ii);
    } else if (tmp_eps[edge2nodes[2*ii]] > tmp_eps[edge2nodes[2*ii +1]]) {
      flux.push_back (-1);
      edges.push_back (ii);
    }
  }

  return make_tuple (tmp_phi, tmp_eps, edges, flux);
}

double wha (double eps1, double eps2, double frac)
{
  return 1.0/ (frac/eps1 + (1-frac)/eps2);
}

double flux_dir (double eps1, double eps2)
{
  return eps1 < eps2 ? 1 : -1;
}

double phi0 (double eps1, double eps2,
             double phi1, double phi2, double frac)
{
  return phi1 + frac*eps2* (phi2-phi1)/ (eps2*frac + eps1* (1-frac));
}

double areaTriangle (const std::array<std::array<double,3>,3> &triangle)
{

  double area;
  std::array<double,3> ab;
  std::array<double,3> ac;

  for (int i = 0; i < 3; ++i) {
    ab[i] = triangle[0][i] - triangle[1][i];
    ac[i] = triangle[0][i] - triangle[2][i];
  }

  area = 0.5* std::hypot (ab[1]*ac[2] - ab[2]*ac[1],
                          ab[0]*ac[2] - ab[2]*ac[0],
                          ab[1]*ac[0] - ab[0]*ac[1]);
  return area;
}

double SphercalAreaTriangle (const std::array<std::array<double,3>,3> &triangle)
{

  double area;
  double rs =2;
  double a,b,c, aa, bb, cc, E;
  double nA,nB,nC;
  std::array<double,3> A = triangle[0];
  std::array<double,3> B = triangle[1];
  std::array<double,3> C = triangle[2];

  nA = std::hypot (A[0], A[1], A[2]);
  nB = std::hypot (B[0], B[1], B[2]);
  nC = std::hypot (C[0], C[1], C[2]);

  a = std::inner_product (std::begin (A), std::end (A), std::begin (B), 0.0);
  b = std::inner_product (std::begin (A), std::end (A), std::begin (C), 0.0);
  c = std::inner_product (std::begin (B), std::end (B), std::begin (C), 0.0);

  a = a/ (nA*nB);
  b = b/ (nA*nC);
  c = c/ (nC*nB);

  a = std::acos (a);
  b = std::acos (b);
  c = std::acos (c);

  if (a>pi/2.0)
    a = pi - a;

  if (b>pi/2.0)
    b = pi - b;

  if (c>pi/2.0)
    c = pi - c;

  aa = std::acos ( (cos (a)-cos (b)*cos (c))/ (sin (b)*sin (c)));
  bb = std::acos ( (cos (b)-cos (c)*cos (a))/ (sin (a)*sin (c)));
  cc = std::acos ( (cos (c)-cos (b)*cos (a))/ (sin (b)*sin (a)));

  E = aa + bb + cc - pi;

  area = rs*rs * E;

  return area;
}

int
poisson_boltzmann::getTriangles (int cubeindex,
                                 std::array<std::array<int,3>,5> &triangles)
{
  int ntriang = 0;
  int i;
  triangles.fill ({}); //set matrix to zero

  for (i=0; triTable[cubeindex][i]!=-1; i+=3) {
    // save all the assigned indexes
    triangles[ntriang][0] = triTable[cubeindex][i ];
    triangles[ntriang][1] = triTable[cubeindex][i+1];
    triangles[ntriang][2] = triTable[cubeindex][i+2];
    ntriang++;
  }

  return ntriang;
}

// ============================================================
//  Electrostatic energy and per-atom potential/field decomposition
//
//  Boundary-integral representation on the molecular surface, from phi and
//  its normal flux. Three summable terms at each atom (self term excluded):
//    level 1: c = direct Coulomb, p = polarization
//    level 2: i = ionic
//  The ionic term is the only expensive one (atoms x surface triangles).
//
//  Options:
//    calc_potential_terms (P), calc_field_terms (F): pot_field.dat columns.
//      The field of a term gives its potential for free, not the other way
//      round, so potentials are written up to level max(P,F) and fields up
//      to level F. E_c is always written together with phi_c.
//    calc_energy (E): energy terms; the surface loop runs once, up to level
//      max(P,F,E). Without salt the ionic term vanishes: level 2 -> 1.
//  pot_field.dat is written only if max(P,F) > 0; otherwise only the charged
//  atoms enter the loops.
// ============================================================
void
poisson_boltzmann::energy_pot_field (ray_cache_t & ray_cache)
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);

  const double inv_4pi = 1.0 / (4.0 * pi);
  const double eps_in = 4.0 * pi * e_0 * e_in * kb * T * Angs / (e * e);
  const double eps_out = 4.0 * pi * e_0 * e_out * kb * T * Angs / (e * e);

  const double C0 = 1.0e3 * N_av * ionic_strength; // [mol/m^3]
  const double k2 = 2.0 * C0 * Angs * Angs * e * e / (e_0 * e_out * kb * T);
  const bool salt = std::sqrt (k2) > 1.e-5;

  const double den_in = 1.0 / eps_in;
  const double constant_pol = (1.0 / eps_out - 1.0 / eps_in) * inv_4pi;
  const double constant_react = (1.0 / eps_out) * inv_4pi;

  // --- levels ---
  const int lev_pot = std::max (calc_potential_term, calc_field_term);
  const int lev_field = calc_field_term;
  const int lev_req = std::max (lev_pot, calc_energy);
  const bool write_file = lev_pot > 0;
  const bool ionic = lev_req >= 2 && salt;
  const bool pot_i = lev_pot >= 2 && salt;
  const bool field_p = lev_field >= 1;
  const bool field_i = lev_field >= 2 && salt;

  // Uniform mesh at the surface: no hanging nodes on border quadrants.
  const bool refined = (loc_refinement == 1 || mesh_shape > 2
                        || (mesh_shape == 2 && refine_box == 1));

  if (rank == 0) {
    std::cout << "\n================ [ Electrostatic Energy ] =================\n";

    if (calc_potential_term < calc_field_term)
      std::cout << "  [INFO] calc_potential_terms raised to " << lev_pot
                << " (= calc_field_terms)\n";

    if (lev_req >= 2 && !salt)
      std::cout << "  [INFO] No salt: the ionic term is zero, computed up to level 1\n";
  }

  // --- atoms entering the loops ---
  // pot_field.dat: all atoms, used in place. Energy only: just the charged
  // atoms, and the full vectors are freed (nothing needs them afterwards).
  std::vector<double> q_chg;
  std::vector<std::array<double, 3>> r_chg;

  if (!write_file) {
    for (size_t ii = 0; ii < charge_atoms.size (); ++ii)
      if (std::fabs (charge_atoms[ii]) > 1.e-5) {
        q_chg.push_back (charge_atoms[ii]);
        r_chg.push_back (pos_atoms[ii]);
      }

    std::vector<double>().swap (charge_atoms);
    std::vector<std::array<double, 3>>().swap (pos_atoms);
  }

  const std::vector<double> &q_at = write_file ? charge_atoms : q_chg;
  const std::vector<std::array<double, 3>> &r_at = write_file ? pos_atoms : r_chg;
  const size_t num_atoms = q_at.size ();

  std::vector<double> phi_c, phi_p, phi_i;
  std::vector<double> field_cx, field_cy, field_cz;
  std::vector<double> field_px, field_py, field_pz;
  std::vector<double> field_ix, field_iy, field_iz;

  if (write_file) {
    phi_c.assign (num_atoms, 0.0);
    phi_p.assign (num_atoms, 0.0);
    field_cx.assign (num_atoms, 0.0);
    field_cy.assign (num_atoms, 0.0);
    field_cz.assign (num_atoms, 0.0);
  }

  if (pot_i)
    phi_i.assign (num_atoms, 0.0);

  if (field_p) {
    field_px.assign (num_atoms, 0.0);
    field_py.assign (num_atoms, 0.0);
    field_pz.assign (num_atoms, 0.0);
  }

  if (field_i) {
    field_ix.assign (num_atoms, 0.0);
    field_iy.assign (num_atoms, 0.0);
    field_iz.assign (num_atoms, 0.0);
  }

  // --- direct Coulomb term (replicated on every rank) ---
  // phi_c is needed for pot_field.dat, so there the energy comes for free and
  // is reported even with calc_coulombic = 0.
  const bool coul_on = write_file || calc_coulombic == 1;
  this->coul_energy = 0.0;

  if (coul_on) {
    double coul = 0.0;

    for (size_t i = 0; i < num_atoms; ++i) {
      const std::array<double,3> &ri = r_at[i];
      const double qi = q_at[i];

      for (size_t j = i + 1; j < num_atoms; ++j) {
        const std::array<double,3> &rj = r_at[j];
        const double qj = q_at[j];
        const double dx = ri[0] - rj[0];
        const double dy = ri[1] - rj[1];
        const double dz = ri[2] - rj[2];
        const double r2 = dx * dx + dy * dy + dz * dz;
        const double r = std::sqrt (r2);

        coul += qi * qj / r;

        if (write_file) {
          const double inv_r3 = den_in / (r2 * r);
          phi_c[i] += qj / r * den_in;
          phi_c[j] += qi / r * den_in;
          field_cx[i] += dx * inv_r3 * qj;
          field_cy[i] += dy * inv_r3 * qj;
          field_cz[i] += dz * inv_r3 * qj;
          field_cx[j] -= dx * inv_r3 * qi;
          field_cy[j] -= dy * inv_r3 * qi;
          field_cz[j] -= dz * inv_r3 * qi;
        }
      }
    }

    this->coul_energy = coul * den_in;
  }

  // --- surface loop: polarization (fluxes) and ionic (triangles) ---
  double charge_pol = 0.0, first_int = 0.0, second_int = 0.0;

  std::array<double,3> h{0}, area_h{0};
  std::array<double,3> V, N;
  std::array<double,8> tmp_eps, tmp_phi;
  std::vector<int> edg, fl_dir;
  std::array<std::array<double,3>,3> vert_triangles, norms_vert;
  std::array<double,3> phi_sup;

  // Per-atom kernels. The optional outputs are compile-time switches
  // (if constexpr): runtime branches in these loops stop vectorization and
  // made the energy-only case ~3x slower.
  const std::true_type yes;
  const std::false_type no;

  // Polarization: one flux element tmp_flux at V. Returns sum_a q_a tmp_flux / r_a.
  auto flux_kernel = [&] (auto pot, auto field, const std::array<double,3> &V,
                          double tmp_flux) {
    double acc = 0.0;

    for (size_t ia = 0; ia < num_atoms; ++ia) {
      const std::array<double,3> &ra = r_at[ia];
      const double dx = ra[0] - V[0];
      const double dy = ra[1] - V[1];
      const double dz = ra[2] - V[2];
      const double r = std::sqrt (dx * dx + dy * dy + dz * dz);
      const double qflux = tmp_flux / r;

      acc += q_at[ia] * qflux;

      if constexpr (decltype (pot)::value)
        phi_p[ia] += qflux * constant_pol;

      if constexpr (decltype (field)::value) {
        const double c = tmp_flux * constant_pol / (r * r * r);
        field_px[ia] += dx * c;
        field_py[ia] += dy * c;
        field_pz[ia] += dz * c;
      }
    }

    return acc;
  };

  // Ionic: one surface triangle (vertex quadrature). Returns sum_a q_a dphi_a.
  auto tri_kernel = [&] (auto pot, auto field,
                         const std::array<std::array<double,3>,3> &vt,
                         const std::array<std::array<double,3>,3> &nt,
                         const std::array<double,3> &ps, double area) {
    double acc = 0.0;

    for (size_t ia = 0; ia < num_atoms; ++ia) {
      const std::array<double,3> &ra = r_at[ia];

      for (int kk = 0; kk < 3; ++kk) {
        const std::array<double,3> dv = {vt[kk][0] - ra[0],
                                         vt[kk][1] - ra[1],
                                         vt[kk][2] - ra[2]
                                        };
        const double r2 = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];
        const double r = std::sqrt (r2);
        const double inv_r3 = 1.0 / (r2 * r);
        const std::array<double,3> &nv = nt[kk];
        const double dot = dv[0]*nv[0] + dv[1]*nv[1] + dv[2]*nv[2];
        const double factor = ps[kk] * inv_4pi * area / 3.0;
        const double dphi = factor * dot * inv_r3;

        acc += q_at[ia] * dphi;

        if constexpr (decltype (pot)::value)
          phi_i[ia] += dphi;

        if constexpr (decltype (field)::value) {
          const double inv_r5 = inv_r3 / r2;
          field_ix[ia] += factor * (-3 * dv[0] * inv_r5 * dot + nv[0] * inv_r3);
          field_iy[ia] += factor * (-3 * dv[1] * inv_r5 * dot + nv[1] * inv_r3);
          field_iz[ia] += factor * (-3 * dv[2] * inv_r5 * dot + nv[2] * inv_r3);
        }
      }
    }

    return acc;
  };

  auto quadrant = this->tmsh.begin_quadrant_sweep ();

  auto set_h = [&] () {
    for (int d = 0; d < 3; ++d)
      h[d] = quadrant->p (d, 7) - quadrant->p (d, 0);

    area_h = {h[1]*h[2]/h[0]*0.25, h[0]*h[2]/h[1]*0.25, h[0]*h[1]/h[2]*0.25};
  };

  // Uniform mesh: all border quadrants have the same size, set it once.
  if (!refined && !border_quad.empty ()) {
    quadrant[border_quad[0]];
    set_h ();
  }

  for (const int ii : border_quad) {
    quadrant[ii];

    if (refined)
      set_h ();

    std::tie (tmp_phi, tmp_eps, edg, fl_dir) = refined
        ? classifyCube_flux (quadrant, tmp_phi, tmp_eps)
        : classifyCube_flux_fast (quadrant, tmp_phi, tmp_eps);

    // --- fluxes (polarization)
    for (int ip = 0; ip < edg.size (); ++ip) {
      const int edge = edg[ip];
      const int axis = edge_axis[edge];
      const int i1 = edge2nodes[2 * edge];
      const int i2 = edge2nodes[2 * edge + 1];

      double fract = 0.0;
      normal_intersection (quadrant, ray_cache, edge, N, fract);

      V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
      V[axis] += fract * h[axis];

      const double tmp_flux =
        - (tmp_phi[i2] - tmp_phi[i1]) * wha (tmp_eps[i1], tmp_eps[i2], fract)
        * fl_dir[ip] * area_h[axis];

      charge_pol += tmp_flux;

      if (!write_file)
        first_int += flux_kernel (no, no, V, tmp_flux);
      else if (!field_p)
        first_int += flux_kernel (yes, no, V, tmp_flux);
      else
        first_int += flux_kernel (yes, yes, V, tmp_flux);
    }

    if (!ionic)
      continue;

    // --- triangles (ionic)
    const int cubeindex = refined ? classifyCube (quadrant, eps_out)
                                  : classifyCube_fast (quadrant, eps_out);
    const int ntriang = cubeindex < 0 ? 0 : getTriangles (cubeindex, triangles);

    for (int itri = 0; itri < ntriang; ++itri) {
      for (int jj = 0; jj < 3; ++jj) {
        const int edge = triangles[itri][jj];
        const int axis = edge_axis[edge];
        const int i1 = edge2nodes[2 * edge];
        const int i2 = edge2nodes[2 * edge + 1];

        double fract = 0.0;
        normal_intersection (quadrant, ray_cache, edge, N, fract);

        V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
        V[axis] += fract * h[axis];

        vert_triangles[jj] = V;
        norms_vert[jj] = N;

        phi_sup[jj] = phi0 (tmp_eps[i1], tmp_eps[i2], tmp_phi[i1], tmp_phi[i2], fract);
      }

      const double area = areaTriangle (vert_triangles);

      if (!pot_i)
        second_int += tri_kernel (no, no, vert_triangles, norms_vert, phi_sup, area);
      else if (!field_i)
        second_int += tri_kernel (yes, no, vert_triangles, norms_vert, phi_sup, area);
      else
        second_int += tri_kernel (yes, yes, vert_triangles, norms_vert, phi_sup, area);
    }
  }

  this->energy_pol = 0.5 * constant_pol * first_int;
  this->energy_react = ionic ? 0.5 * (second_int - first_int * constant_react) : 0.0;

  auto reduce_double = [&] (double &x) {
    MPI_Reduce (rank == 0 ? MPI_IN_PLACE : &x, &x, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
  };

  auto reduce_vec = [&] (std::vector<double> &v) {
    MPI_Reduce (rank == 0 ? MPI_IN_PLACE : v.data (),
                rank == 0 ? v.data () : nullptr,
                (int) v.size (), MPI_DOUBLE, MPI_SUM, 0, mpicomm);
  };

  reduce_double (charge_pol);
  reduce_double (energy_pol);
  reduce_double (energy_react);

  for (std::vector<double> *v : {&phi_p, &phi_i, &field_px, &field_py, &field_pz,
                                 &field_ix, &field_iy, &field_iz})
    reduce_vec (*v);

  if (rank != 0)
    return;

  constexpr int label_width = 50;
  constexpr int precision = 16;

  std::cout << std::left << std::setw (label_width) << "  Net charge [e]:"
            << std::setprecision (precision) << net_charge << "\n";

  std::cout << std::left << std::setw (label_width) << "  Flux charge [e]:"
            << std::setprecision (precision) << charge_pol / (4.0 * pi) << "\n";

  std::cout << std::left << std::setw (label_width) << "  Polarization energy [kT]:"
            << std::setprecision (precision) << energy_pol << "\n";

  if (lev_req >= 2) {
    std::cout << std::left << std::setw (label_width) << "  Direct ionic energy [kT]:"
              << std::setprecision (precision) << energy_react << "\n";
  }

  if (coul_on) {
    std::cout << std::left << std::setw (label_width) << "  Coulombic energy [kT]:"
              << std::setprecision (precision) << coul_energy << "\n";
  }

  std::cout << std::left << std::setw (label_width) << "  Sum of electrostatic energy contributions [kT]:"
            << std::setprecision (precision)
            << (energy_pol + energy_react + coul_energy) << "\n";

  std::cout << "===========================================================\n";

  if (!write_file)
    return;

  // The ionic surface integral also contains the reaction field of the
  // polarization charge seen from the solvent: remove it.
  const double p2i = constant_react / constant_pol;

  std::ofstream fout ("pot_field.dat");
  fout << "# index    x    y    z    phi_c    phi_p    ";

  if (pot_i)
    fout << "phi_i    ";

  fout << "Ex_c    Ey_c    Ez_c";

  if (field_p)
    fout << "   Ex_p    Ey_p    Ez_p";

  if (field_i)
    fout << "   Ex_i    Ey_i    Ez_i";

  fout << "\n";

  for (size_t i = 0; i < num_atoms; ++i) {
    fout << std::setw (5) << i + 1 << "  "
         << std::setw (8) << r_at[i][0] << "  "
         << std::setw (8) << r_at[i][1] << "  "
         << std::setw (8) << r_at[i][2] << "  "
         << std::setw (8) << phi_c[i] << "  "
         << std::setw (8) << phi_p[i] << "  ";

    if (pot_i)
      fout << std::setw (8) << phi_i[i] - phi_p[i] * p2i << "  ";

    fout << std::setw (8) << field_cx[i] << "  "
         << std::setw (8) << field_cy[i] << "  "
         << std::setw (8) << field_cz[i] << "  ";

    if (field_p)
      fout << std::setw (8) << field_px[i] << "  "
           << std::setw (8) << field_py[i] << "  "
           << std::setw (8) << field_pz[i] << "  ";

    if (field_i)
      fout << std::setw (8) << field_ix[i] - field_px[i] * p2i << "  "
           << std::setw (8) << field_iy[i] - field_py[i] * p2i << "  "
           << std::setw (8) << field_iz[i] - field_pz[i] * p2i << "  ";

    fout << "\n";
  }

  fout.close ();
  std::cout << "Atom potentials and fields written to 'pot_field.dat'\n";
}

// ============================================================
//  Nonlinear excess ionic free energy (1:1 salt, sinh or steric model)
//
//  For the nonlinear PBE the surface-integral partition computed by
//  energy_pot_field() still yields
//      G_coul + G_pol + G_ion_dir = 1/2 sum_i q_i phi(r_i)
//  (the Green identity behind it only needs -div(eps grad phi) = rho_s in
//  the solvent, whatever rho_s(phi) is). What is missing is the excess term
//      G_exc = -int_{solvent} [ 1/2 rho_s phi + (P - P0) ] dV
//            = 2 kT n_b int_{solvent} f_exc(psi) dV,
//  with f_exc = psi/2 g(psi) - osm(psi) from ion_model_t:
//      sinh  : f_exc = psi/2 sinh(psi) - cosh(psi) + 1              (>= 0)
//      steric: f_exc = psi/2 sinh(psi)/D(psi) - log(D(psi))/nu,
//              D = 1 + nu (cosh(psi) - 1)   (P - P0 = kT/a^3 log D)
//  It vanishes identically in the linear limit (O(psi^4)).
//
//  In code units (psi dimensionless, lengths in Angstrom, charges in e) the
//  nodal reaction coefficient is C = reaction_nodes = 8 pi n_b A^3 (0 inside
//  the molecule), hence  G_exc/kT = (1/4pi) int C(r) f_exc(psi) dV.
//  The volume integral uses bim3a_rhs_frac with the SAME cut-cell quadrature
//  as the Newton reaction term, so energy and residual are consistent.
//  Requires reaction_nodes / ones / phi still allocated (true after
//  newton_solve; assemple_system_matrix would free them).
// ============================================================
void
poisson_boltzmann::energy_excess_nonlinear (ray_cache_t & ray_cache)
{
  int size, rank;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);

  if (!reaction_nodes || !ones || !phi) {
    if (rank == 0)
      std::cerr << "  [energy_exc] ERROR: reaction_nodes/ones/phi not allocated. "
                   "energy_excess_nonlinear() must run after newton_solve().\n";
    return;
  }

  const double inv_4pi = 1.0 / (4.0 * pi);
  const std::size_t n = tmsh.num_owned_nodes ();

  auto func_frac = [&] (tmesh_3d::quadrant_iterator& quadrant) {
    return cube_fraction_intersection (quadrant, ray_cache);
  };

  // Nodal quadrature weights w_i = solvent-fraction patch volume of node i.
  // bim3a_rhs_frac is linear in its nodal argument (hanging nodes spread
  // their volume evenly over the parents), so int_{solvent} g dV = sum_i w_i g_i
  // with g evaluated at owned nodes only. Computing w once avoids one
  // cut-cell sweep per integrand.
  distributed_vector w (n, mpicomm);
  w.get_owned_data ().assign (n, 0.0);
  bim3a_rhs_frac (tmsh, *ones, *ones, w, func_frac);
  if (size > 1)
    w.assemble ();

  // Integrands are all multiplied by C and set to zero where C == 0, to
  // avoid 0*inf near the point charges (exactly as in assemble_newton_system):
  //   C f_exc(psi)      -> G_exc
  //   C psi/2 g(psi)    -> -1/2 int rho_s phi
  //   C osm(psi)        -> int (P - P0)
  //   C g(psi)          -> -4pi * mobile ion charge
  //   C (n_+/n_b - 1)   -> 8pi * cation excess in the solvent (n_+ = n_b e^{-psi}/D)
  //   C (n_-/n_b - 1)   -> 8pi * anion excess in the solvent  (n_- = n_b e^{+psi}/D)
  //   w [C != 0]        -> ion-accessible volume
  // G_exc is NOT recomputed as the difference of the two addends: for the
  // steric model they are both O(psi/nu) at large psi and would cancel.
  // n_pm/n_b - 1 = (expm1(-+psi) - nu coshm1(psi)) / D, without cancellation.
  double loc[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  {
    auto & Cd = reaction_nodes->get_owned_data ();
    auto & ud = phi->get_owned_data ();
    auto & wd = w.get_owned_data ();

    for (std::size_t i = 0; i < n; ++i) {
      const double Ci = Cd[i];
      if (Ci == 0.0)
        continue;
      const double ui = ud[i];
      const double g  = ion_model.g (ui);
      const double cw = Ci * wd[i];
      loc[0] += cw * ion_model.f_exc (ui);
      loc[1] += cw * 0.5 * ui * g;
      loc[2] += cw * ion_model.osm (ui);
      loc[3] += cw * g;
      const double d  = ion_model.D (ui);
      const double sc = ion_model.nu * ion_model_t::coshm1 (ui);
      loc[4] += cw * (std::expm1 (-ui) - sc) / d;
      loc[5] += cw * (std::expm1 (ui) - sc) / d;
      loc[6] += wd[i];
    }
  }

  // Box volume, for the preferential interaction coefficients
  // Gamma_pm = int_solvent (n_pm - n_b) dV - n_b V_excl,  V_excl = V_box - V_acc,
  // i.e. the ion excess over the whole box relative to bulk solution (the
  // definition of ion counting experiments, Bai et al. 2007).
  double vbox_loc = 0.0;
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    vbox_loc += (quadrant->p (0, 7) - quadrant->p (0, 0))
                * (quadrant->p (1, 7) - quadrant->p (1, 0))
                * (quadrant->p (2, 7) - quadrant->p (2, 0));

  double glob[7];
  std::copy (loc, loc + 7, glob);
  double vbox = vbox_loc;
  if (size > 1) {
    MPI_Allreduce (loc, glob, 7, MPI_DOUBLE, MPI_SUM, mpicomm);
    MPI_Allreduce (&vbox_loc, &vbox, 1, MPI_DOUBLE, MPI_SUM, mpicomm);
  }

  energy_exc = inv_4pi * glob[0];
  const double energy_half = inv_4pi * glob[1];
  const double energy_osm  = inv_4pi * glob[2];
  const double charge_ion  = -inv_4pi * glob[3];
  const double n_b    = 1.0e3 * N_av * ionic_strength * Angs * Angs * Angs;  // 1/A^3
  const double v_acc  = glob[6];
  const double v_excl = vbox - v_acc;
  const double gamma_p = 0.5 * inv_4pi * glob[4] - n_b * v_excl;
  const double gamma_m = 0.5 * inv_4pi * glob[5] - n_b * v_excl;

  if (rank == 0) {
    constexpr int label_width = 50;
    constexpr int precision = 16;

    std::cout << "\n============ [ Nonlinear excess ionic energy ] ============\n";
    std::cout << std::left << std::setw (label_width) << "  Ion model:"
              << ion_model.name ();
    if (ion_model.nu > 0.0)
      std::cout << "  (a = " << ion_size << " A, nu = " << ion_model.nu << ")";
    std::cout << "\n";
    std::cout << std::left << std::setw (label_width) << "  -1/2 int rho_s phi [kT]:"
              << std::setprecision (precision) << energy_half << "\n";
    std::cout << std::left << std::setw (label_width) << "  Osmotic term int (P-P0) [kT]:"
              << std::setprecision (precision) << energy_osm << "\n";
    std::cout << std::left << std::setw (label_width) << "  Excess ionic energy [kT]:"
              << std::setprecision (precision) << energy_exc << "\n";
    std::cout << std::left << std::setw (label_width) << "  Mobile ion charge [e]:"
              << std::setprecision (precision) << charge_ion << "\n";
    std::cout << std::left << std::setw (label_width) << "  Ion-excluded volume [A^3]:"
              << std::setprecision (precision) << v_excl << "\n";
    std::cout << std::left << std::setw (label_width) << "  Cation excess Gamma+ [ions]:"
              << std::setprecision (precision) << gamma_p << "\n";
    std::cout << std::left << std::setw (label_width) << "  Anion excess Gamma- [ions]:"
              << std::setprecision (precision) << gamma_m << "\n";
    std::cout << std::left << std::setw (label_width) << "  Total nonlinear free energy [kT]:"
              << std::setprecision (precision)
              << (energy_pol + energy_react + coul_energy + energy_exc) << "\n";
    std::cout << "===========================================================\n";
  }
}

void
poisson_boltzmann::write_potential_on_surface (ray_cache_t & ray_cache)
{
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  const double eps_in = 4.0 * pi * e_0 * e_in * kb * T * Angs / (e * e);
  const double eps_out = 4.0 * pi * e_0 * e_out * kb * T * Angs / (e * e);
  const double C_0 = 1.0e3 * N_av * ionic_strength;
  const double k2 = 2.0 * C_0 * Angs * Angs * e * e / (e_0 * e_out * kb * T);
  const double k = std::sqrt (k2);

  std::array<double, 3> V, N, h, area_h;
  std::array<double, 8> tmp_eps, tmp_phi;
  std::vector<int> edg, fl_dir;
  std::array<std::array<double, 3>, 3> vert_triangles, norms_vert;
  std::array<double, 3> phi_sup;

  int cubeindex = -1;
  int edge, i1 = 0, i2 = 0, ntriang = 0;
  double fract;
  double tmp_phi_1 = 0.0, tmp_phi_2 = 0.0;
  double tmp_eps_1 = 0.0, tmp_eps_2 = 0.0;

  auto quadrant = this->tmsh.begin_quadrant_sweep();

  // File setup
  // std::string filename_nodes = "phi_nodes_" + pqrfilename + ".txt";
  // std::string filename_surf  = "phi_surf_" + pqrfilename + ".txt";

  // std::ofstream phi_nodes_txt(filename_nodes);
  // std::ofstream phi_surf_txt(filename_surf);
  // FILE* phi_nod_delphi = std::fopen("filename_nodes_delphi.txt", "w");
  // FILE* phi_sup_delphi = std::fopen("filename_sup_delphi.txt", "w");
  std::vector<std::string> phi_nodes_local;
  std::vector<std::string> phi_surf_local;
  // std::vector<std::string> phi_nodes_delphi_local;
  // std::vector<std::string> phi_sup_delphi_local;

  for (const int ii : border_quad) {
    quadrant[ii];
    cubeindex = classifyCube (quadrant, eps_out);

    // Compute edge lengths and area scale factors
    h[0] = quadrant->p (0, 7) - quadrant->p (0, 0);
    h[1] = quadrant->p (1, 7) - quadrant->p (1, 0);
    h[2] = quadrant->p (2, 7) - quadrant->p (2, 0);
    area_h[0] = h[1] * h[2] / h[0] * 0.25;
    area_h[1] = h[0] * h[2] / h[1] * 0.25;
    area_h[2] = h[0] * h[1] / h[2] * 0.25;

    std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux (quadrant, tmp_phi, tmp_eps);
    ntriang = getTriangles (cubeindex, triangles);

    for (int t = 0; t < ntriang; ++t) {
      for (int j = 0; j < 3; ++j) {
        edge = triangles[t][j];
        i1 = edge2nodes[2 * edge];
        i2 = edge2nodes[2 * edge + 1];

        // Intersection and vertex
        V[0] = quadrant->p (0, i1);
        V[1] = quadrant->p (1, i1);
        V[2] = quadrant->p (2, i1);

        normal_intersection (quadrant, ray_cache, edge, N, fract);
        V[edge_axis[edge]] += fract * h[edge_axis[edge]];

        vert_triangles[j] = V;
        norms_vert[j] = N;

        // Interpolate phi/epsilon at i1
        if (!quadrant->is_hanging (i1)) {
          tmp_phi_1 = (*phi)[quadrant->gt (i1)];
          tmp_eps_1 = (*epsilon_nodes)[quadrant->gt (i1)];
        } else {
          tmp_phi_1 = tmp_eps_1 = 0.0;
          int np = quadrant->num_parents (i1);

          for (int k = 0; k < np; ++k) {
            tmp_phi_1 += (*phi)[quadrant->gparent (k, i1)] / np;
            tmp_eps_1 += (*epsilon_nodes)[quadrant->gparent (k, i1)] / np;
          }
        }

        // Interpolate phi/epsilon at i2
        if (!quadrant->is_hanging (i2)) {
          tmp_phi_2 = (*phi)[quadrant->gt (i2)];
          tmp_eps_2 = (*epsilon_nodes)[quadrant->gt (i2)];
        } else {
          tmp_phi_2 = tmp_eps_2 = 0.0;
          int np = quadrant->num_parents (i2);

          for (int k = 0; k < np; ++k) {
            tmp_phi_2 += (*phi)[quadrant->gparent (k, i2)] / np;
            tmp_eps_2 += (*epsilon_nodes)[quadrant->gparent (k, i2)] / np;
          }
        }

        // Interpolate phi on surface
        phi_sup[j] = phi0 (tmp_eps_1, tmp_eps_2, tmp_phi_1, tmp_phi_2, fract);

        std::ostringstream oss;
        oss << std::scientific << std::setprecision (5)
            << quadrant->p (0, i1) << " "
            << quadrant->p (1, i1) << " "
            << quadrant->p (2, i1) << " "
            << tmp_phi_1 << "\n"
            << quadrant->p (0, i2) << " "
            << quadrant->p (1, i2) << " "
            << quadrant->p (2, i2) << " "
            << tmp_phi_2;
        phi_nodes_local.push_back (oss.str());
        oss.str ("");
        oss.clear();
        oss << std::scientific << std::setprecision (5)
            << V[0] << " "
            << V[1] << " "
            << V[2] << " "
            << phi_sup[j];
        phi_surf_local.push_back (oss.str());
        oss.str ("");
        oss.clear();
        // // Write to ASCII
        // phi_nodes_txt << quadrant->p(0, i1) << "  " << quadrant->p(1, i1) << "  " << quadrant->p(2, i1) << "  " << tmp_phi_1 << "\n";
        // phi_nodes_txt << quadrant->p(0, i2) << "  " << quadrant->p(1, i2) << "  " << quadrant->p(2, i2) << "  " << tmp_phi_2 << "\n";

        // phi_surf_txt  << V[0] << "  " << V[1] << "  " << V[2] << "  " << phi_sup[j] << "\n";

        // // Write to Delphi PDB-like format
        // std::fprintf(phi_nod_delphi,
        // "\nATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%8.4f%8.4f",
        // 1, "X", "XXX", " ", 0,
        // quadrant->p(0, i1), quadrant->p(1, i1), quadrant->p(2, i1), tmp_phi_1, tmp_phi_2);

        // std::fprintf(phi_nod_delphi,
        // "\nATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%8.4f%8.4f",
        // 1, "X", "XXX", " ", 0,
        // quadrant->p(0, i2), quadrant->p(1, i2), quadrant->p(2, i2), tmp_phi_1, tmp_phi_2);

        // std::fprintf(phi_sup_delphi,
        // "\nATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%8.4f%8.4f",
        // 1, "X", "XXX", " ", 0,
        // V[0], V[1], V[2], phi_sup[j], 0.0);
      }
    }
  }

  auto gather_and_write = [&] (const std::string& filename,
  const std::vector<std::string>& local_lines) {
    if (rank == 0) {
      std::ofstream ofs (filename);

      for (const auto& line : local_lines)
        ofs << line << "\n";

      for (int r = 1; r < size; ++r) {
        int n_lines;
        MPI_Recv (&n_lines, 1, MPI_INT, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int i = 0; i < n_lines; ++i) {
          char buf[512];
          MPI_Recv (buf, 512, MPI_CHAR, r, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          ofs << buf << "\n";
        }
      }

      ofs.close();
    } else {
      int n_lines = static_cast<int> (local_lines.size());
      MPI_Send (&n_lines, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

      for (const auto& line : local_lines) {
        MPI_Send (line.c_str(), static_cast<int> (line.size()) + 1,
                  MPI_CHAR, 0, 1, MPI_COMM_WORLD);
      }
    }
  };

  gather_and_write ("phi_nodes.txt", phi_nodes_local);
  gather_and_write ("phi_surf.txt", phi_surf_local);
  // Close files
  // phi_nodes_txt.close();
  // phi_surf_txt.close();
  // std::fclose(phi_nod_delphi);
  // std::fclose(phi_sup_delphi);
}


double
poisson_boltzmann::coulomb_boundary_conditions (double x, double y, double z)
{
  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/ (e*e); //adim e_out
  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/ (e_0*e_out*kb*T);

  double dist = 0.0;
  double pot = 0.0;
  double k = std::sqrt (k2);

  
  if (pos_atoms.size () != charge_atoms.size ()) {
    return 0.0;
  }

  for (std::size_t i = 0; i < charge_atoms.size (); ++i) {
    dist = std::hypot (pos_atoms[i][0] - x,
                      pos_atoms[i][1] - y,
                      pos_atoms[i][2] - z);

    if (dist > 1.0e-12) {
      pot += charge_atoms[i] * std::exp (-k * dist)
            / (dist * eps_out);
    }
  }
  return pot;
}

double
poisson_boltzmann::analytic_solution (double x, double y, double z)
{
  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/ (e*e); //adim e_out
  double eps_in = 4.0*pi*e_0*e_in*kb*T*Angs/ (e*e); //adim e_in
  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/ (e_0*e_out*kb*T);
  double k = std::sqrt (k2);
  double rs = 2.0;
  double pippo_out = eps_out* (1.0+k*rs);
  double pippo_in = eps_in*rs;
  double dist= 0.0;
  double pot = 0.0;

  dist = std::hypot (x, y, z);

  if (dist <= rs)
    pot = 1.0/ (pippo_out*rs) + (rs-dist)/ (pippo_in*dist);
  else
    pot = exp (k* (rs-dist))/ (dist*pippo_out);

  return pot;
}



void
poisson_boltzmann::analitic_potential ()
{
  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/ (e*e); //adim e_out
  double eps_in = 4.0*pi*e_0*e_in*kb*T*Angs/ (e*e); //adim e_in
  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/ (e_0*e_out*kb*T);
  double k = std::sqrt (k2);
  double rs = 2.0;
  double pippo_out = eps_out* (1.0+k*rs);
  double pippo_in = eps_in*rs;
  double dist= 0.0;


  distributed_vector phi_an (tmsh.num_owned_nodes (), mpicomm);
  // distributed_vector field_an (tmsh.num_owned_nodes (), mpicomm);
  bim3a_solution_with_ghosts (tmsh, phi_an, replace_op);
  // bim3a_solution_with_ghosts (tmsh, field_an, replace_op);

  for (auto quadrant = this->tmsh.begin_quadrant_sweep ();
       quadrant != this->tmsh.end_quadrant_sweep ();
       ++quadrant) {
    for (const NS::Atom& i : atoms) {
      if (std::fabs (i.charge) > std::numeric_limits<double>::epsilon ()) {
        for (int ii = 0; ii < 8; ++ii) {
          if (! quadrant->is_hanging (ii)) {
            dist = std::hypot ( (i.pos[0] - quadrant->p (0,ii)),
                                (i.pos[1] - quadrant->p (1,ii)),
                                (i.pos[2] - quadrant->p (2,ii)));

            if (dist <= rs) {
              phi_an[quadrant->gt (ii)] = 1.0/ (pippo_out*rs) +
                                          i.charge* (rs-dist)/ (pippo_in*dist);
              // field_an[quadrant->gt (ii)] = (1.0 + (10.0-dist)/dist)*
              // i.charge/(pippo_in*dist);
            } else
              phi_an[quadrant->gt (ii)] = i.charge*exp (k* (rs-dist))/ (dist*pippo_out);

            // field_an[quadrant->gt (ii)] = phi_an[quadrant->gt (ii)]*(k +1.0/dist);
          }

          // else{
          // for (int jj = 0; jj < quadrant->num_parents (ii); ++jj){
          // phi_an[quadrant->gparent (jj, ii)] += 0.0;
          // field_an[quadrant->gparent (jj, ii)] += 0.0;
          // }
          // }
        }
      }
    }
  }

  phi_an.assemble (replace_op);
  tmsh.octbin_export ("phi_an_0", phi_an);
  // field_an.assemble (replace_op);
  // tmsh.octbin_export ("field_an_0", field_an);
}

bool
poisson_boltzmann::controlla_coordinate (int i, const p8est_quadrant_t *quadrant)
{

  double tol = p4esttol * (rr[0]-ll[0]);

  if (mesh_shape == 2)
    tol = p4esttol * (r_c[0]-l_c[0]);


  bool retval = false;

  p8est_quadrant_t node;
  int ii, jj;
  double vxyz [3 * 8] = {0,0,0, 0,0,0, 0,0,0, 0,0,0,
                         0,0,0, 0,0,0, 0,0,0, 0,0,0
                        };

  for (ii = 0; ii < 8; ++ii) {
    p8est_quadrant_corner_node (quadrant, ii, &node);
    p8est_qcoord_to_vertex (tmsh.p8est->connectivity, 0,
                            node.x, node.y, node.z, & (vxyz[3 * ii]));
  }

  double l, r, t, b, f, bk;

  l = vxyz[0];
  r = vxyz[3*7];

  f = vxyz[1];
  bk = vxyz[3*7 +1];

  b = vxyz[2];
  t = vxyz[3*7 +2];

  // retval = (atoms[i].pos[0] > l- tol) && (atoms[i].pos[0] <= r- tol); //make sure that the charge is assigned only once
  // retval = retval && (atoms[i].pos[1] > f- tol) && (atoms[i].pos[1] <= bk- tol);
  // retval = retval && (atoms[i].pos[2] > b- tol) && (atoms[i].pos[2] <= t- tol);
  retval = (pos_atoms[i][0] > l- tol) && (pos_atoms[i][0] <= r- tol); //make sure that the charge is assigned only once
  retval = retval && (pos_atoms[i][1] > f- tol) && (pos_atoms[i][1] <= bk- tol);
  retval = retval && (pos_atoms[i][2] > b- tol) && (pos_atoms[i][2] <= t- tol);

  return retval;
}


int
poisson_boltzmann::cerca_atomo (p8est_t * p4est,
                                p4est_topidx_t which_tree,
                                p8est_quadrant_t * quadrant,
                                p4est_locidx_t local_num,
                                void *point)
{
  int *pt = (int *) point;
  // std::cout << "\n pt: " <<*pt <<"\n"<< std::endl;
  bool tf = controlla_coordinate (*pt, quadrant);

  if (tf) {
    if (local_num >= 0) {
      auto quadrant = this->tmsh.begin_quadrant_sweep ();
      quadrant[local_num];
      tmesh_3d::quadrant_t qi = tmsh.current_quadrant;
      lookup_table.emplace (*pt, qi);

    }

    return 1;
  }

  return 0;

}

void 
poisson_boltzmann::search_points()
{
    size_t count = charge_atoms.size();
    auto base = std::make_unique<int[]>(count);

    for (size_t ii = 0; ii < count; ++ii)
        base[ii] = ii;

    sc_array_t* points = sc_array_new_data(
        base.get(), sizeof(int), count);

    // Set the global pointer for the callback
    pb_global_wrapper = this;

    // call the search function
    p8est_search_local(
        tmsh.p8est,
        0,
        NULL,
        cerca_atomo_wrapper,
        points
    );
}

void
poisson_boltzmann::write_dataset (ray_cache_t & ray_cache)
{
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0)
    std::cout << "\n================ [ Computing vertex quantities ] =================\n";

  const double eps_in  = 4.0 * pi * e_0 * e_in  * kb * T * Angs / (e * e);
  const double eps_out = 4.0 * pi * e_0 * e_out * kb * T * Angs / (e * e);

  std::array<double, 3> N;
  std::array<double, 3> V;
  std::array<double, 3> h;
  std::array<double, 8> tmp_eps;
  std::array<double, 8> tmp_phi;
  std::vector<int> edg;
  std::vector<int> fl_dir;

  int cubeindex = -1;
  std::unordered_map<EdgeKey, VertexData, EdgeHash> edgeMap;

  auto quadrant = this->tmsh.begin_quadrant_sweep ();

  if (!border_quad.empty ()) {
    quadrant[border_quad[0]];
    h[0] = quadrant->p (0, 7) - quadrant->p (0, 0);
    h[1] = quadrant->p (1, 7) - quadrant->p (1, 0);
    h[2] = quadrant->p (2, 7) - quadrant->p (2, 0);
  }

  for (const int ii : border_quad) {
    quadrant[ii];
    cubeindex = classifyCube_fast (quadrant, eps_out);
    std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux_fast (quadrant, tmp_phi, tmp_eps);

    for (int ip = 0; ip < (int)edg.size (); ++ip) {
      const int edge = edg[ip];
      const int axis = edge_axis[edge];
      const int i1   = edge2nodes[2 * edge];
      const int i2   = edge2nodes[2 * edge + 1];

      EdgeKey triEdges{
        static_cast<int>(std::min (quadrant->gt (i1), quadrant->gt (i2))),
        static_cast<int>(std::max (quadrant->gt (i1), quadrant->gt (i2)))
      };

      auto & vertexData = edgeMap[triEdges];

      vertexData.axis = axis;

      double fract = 0.0;
      normal_intersection (quadrant, ray_cache, edge, N, fract);

      for (int jj = 0; jj < 3; ++jj) {
        int nu = (axis + jj) % 3;
        vertexData.N[jj] = N[nu];
      }

      vertexData.pos1[0] = quadrant->p (0, i1);
      vertexData.pos1[1] = quadrant->p (1, i1);
      vertexData.pos1[2] = quadrant->p (2, i1);

      vertexData.pos0[0] = vertexData.pos1[0];
      vertexData.pos0[1] = vertexData.pos1[1];
      vertexData.pos0[2] = vertexData.pos1[2];
      vertexData.pos0[edge_axis[edge]] += fract * h[edge_axis[edge]];

      vertexData.phi1 = tmp_phi[i1];
      vertexData.phi2 = tmp_phi[i2];
      vertexData.eps1 = tmp_eps[i1];
      vertexData.eps2 = tmp_eps[i2];

      vertexData.alpha = fract;
      vertexData.phi0  = phi0 (tmp_eps[i1], tmp_eps[i2], tmp_phi[i1], tmp_phi[i2], fract);

      const int i1_nn1   = edge2nodes_nn1[2 * edge];
      const int i1_nn2   = edge2nodes_nn1[2 * edge + 1];
      const int i2_nn1   = edge2nodes_nn2[2 * edge];
      const int i2_nn2   = edge2nodes_nn2[2 * edge + 1];
      const int ind1_nn1 = edge2index_nn1[2 * edge];
      const int ind1_nn2 = edge2index_nn1[2 * edge + 1];
      const int ind2_nn1 = edge2index_nn2[2 * edge];
      const int ind2_nn2 = edge2index_nn2[2 * edge + 1];

      vertexData.phi1_nn[ind1_nn1] = tmp_phi[i1_nn1];
      vertexData.phi1_nn[ind1_nn2] = tmp_phi[i1_nn2];
      vertexData.eps1_nn[ind1_nn1] = tmp_eps[i1_nn1];
      vertexData.eps1_nn[ind1_nn2] = tmp_eps[i1_nn2];
      vertexData.phi2_nn[ind2_nn1] = tmp_phi[i2_nn1];
      vertexData.phi2_nn[ind2_nn2] = tmp_phi[i2_nn2];
      vertexData.eps2_nn[ind2_nn1] = tmp_eps[i2_nn1];
      vertexData.eps2_nn[ind2_nn2] = tmp_eps[i2_nn2];
    }
  }

  {
    std::ostringstream fname;
    fname << "vertexdata_rank" << rank << ".csv";
    std::ofstream ofs (fname.str ());
    if (!ofs) {
      std::cerr << "Error: cannot open " << fname.str () << " for writing\n";
      return;
    }

    ofs << "phi1,eps1,alpha,N_nu,N_nu1,N_nu2,"
        << "phi2,eps2,"
        << "phi_perp_1_p_1,eps_perp_1_p_1,"
        << "phi_perp_1_m_1,eps_perp_1_m_1,"
        << "phi_perp_2_p_1,eps_perp_2_p_1,"
        << "phi_perp_2_m_1,eps_perp_2_m_1,"
        << "phi_perp_1_p_2,eps_perp_1_p_2,"
        << "phi_perp_1_m_2,eps_perp_1_m_2,"
        << "phi_perp_2_p_2,eps_perp_2_p_2,"
        << "phi_perp_2_m_2,eps_perp_2_m_2,"
        << "x0,y0,z0,phi0,x1,y1,z1,axis\n";

    ofs << std::scientific << std::setprecision (8);

    for (const auto & [k, vd] : edgeMap) {
      ofs << vd.phi1 << "," << vd.eps1 << "," << vd.alpha << ","
          << vd.N[0] << "," << vd.N[1] << "," << vd.N[2] << ","
          << vd.phi2 << "," << vd.eps2 << ",";

      ofs << vd.phi1_nn[0] << "," << vd.eps1_nn[0] << ","
          << vd.phi1_nn[1] << "," << vd.eps1_nn[1] << ","
          << vd.phi1_nn[2] << "," << vd.eps1_nn[2] << ","
          << vd.phi1_nn[3] << "," << vd.eps1_nn[3] << ",";

      ofs << vd.phi2_nn[0] << "," << vd.eps2_nn[0] << ","
          << vd.phi2_nn[1] << "," << vd.eps2_nn[1] << ","
          << vd.phi2_nn[2] << "," << vd.eps2_nn[2] << ","
          << vd.phi2_nn[3] << "," << vd.eps2_nn[3] << ",";

      ofs << vd.pos0[0] << "," << vd.pos0[1] << "," << vd.pos0[2] << "," << vd.phi0 << ",";
      ofs << vd.pos1[0] << "," << vd.pos1[1] << "," << vd.pos1[2] << "," << vd.axis << "\n";
    }
  }

  MPI_Barrier (MPI_COMM_WORLD);

  if (rank == 0) {
    std::ofstream final ("vertexdata.csv");
    if (!final) {
      std::cerr << "Error: cannot open vertexdata.csv for writing\n";
      return;
    }

    final << "phi1,eps1,alpha,N_nu,N_nu1,N_nu2,"
          << "phi2,eps2,"
          << "phi_perp_1_p_1,eps_perp_1_p_1,"
          << "phi_perp_1_m_1,eps_perp_1_m_1,"
          << "phi_perp_2_p_1,eps_perp_2_p_1,"
          << "phi_perp_2_m_1,eps_perp_2_m_1,"
          << "phi_perp_1_p_2,eps_perp_1_p_2,"
          << "phi_perp_1_m_2,eps_perp_1_m_2,"
          << "phi_perp_2_p_2,eps_perp_2_p_2,"
          << "phi_perp_2_m_2,eps_perp_2_m_2,"
          << "x0,y0,z0,phi0,x1,y1,z1,axis\n";

    final << std::scientific << std::setprecision (8);

    for (int r = 0; r < size; ++r) {
      std::ostringstream fname;
      fname << "vertexdata_rank" << r << ".csv";
      std::ifstream ifs (fname.str ());
      if (!ifs) continue;

      std::string line;
      bool first_line = true;
      while (std::getline (ifs, line)) {
        if (first_line) { first_line = false; continue; }
        final << line << "\n";
      }
      ifs.close ();
      std::filesystem::remove (fname.str ());
    }

    std::cout << "\nAll data merged into vertexdata.csv\n";
  }
}
