/*
 *  Copyright (C) 2021-2026 Vincenzo Di Florio
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

#include "pb_membrane.h"
#include "readpdb.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <mpi.h>

// ─── Lipid I/O ───────────────────────────────────────────────────────────────

void
read_lipids (poisson_boltzmann &pb)
{
  pb.lipid_atoms.clear ();
  pb.pos_lipid_atoms.clear ();
  pb.charge_lipid_atoms.clear ();
  pb.r_lipid_atoms.clear ();

  if (pb.lipid_filetype == "pqr") {
    std::ifstream f (pb.lipid_file);

    if (!f)
      throw std::runtime_error ("read_lipids: cannot open " + pb.lipid_file);

    NS::Atom a;

    while (f >> a)
      pb.lipid_atoms.push_back (a);

  } else if (pb.lipid_filetype == "pdb") {
    read_pdb (pb.lipid_file, pb.lipid_atoms);
    load_radii (pb.radiusfilename, pb.lipid_atoms);
    load_charges (pb.chargefilename, pb.lipid_atoms);

  } else {
    throw std::runtime_error ("read_lipids: unknown filetype '" + pb.lipid_filetype + "'");
  }

  for (const NS::Atom& a : pb.lipid_atoms) {
    pb.pos_lipid_atoms.push_back ({a.pos[0], a.pos[1], a.pos[2]});
    pb.charge_lipid_atoms.push_back (a.charge);
    pb.r_lipid_atoms.push_back (a.radius);
  }

  int rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  if (rank == 0)
    std::cout << "  Lipid atoms read     : " << pb.lipid_atoms.size () << '\n';
}

void
broadcast_lipid_vectors (poisson_boltzmann &pb)
{
  int n = static_cast<int> (pb.charge_lipid_atoms.size ());
  MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  pb.pos_lipid_atoms.resize (n);
  pb.charge_lipid_atoms.resize (n);
  pb.r_lipid_atoms.resize (n);

  MPI_Bcast (pb.pos_lipid_atoms.data (), n * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (pb.charge_lipid_atoms.data (), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (pb.r_lipid_atoms.data (), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Rebuild lipid_atoms on non-zero ranks from the broadcast flat arrays
  pb.lipid_atoms.resize (n);

  for (int i = 0; i < n; ++i) {
    pb.lipid_atoms[i].pos[0] = pb.pos_lipid_atoms[i][0];
    pb.lipid_atoms[i].pos[1] = pb.pos_lipid_atoms[i][1];
    pb.lipid_atoms[i].pos[2] = pb.pos_lipid_atoms[i][2];
    pb.lipid_atoms[i].charge = pb.charge_lipid_atoms[i];
    pb.lipid_atoms[i].radius = pb.r_lipid_atoms[i];
  }
}

// ─── NanoShaper supercell ────────────────────────────────────────────────────

std::vector<NS::Atom>
build_ns_supercell (const poisson_boltzmann &pb)
{
  // Start with the protein atoms
  std::vector<NS::Atom> supercell = pb.atoms;

  // Append lipid atoms of the central cell only.
  // Periodicity in xy is handled by applying periodic BCs on the potential,
  // not by replicating the geometry into NanoShaper.
  for (const NS::Atom& a : pb.lipid_atoms)
    supercell.push_back (a);

  return supercell;
}

// ─── Post-NanoShaper lipid trimming ──────────────────────────────────────────

void
trim_lipid_atoms (poisson_boltzmann &pb)
{
  // Remove lipid atoms that lie outside the computational domain in xy.
  // These atoms were included in the NanoShaper supercell so that the
  // membrane surface closes at the periodic boundary, but they must not
  // contribute to the FEM assembly.
  const double x_lo = pb.l_cr[0], x_hi = pb.r_cr[0];
  const double y_lo = pb.l_cr[1], y_hi = pb.r_cr[1];

  std::vector<NS::Atom> kept;
  kept.reserve (pb.lipid_atoms.size ());
  std::size_t out_x = 0, out_y = 0, out_xy = 0;

  for (const NS::Atom& a : pb.lipid_atoms) {
    bool in_x = (a.pos[0] >= x_lo && a.pos[0] <= x_hi);
    bool in_y = (a.pos[1] >= y_lo && a.pos[1] <= y_hi);

    if (in_x && in_y) {
      kept.push_back (a);
    } else {
      if (!in_x && !in_y) ++out_xy;
      else if (!in_x) ++out_x;
      else ++out_y;
    }
  }

  int rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    std::cout << "\n=== [ Lipid atom xy-trim debug ] ===\n";
    std::cout << "  FEM domain  x : [" << x_lo << ", " << x_hi << "] Å\n";
    std::cout << "  FEM domain  y : [" << y_lo << ", " << y_hi << "] Å\n";
    std::cout << "  Lipid atoms before trim : " << pb.lipid_atoms.size () << '\n';
    std::cout << "  Removed (outside x only): " << out_x << '\n';
    std::cout << "  Removed (outside y only): " << out_y << '\n';
    std::cout << "  Removed (outside x & y) : " << out_xy << '\n';
    std::cout << "  Lipid atoms after trim  : " << kept.size () << '\n';
    std::cout << "====================================\n\n";
  }

  pb.lipid_atoms = std::move (kept);

  // Rebuild flat vectors to stay in sync with lipid_atoms
  pb.pos_lipid_atoms.clear ();
  pb.charge_lipid_atoms.clear ();
  pb.r_lipid_atoms.clear ();
  pb.pos_lipid_atoms.reserve (pb.lipid_atoms.size ());
  pb.charge_lipid_atoms.reserve (pb.lipid_atoms.size ());
  pb.r_lipid_atoms.reserve (pb.lipid_atoms.size ());

  for (const NS::Atom& a : pb.lipid_atoms) {
    pb.pos_lipid_atoms.push_back ({a.pos[0], a.pos[1], a.pos[2]});
    pb.charge_lipid_atoms.push_back (a.charge);
    pb.r_lipid_atoms.push_back (a.radius);
  }
}

void
zero_boundary_residue_charges (poisson_boltzmann &pb)
{
  // Zero the charge of every atom belonging to a membrane residue that has
  // at least one atom within (max_radius + probe_radius) of the xy domain
  // boundary.  This prevents singularities near the periodic faces and
  // ensures postprocessing consistency.
  if (pb.lipid_atoms.empty ()) return;

  const double x_lo = pb.l_cr[0], x_hi = pb.r_cr[0];
  const double y_lo = pb.l_cr[1], y_hi = pb.r_cr[1];

  double max_r = 0.0;

  for (const NS::Atom& a : pb.lipid_atoms)
    max_r = std::max (max_r, static_cast<double> (a.radius));

  const double cutoff = max_r + pb.surf_param;

  // Collect residue IDs (chain string + resNum) that straddle the boundary.
  // Also store (resName, chain, resNum) for the diagnostic print.
  using ResID = std::pair<std::string, int>;
  std::set<ResID> boundary_res;
  // Map ResID → resName for the diagnostic listing
  std::map<ResID, std::string> res_name_map;

  for (const NS::Atom& a : pb.lipid_atoms) {
    if (a.pos[0] < x_lo + cutoff || a.pos[0] > x_hi - cutoff ||
        a.pos[1] < y_lo + cutoff || a.pos[1] > y_hi - cutoff) {
      ResID rid {a.ai.chain, a.ai.resNum};
      boundary_res.insert (rid);
      res_name_map.emplace (rid, a.ai.resName);
    }
  }

  // Zero charges of all atoms in those residues
  std::size_t n_zeroed = 0;

  for (std::size_t i = 0; i < pb.lipid_atoms.size (); ++i) {
    NS::Atom &a = pb.lipid_atoms[i];

    if (boundary_res.count ({a.ai.chain, a.ai.resNum})) {
      a.charge = 0.0;
      pb.charge_lipid_atoms[i] = 0.0;
      ++n_zeroed;
    }
  }

  int rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    std::cout << "\n=== [ Boundary residue charge-zeroing debug ] ===\n";
    std::cout << "  Cutoff (max_r + probe_r): " << cutoff << " Å\n";
    std::cout << "  FEM domain  x : [" << x_lo << ", " << x_hi << "] Å\n";
    std::cout << "  FEM domain  y : [" << y_lo << ", " << y_hi << "] Å\n";
    std::cout << "  Zeroed residues (" << boundary_res.size () << "):\n";

    for (const auto& rid : boundary_res) {
      const std::string &rname = res_name_map.at (rid);
      std::cout << "    " << rname
                << "  chain=" << (rid.first.empty () ? "-" : rid.first)
                << "  resNum=" << rid.second << '\n';
    }

    std::cout << "  Total atoms zeroed: " << n_zeroed << '\n';
    std::cout << "=================================================\n\n";
  }
}

// ─── Mortar assembly ─────────────────────────────────────────────────────────

void
assemble_mortar_block (distributed_sparse_matrix        &A,
                       tmesh_3d::quadrant_iterator       q,
                       int                               face_idx,
                       const MortarFacePair             &pair,
                       double                            sign,
                       int                               nm)
{
  const auto &fnodes = q->ef (face_idx);
  if (fnodes.empty ()) return;

  const int d0 = pair.d0, d1 = pair.d1;

  // Step 1 — tangential coordinates of the 4 face nodes
  // fnodes layout (same for all faces in d0/d1):
  //   [0]: (d0_min, d1_min)  [1]: (d0_max, d1_min)
  //   [2]: (d0_min, d1_max)  [3]: (d0_max, d1_max)
  double c0[2], c1[2];
  c0[0] = q->p (d0, fnodes[0]);  // d0_min
  c0[1] = q->p (d0, fnodes[1]);  // d0_max
  c1[0] = q->p (d1, fnodes[0]);  // d1_min
  c1[1] = q->p (d1, fnodes[2]);  // d1_max

  // Step 2 — Dirichlet filter: skip nodes at z boundaries (d1==2 for both pairs)
  const double z_lo = pair.d1_min;
  const double z_hi = pair.d1_min + pair.L1;
  const double tol  = pair.L1 * 1e-10;
  bool keep[4];
  for (int fi = 0; fi < 4; ++fi) {
    double zi = q->p (2, fnodes[fi]);
    keep[fi] = (zi > z_lo + tol) && (zi < z_hi - tol);
  }

  // Step 3 — mortar cell index (physical offset accounted for)
  double hm0 = pair.L0 / nm;
  double hm1 = pair.L1 / nm;
  int j0 = static_cast<int> ((c0[0] - pair.d0_min) / hm0);
  int j1 = static_cast<int> ((c1[0] - pair.d1_min) / hm1);
  j0 = std::min (j0, nm - 1);
  j1 = std::min (j1, nm - 1);

  // Step 4 — 1D mixed hat-function mass matrices (2×2 each)
  // M[i][jl] = ∫ φ_i(t) ψ_jl(t) dt
  // φ_i : FEM hat on [c_min, c_max], values (1,0) for i=0, (0,1) for i=1
  // ψ_jl: mortar hat, evaluated at c_min and c_max → psi[jl][0], psi[jl][1]
  double orig0 = pair.d0_min + j0 * hm0;
  double orig1 = pair.d1_min + j1 * hm1;

  double psi_0[2][2], psi_1[2][2];
  psi_0[0][0] = (hm0 - (c0[0] - orig0)) / hm0;
  psi_0[0][1] = (hm0 - (c0[1] - orig0)) / hm0;
  psi_0[1][0] =         (c0[0] - orig0)  / hm0;
  psi_0[1][1] =         (c0[1] - orig0)  / hm0;

  psi_1[0][0] = (hm1 - (c1[0] - orig1)) / hm1;
  psi_1[0][1] = (hm1 - (c1[1] - orig1)) / hm1;
  psi_1[1][0] =         (c1[0] - orig1)  / hm1;
  psi_1[1][1] =         (c1[1] - orig1)  / hm1;

  double h0 = c0[1] - c0[0], h1 = c1[1] - c1[0];
  double M_0[2][2], M_1[2][2];

  // Exact: ∫_a^b f·g dt = (b-a)/6·(2·f_a·g_a + f_a·g_b + f_b·g_a + 2·f_b·g_b)
  for (int jl = 0; jl < 2; ++jl) {
    M_0[0][jl] = h0 / 6.0 * (2.0 * psi_0[jl][0] +       psi_0[jl][1]);
    M_0[1][jl] = h0 / 6.0 * (      psi_0[jl][0] + 2.0 * psi_0[jl][1]);
    M_1[0][jl] = h1 / 6.0 * (2.0 * psi_1[jl][0] +       psi_1[jl][1]);
    M_1[1][jl] = h1 / 6.0 * (      psi_1[jl][0] + 2.0 * psi_1[jl][1]);
  }
  // Trapezi (per confronto, O(h^2) di errore di quadratura):
  // for (int jl = 0; jl < 2; ++jl) {
  //   M_0[0][jl] = h0 / 2.0 * psi_0[jl][0];   // φ_0 vale 1 in c_min, 0 in c_max
  //   M_0[1][jl] = h0 / 2.0 * psi_0[jl][1];   // φ_1 vale 0 in c_min, 1 in c_max
  //   M_1[0][jl] = h1 / 2.0 * psi_1[jl][0];
  //   M_1[1][jl] = h1 / 2.0 * psi_1[jl][1];
  // }

  // Step 5 — accumulate C and C^T into A (tensor product, += element assembly)
  // d0/d1 local sub-index of each face node (same pattern for all faces)
  static const int ii0[4] = {0, 1, 0, 1};
  static const int ii1[4] = {0, 0, 1, 1};

  for (int fi = 0; fi < 4; ++fi) {
    if (!keep[fi]) continue;
    size_t row = static_cast<size_t> (q->gt (fnodes[fi]));

    for (int dj0 = 0; dj0 < 2; ++dj0)
      for (int dj1 = 0; dj1 < 2; ++dj1) {
        size_t col = pair.mortar_offset
                     + static_cast<size_t> ((j1 + dj1) * (nm + 1) + (j0 + dj0));
        double val = sign * M_0[ii0[fi]][dj0] * M_1[ii1[fi]][dj1];
        A[row][col] += val;   // C block:   FEM row,    mortar col
        A[col][row] += val;   // C^T block: mortar row, FEM col
      }
  }
}


