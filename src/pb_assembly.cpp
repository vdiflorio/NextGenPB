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
#include "pb_membrane.h"
#include <bim_distributed_vector.h>
#include <quad_operators_3d.h>
#include <mumps_class.h>
#include <lis_class.h>
#include <cmath>

void
poisson_boltzmann::create_density_map (ray_cache_t &ray_cache)
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
  this->rho_fixed = std::make_unique<distributed_vector> (tmsh.num_owned_nodes(), mpicomm);
  this->rho_fixed->get_owned_data().assign (tmsh.num_owned_nodes(), 0.0);

  // ------------------------------------------------------------
  // Vector of ones at nodes
  // ------------------------------------------------------------
  this->ones = std::make_unique<distributed_vector> (tmsh.num_owned_nodes(), mpicomm);
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
    bim3a_solution_with_ghosts (tmsh, *ones, replace_op);

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
  // Sync ghost nodes of rho_fixed in parallel runs
  // ------------------------------------------------------------
  if (size > 1)
    bim3a_solution_with_ghosts (tmsh, *rho_fixed);
}

void
poisson_boltzmann::assemple_system_matrix (ray_cache_t &ray_cache)
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
  auto func_frac = [&] (tmesh_3d::quadrant_iterator &quadrant) {
    return cube_fraction_intersection (quadrant, ray_cache);
  };

  if (rank == 0) {
    std::cout << "\n========== [ System Assembly ] ==========\n";
    if (periodic_x || periodic_y) {
      std::cout << "  PBC active:  x=" << periodic_x << "  y=" << periodic_y << '\n';
      std::cout << "  Mortar level (minlevel) : " << minlevel << '\n';
      std::cout << "  nm = 2^minlevel         : " << (1 << minlevel) << " cells/direction\n";
    } else {
      std::cout << "  PBC: disabled\n";
    }
  }

  check_pbc_face_conformity ();
  if (periodic_x || periodic_y)
    build_pbc_node_map ();

  // Allocate sparse matrix A and RHS vector
  A = std::make_unique<distributed_sparse_matrix> (mpicomm);
  A->set_ranges (tmsh.num_owned_nodes());

  // === [MORTAR BLOCK — DISABLED FOR STRONG PBC] Zone A: augmented matrix resize ===
  // if (periodic_x || periodic_y) {
  //   int nm     = (1 << minlevel);
  //   ndofm      = (nm + 1) * (nm + 1);
  //   int n_pairs = (periodic_x ? 1 : 0) + (periodic_y ? 1 : 0);
  //   size_t N_phys = A->rows ();
  //   A->resize (N_phys + n_pairs * ndofm);
  //   A->m = static_cast<int> (tmsh.num_global_nodes ()) + n_pairs * ndofm;
  //   if (rank == 0) {
  //     std::cout << "  Physical DOFs (N_phys)  : " << N_phys << '\n';
  //     std::cout << "  Active mortar pairs     : " << n_pairs << '\n';
  //     std::cout << "  Mortar DOFs per pair    : " << ndofm << "  (= (nm+1)^2)\n";
  //     std::cout << "  Augmented system size   : " << N_phys + n_pairs * ndofm << '\n';
  //   }
  // } else {
  //   if (rank == 0)
  //     std::cout << "  Physical DOFs           : " << A->rows () << '\n';
  // }
  // === [END MORTAR BLOCK Zone A] ===

  if (rank == 0)
    std::cout << "  Physical DOFs           : " << A->rows () << '\n';

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

  // Add reaction term (Stern layer present or fractional formulation)
  if (stern_layer_surf == 1) {
    // Standard Stern-layer reaction
    bim3a_reaction (tmsh, reaction, *ones, *A);
  } else {
    // Fractional reaction for intersected cells
    bim3a_reaction_frac (tmsh, (*reaction_nodes), *ones, *A, func_frac);
  }

  // --- Membrane applied potential: RHS source +c·φ̄ ---
  // The membrane LPBE reaction term is c(φ − φ̄): the operator c·φ is already in
  // A (above), while the known source c·φ̄ belongs on the RHS. The box-method
  // reaction is lumped (diagonal M_c), so the source equals M_c·φ̄ exactly when
  // assembled with the SAME operator and fractional weights as the matrix term.
  // Accumulated into a zeroed temporary, then added to *rhs before the single
  // rhs->assemble() below. Skipped entirely when phi_bar is null (V̄ = 0).
  if (phi_bar) {
    distributed_vector rhs_src (tmsh.num_owned_nodes (), mpicomm);
    rhs_src.get_owned_data ().assign (tmsh.num_owned_nodes (), 0.0);
    if (stern_layer_surf == 1)
      bim3a_rhs (tmsh, reaction, *phi_bar, rhs_src);
    else
      bim3a_rhs_frac (tmsh, *reaction_nodes, *phi_bar, rhs_src, func_frac);

    auto &rd = rhs->get_owned_data ();
    auto &sd = rhs_src.get_owned_data ();
    for (size_t i = 0; i < rd.size (); ++i)
      rd[i] += sd[i];

    phi_bar.reset ();
  }

  // Reaction-related vectors are no longer required
  reaction_nodes.reset();
  ones.reset();
  std::vector<double>().swap (reaction);
  std::vector<double>().swap (marker);

  // ------------------------------------------------------------
  // 3) MORTAR ASSEMBLY (PBC, before Dirichlet)
  // ------------------------------------------------------------
  // === [MORTAR BLOCK — DISABLED FOR STRONG PBC] Zone B: mortar assembly ===
  // if (periodic_x || periodic_y) {
  //   int    nm      = (1 << minlevel);
  //   int    n_pairs = (periodic_x ? 1 : 0) + (periodic_y ? 1 : 0);
  //   size_t N_own   = static_cast<size_t> (tmsh.num_owned_nodes ());
  //   size_t N_global = static_cast<size_t> (tmsh.num_global_nodes ());
  //   double Lx = r_cr[0] - l_cr[0];
  //   double Ly = r_cr[1] - l_cr[1];
  //   double Lz = r_cr[2] - l_cr[2];
  //
  //   std::vector<MortarFacePair> pairs;
  //   size_t cumulative = 0;
  //   if (periodic_x) {
  //     pairs.push_back ({ 0, 1,  1, 2,  l_cr[1], Ly,  l_cr[2], Lz,
  //                        N_own    + cumulative,
  //                        N_global + cumulative,
  //                        cumulative });
  //     cumulative += ndofm;
  //   }
  //   if (periodic_y) {
  //     pairs.push_back ({ 2, 3,  0, 2,  l_cr[0], Lx,  l_cr[2], Lz,
  //                        N_own    + cumulative,
  //                        N_global + cumulative,
  //                        cumulative });
  //     cumulative += ndofm;
  //   }
  //
  //   mortar_C.assign (static_cast<size_t> (n_pairs) * ndofm * N_global, 0.0);
  //   bool fill_mortar_rows = (linear_solver_name == "lis") && (size == 1);
  //
  //   for (auto& pair : pairs) {
  //     const char* dir = (pair.face_left < 2) ? "x" : "y";
  //     if (rank == 0)
  //       std::cout << "  Mortar assembly: face pair " << dir << "± ...\n";
  //     for (auto q = tmsh.begin_quadrant_sweep ();
  //          q != tmsh.end_quadrant_sweep (); ++q) {
  //       if (!q->ef (pair.face_left).empty ())
  //         assemble_mortar_block (*A, mortar_C, N_global, q,
  //                                pair.face_left,  pair, +1.0, nm, fill_mortar_rows);
  //       if (!q->ef (pair.face_right).empty ())
  //         assemble_mortar_block (*A, mortar_C, N_global, q,
  //                                pair.face_right, pair, -1.0, nm, fill_mortar_rows);
  //     }
  //     if (rank == 0)
  //       std::cout << "  Mortar assembly: face pair " << dir << "± done\n";
  //   }
  //
  //   if (size > 1) {
  //     MPI_Reduce (rank == 0 ? MPI_IN_PLACE : mortar_C.data (),
  //                 mortar_C.data (),
  //                 static_cast<int> (mortar_C.size ()), MPI_DOUBLE,
  //                 MPI_SUM, 0, mpicomm);
  //   }
  //
  //   if (linear_solver_name == "lis") {
  //     double local_max = 0.0;
  //     for (int ii = A->range_start (); ii < A->range_end (); ++ii)
  //       if ((*A)[ii].count (ii))
  //         local_max = std::max (local_max, std::abs ((*A)[ii].at (ii)));
  //     if (size > 1)
  //       MPI_Allreduce (&local_max, &mortar_eps_reg, 1, MPI_DOUBLE, MPI_MAX, mpicomm);
  //     else
  //       mortar_eps_reg = local_max;
  //     if (mortar_eps_reg <= 0.0) mortar_eps_reg = 1.0;
  //     if (size == 1) {
  //       for (auto& pair : pairs)
  //         for (size_t kk = 0; kk < static_cast<size_t> (ndofm); ++kk)
  //           (*A)[pair.mortar_row_offset + kk][pair.mortar_col_offset + kk] += mortar_eps_reg;
  //     }
  //     if (rank == 0)
  //       std::cout << "  Mortar diagonal regularization: eps = " << mortar_eps_reg << '\n';
  //   }
  //
  //   int n_bd = 0;
  //   for (auto& pair : pairs) {
  //     for (int j0 = 0; j0 <= nm; ++j0) {
  //       for (int j1_bd : {0, nm}) {
  //         size_t k    = static_cast<size_t> (j1_bd) * (nm + 1) + j0;
  //         size_t row  = pair.mortar_row_offset  + k;
  //         size_t gcol = pair.mortar_col_offset  + k;
  //         size_t drow = pair.mortar_dense_offset + k;
  //         if (fill_mortar_rows) {
  //           (*A)[row].clear ();
  //           (*A)[row][gcol] = 1.0;
  //         }
  //         if (rank == 0)
  //           for (size_t c = 0; c < N_global; ++c)
  //             mortar_C[drow * N_global + c] = 0.0;
  //         ++n_bd;
  //       }
  //     }
  //   }
  //   if (rank == 0)
  //     std::cout << "  Mortar z-boundary constraints: " << n_bd << " DOFs fixed to zero\n";
  // }
  // === [END MORTAR BLOCK Zone B] ===

  // ------------------------------------------------------------
  // 4) APPLY BOUNDARY CONDITIONS
  // ------------------------------------------------------------
  // Build Dirichlet boundary list
  dirichlet_bcs3 bcs;

  if (bc == 1) { // (Possibly inhomogeneous) Dirichlet BC
    if (std::fabs (pot_bc) > 1.e-5 && rank == 0)
      std::cerr << "[WARNING] Boundary conditions may be inaccurate!!\n";

    // Membrane voltage: the lower electrode (z-, face 4) stays grounded while the
    // upper electrode (z+, face 5) is held at applied_potential. Without membrane
    // voltage every Dirichlet face is homogeneous (legacy behaviour).
    const bool membrane_voltage =
      membrane_enabled && std::fabs (applied_potential) > 1.e-12;
    if (membrane_voltage && rank == 0)
      std::cout << "  Membrane voltage: z- = 0, z+ = " << applied_potential
                << " kT/e\n";

    for (auto const& ibc : bcells) {
      if (periodic_x && (ibc.second == 0 || ibc.second == 1)) continue;
      if (periodic_y && (ibc.second == 2 || ibc.second == 3)) continue;
      const double bc_val =
        (membrane_voltage && ibc.second == 5) ? applied_potential : 0.0;
      bcs.emplace_back (ibc.first, ibc.second,
      [bc_val] (double, double, double) {
        return bc_val;
      });
    }

    bim3a_dirichlet_bc (tmsh, bcs, *A, *rhs);
  }

  if (bc == 2) { // Coulombic Dirichlet BC
    for (auto const& ibc : bcells) {
      if (periodic_x && (ibc.second == 0 || ibc.second == 1)) continue;
      if (periodic_y && (ibc.second == 2 || ibc.second == 3)) continue;
      bcs.emplace_back (ibc.first, ibc.second,
      [&] (double x, double y, double z) {
        return coulomb_boundary_conditions (x, y, z);
      });
    }

    bim3a_dirichlet_bc (tmsh, bcs, *A, *rhs);
  }

  if (bc == 3) { // Analytic Dirichlet BC (sphere test case)
    for (auto const& ibc : bcells) {
      if (periodic_x && (ibc.second == 0 || ibc.second == 1)) continue;
      if (periodic_y && (ibc.second == 2 || ibc.second == 3)) continue;
      bcs.emplace_back (ibc.first, ibc.second,
      [&] (double x, double y, double z) {
        return analytic_solution (x, y, z);
      });
    }

    bim3a_dirichlet_bc (tmsh, bcs, *A, *rhs);
  }

  if (rank == 0 && (periodic_x || periodic_y)) {
    static const char* fname[6] = {"x-","x+","y-","y+","z-","z+"};
    int cnt[6] = {};
    for (auto const& ibc : bcells) cnt[ibc.second]++;
    std::cout << "  Dirichlet BC nodes per face (0=skipped for PBC):\n";
    for (int f = 0; f < 6; ++f)
      std::cout << "    face " << fname[f] << " : "
                << (((f<2)&&periodic_x)||((f>=2&&f<4)&&periodic_y) ? 0 : cnt[f])
                << '\n';
  }

  // ------------------------------------------------------------
  // 5) EXTEND RHS FOR MORTAR DOFs
  // ------------------------------------------------------------
  // === [MORTAR BLOCK — DISABLED FOR STRONG PBC] Zone C: RHS extension ===
  // if (periodic_x || periodic_y) {
  //   int n_pairs = (periodic_x ? 1 : 0) + (periodic_y ? 1 : 0);
  //   if (linear_solver_name != "lis" || size == 1 || rank == size - 1)
  //     rhs->get_owned_data ().resize (
  //       rhs->get_owned_data ().size () + n_pairs * ndofm, 0.0);
  // }
  // === [END MORTAR BLOCK Zone C] ===

  // ------------------------------------------------------------
  // 6) PARALLEL ASSEMBLY (MPI)
  // ------------------------------------------------------------

  if (size > 1) {
    A->assemble();
    rhs->assemble();
  }
}

void
poisson_boltzmann::mumps_compute_electric_potential (ray_cache_t &ray_cache)
{
  int rank, size;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);


  mumps mumps_solver;

  std::vector<double> vals;
  std::vector<int> irow, jcol;

  int n_pairs  = (periodic_x ? 1 : 0) + (periodic_y ? 1 : 0);
  int n_system = static_cast<int> (tmsh.num_global_nodes ());

  // -----------------------------------------------------------------------
  // Build distributed COO for MUMPS.
  //
  // Physical rows: each rank exports its owned rows [is, ie) via aij(flag=false).
  // Mortar column indices in physical rows are already global (pair.mortar_col_offset + k
  // = N_global + cumulative + k), so no correction is needed.
  //
  // Mortar rows: rank 0 converts mortar_C (dense C^T, already MPI_Reduce'd) to COO,
  // and appends identity entries for z-boundary DOFs.  Other ranks contribute nothing
  // to the mortar block.
  // -----------------------------------------------------------------------
  (*A).aij (vals, irow, jcol, mumps_solver.get_index_base ());   // flag=false, physical rows

  // === [MORTAR BLOCK — DISABLED] Zone D: mortar COO rows ===
  // if (periodic_x || periodic_y) {
  //   int    base     = mumps_solver.get_index_base ();
  //   int    nm       = (1 << minlevel);
  //   size_t N_global = static_cast<size_t> (tmsh.num_global_nodes ());
  //
  //   if (rank == 0) {
  //     size_t n_mortar_rows = static_cast<size_t> (n_pairs * ndofm);
  //     for (size_t drow = 0; drow < n_mortar_rows; ++drow) {
  //       int global_row = static_cast<int> (N_global + drow) + base;
  //
  //       // C^T entries from the dense accumulator (z-boundary rows were zeroed in assembly)
  //       for (size_t c = 0; c < N_global; ++c) {
  //         double v = mortar_C[drow * N_global + c];
  //         if (v != 0.0) {
  //           irow.push_back (global_row);
  //           jcol.push_back (static_cast<int> (c) + base);
  //           vals.push_back (v);
  //         }
  //       }
  //
  //       // Identity diagonal for z-boundary mortar DOFs (j1 = 0 or nm)
  //       int j1 = static_cast<int> (drow % static_cast<size_t> (ndofm)) / (nm + 1);
  //       if (j1 == 0 || j1 == nm) {
  //         irow.push_back (global_row);
  //         jcol.push_back (global_row);
  //         vals.push_back (1.0);
  //       }
  //     }
  //     std::cout << "  [Mortar] COO entries appended by rank 0: "
  //               << irow.size () << " total\n";
  //   }
  // }
  // === [END MORTAR BLOCK Zone D] ===

  // -----------------------------------------------------------------------
  // Step 3 — Gather the full RHS on rank 0 (PBC only).
  // Moved before the elimination block so pbc_rhs_buf is available there.
  // rhs_counts / rhs_displs are reused in step 7 for the scatter.
  // -----------------------------------------------------------------------
  std::vector<double> pbc_rhs_buf;
  std::vector<int>    rhs_counts (size, 0), rhs_displs (size, 0);

  if (periodic_x || periodic_y) {
    int local_n = static_cast<int> (tmsh.num_owned_nodes ());
    MPI_Gather (&local_n, 1, MPI_INT,
                rhs_counts.data (), 1, MPI_INT, 0, mpicomm);
    if (rank == 0) {
      for (int r = 1; r < size; ++r)
        rhs_displs[r] = rhs_displs[r - 1] + rhs_counts[r - 1];
      pbc_rhs_buf.assign (n_system, 0.0);
    }
    MPI_Gatherv (rhs->get_owned_data ().data (), local_n, MPI_DOUBLE,
                 pbc_rhs_buf.data (), rhs_counts.data (), rhs_displs.data (),
                 MPI_DOUBLE, 0, mpicomm);
  }

  // -----------------------------------------------------------------------
  // PBC strong enforcement: eliminate right-face nodes (steps 1–4).
  //
  // Steps 1+2: remap right→left in irow/jcol (all ranks) and merge RHS
  //   entries (rank 0 only, pbc_rhs_buf is rank-0-only).
  //   y-pairs first, x-pairs second: chain xR_yR→(y)→xR_yL→(x)→xL_yL.
  // Steps 3+4: build old_to_new (all ranks, deterministic), renumber
  //   irow/jcol, extract rhs_new[N_new] (rank 0 only).
  // -----------------------------------------------------------------------
  const size_t        N_phys_saved = static_cast<size_t> (n_system);
  std::vector<int>    old_to_new;
  std::vector<double> rhs_new;

  if (periodic_x || periodic_y) {
    const int base = mumps_solver.get_index_base ();

    // Step 1 — y-pairs: remap in irow/jcol (all ranks)
    for (size_t k = 0; k < irow.size (); ++k) {
      auto ity = pbc_y_right_to_left.find (static_cast<size_t> (irow[k] - base));
      if (ity != pbc_y_right_to_left.end ())
        irow[k] = static_cast<int> (ity->second) + base;
      ity = pbc_y_right_to_left.find (static_cast<size_t> (jcol[k] - base));
      if (ity != pbc_y_right_to_left.end ())
        jcol[k] = static_cast<int> (ity->second) + base;
    }
    if (rank == 0)
      for (auto it = pbc_y_right_to_left.begin ();
           it != pbc_y_right_to_left.end (); ++it) {
        pbc_rhs_buf[it->second] += pbc_rhs_buf[it->first];
        pbc_rhs_buf[it->first]   = 0.0;
      }

    // Step 2 — x-pairs: remap in irow/jcol (all ranks)
    // xR_yR was already replaced by xR_yL in step 1; the x-pass now
    // maps xR_yL → xL_yL, completing the corner chain.
    for (size_t k = 0; k < irow.size (); ++k) {
      auto itx = pbc_x_right_to_left.find (static_cast<size_t> (irow[k] - base));
      if (itx != pbc_x_right_to_left.end ())
        irow[k] = static_cast<int> (itx->second) + base;
      itx = pbc_x_right_to_left.find (static_cast<size_t> (jcol[k] - base));
      if (itx != pbc_x_right_to_left.end ())
        jcol[k] = static_cast<int> (itx->second) + base;
    }
    if (rank == 0)
      for (auto it = pbc_x_right_to_left.begin ();
           it != pbc_x_right_to_left.end (); ++it) {
        pbc_rhs_buf[it->second] += pbc_rhs_buf[it->first];
        pbc_rhs_buf[it->first]   = 0.0;
      }

    // Step 3 — build old_to_new (all ranks, same deterministic computation)
    std::vector<bool> is_right (N_phys_saved, false);
    for (auto it = pbc_y_right_to_left.begin ();
         it != pbc_y_right_to_left.end (); ++it)
      is_right[it->first] = true;
    for (auto it = pbc_x_right_to_left.begin ();
         it != pbc_x_right_to_left.end (); ++it)
      is_right[it->first] = true;

    old_to_new.assign (N_phys_saved, -1);
    int new_idx = 0;
    for (size_t i = 0; i < N_phys_saved; ++i)
      if (!is_right[i])
        old_to_new[i] = new_idx++;
    n_system = new_idx;

    if (rank == 0)
      std::cout << "  [PBC strong] N_new = " << n_system
                << "  (eliminated " << (N_phys_saved - static_cast<size_t> (n_system))
                << " right-face nodes)\n";

    // Step 4 — renumber irow/jcol (all ranks), extract rhs_new (rank 0)
    for (size_t k = 0; k < irow.size (); ++k) {
      irow[k] = old_to_new[irow[k] - base] + base;
      jcol[k] = old_to_new[jcol[k] - base] + base;
    }
    if (rank == 0) {
      rhs_new.resize (n_system);
      for (size_t i = 0; i < N_phys_saved; ++i)
        if (old_to_new[i] >= 0)
          rhs_new[old_to_new[i]] = pbc_rhs_buf[i];
    }
  }

  mumps_solver.set_lhs_distributed ();
  mumps_solver.set_distributed_lhs_structure (n_system, irow, jcol);
  mumps_solver.set_distributed_lhs_data (vals);

  if (periodic_x || periodic_y) {
    if (rank == 0)
      mumps_solver.set_rhs (rhs_new);
  } else {
    mumps_solver.set_rhs_distributed (*rhs);
  }

  rhs.reset ();
  A.reset ();

  if (rank == 0)
    std::cout << "  System size for MUMPS: " << n_system << '\n';

  std::cout << "mumps_solver.analyze () = "
            << mumps_solver.analyze ()
            << std::endl;
  std::cout << "mumps_solver.factorize () = "
            << mumps_solver.factorize ()
            << std::endl;
  std::cout << "mumps_solver.solve () = "
            << mumps_solver.solve ()
            << std::endl;

  phi = std::make_unique<distributed_vector> (tmsh.num_owned_nodes ());

  if (periodic_x || periodic_y) {
    // -----------------------------------------------------------------------
    // Step 6 — Expand sol_new[N_new] → phi_full[N_phys] on rank 0.
    //
    // x-pairs filled before y-pairs so corner chain resolves correctly:
    //   x-fill: phi[xR_yL] = phi[xL_yL]
    //   y-fill: phi[xR_yR] = phi[xR_yL] = phi[xL_yL]
    // -----------------------------------------------------------------------
    if (rank == 0) {
      pbc_rhs_buf.assign (N_phys_saved, 0.0);
      for (size_t i = 0; i < N_phys_saved; ++i)
        if (old_to_new[i] >= 0)
          pbc_rhs_buf[i] = rhs_new[old_to_new[i]];
      for (auto it = pbc_x_right_to_left.begin ();
           it != pbc_x_right_to_left.end (); ++it)
        pbc_rhs_buf[it->first] = pbc_rhs_buf[it->second];
      for (auto it = pbc_y_right_to_left.begin ();
           it != pbc_y_right_to_left.end (); ++it)
        pbc_rhs_buf[it->first] = pbc_rhs_buf[it->second];
    }

    // Step 7 — Scatter phi_full[N_phys] to all ranks.
    int local_n = static_cast<int> (tmsh.num_owned_nodes ());
    MPI_Scatterv (pbc_rhs_buf.data (), rhs_counts.data (), rhs_displs.data (),
                  MPI_DOUBLE,
                  phi->get_owned_data ().data (), local_n,
                  MPI_DOUBLE, 0, mpicomm);
  } else {
    (*phi) = mumps_solver.get_distributed_solution ();
  }

  if (size > 1)
    bim3a_solution_with_ghosts (tmsh, *phi, replace_op);

  // Collective solution debug: each rank contributes its owned face nodes,
  // then MPI_Gatherv assembles the complete picture on rank 0.
  {
    static const char* fname[6] = {"x-","x+","y-","y+","z-","z+"};
    const int phi_is = static_cast<int> (phi->get_range_start ());
    const int phi_ie = phi_is + static_cast<int> (phi->get_owned_data ().size ());

    if (rank == 0)
      std::cout << "\n=== [ Solution debug ] ===\n";

    for (int face : {0, 1, 2, 3, 4, 5}) {
      bool is_pbc = ((face < 2) && periodic_x) || ((face >= 2 && face < 4) && periodic_y);

      // Collect owned face nodes (map deduplicates automatically)
      std::map<int, double> face_map;
      for (auto q = tmsh.begin_quadrant_sweep ();
           q != tmsh.end_quadrant_sweep (); ++q)
        for (auto local_idx : q->ef (face)) {
          int g = q->gt (local_idx);
          if (g >= phi_is && g < phi_ie)
            face_map[g] = phi->get_owned_data ()[g - phi_is];
        }

      std::vector<int>    local_gidx;
      std::vector<double> local_vals;
      local_gidx.reserve (face_map.size ());
      local_vals.reserve (face_map.size ());
      for (const auto& kv : face_map) {
        local_gidx.push_back (kv.first);
        local_vals.push_back (kv.second);
      }

      int local_n = static_cast<int> (local_gidx.size ());
      std::vector<int> all_ns (size, 0), displs (size, 0);
      MPI_Gather (&local_n, 1, MPI_INT, all_ns.data (), 1, MPI_INT, 0, mpicomm);

      int total_n = 0;
      std::vector<int>    all_gidx;
      std::vector<double> all_vals;
      if (rank == 0) {
        for (int r = 1; r < size; ++r)
          displs[r] = displs[r - 1] + all_ns[r - 1];
        total_n = displs[size - 1] + all_ns[size - 1];
        all_gidx.resize (total_n);
        all_vals.resize (total_n);
      }

      MPI_Gatherv (local_gidx.data (), local_n, MPI_INT,
                   all_gidx.data (), all_ns.data (), displs.data (), MPI_INT, 0, mpicomm);
      MPI_Gatherv (local_vals.data (), local_n, MPI_DOUBLE,
                   all_vals.data (), all_ns.data (), displs.data (), MPI_DOUBLE, 0, mpicomm);

      if (rank == 0 && total_n > 0) {
        double fmin = all_vals[0], fmax = all_vals[0];
        for (int i = 1; i < total_n; ++i) {
          fmin = std::min (fmin, all_vals[i]);
          fmax = std::max (fmax, all_vals[i]);
        }
        std::cout << "  face " << fname[face]
                  << (is_pbc ? " [PBC]" : " [Dir]")
                  << "  nodes=" << total_n
                  << "  phi=[" << fmin << ", " << fmax << "]\n";
      }
    }

    if (rank == 0)
      std::cout << "=========================\n";
  }

  mumps_solver.cleanup ();
}



void
poisson_boltzmann::lis_compute_electric_potential (ray_cache_t &ray_cache)
{
  int rank, size;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);

  bool pbc = (periodic_x || periodic_y);

  std::vector<double> vals;
  std::vector<int>    irow, jcol;

  LIS_INT is, ie;

  if (pbc && size > 1) {
    // ── Multi-rank PBC ─────────────────────────────────────────────────────
    // Physical block: each rank exports only its owned rows (flag=false).
    // Mortar rows: rank size-1 appends them manually from mortar_C.
    // mortar_C was MPI_Reduce'd to rank 0 in assemple_system_matrix; broadcast
    // here so rank size-1 can build the mortar CSR block.
    // ───────────────────────────────────────────────────────────────────────

    int    n_pairs    = (periodic_x ? 1 : 0) + (periodic_y ? 1 : 0);
    int    nm         = (1 << minlevel);
    int    N_global   = static_cast<int> (tmsh.num_global_nodes ());
    int    mortar_C_n = n_pairs * ndofm;

    (*A).csr (vals, jcol, irow, 0, false);
    A.reset ();

    MPI_Bcast (mortar_C.data (), mortar_C_n * N_global, MPI_DOUBLE, 0, mpicomm);

    if (rank == 0)
      std::cout << "  [PBC multi-rank LIS] eps_reg = " << mortar_eps_reg << '\n';

    if (rank == size - 1) {
      for (int drow = 0; drow < mortar_C_n; ++drow) {
        int  j1      = (drow % ndofm) / (nm + 1);
        bool is_zbdy = (j1 == 0 || j1 == nm);
        int  gcol    = N_global + drow;

        if (is_zbdy) {
          jcol.push_back (gcol);
          vals.push_back (1.0);
        } else {
          for (int c = 0; c < N_global; ++c) {
            double v = mortar_C[static_cast<size_t> (drow) * N_global + c];
            if (v != 0.0) {
              jcol.push_back (c);
              vals.push_back (v);
            }
          }
          jcol.push_back (gcol);
          vals.push_back (mortar_eps_reg);
        }
        irow.push_back (static_cast<int> (jcol.size ()));
      }
    }

    LIS_INT n_own = static_cast<LIS_INT> (tmsh.num_owned_nodes ());
    LIS_INT ln    = n_own + (rank == size - 1
                              ? static_cast<LIS_INT> (mortar_C_n) : 0);
    LIS_INT nnz   = static_cast<LIS_INT> (vals.size ());

    LIS_VECTOR rhs_lis;
    lis_vector_create (mpicomm, &rhs_lis);
    lis_vector_set_size (rhs_lis, ln, 0);
    lis_vector_get_range (rhs_lis, &is, &ie);
    for (LIS_INT i = is; i < ie; i++)
      lis_vector_set_value (LIS_INS_VALUE, i,
                            rhs->get_owned_data ()[i - is], rhs_lis);
    rhs.reset ();

    LIS_VECTOR phi_lis;
    lis_vector_create (mpicomm, &phi_lis);
    lis_vector_set_size (phi_lis, ln, 0);

    LIS_MATRIX A_lis;
    lis_matrix_create (mpicomm, &A_lis);
    lis_matrix_set_size (A_lis, ln, 0);
    lis_matrix_set_csr (nnz, irow.data (), jcol.data (), vals.data (), A_lis);
    lis_matrix_assemble (A_lis);

    LIS_SOLVER solver;
    lis_solver_create (&solver);
    std::string opts = linear_solver_options;
    lis_solver_set_option (&opts[0], solver);
    lis_solve (A_lis, rhs_lis, phi_lis, solver);
    lis_solver_destroy (solver);
    lis_vector_destroy (rhs_lis);

    phi = std::make_unique<distributed_vector> (tmsh.num_owned_nodes (), mpicomm);
    lis_vector_get_values (phi_lis, is, n_own, phi->get_owned_data ().data ());
    lis_vector_destroy (phi_lis);

  } else {
    // ── Non-PBC or single-rank ─────────────────────────────────────────────
    // flag=pbc: single-rank PBC exports mortar rows already in A (fill_mortar_rows=true).
    // flag=false for non-PBC or multi-rank without PBC.
    // ───────────────────────────────────────────────────────────────────────
    (*A).csr (vals, jcol, irow, 0, false);
    A.reset ();

    if (pbc) {
      // Single-rank PBC: strong enforcement via CSR elimination.
      //
      // Steps 1+2: remap jcol (y-pass then x-pass, 0-based indices).
      // row_dest[]: canonical destination row for each source row.
      // Steps 3+4: old_to_new, N_new, two-pass CSR build + column sort.
      // RHS merge identical to MUMPS path.
      const size_t N_phys = static_cast<size_t> (tmsh.num_global_nodes ());

      // Step 1: y-pairs column remap
      for (size_t k = 0; k < jcol.size (); ++k) {
        auto ity = pbc_y_right_to_left.find (static_cast<size_t> (jcol[k]));
        if (ity != pbc_y_right_to_left.end ())
          jcol[k] = static_cast<int> (ity->second);
      }
      // Step 2: x-pairs column remap
      for (size_t k = 0; k < jcol.size (); ++k) {
        auto itx = pbc_x_right_to_left.find (static_cast<size_t> (jcol[k]));
        if (itx != pbc_x_right_to_left.end ())
          jcol[k] = static_cast<int> (itx->second);
      }

      // row_dest[r]: canonical left-face destination for row r.
      // y-pass sets yR->yL; x-pass follows chain so xR_yR->(y)xL_yR->(x)xL_yL.
      std::vector<size_t> row_dest (N_phys);
      for (size_t i = 0; i < N_phys; ++i) row_dest[i] = i;
      for (auto it = pbc_y_right_to_left.begin ();
           it != pbc_y_right_to_left.end (); ++it)
        row_dest[it->first] = it->second;
      for (auto it = pbc_x_right_to_left.begin ();
           it != pbc_x_right_to_left.end (); ++it)
        row_dest[it->first] = row_dest[it->second];

      // Step 3: is_right[], old_to_new[], N_new
      std::vector<bool> is_right (N_phys, false);
      for (auto it = pbc_y_right_to_left.begin ();
           it != pbc_y_right_to_left.end (); ++it)
        is_right[it->first] = true;
      for (auto it = pbc_x_right_to_left.begin ();
           it != pbc_x_right_to_left.end (); ++it)
        is_right[it->first] = true;

      std::vector<int> old_to_new (N_phys, -1);
      int new_idx = 0;
      for (size_t i = 0; i < N_phys; ++i)
        if (!is_right[i]) old_to_new[i] = new_idx++;
      int N_new = new_idx;

      if (rank == 0)
        std::cout << "  [PBC strong LIS] N_new = " << N_new
                  << "  (eliminated "
                  << (N_phys - static_cast<size_t> (N_new))
                  << " right-face nodes)\n";

      // Step 4: build new N_new x N_new CSR (two-pass: count then fill)
      std::vector<int> new_irow (N_new + 1, 0);
      for (size_t r = 0; r < N_phys; ++r) {
        int nd = old_to_new[row_dest[r]];
        for (int k = irow[r]; k < irow[r + 1]; ++k)
          if (old_to_new[jcol[k]] >= 0)
            new_irow[nd + 1]++;
      }
      for (int r = 0; r < N_new; ++r) new_irow[r + 1] += new_irow[r];
      int new_nnz = new_irow[N_new];

      std::vector<int>    new_jcol (new_nnz);
      std::vector<double> new_vals (new_nnz);
      {
        std::vector<int> pos (new_irow.begin (), new_irow.begin () + N_new);
        for (size_t r = 0; r < N_phys; ++r) {
          int nd = old_to_new[row_dest[r]];
          for (int k = irow[r]; k < irow[r + 1]; ++k) {
            int nc = old_to_new[jcol[k]];
            if (nc >= 0) {
              new_jcol[pos[nd]] = nc;
              new_vals[pos[nd]] = vals[k];
              pos[nd]++;
            }
          }
        }
      }
      // Sort column indices per row, then deduplicate (row merging produces
      // duplicate (row,col) entries that must be summed, not overwritten).
      {
        std::vector<int>    perm, tmp_j;
        std::vector<double> tmp_v;
        for (int r = 0; r < N_new; ++r) {
          int rb = new_irow[r], re = new_irow[r + 1], n = re - rb;
          if (n <= 1) continue;
          perm.resize (n); tmp_j.resize (n); tmp_v.resize (n);
          for (int i = 0; i < n; ++i) perm[i] = i;
          std::sort (perm.begin (), perm.end (),
                     [&] (int a, int b)
                     { return new_jcol[rb + a] < new_jcol[rb + b]; });
          for (int i = 0; i < n; ++i) {
            tmp_j[i] = new_jcol[rb + perm[i]];
            tmp_v[i] = new_vals[rb + perm[i]];
          }
          std::copy (tmp_j.begin (), tmp_j.end (), new_jcol.begin () + rb);
          std::copy (tmp_v.begin (), tmp_v.end (), new_vals.begin () + rb);
        }
      }
      // Deduplicate: adjacent equal-column entries per row are summed.
      std::vector<int>    ded_irow (N_new + 1, 0);
      std::vector<int>    ded_jcol;
      std::vector<double> ded_vals;
      ded_jcol.reserve (new_nnz);
      ded_vals.reserve (new_nnz);
      for (int r = 0; r < N_new; ++r) {
        ded_irow[r] = static_cast<int> (ded_jcol.size ());
        const int rb = new_irow[r], re = new_irow[r + 1];
        for (int k = rb; k < re; ++k) {
          if (static_cast<int> (ded_jcol.size ()) > ded_irow[r] &&
              ded_jcol.back () == new_jcol[k])
            ded_vals.back () += new_vals[k];
          else {
            ded_jcol.push_back (new_jcol[k]);
            ded_vals.push_back (new_vals[k]);
          }
        }
      }
      ded_irow[N_new] = static_cast<int> (ded_jcol.size ());
      int ded_nnz = static_cast<int> (ded_jcol.size ());

      // Merge RHS (y-pass then x-pass; single rank owns all N_phys entries)
      std::vector<double> rhs_vec (rhs->get_owned_data ().begin (),
                                   rhs->get_owned_data ().end ());
      rhs.reset ();
      for (auto it = pbc_y_right_to_left.begin ();
           it != pbc_y_right_to_left.end (); ++it) {
        rhs_vec[it->second] += rhs_vec[it->first];
        rhs_vec[it->first]   = 0.0;
      }
      for (auto it = pbc_x_right_to_left.begin ();
           it != pbc_x_right_to_left.end (); ++it) {
        rhs_vec[it->second] += rhs_vec[it->first];
        rhs_vec[it->first]   = 0.0;
      }
      std::vector<double> rhs_new (N_new);
      for (size_t i = 0; i < N_phys; ++i)
        if (old_to_new[i] >= 0) rhs_new[old_to_new[i]] = rhs_vec[i];

      // LIS solve on N_new x N_new
      LIS_INT lis_n   = static_cast<LIS_INT> (N_new);
      LIS_INT lis_nnz = static_cast<LIS_INT> (ded_nnz);

      LIS_VECTOR rhs_lis;
      lis_vector_create (mpicomm, &rhs_lis);
      lis_vector_set_size (rhs_lis, lis_n, 0);
      lis_vector_get_range (rhs_lis, &is, &ie);
      for (LIS_INT i = is; i < ie; ++i)
        lis_vector_set_value (LIS_INS_VALUE, i, rhs_new[i], rhs_lis);

      LIS_VECTOR phi_lis;
      lis_vector_create (mpicomm, &phi_lis);
      lis_vector_set_size (phi_lis, lis_n, 0);
      lis_vector_get_range (phi_lis, &is, &ie);

      LIS_MATRIX A_lis;
      lis_matrix_create (mpicomm, &A_lis);
      lis_matrix_set_size (A_lis, lis_n, 0);
      lis_matrix_set_csr (lis_nnz, ded_irow.data (), ded_jcol.data (),
                          ded_vals.data (), A_lis);
      lis_matrix_assemble (A_lis);

      LIS_SOLVER solver;
      lis_solver_create (&solver);
      std::string opts = linear_solver_options;
      lis_solver_set_option (&opts[0], solver);
      lis_solve (A_lis, rhs_lis, phi_lis, solver);
      lis_solver_destroy (solver);
      lis_vector_destroy (rhs_lis);

      // Expand sol_new[N_new] -> phi_full[N_phys]; x-fill before y-fill
      std::vector<double> sol_new (N_new);
      lis_vector_get_values (phi_lis, is, lis_n, sol_new.data ());
      lis_vector_destroy (phi_lis);

      phi = std::make_unique<distributed_vector> (tmsh.num_owned_nodes (), mpicomm);
      auto& phi_data = phi->get_owned_data ();
      phi_data.assign (N_phys, 0.0);
      for (size_t i = 0; i < N_phys; ++i)
        if (old_to_new[i] >= 0) phi_data[i] = sol_new[old_to_new[i]];
      for (auto it = pbc_x_right_to_left.begin ();
           it != pbc_x_right_to_left.end (); ++it)
        phi_data[it->first] = phi_data[it->second];
      for (auto it = pbc_y_right_to_left.begin ();
           it != pbc_y_right_to_left.end (); ++it)
        phi_data[it->first] = phi_data[it->second];

    } else {
      // Non-PBC
      LIS_INT i, ln;
      LIS_VECTOR rhs_lis;
      ln = static_cast<LIS_INT> (rhs->get_owned_data ().size ());

      lis_vector_create (mpicomm, &rhs_lis);
      lis_vector_set_size (rhs_lis, ln, 0);
      lis_vector_get_range (rhs_lis, &is, &ie);
      for (i = is; i < ie; i++)
        lis_vector_set_value (LIS_INS_VALUE, i,
                              rhs->get_owned_data ()[i - is], rhs_lis);
      rhs.reset ();

      LIS_VECTOR phi_lis;
      lis_vector_create (mpicomm, &phi_lis);
      lis_vector_set_size (phi_lis, ln, 0);
      lis_vector_get_range (phi_lis, &is, &ie);

      LIS_INT n   = ln;
      LIS_INT nnz = static_cast<LIS_INT> (vals.size ());

      LIS_MATRIX A_lis;
      lis_matrix_create (mpicomm, &A_lis);
      lis_matrix_set_size (A_lis, n, 0);
      lis_matrix_set_csr (nnz, &irow[0], &jcol[0], &vals[0], A_lis);
      lis_matrix_assemble (A_lis);

      LIS_SOLVER solver;
      lis_solver_create (&solver);
      std::string opts = linear_solver_options;
      lis_solver_set_option (&opts[0], solver);
      lis_solve (A_lis, rhs_lis, phi_lis, solver);
      lis_solver_destroy (solver);
      lis_vector_destroy (rhs_lis);

      phi = std::make_unique<distributed_vector> (tmsh.num_owned_nodes (), mpicomm);
      LIS_INT n_phys = static_cast<LIS_INT> (tmsh.num_owned_nodes ());
      lis_vector_get_values (phi_lis, is, n_phys, phi->get_owned_data ().data ());
      lis_vector_destroy (phi_lis);
    }
  }

  if (size > 1)
    bim3a_solution_with_ghosts (tmsh, *phi, replace_op);

  // Collective debug — same pattern as MUMPS: each rank contributes its owned
  // face nodes via MPI_Gatherv, rank 0 prints the summary.
  {
    static const char* fname[6] = {"x-","x+","y-","y+","z-","z+"};
    const int phi_is = static_cast<int> (phi->get_range_start ());
    const int phi_ie = phi_is + static_cast<int> (phi->get_owned_data ().size ());

    if (rank == 0)
      std::cout << "\n=== [ Solution debug (LIS) ] ===\n";

    for (int face : {0, 1, 2, 3, 4, 5}) {
      bool is_pbc = ((face < 2) && periodic_x) || ((face >= 2 && face < 4) && periodic_y);

      std::map<int, double> face_map;
      for (auto q = tmsh.begin_quadrant_sweep ();
           q != tmsh.end_quadrant_sweep (); ++q)
        for (auto local_idx : q->ef (face)) {
          int g = q->gt (local_idx);
          if (g >= phi_is && g < phi_ie)
            face_map[g] = phi->get_owned_data ()[g - phi_is];
        }

      std::vector<int>    local_gidx;
      std::vector<double> local_vals;
      local_gidx.reserve (face_map.size ());
      local_vals.reserve (face_map.size ());
      for (const auto& kv : face_map) {
        local_gidx.push_back (kv.first);
        local_vals.push_back (kv.second);
      }

      int local_n = static_cast<int> (local_gidx.size ());
      std::vector<int> all_ns (size, 0), displs (size, 0);
      MPI_Gather (&local_n, 1, MPI_INT, all_ns.data (), 1, MPI_INT, 0, mpicomm);

      int total_n = 0;
      std::vector<int>    all_gidx;
      std::vector<double> all_vals;
      if (rank == 0) {
        for (int r = 1; r < size; ++r)
          displs[r] = displs[r - 1] + all_ns[r - 1];
        total_n = displs[size - 1] + all_ns[size - 1];
        all_gidx.resize (total_n);
        all_vals.resize (total_n);
      }

      MPI_Gatherv (local_gidx.data (), local_n, MPI_INT,
                   all_gidx.data (), all_ns.data (), displs.data (), MPI_INT,
                   0, mpicomm);
      MPI_Gatherv (local_vals.data (), local_n, MPI_DOUBLE,
                   all_vals.data (), all_ns.data (), displs.data (), MPI_DOUBLE,
                   0, mpicomm);

      if (rank == 0 && total_n > 0) {
        double fmin = all_vals[0], fmax = all_vals[0];
        for (int ii = 1; ii < total_n; ++ii) {
          fmin = std::min (fmin, all_vals[ii]);
          fmax = std::max (fmax, all_vals[ii]);
        }
        std::cout << "  face " << fname[face]
                  << (is_pbc ? " [PBC]" : " [Dir]")
                  << "  nodes=" << total_n
                  << "  phi=[" << fmin << ", " << fmax << "]\n";
      }
    }

    if (rank == 0)
      std::cout << "================================\n";
  }
}


double
poisson_boltzmann::coulomb_boundary_conditions (double x, double y, double z)
{
  double eps_out = 4.0 * pi * e_0 * e_out * kb * T * Angs / (e *e); //adim e_out
  double C_0 = 1.0e3 * N_av * ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0 * C_0 * Angs * Angs * e * e / (e_0 *e_out *kb *T);

  double dist = 0.0;
  double pot = 0.0;
  double k = std::sqrt (k2);

  for (const NS::Atom& i : atoms) {
    dist = std::hypot ( (i.pos[0] - x), (i.pos[1] - y), (i.pos[2] - z));
    pot += i.charge * exp (-k *dist) / (dist *eps_out);
  }

  return pot;
}

