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

  // Allocate sparse matrix A and RHS vector
  A = std::make_unique<distributed_sparse_matrix> (mpicomm);
  A->set_ranges (tmsh.num_owned_nodes());

  if (periodic_x || periodic_y) {
    int nm     = (1 << minlevel);
    ndofm      = (nm + 1) * (nm + 1);
    int n_pairs = (periodic_x ? 1 : 0) + (periodic_y ? 1 : 0);
    size_t N_phys = A->rows ();
    A->resize (N_phys + n_pairs * ndofm);
    A->m = N_phys + n_pairs * ndofm;
    if (rank == 0) {
      std::cout << "  Physical DOFs (N_phys)  : " << N_phys << '\n';
      std::cout << "  Active mortar pairs     : " << n_pairs << '\n';
      std::cout << "  Mortar DOFs per pair    : " << ndofm << "  (= (nm+1)^2)\n";
      std::cout << "  Augmented system size   : " << N_phys + n_pairs * ndofm << '\n';
    }
  } else {
    if (rank == 0)
      std::cout << "  Physical DOFs           : " << A->rows () << '\n';
  }

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

  // Reaction-related vectors are no longer required
  reaction_nodes.reset();
  ones.reset();
  std::vector<double>().swap (reaction);
  std::vector<double>().swap (marker);

  // ------------------------------------------------------------
  // 3) MORTAR ASSEMBLY (PBC, before Dirichlet)
  // ------------------------------------------------------------

  if (periodic_x || periodic_y) {
    int    nm     = (1 << minlevel);
    size_t N_phys = static_cast<size_t> (tmsh.num_owned_nodes ());
    double Lx = r_cr[0] - l_cr[0];
    double Ly = r_cr[1] - l_cr[1];
    double Lz = r_cr[2] - l_cr[2];

    // Build only the active pairs; mortar offsets are assigned sequentially
    // so that inactive directions leave no zero rows in the augmented matrix.
    std::vector<MortarFacePair> pairs;
    size_t mortar_offset = N_phys;
    if (periodic_x) {
      pairs.push_back ({ 0, 1,  1, 2,  l_cr[1], Ly,  l_cr[2], Lz,  mortar_offset });
      mortar_offset += ndofm;
    }
    if (periodic_y) {
      pairs.push_back ({ 2, 3,  0, 2,  l_cr[0], Lx,  l_cr[2], Lz,  mortar_offset });
      mortar_offset += ndofm;
    }

    for (auto& pair : pairs) {
      const char* dir = (pair.face_left < 2) ? "x" : "y";
      if (rank == 0)
        std::cout << "  Mortar assembly: face pair " << dir << "± ...\n";

      for (auto q = tmsh.begin_quadrant_sweep ();
           q != tmsh.end_quadrant_sweep (); ++q) {
        if (!q->ef (pair.face_left).empty ())
          assemble_mortar_block (*A, q, pair.face_left,  pair, +1.0, nm);
        if (!q->ef (pair.face_right).empty ())
          assemble_mortar_block (*A, q, pair.face_right, pair, -1.0, nm);
      }

      if (rank == 0)
        std::cout << "  Mortar assembly: face pair " << dir << "± done\n";
    }

    // Mortar DOFs at j1=0 (z_lo) and j1=nm (z_hi) touch only physical nodes
    // on the Dirichlet z-boundaries, which are skipped by the keep[] filter in
    // assemble_mortar_block.  Those mortar rows are structurally zero → singular.
    // Enforce F=0 there with an identity row (RHS is already zero).
    int n_bd = 0;
    for (auto& pair : pairs) {
      for (int j0 = 0; j0 <= nm; ++j0) {
        for (int j1_bd : {0, nm}) {
          size_t col = pair.mortar_offset + static_cast<size_t> (j1_bd) * (nm + 1) + j0;
          (*A)[col].clear ();
          (*A)[col][col] = 1.0;
          ++n_bd;
        }
      }
    }
    if (rank == 0)
      std::cout << "  Mortar z-boundary constraints: " << n_bd << " DOFs fixed to zero\n";
  }

  // ------------------------------------------------------------
  // 4) APPLY BOUNDARY CONDITIONS
  // ------------------------------------------------------------
  // Build Dirichlet boundary list
  dirichlet_bcs3 bcs;

  if (bc == 1) { // Homogeneous Dirichlet BC
    if (std::fabs (pot_bc) > 1.e-5 && rank == 0)
      std::cerr << "[WARNING] Boundary conditions may be inaccurate!!\n";

    for (auto const& ibc : bcells) {
      if (periodic_x && (ibc.second == 0 || ibc.second == 1)) continue;
      if (periodic_y && (ibc.second == 2 || ibc.second == 3)) continue;
      bcs.emplace_back (ibc.first, ibc.second,
      [] (double, double, double) {
        return 0.0;
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
  // The mortar rows of the K_static rhs are zero by construction (C^T u = 0).
  // We append n_pairs*ndofm zeros to the local owned data so that the buffer
  // passed to MUMPS covers the full augmented system.
  // Note: rhs->size() (from distributed_vector::ranges) still reports N_phys.

  if (periodic_x || periodic_y) {
    int n_pairs = (periodic_x ? 1 : 0) + (periodic_y ? 1 : 0);
    rhs->get_owned_data ().resize (
      rhs->get_owned_data ().size () + n_pairs * ndofm, 0.0);
  }

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
  int n_system = (periodic_x || periodic_y)
                 ? tmsh.num_global_nodes () + n_pairs * ndofm
                 : tmsh.num_global_nodes ();

  if (periodic_x || periodic_y) {
    // distributed_sparse_matrix::aij() only exports [is, ie) rows; mortar rows
    // (N_phys..N_phys+n_pairs*ndofm-1) are outside that window.
    // flag=true calls sparse_matrix::aij() which iterates ALL rows in the vector,
    // so mortar rows (C^T and identity) are included.
    // Use centralized MUMPS input (ICNTL(18)=0, rank-0 only) — no set_lhs_distributed().
    (*A).aij (vals, irow, jcol, mumps_solver.get_index_base (), true);
    mumps_solver.set_lhs_structure (n_system, irow, jcol);
    mumps_solver.set_lhs_data (vals);
  } else {
    (*A).aij (vals, irow, jcol, mumps_solver.get_index_base ());
    mumps_solver.set_lhs_distributed ();
    mumps_solver.set_distributed_lhs_structure (n_system, irow, jcol);
    mumps_solver.set_distributed_lhs_data (vals);
  }

  // For the PBC path, set_rhs() stores a raw pointer into rhs->get_owned_data().
  // If we reset rhs before solve(), that pointer dangles and MUMPS writes to
  // freed memory.  Copy the buffer first so it outlives rhs.reset().
  // After solve() MUMPS writes the solution back in-place into this buffer.
  std::vector<double> pbc_rhs_buf;

  if (periodic_x || periodic_y) {
    if (rank == 0) {
      pbc_rhs_buf = rhs->get_owned_data ();   // copy before rhs is destroyed
      mumps_solver.set_rhs (pbc_rhs_buf);     // id.rhs → pbc_rhs_buf (persists)
    }
  } else {
    mumps_solver.set_rhs_distributed (*rhs);
  }

  rhs.reset ();
  A.reset ();

  if (rank == 0) {
    std::cout << "  System size: " << n_system;
    if (periodic_x || periodic_y)
      std::cout << "  (N_phys=" << tmsh.num_global_nodes ()
                << " + " << n_pairs << "*ndofm=" << n_pairs * ndofm << ")";
    std::cout << '\n';
  }

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
    // MUMPS wrote the solution in-place into pbc_rhs_buf.
    // Copy only the physical DOFs; the mortar multipliers (tail of pbc_rhs_buf)
    // are discarded — post-processing only needs the physical solution.
    if (rank == 0) {
      size_t n_own = tmsh.num_owned_nodes ();
      phi->get_owned_data ().assign (pbc_rhs_buf.begin (),
                                     pbc_rhs_buf.begin () + n_own);
    }
  } else {
    (*phi) = mumps_solver.get_distributed_solution ();
  }

  if (size > 1)
    bim3a_solution_with_ghosts (tmsh, *phi, replace_op);

  if (rank == 0) {
    const auto& sol = phi->get_owned_data ();
    std::cout << "\n=== [ Solution debug ] ===\n";
    static const char* fname[6] = {"x-","x+","y-","y+","z-","z+"};
    for (int face : {0, 1, 2, 3, 4, 5}) {
      bool is_pbc = ((face < 2) && periodic_x) || ((face >= 2 && face < 4) && periodic_y);
      std::vector<std::pair<size_t, double>> fnodes;
      for (auto q = tmsh.begin_quadrant_sweep ();
           q != tmsh.end_quadrant_sweep (); ++q) {
        for (auto local_idx : q->ef (face)) {
          // ef() returns LOCAL node indices; use gt() to convert to global DOF
          size_t global = static_cast<size_t> (q->gt (local_idx));
          if (global < sol.size ())
            fnodes.emplace_back (global, sol[global]);
        }
      }
      std::sort (fnodes.begin (), fnodes.end ());
      fnodes.erase (std::unique (fnodes.begin (), fnodes.end (),
        [] (auto& a, auto& b) { return a.first == b.first; }), fnodes.end ());
      if (fnodes.empty ()) continue;
      double fmin = fnodes[0].second, fmax = fnodes[0].second;
      for (auto& p : fnodes) { fmin = std::min (fmin, p.second); fmax = std::max (fmax, p.second); }
      std::cout << "  face " << fname[face]
                << (is_pbc ? " [PBC]" : " [Dir]")
                << "  nodes=" << fnodes.size ()
                << "  phi=[" << fmin << ", " << fmax << "]\n";
    }
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

  //CSR
  std::vector<double> vals;
  std::vector<int> irow, jcol;


  (*A).csr (vals, jcol, irow);

  // lis RHS
  // ln = local size including mortar DOFs appended in assemple_system_matrix
  LIS_INT i, is, ie, n_rhs, ln;
  LIS_VECTOR rhs_lis;
  ln = rhs->get_owned_data ().size ();  // N_phys_owned [+ n_pairs*ndofm for PBC]

  lis_vector_create (mpicomm, &rhs_lis);
  lis_vector_set_size (rhs_lis, ln, 0);
  lis_vector_get_range (rhs_lis, &is, &ie);

  for (i = is; i < ie; i++)
    lis_vector_set_value (LIS_INS_VALUE, i, rhs->get_owned_data ()[i - is], rhs_lis);

  //cleaning of rhs
  rhs.reset ();

  //lis_vector_print(rhs_lis);
  // lis PHI
  LIS_VECTOR phi_lis;

  lis_vector_create (mpicomm, &phi_lis);
  lis_vector_set_size (phi_lis, ln, 0);
  lis_vector_get_range (phi_lis, &is, &ie);

  // lis MATRIX
  LIS_INT n, nnz; //n: matrix dim ; nnz: numb of non zero elems
  LIS_INT *index; //array of integer containing the col index of non zero elems
  LIS_INT *ptr; //array of integer with starting points of rows
  LIS_SCALAR *value; //array of double stores non-zero elements of matrix A along the row
  LIS_MATRIX A_lis; //array of integer containing the col index of non zero elems

  nnz = vals.size ();  // total nnz from csr(), includes mortar rows for PBC
  n   = ln;            // augmented local size (= num_owned_nodes + n_pairs*ndofm for PBC)

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

  phi = std::make_unique<distributed_vector> (tmsh.num_owned_nodes (), mpicomm);

  // For PBC ln > num_owned_nodes(): extract only the physical DOFs.
  LIS_INT n_phys = tmsh.num_owned_nodes ();
  lis_vector_get_values (phi_lis, is, n_phys, phi->get_owned_data ().data ());

  lis_vector_destroy (phi_lis);

  if (size > 1)
    bim3a_solution_with_ghosts (tmsh, *phi, replace_op);

  if (rank == 0) {
    const auto& sol = phi->get_owned_data ();
    std::cout << "\n=== [ Solution debug (LIS) ] ===\n";
    static const char* fname[6] = {"x-","x+","y-","y+","z-","z+"};
    for (int face : {0, 1, 2, 3, 4, 5}) {
      bool is_pbc = ((face < 2) && periodic_x) || ((face >= 2 && face < 4) && periodic_y);
      std::vector<std::pair<size_t, double>> fnodes;
      for (auto q = tmsh.begin_quadrant_sweep ();
           q != tmsh.end_quadrant_sweep (); ++q) {
        for (auto local_idx : q->ef (face)) {
          size_t global = static_cast<size_t> (q->gt (local_idx));
          if (global < sol.size ())
            fnodes.emplace_back (global, sol[global]);
        }
      }
      std::sort (fnodes.begin (), fnodes.end ());
      fnodes.erase (std::unique (fnodes.begin (), fnodes.end (),
        [] (auto& a, auto& b) { return a.first == b.first; }), fnodes.end ());
      if (fnodes.empty ()) continue;
      double fmin = fnodes[0].second, fmax = fnodes[0].second;
      for (auto& p : fnodes) { fmin = std::min (fmin, p.second); fmax = std::max (fmax, p.second); }
      std::cout << "  face " << fname[face]
                << (is_pbc ? " [PBC]" : " [Dir]")
                << "  nodes=" << fnodes.size ()
                << "  phi=[" << fmin << ", " << fmax << "]\n";
    }
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

