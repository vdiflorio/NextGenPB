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
#include <bim_distributed_vector.h>
#include <quad_operators_3d.h>
#include <mumps_class.h>
#include <lis_class.h>
#include <cmath>

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

  phi = std::make_unique<distributed_vector> (tmsh.num_owned_nodes ());
  (*phi) = mumps_solver.get_distributed_solution ();

  if (size > 1)
    bim3a_solution_with_ghosts (tmsh, *phi, replace_op);

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

  for (i=is; i<ie; i++)
    lis_vector_set_value (LIS_INS_VALUE, i, rhs->get_owned_data ()[i-is], rhs_lis);

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

  phi = std::make_unique<distributed_vector> (tmsh.num_owned_nodes (), mpicomm);

  lis_vector_get_values (phi_lis, is, ln, phi->get_owned_data ().data ());

  lis_vector_destroy (phi_lis);

  if (size > 1)
    bim3a_solution_with_ghosts (tmsh, *phi, replace_op);

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

  for (const NS::Atom& i : atoms) {
    dist = std::hypot ( (i.pos[0] - x), (i.pos[1] - y), (i.pos[2] - z));
    pot += i.charge*exp (-k*dist)/ (dist*eps_out);
  }

  return pot;
}

