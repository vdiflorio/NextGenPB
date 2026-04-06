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

/**
 * @file pb_class.h
 * @brief Main header for the NextGenPB Poisson-Boltzmann solver.
 *
 * Defines the `poisson_boltzmann` struct, which encapsulates all data and
 * methods needed to set up, solve, and post-process the linearized or
 * non-linearized Poisson-Boltzmann equation on an adaptive octree mesh
 * (p4est/p8est) using finite elements (BIM++) in a distributed-memory
 * MPI environment.
 *
 * The solver pipeline:
 * 1. **I/O**:        parse_options(), read_atoms_from_*(), broadcast_vectors()
 * 2. **Mesh**:       create_mesh*(), init_tmesh*()
 * 3. **Surface**:    create_markers(), refine_surface()
 * 4. **Assembly**:   assemple_system_matrix()
 * 5. **Solve**:      mumps_compute_electric_potential() or lis_compute_electric_potential()
 * 6. **Post-proc**:  energy(), pot_field*(), write_potential_on_*()
 */
#ifndef HAVE_PB_CLASS_H
#define HAVE_PB_CLASS_H

#include <bim_timing.h>
#include <tmesh_3d.h>

#define TIC() MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
#define TOC(S) MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc (S); }
#include <cmath>

const double p4esttol = 1 / std::pow (2, P8EST_QMAXLEVEL);

#include <array>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <memory>

#include "raytracer.h"
#include <nanoshaper.h>


#include "wrapper_search.h"


/// @defgroup constants Physical Constants
/// @{
constexpr double e_0  = 8.85418781762e-12; ///< Permittivity of free space [F/m]
constexpr double kb   = 1.380649e-23;      ///< Boltzmann constant [J/K]
constexpr double e    = 1.602176634e-19;   ///< Elementary charge [C]
constexpr double N_av = 6.022e23;          ///< Avogadro's number [mol^-1]
constexpr double Angs = 1e-10;             ///< Angstrom to metre conversion [m]
constexpr double pi   = 3.14159265358979323846; ///< Mathematical constant π
/// @}


/**
 * @brief Main struct for the NextGenPB Poisson-Boltzmann solver.
 *
 * Holds all configuration parameters, mesh data, distributed linear algebra
 * objects, and solver methods. Instantiated once per run; all MPI ranks share
 * the same struct layout, with data partitioned via p4est and BIM++.
 */
struct
  poisson_boltzmann {

  /// @name Octree mesh connectivity (p4est)
  /// @{
  p4est_topidx_t simple_conn_num_vertices; ///< Number of mesh vertices
  p4est_topidx_t simple_conn_num_trees;    ///< Number of trees in the octree
  std::unique_ptr<double[]> simple_conn_p;          ///< Vertex coordinate array
  std::unique_ptr<p4est_topidx_t[]> simple_conn_t;  ///< Tree-to-vertex connectivity
  std::vector<std::pair<p4est_topidx_t, p4est_topidx_t>> bcells; ///< Boundary cell pairs
  /// @}

  /// @name Atomic data
  /// @{
  std::vector<NS::Atom> atoms;                  ///< Full atom list (NanoShaper format)
  std::vector<std::array<double,3>> pos_atoms;  ///< Atomic positions [Å]
  std::vector<double> charge_atoms;             ///< Partial charges [e]
  std::vector<double> r_atoms;                  ///< Van der Waals radii [Å]
  std::vector<int> index_atoms;                 ///< Atom index mapping
  /// @}

  /// @name Input / file options
  /// @{
  std::string filetype;         ///< Input file type ("pqr", "pdb", …)
  std::string radiusfilename;   ///< File with van der Waals radii
  std::string chargefilename;   ///< File with partial charges
  std::string name_pqr;         ///< Base name for PQR output
  int write_pqr;                ///< Write PQR file if non-zero
  /// @}

  /// @name Domain bounding boxes
  /// @{
  double cc[3];      ///< Center of the bounding box [Å]
  double ll[3];      ///< Cubic mesh: min (x,y,z) corner [Å]
  double rr[3];      ///< Cubic mesh: max (x,y,z) corner [Å]
  double l_c[3];     ///< Stretched mesh: min (x,y,z) corner [Å]
  double r_c[3];     ///< Stretched mesh: max (x,y,z) corner [Å]
  double l_cr[3];    ///< Refined sub-box: min corner [Å]
  double r_cr[3];    ///< Refined sub-box: max corner [Å]
  double l_box[3];   ///< Focusing sub-box: min corner [Å]
  double r_box[3];   ///< Focusing sub-box: max corner [Å]
  double pot_bc = 0.0; ///< Potential value at outer boundary [kT/e]
  p4est_topidx_t num_trees[3]; ///< Number of octree trees per axis
  double len;        ///< Side length of the coarsest tree [Å]
  double cc_focusing[3]; ///< Center for focusing mesh [Å]
  int n_grid;        ///< Grid resolution for focusing region
  /// @}

  /// @name Mesh control parameters
  /// @{
  int maxlevel;           ///< Maximum octree refinement level
  int minlevel;           ///< Minimum octree refinement level
  int unilevel;           ///< Uniform (base) refinement level
  int outlevel;           ///< Refinement level in outer region
  int loc_refinement;     ///< Enable local refinement near surface
  int mesh_shape;         ///< 0 = cubic, 1 = stretched, 2 = custom
  int refine_box;         ///< Enable box-based local refinement
  int rand_center;        ///< Randomise center (avoid symmetry artifacts)
  int scale_level;        ///< Scale-based refinement level
  int scale_level_min_box;///< Min scale level inside the box
  double scale, scale_min, scale_max; ///< Scaling parameters for stretched mesh
  double perfil1, perfil2;            ///< Stretching profile control points
  int loc_ref = 0;   ///< Local refinement flag
  int aligned = 0;   ///< Align mesh to axes flag
  /// @}

  /// @name Physical model parameters
  /// @{
  int linearized;          ///< 1 = linearized PBE, 0 = non-linear
  int bc;                  ///< Boundary condition type (0 = Coulomb, …)
  double e_in;             ///< Dielectric constant inside the solute
  double e_out;            ///< Dielectric constant of the solvent
  double ionic_strength;   ///< Ionic strength [mol/L]
  double T;                ///< Temperature [K]
  int calc_energy;         ///< Compute solvation energy if non-zero
  double energy_pol  = 0.0;  ///< Polarisation energy contribution [kT]
  double energy_react = 0.0; ///< Reaction field energy contribution [kT]
  double coul_energy  = 0.0; ///< Coulombic energy contribution [kT]
  int calc_coulombic;      ///< Include Coulombic energy term
  int calc_potential_term; ///< Compute potential on atoms
  int calc_field_term;     ///< Compute electric field on atoms
  /// @}

  /// @name Molecular surface parameters
  /// @{
  NS::surface_type surf_type;    ///< Surface type (SES, SAS, …)
  int surf_type_num = 0;         ///< Numeric code for surface type
  double surf_param;             ///< Surface smoothing parameter
  double prb_radius = 1.4;       ///< Probe radius for SES [Å]
  int stern_layer_surf;          ///< Enable Stern layer on solute surface
  double stern_layer;            ///< Stern layer thickness [Å]
  unsigned num_threads;          ///< Number of threads for NanoShaper
  /// @}

  /// @name Solver options
  /// @{
  std::string linear_solver_name;    ///< "mumps" or "lis"
  std::string linear_solver_options; ///< Options string passed to LIS
  /// @}

  /// @name MPI and mesh objects
  /// @{
  MPI_Comm mpicomm; ///< MPI communicator
  tmesh_3d tmsh;    ///< Adaptive octree mesh (BIM++ wrapper around p8est)
  /// @}

  /// @name File name configuration
  /// @{
  std::string optionsfilename;  ///< Path to options/config file
  std::string pqrfilename;      ///< Path to input PQR file
  std::string p4estfilename;    ///< Path to p4est checkpoint file
  std::string surffilename;     ///< Path to NanoShaper surface file
  std::string markerfilename;   ///< Path to saved marker file
  std::string pqr_atoms;        ///< Inline PQR data (alternative to file)
  /// @}

  /// @name Post-processing flags
  /// @{
  int atoms_write;     ///< Write potential on atoms if non-zero
  int surf_write;      ///< Write potential on surface if non-zero
  std::string map_type; ///< Type of potential map to export
  int potential_map;   ///< Enable potential map output
  int eps_map;         ///< Enable dielectric map output
  std::map<int, tmesh_3d::quadrant_t> lookup_table; ///< Quadrant lookup cache
  /// @}

  /// @name Per-quadrant field arrays (cell-centred)
  /// @{
  std::vector<double> marker;      ///< Quadrant classification: 0=inside, 0.5=boundary, 1=outside
  std::vector<double> marker_k;    ///< Quadrant classification including Stern layer
  std::vector<double> epsilon;     ///< Dielectric constant per quadrant
  std::vector<double> epsilon_in;  ///< Interior dielectric per quadrant
  std::vector<double> epsilon_out; ///< Exterior dielectric per quadrant
  std::vector<double> reaction;    ///< Reaction (κ²) parameter per quadrant
  std::vector<int> border_quad;    ///< Indices of boundary quadrants
  std::vector<double> const_ones;  ///< Constant-1 auxiliary vector
  /// @}

  /// @name Distributed node-centred fields (BIM++)
  /// @{
  std::unique_ptr<distributed_vector> markn;          ///< Node-level marker
  std::unique_ptr<distributed_vector> epsilon_nodes;  ///< Node-level dielectric
  std::unique_ptr<distributed_vector> reaction_nodes; ///< Node-level κ² field
  double net_charge; ///< Total charge of the solute [e]
  std::unique_ptr<distributed_vector> phi;       ///< Electric potential solution [kT/e]
  std::unique_ptr<distributed_vector> rho_fixed; ///< Fixed charge density RHS
  std::unique_ptr<distributed_vector> ones;      ///< Auxiliary all-ones vector
  std::unique_ptr<distributed_vector> rhs;       ///< Right-hand side vector
  std::unique_ptr<distributed_sparse_matrix> A;  ///< System matrix
  /// @}

  /// @name Membrane and periodic boundary conditions
  /// @{

  bool membrane_enabled = false; ///< Enable membrane mode

  /// @name Lipid atom data (read from a separate PQR/PDB file)
  /// @{
  std::string lipid_file;                             ///< Path to lipid PQR/PDB file
  std::string lipid_filetype;                         ///< "pqr" or "pdb"
  std::vector<NS::Atom>             lipid_atoms;      ///< Lipid atom list
  std::vector<std::array<double,3>> pos_lipid_atoms;  ///< Lipid positions [Å]
  std::vector<double>               charge_lipid_atoms; ///< Lipid partial charges [e]
  std::vector<double>               r_lipid_atoms;    ///< Lipid van der Waals radii [Å]
  /// @}

  /// @name Periodic cell geometry
  /// @{
  bool   periodic_x    = false; ///< Periodic boundary condition in x
  bool   periodic_y    = false; ///< Periodic boundary condition in y
  double cell_length_x = 0.0;  ///< Unit cell length in x [Å]
  double cell_length_y = 0.0;  ///< Unit cell length in y [Å]
  // [PLACEHOLDER] NanoShaper replication parameters (TBD)
  /// @}

  /// @name Membrane dielectric and ionic profile
  /// @{
  double e_mem          = 2.0;  ///< Membrane dielectric constant
  double kappa_in       = 0.0;  ///< Ionic strength — intracellular side [mol/L]
  double kappa_out      = 0.0;  ///< Ionic strength — extracellular side [mol/L]
  bool   stern_membrane = false; ///< Enable Stern layer on membrane surface
  double stern_membrane_d = 0.0; ///< Stern layer thickness on membrane [Å]
  /// @}

  /// @}

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

  /**
   * @brief Determines whether a point is inside the Stern layer along a specified direction.
   *
   * This function evaluates whether a point defined by coordinates `(x, y, z)` lies
   * inside, outside, or on the boundary of the Stern layer for a given direction.
   * The evaluation uses precomputed ray intersections stored in a ray cache.
   *
   * @param ray_cache A reference to the ray tracing cache that stores intersection
   *                  data and handles ray operations.
   * @param x The x-coordinate of the point.
   * @param y The y-coordinate of the point.
   * @param z The z-coordinate of the point.
   * @param dir The direction of the evaluation:
   *            - `0`: Evaluate in the yz-plane.
   *            - `1`: Evaluate in the xz-plane.
   *            - `2`: Evaluate in the xy-plane.
   *
   * @return A value indicating the position of the point relative to the Stern layer:
   *         - `0.0`: The point is outside the Stern layer.
   *         - `1.0`: The point is inside the Stern layer.
   *         - `-1.0`: The point requires additional ray tracing or data is unavailable.
   *
   * ### Algorithm Details
   * - Coordinates `(x, y, z)` are reordered based on the evaluation direction (`dir`).
   * - Intersections along the specified direction are retrieved from the ray cache.
   * - If the point lies outside all intersections, it is marked as outside.
   * - Iteratively evaluates whether the point alternates between inside and outside
   *   based on the intersection list, accounting for the Stern layer thickness.
   *
   * ### Notes
   * - Requires the `ray_cache` to be properly initialized and populated with
   *   intersection data.
   * - Assumes a uniform thickness for the Stern layer, defined as `stern_layer`.
   * - The `sign` variable alternates to evaluate the nesting of intersections.
   */
  double
  is_in_ns_surf_stern (ray_cache_t & ray_cache, double x, double y, double z, int dir);

  /// @brief p4est refinement callback that always returns 1 (refine every quadrant).
  static int
  uniform_refinement (tmesh_3d::quadrant_iterator quadrant)
  {
    return 1;
  }

  /// @name Mesh generation
  /// @{

  /// @brief Build the adaptive octree mesh from atomic positions (standard mode).
  void
  create_mesh ();

  /// @brief Build the adaptive octree mesh using a NanoShaper-defined domain.
  void
  create_mesh_ns ();

  /// @brief Build the adaptive octree mesh with scale-based stretching.
  void
  create_mesh_scale ();

  /// @}

  /// @name Option parsing and I/O
  /// @{

  /// @brief Parse command-line arguments and populate all configuration fields.
  /// @param argc Argument count from main().
  /// @param argv Argument values from main().
  /// @return 0 on success, non-zero on error.
  int
  parse_options (int argc, char **argv);

  /// @brief Print the current configuration to stdout (rank 0 only).
  void
  print_options ();

  /// @brief Read atoms from a PQR-format input stream.
  /// @param inputfile Opened input stream for a .pqr file.
  void
  read_atoms_from_pqr (std::basic_istream<char> &inputfile);

  /// @brief Read atoms from a PDB-format input stream (charges and radii loaded separately).
  /// @param inputfile Opened input stream for a .pdb file.
  void
  read_atoms_from_pdb (std::basic_istream<char> &inputfile);

  /// @brief Read atoms from the inline `pqr_atoms` string field.
  void
  read_atoms_from_class ();

  /// @brief Broadcast atomic data (positions, charges, radii) from rank 0 to all MPI ranks.
  void
  broadcast_vectors ();

  /// @brief Write atomic data to a PQR-format output stream.
  /// @param outputfile Opened output stream for a .pqr file.
  void
  write_atoms_to_pqr (std::basic_ostream<char> &outputfile);

  friend std::basic_istream<char>&
  operator>> (std::basic_istream<char>& inputfile, NS::Atom &a);

  friend std::basic_istream<char>&
  operator>> (std::basic_istream<char>& inputfile, std::array<float,5> &a);

  /// @brief Write the electric potential evaluated at each atom position to file.
  void
  write_potential_on_atoms ();

  /// @brief Write the electric potential at atom positions using fast interpolation.
  void
  write_potential_on_atoms_fast ();

  /// @}

  /// @name Mesh initialisation and refinement
  /// @{

  /// @brief Initialise the p4est forest with uniform refinement (level = unilevel).
  void
  init_tmesh ();

  /// @brief Initialise the p4est forest with scale-based refinement.
  void
  init_tmesh_with_refine_scale ();

  /// @brief Initialise the p4est forest with scale-based refinement inside a box.
  void
  init_tmesh_with_refine_box_scale ();

  /// @}

  /// @name Surface detection and marker creation
  /// @{

  /// @brief Test whether atom @p i is geometrically inside quadrant @p q.
  /// @param i Atom to test.
  /// @param q Octree quadrant to test against.
  /// @return True if the atom centre lies inside the quadrant.
  bool
  is_in (const NS::Atom& i, tmesh_3d::quadrant_iterator q);

  /// @brief Test whether atom @p i is geometrically inside a refined quadrant @p q.
  bool
  is_in_ref (const NS::Atom& i, tmesh_3d::quadrant_iterator q);

  /// @brief Refine quadrants near the molecular surface until maxlevel is reached.
  /// @param ray_cache Populated ray cache for surface intersection queries.
  void
  refine_surface (ray_cache_t & ray_cache);

  /// @brief Refine only the quadrants that straddle the molecular surface (no interior refinement).
  /// @param ray_cache Populated ray cache for surface intersection queries.
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
   * - Updates reaction and dielectric properties for nodes inside the molecule or
   *   Stern layer.
   * - Handles MPI-based parallelism, including barrier synchronization and data
   *   exchange.
   * - Ensures rays are calculated and cached for points near molecular boundaries.
   *
   * ### Stern Layer Handling
   * If `stern_layer_surf` is set to `1`, the function also processes Stern layer
   * interactions, updating markers and reaction values accordingly.
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

  /// @}

  /// @name Mesh and field export
  /// @{

  /// @brief Export the marked mesh to VTK format for visualisation.
  /// @param ray_cache Ray cache used to evaluate the surface field.
  void
  export_tmesh (ray_cache_t & ray_cache);

  /// @brief Export the electric potential field to a binary map file.
  /// @param ray_cache Ray cache used during export.
  void
  export_potential_map (ray_cache_t & ray_cache);

  /// @brief Export the mesh with marker values to VTK format.
  void
  export_marked_tmesh ();

  /// @brief Save the current p4est forest to a checkpoint file.
  void
  export_p4est ();

  /// @}

  /// @name System assembly and solvers
  /// @{

  /// @brief Assemble the FEM system matrix A and right-hand side rhs.
  ///
  /// Iterates over all local quadrants, integrates the bilinear form
  /// (ε ∇φ, ∇v) + (κ² φ, v) using quadrature, and assembles the
  /// distributed sparse matrix A and vector rhs with boundary conditions.
  ///
  /// @param ray_cache Ray cache for surface intersection queries.
  void
  assemple_system_matrix (ray_cache_t & ray_cache);

  /// @brief Compute the fixed-charge density map on the octree nodes.
  /// @param ray_cache Ray cache for surface intersection queries.
  void
  create_density_map (ray_cache_t & ray_cache);

  /// @brief Solve the system Aφ = rhs using the MUMPS direct solver.
  /// @param ray_cache Ray cache (used for ghost synchronisation after solve).
  void
  mumps_compute_electric_potential (ray_cache_t & ray_cache);

  /// @brief Solve the system Aφ = rhs using the LIS iterative solver.
  /// @param ray_cache Ray cache (used for ghost synchronisation after solve).
  void
  lis_compute_electric_potential (ray_cache_t & ray_cache);

  /// @}

  /// @name Marching cubes helpers
  /// @{

  /// @brief Classify a quadrant by the standard marching-cubes algorithm.
  /// @param quadrant Quadrant to classify.
  /// @param isolevel Iso-value for surface extraction.
  /// @return Cube index (0–255) for the marching-cubes lookup table.
  int
  classifyCube (tmesh_3d::quadrant_iterator& quadrant,double isolevel);

  /// @brief Classify a quadrant using the fast (precomputed) marching-cubes path.
  /// @param quadrant Quadrant to classify.
  /// @param isolevel Iso-value for surface extraction.
  /// @return Cube index (0–255).
  int
  classifyCube_fast (tmesh_3d::quadrant_iterator& quadrant,double isolevel);

  /// @brief Classify a quadrant and return potential/dielectric values on edges with flux information.
  /// @param quadrant  Quadrant to classify.
  /// @param tmp_phi   Output: potential values at the 8 nodes.
  /// @param tmp_eps   Output: dielectric values at the 8 nodes.
  /// @return Tuple (phi, eps, edge_indices, flux_directions).
  std::tuple<std::array<double,8>, std::array<double,8>, std::vector<int>,std::vector<int>>
  classifyCube_flux (tmesh_3d::quadrant_iterator& quadrant,
                     std::array<double,8>& tmp_phi,
                     std::array<double,8>& tmp_eps);

  /// @brief Fast version of classifyCube_flux using precomputed fields.
  std::tuple<std::array<double,8>, std::array<double,8>, std::vector<int>,std::vector<int> >
  classifyCube_flux_fast (tmesh_3d::quadrant_iterator& quadrant,
                          std::array<double,8>& tmp_phi,
                          std::array<double,8>& tmp_eps);

  /// @}

  /// @name Post-processing: energy and fields
  /// @{

  /// @brief Compute and print the solvation energy components (standard path).
  /// @param ray_cache Ray cache for surface intersection queries.
  void
  energy (ray_cache_t & ray_cache);

  /// @brief Compute and print the solvation energy components using fast interpolation.
  /// @param ray_cache Ray cache for surface intersection queries.
  void
  energy_fast (ray_cache_t & ray_cache);

  /// @brief Write the electric potential sampled on the molecular surface to file.
  /// @param ray_cache Ray cache for surface intersection queries.
  void
  write_potential_on_surface (ray_cache_t & ray_cache);

  /// @brief Compute the Coulomb (Debye-Hückel) boundary potential at a point.
  /// @param x X-coordinate [Å].
  /// @param y Y-coordinate [Å].
  /// @param z Z-coordinate [Å].
  /// @return Boundary potential value [kT/e].
  double
  coulomb_boundary_conditions (double x, double y, double z);

  /// @brief Evaluate the analytical Debye-Hückel solution at a single point.
  /// @param x X-coordinate [Å].
  /// @param y Y-coordinate [Å].
  /// @param z Z-coordinate [Å].
  /// @return Analytical potential value [kT/e].
  double
  analytic_solution (double x, double y, double z);

  /// @brief Compute the analytical Debye-Hückel potential on the full mesh and export.
  void
  analitic_potential ();

  /// @brief Replace each nodal value in @p phi with its absolute value (in-place).
  /// @param phi Distributed vector to modify.
  void
  abs_value_field (distributed_vector &phi);

  /// @brief Compute the fraction of each quadrant edge that lies inside the molecular surface.
  /// @param quadrant  Quadrant to analyse.
  /// @param ray_cache Ray cache for surface queries.
  /// @return Array of 12 fractional intersection lengths (one per edge).
  std::array<double,12>
  cube_fraction_intersection (tmesh_3d::quadrant_iterator& quadrant,
                              const ray_cache_t & ray_cache);

  /// @brief Compute the surface normal and edge fraction at a boundary edge intersection.
  /// @param quadrant  Quadrant containing the edge.
  /// @param ray_cache Ray cache for surface queries.
  /// @param edge      Edge index (0–11).
  /// @param norm      Output: outward unit normal at the intersection point.
  /// @param frac      Output: fractional distance from the first node to the intersection.
  void
  normal_intersection (tmesh_3d::quadrant_iterator& quadrant,
                       const ray_cache_t & ray_cache,
                       int edge, std::array<double,3> &norm,double &frac);

  /// @brief Look up the triangle connectivity for a marching-cubes cube index.
  /// @param cubeindex   Cube classification index (0–255).
  /// @param triangles   Output: up to 5 triangles, each with 3 vertex indices.
  /// @return Number of triangles generated.
  int
  getTriangles (int cubeindex,
                std::array<std::array<int,3>,5> &triangles);

  /// @}

  /// @name p4est search callbacks
  /// @{

  /// @brief p4est point-search callback: checks if atom @p i is inside @p quadrant.
  /// @param i        Atom index.
  /// @param quadrant p4est quadrant to test.
  /// @return True if the atom lies inside the quadrant.
  bool
  controlla_coordinate (int i, const p8est_quadrant_t *quadrant);

  /// @brief p4est search callback used during atom-to-quadrant mapping.
  int
  cerca_atomo (p8est_t * p4est,
               p4est_topidx_t which_tree,
               p8est_quadrant_t * quadrant,
               p4est_locidx_t local_num,
               void *point);

  /// @brief Perform p4est point search to map all atoms to their owning quadrants.
  void
  search_points ();

  /// @}

  /// @name Electric potential field output
  /// @{

  /// @brief Compute and export the electric potential and field on atom positions (fast path).
  /// @param ray_cache Ray cache for surface intersection queries.
  void
  pot_field_fast (ray_cache_t & ray_cache);

  /// @brief Compute and export the electric potential and field on atom positions (standard path).
  /// @param ray_cache Ray cache for surface intersection queries.
  void
  pot_field (ray_cache_t & ray_cache);

  /// @}
};

std::basic_istream<char>&
operator>> (std::basic_istream<char>& inputfile, NS::Atom &a);
std::basic_istream<char>&
operator>> (std::basic_istream<char>& inputfile, std::array<float,5> &a);
#endif
