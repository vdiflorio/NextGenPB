//Prova del calcolo dell'energia

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

// Problem parameters
constexpr double e_0 = 8.85418781762e-12;	//Dielectric void const [F/m]
constexpr double kb = 1.380649e-23;		//Boltzmann constant [J/K]
constexpr double T = 273.15 + 25;		//Temperature [K]
//constexpr double T = 297.33421190000001;
constexpr double e = 1.602176634e-19;  	//Charge of an electron [C]
constexpr double N_av = 6.022e23;      	//Avogadro Number [mol^-1]
constexpr double Angs = 1e-10;         	//Angstrom [m]
constexpr double pi = 3.14159265358979323846; 

struct
poisson_boltzmann
{

  p4est_topidx_t simple_conn_num_vertices;
  p4est_topidx_t simple_conn_num_trees;
  //std::unique_ptr<double> simple_conn_p;
  //std::unique_ptr<double> simple_conn_p;
  double *simple_conn_p;
  p4est_topidx_t *simple_conn_t;
  std::vector<std::pair<p4est_topidx_t, p4est_topidx_t>> bcells;
  
  std::vector<NS::Atom> atoms;

  //Cubic mesh:
  double ll; //min value between all the coordinates 
  double rr; //max value between all the coordinates
  
  //Stretched mesh:
  double l_c[3]; //min x, y, z value
  double r_c[3]; //max x, y, z value
  double l_cr[3]; //refine box min x, y, z value
  double r_cr[3]; //refine box max x, y, z value
  
  //number of trees
  p4est_topidx_t num_trees[3];
  
  //mesh:
  int maxlevel;
  int minlevel;
  int unilevel;
  int outlevel;
  int mesh_shape;

  int maxlevel1 = 10, minlevel1 = 1;
  
  //model:
  int linearized;
  int bc;
  double e_in, e_out, ionic_strength; //[M]
  
  //surface:
  NS::surface_type surf_type;
  double surf_param;
  double stern_layer;
  unsigned num_threads;
  
  //algorithm:
  std::string linear_solver_name;
  std::string linear_solver_options;
  
  MPI_Comm mpicomm;
  tmesh_3d tmsh;

  std::string optionsfilename;
  std::string pqrfilename;
  std::string p4estfilename;
  std::string surffilename; 
  std::string markerfilename;

  std::vector<double> marker; 
  std::vector<double> epsilon; 
  std::vector<double> epsilon_in; 
  std::vector<double> epsilon_out; 
  // std::vector<double> rho_fixed; 
  std::vector<double> reaction; 
  std::vector<double> ones_in;

  std::unique_ptr<distributed_vector> phi;
  std::unique_ptr<distributed_vector> rho_fixed;
  
  std::unique_ptr<distributed_vector> phi_energy;
  std::unique_ptr<distributed_vector> rho_fixed_energy;


  poisson_boltzmann (int maxlevel_ = 9, int minlevel_ = 3, int unilevel_ = 5, int mesh_shape_ = 1,
                     int bc_ = 1, int linearized_ = 1, 
                     double e_in_ = 2.0, double e_out_ = 80.0, double ionic_strength_ = 0.145,
                     std::string linear_solver_name_ = "mumps", std::string linear_solver_options_ = "",
                     MPI_Comm mpicomm_ = MPI_COMM_WORLD)
    : maxlevel(maxlevel_),
      minlevel(minlevel_),
      unilevel(unilevel_),
      mesh_shape(mesh_shape_),
      bc(bc_),
      linearized(linearized_),
      e_in(e_in_),
      e_out(e_out_),
      ionic_strength(ionic_strength_),
      linear_solver_name(linear_solver_name_),
      linear_solver_options(linear_solver_options_),
      mpicomm(mpicomm_),
      tmsh(mpicomm)
  {  };

  double
  levelsetfun (double x, double y, double z);
  
  double
  is_in_ns_surf (ray_cache_t & ray_cache, double x, double y, double z);

  static int
  uniform_refinement (tmesh_3d::quadrant_iterator quadrant)
  { return 1; }
  
  
  void
  create_mesh ();
  
  int
  parse_options (int argc, char **argv);
  
  void 
  print_options ();
  
  void
  read_atoms_from_pqr (std::basic_istream<char> &inputfile);
  
  void
  write_atoms_to_pqr (std::basic_ostream<char> &outputfile);
  
  friend std::basic_istream<char>& 
  operator>> (std::basic_istream<char>& inputfile, NS::Atom &a);

  void
  init_tmesh ();
  
  void
  init_tmesh_with_refine_box ();

  bool
  is_in (const NS::Atom& i, tmesh_3d::quadrant_iterator q); 
  
  void
  refine_surface (ray_cache_t & ray_cache);

  void
  refine_only_surface (ray_cache_t & ray_cache);
  
  void
  create_markers (ray_cache_t & ray_cache);
  
  void
  export_tmesh (ray_cache_t & ray_cache);

  void
  export_marked_tmesh ();

  void
  export_p4est ();

  void
  mumps_compute_electric_potential ();
  
  void
  lis_compute_electric_potential ();
  
  
  void 
  surface_integrals_energy();
  
  // void 
  // surface_integrals_energy(ray_cache_t & ray_cache);
  
  void
  grid_energy(distributed_vector &phi, distributed_vector &rho_fixed, std::string int_rule, double &energy);
  
  void
  grid_energy_electric_field(distributed_vector &phi, std::vector<double> &epsilon, std::string int_rule, double &energy);
  
  double
  coulomb_boundary_conditions(double x, double y, double z);
  
  void
  analitic_potential();
  
  void
  abs_value_field(distributed_vector &phi);
  
};

std::basic_istream<char>& 
operator>>(std::basic_istream<char>& inputfile, NS::Atom &a);

#endif
