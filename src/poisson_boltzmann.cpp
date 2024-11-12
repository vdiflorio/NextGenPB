#include <mpi.h>

#include "pb_class.h"

#include "wrapper_search.h"

int 
cerca_atomo_wrapper(p8est_t *p4est, 
                    p4est_topidx_t which_tree, 
                    p8est_quadrant_t *quadrant, 
                    p4est_locidx_t local_num, 
                    void *point){

  poisson_boltzmann * pb_wrapper = static_cast<poisson_boltzmann*> (pb_global);
  return pb_wrapper->cerca_atomo(p4est, which_tree, quadrant, local_num,point);
}

// static char filename[255];

double crossings_t::start[3] = {0., 0., 0.};
double crossings_t::end[3] = {0., 0., 0.};

int_coord_t ray_cache_t::count_cache = 0;
int_coord_t ray_cache_t::count_new = 0;


void
print_map (const std::array<std::map<std::array<double, 2>, crossings_t, map_compare>, 3>& r);

void
print_point (const std::array<std::vector<std::array<double, 2>>,3>& r);
void
save_ray_cache (nlohmann::json& j, const std::map<std::array<double, 2>, crossings_t, map_compare>& r);


int
main (int argc, char **argv)
{

  MPI_Init (&argc, &argv);

  // int recursive, partforcoarsen, balance;
  MPI_Comm mpicomm = MPI_COMM_WORLD;
  int rank, size;
  tmesh_3d tmsh;

  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  poisson_boltzmann pb;
  ray_cache_t ray_cache;

  pb_global = (void *) (&pb);

  if (pb.parse_options (argc, argv))
    return 1;

  if (rank == 0){
    std::ifstream inputfile (pb.pqrfilename);
    pb.read_atoms_from_pqr (inputfile);
    inputfile.close ();
    pb.read_atoms_from_class ();
  }

  MPI_Barrier (mpicomm);
  if (size > 1){
    pb.broadcast_vectors();
  }
  // if (rank == 0) {
  //   std::cout << "Atom : " << std::endl;
  //   pb.write_atoms_to_pqr (std::cout);
  //   pb.print_options ();
  // }

  MPI_Barrier (mpicomm);

  TIC ();
  // pb.create_mesh_ns ();
  pb.create_mesh ();
  TOC ("create_mesh");

  std::vector<double>().swap(pb.r_atoms);

  TIC ();
  if (rank == 0){
    ray_cache.init_analytical_surf_ns (pb.atoms, pb.surf_type, pb.surf_param, pb.stern_layer, pb.num_threads, pb.l_cr, pb.r_cr, pb.scale);
    ray_cache.ns->clean();
    std::vector<NS::Atom>().swap(pb.atoms);
  }

  TOC ("init analytical surf");
  MPI_Barrier (mpicomm);

  TIC ();
  if ( pb.mesh_shape == 0 || pb.refine_box == 1 ||pb.mesh_shape == 4 )
    pb.init_tmesh_with_refine_scale ();
  else if (pb.mesh_shape == 3)
    pb.init_tmesh_with_refine_box_scale ();
  else
    pb.init_tmesh ();

  TOC ("init_tmesh");

  // TIC ();

  // if ( rank == 0) {
  //   ray_cache.init_analytical_surf_ns (pb.atoms, pb.surf_type, pb.surf_param, pb.stern_layer, pb.num_threads, pb.l_cr, pb.r_cr, pb.scale);
  // }

  // TOC ("init analytical surf");
  // MPI_Barrier (mpicomm);
  TIC ();

  if (pb.loc_refinement == 1 || pb.mesh_shape == 4)
    pb.refine_surface (ray_cache);

  TOC ("refine the box");
  
  // if ( rank == 0)
  //   ray_cache.ns->clean();
  


  TIC ();
  pb.create_markers (ray_cache);
  // pb.create_markers_prova (ray_cache);
  TOC ("create element markers");


  
  
  TIC ();
  if (pb.linear_solver_name == "mumps")
    pb.mumps_compute_electric_potential (ray_cache);
  else if (pb.linear_solver_name == "lis")
    pb.lis_compute_electric_potential (ray_cache);
  else {
    std::cerr << "Invalid linear solver selected" << std::endl;
    return 1;
  }
  TOC ("compute potential");

  if (pb.atoms_write == 1) {
    TIC ();
      pb.write_potential_on_atoms_fast ();
    TOC ("Write potential on atoms")
  }
  
  if (pb.calc_energy > 0) {
    TIC ();
    // pb.energy (ray_cache);
    pb.energy_fast (ray_cache);
    TOC ("compute energy")
  }


  if (pb.potential_map == 1){
    TIC ();
    pb.export_potential_map (ray_cache);
    TOC ("export potential map");
  }

  // TIC ();
  // pb.export_tmesh (ray_cache);
  // TOC ("export tmesh");

  // TIC ();
  // pb.export_marked_tmesh ();
  // TOC ("export marked tmesh");

  
  if (rank == 0) {
    print_timing_report();

    // if (pb.surf_type != 2)
    // {
    // //Save ray_cache:
    // for (int i = 0; i < 3; ++i)
    // {
    // // std::cout << "Direzione: "<< i <<std::endl;
    // // std::cout << std::endl;

    // nlohmann::json j;
    // save_ray_cache (j, ray_cache.rays[i]);

    // // std::cout << "Count cached rays: " << ray_cache.count_cache_dir[i] << std::endl;
    // // std::cout << "Count new rays: " << ray_cache.count_new_dir[i] << std::endl;
    // // std::cout << std::endl;
    // std::ofstream ray_cached_file;
    // std::string filename = "ray_cache_";
    // std::string extension = ".json";
    // filename += std::to_string(i);
    // filename += extension;
    // ray_cached_file.open (filename.c_str ());

    // if (ray_cached_file.is_open ())
    // ray_cached_file << j;

    // ray_cached_file.close ();

    // //Alternative way to save results:
    // print_map (ray_cache.rays);
    // }
    // }
  }

  MPI_Barrier (mpicomm);

  MPI_Finalize ();

  return 0;

}


void
print_point (const std::array<std::vector<std::array<double, 2>>,3>& r)
{
  std::ofstream ray_cached_file;

  // ray_cached_file.open ("ray_cache_ns.txt");
  for (int i = 0; i < 3; ++i) {
    std::string filename = "ray_point_ns_";
    std::string extension = ".txt";
    filename += std::to_string (i);
    filename += extension;
    ray_cached_file.open (filename.c_str ());

    if (ray_cached_file.is_open ()) {
      for (auto it = r[i].begin(); it !=r[i].end(); ++it) {
        ray_cached_file <<std::setprecision (9)<< "[[" << (*it)[0] << ", " << (*it)[1] << "]" << std::endl;
      }
    }

    ray_cached_file.close ();
  }

}


void
print_map (const std::array<std::map<std::array<double, 2>, crossings_t, map_compare>, 3>& r)
{
  std::ofstream ray_cached_file;

  // ray_cached_file.open ("ray_cache_ns.txt");
  for (int i = 0; i < 3; ++i) {
    std::string filename = "ray_cache_ns_";
    std::string extension = ".txt";
    filename += std::to_string (i);
    filename += extension;
    ray_cached_file.open (filename.c_str ());

    if (ray_cached_file.is_open ()) {
      ray_cached_file << "Count cached rays: " << ray_cache_t::count_cache << std::endl;
      ray_cached_file << "Count new rays: " << ray_cache_t::count_new << std::endl;
      int count = 0;

      for (auto it : r[i]) {
        ray_cached_file <<std::setprecision (9)<< "[[" << it.first.at (0) << ", " << it.first.at (1) << "]";

        if (it.second.inters.size ()>0) {
          count++;
        }

        for (int i = 0; i < it.second.inters.size (); i++)
          ray_cached_file << ", " << it.second.inters[i];

        ray_cached_file << "]" << std::endl;
      }

      ray_cached_file << std::endl;
      ray_cached_file << "num raggi inters: " << count << std::endl;
    }

    ray_cached_file.close ();
  }

}
/*
void
save_ray_cache (nlohmann::json& j, const std::map<std::array<double, 2>, crossings_t, map_compare>& r)
{
  std::vector<double> ints;

  for (auto it : r) {
    ints.clear();

    j += nlohmann::json{{"ray", it.first}, {"inters", (it.second.inters)}, {"normals", (it.second.normals)}};

    // if (it.second.inters.size()>0)
    // {
    // for (int i = 0; i < it.second.inters.size(); ++i)
    // {
    // std::cout<< std::hypot(it.first[0],it.first[1], it.second.inters[i]) << "  ";
    // }
    // std::cout<<std::endl;
    // }
  }
}

*/