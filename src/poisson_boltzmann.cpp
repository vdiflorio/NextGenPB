#include <mpi.h>

#include "pb_class.h"

static char filename[255];

double crossings_t::start[3] = {0., 0., 0.};
double crossings_t::end[3] = {0., 0., 0.};

int_coord_t ray_cache_t::count_cache = 0;
int_coord_t ray_cache_t::count_new = 0;


void
print_map(const std::map<std::array<double, 2> , crossings_t, map_compare>& r); 

void
save_ray_cache (nlohmann::json& j, const std::map<std::array<double, 2> , crossings_t, map_compare>& r);


int
main (int argc, char **argv)
{

  MPI_Init (&argc, &argv);
  
  int                   recursive, partforcoarsen, balance;
  MPI_Comm              mpicomm = MPI_COMM_WORLD;  
  int                   rank, size;
  tmesh_3d              tmsh;

  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);
 
  poisson_boltzmann pb;
  ray_cache_t ray_cache;
  
  if (pb.parse_options (argc, argv))
    return 1;
  
  std::ifstream inputfile (pb.pqrfilename);
  pb.read_atoms_from_pqr (inputfile);
  inputfile.close ();
  
  if (rank == 0)
    {
      std::cout << "Atom : " << std::endl;
      pb.write_atoms_to_pqr (std::cout);
      pb.print_options ();
    }
  
  MPI_Barrier (mpicomm);
  
  TIC ();
  pb.create_mesh ();
  TOC ("create_mesh");
  

  TIC ();

  if (pb.refine_box == 1)
    pb.init_tmesh_with_refine_box ();
  else
    pb.init_tmesh ();
  
  TOC ("init_tmesh");
  
  TIC ();
  if (pb.surf_type != 2 && rank == 0)
    ray_cache.init_analytical_surf (pb.atoms, pb.surf_type, pb.surf_param, pb.stern_layer, pb.num_threads);
  TOC ("init analytical surf");
  TIC ();
  pb.refine_surface (ray_cache);
  TOC ("refine the box");
  
  TIC ();
  pb.create_markers (ray_cache);
  TOC ("create element markers");
  
  TIC ();
  if (pb.linear_solver_name == "mumps")
     pb.mumps_compute_electric_potential (ray_cache);
  else if (pb.linear_solver_name == "lis")
     pb.lis_compute_electric_potential (ray_cache);
  else 
    {
       std::cerr << "Invalid linear solver selected" << std::endl;
       return 1;
    }
  TOC ("compute potential");
  
  TIC ();
    pb.energy(ray_cache);
  TOC ("compute energy")
  
  TIC ();
  pb.export_tmesh (ray_cache);
  TOC ("export tmesh");

  TIC ();
  pb.export_marked_tmesh ();
  TOC ("export marked tmesh");

  // pb.analitic_potential();

  if (rank == 0) 
    { 
      print_timing_report(); 
      
      if (pb.surf_type != 2)
        {
          //Save ray_cache:
          for (int i = 0; i < 3; ++i)
          { 
            // std::cout << "Direzione: "<< i <<std::endl;
            // std::cout << std::endl;
            
            nlohmann::json j;
            save_ray_cache (j, ray_cache.rays[i]);
            
            // std::cout << "Count cached rays: " << ray_cache.count_cache_dir[i] << std::endl;
            // std::cout << "Count new rays: " << ray_cache.count_new_dir[i] << std::endl;
            // std::cout << std::endl;
            std::ofstream ray_cached_file;
            std::string filename = "ray_cache_";
            std::string extension = ".json";
            filename += std::to_string(i);
            filename += extension;
            ray_cached_file.open (filename.c_str ());
    
            if (ray_cached_file.is_open ())
              ray_cached_file << j;
      
            ray_cached_file.close ();
            
            //Alternative way to save results:
            //print_map (ray_cache.rays); 
          }
        }
    }
  
  MPI_Barrier (mpicomm);
     
  MPI_Finalize ();
  
  return 0;
  
}


// void
// print_map(const std::map<std::array<double, 2> , crossings_t, map_compare>& r) 
// {
//   std::ofstream ray_cached_file;
//   ray_cached_file.open ("ray_cache_2.txt");
//   if (ray_cached_file.is_open ())
//   {
//     ray_cached_file << "Count cached rays: " << ray_cache_t::count_cache << std::endl;
//     ray_cached_file << "Count new rays: " << ray_cache_t::count_new << std::endl;
//     for (auto it : r)
//     {
//       ray_cached_file << "[[" << it.first.at(0) << ", " << it.first.at(1) << "]";
//       for (int i = 0; i < it.second.inters.size (); i++) 
//         ray_cached_file << ", " << it.second.inters[i].first;
//       ray_cached_file << "]" << std::endl;
//     }
//   }
//   ray_cached_file.close (); 
// }

void
save_ray_cache (nlohmann::json& j, const std::map<std::array<double, 2> , crossings_t, map_compare>& r)
{
  std::vector<double> ints;
  for (auto it : r)
  {
    ints.clear();
    
    j += nlohmann::json{{"ray", it.first}, {"inters", (it.second.inters)}, {"normals", (it.second.normals)}};
    
    // if (it.second.inters.size()>0)
    // {
    //   for (int i = 0; i < it.second.inters.size(); ++i)
    //   {
    //     std::cout<< std::hypot(it.first[0],it.first[1], it.second.inters[i]) << "  ";
    //   }
    //   std::cout<<std::endl;
    // }
  }
}

