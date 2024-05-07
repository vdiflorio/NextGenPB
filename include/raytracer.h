#ifndef RAYTRACER_H
#define RAYTRACER_H

#include <map>
#include <set>
#include <vector>
#include <array>
#include <mpi.h>
#include <nanoshaper.h>
#include <raytracer_datatype.h>

#include "json.hpp"
#include "serialize.h"

using int_coord_t = unsigned long long int;

// struct
//   crossings_t {

//   static double start[3], end[3];

//   unsigned dir;
//   bool init = 0; //1 when the ray is initilized
//   double point[2]; //point x and z coords: the ones that prescribe the ray

//   std::vector<double> inters; //intersections
//   std::vector<double> normals; //normals

// };

// struct
//   map_compare {
//   bool
//   operator () (const std::array<double, 2> & arr1, const std::array<double, 2> & arr2) const;
// };

struct
  ray_cache_t {

  std::array<std::map<std::array<double, 2>, crossings_t, map_compare>, 3> rays; //map that contains all the rays in the 3 direction

  static int_coord_t count_cache;
  static int_coord_t count_new;
  int count_new_dir[3] = {0, 0, 0};
  int count_cache_dir[3] = {0, 0, 0};
  // unsigned dir; //direction of the ray
  NS::NanoShaper ns;

  // int num_req_rays = 0;
  int num_req_rays[3] = {0, 0, 0};
  std::array<std::set<std::array<double, 2>, map_compare>, 3> rays_list; //list of req rays

  crossings_t &
  operator() (double x0, double x1, unsigned dir = 1);

  void
  fill_cache ();

  void
  init_analytical_surf (const std::vector<NS::Atom> & atoms, const NS::surface_type & surf_type,
                        const double & surf_param, const double & stern_layer, const unsigned & num_threads,const std::string* configFile=nullptr);

  void
  init_analytical_surf_ns (const std::vector<NS::Atom> & atoms, const NS::surface_type & surf_type,
                          const double & surf_param, const double & stern_layer, const unsigned & num_threads,const std::string* configFile=nullptr);
  void
  compute_ns_inters (crossings_t & ct);

  std::vector<unsigned char>
  write_map (const std::map<std::array<double, 2>, crossings_t, map_compare>& container);

  void
  read_map (const std::vector<unsigned char>& data,
            std::map<std::array<double, 2>, crossings_t, map_compare>& container);

  std::vector<unsigned char>
  write_ct (const crossings_t& ct);

  void
  read_ct (const std::vector<unsigned char>& data, crossings_t& ct);
};

#endif //RAYTRACER_H
