#ifndef RAYTRACER_H
#define RAYTRACER_H

#include <map>
#include <set>
#include <vector>
#include <array>
#include <mpi.h>
#include <nanoshaper.h>

#include "json.hpp" 
#include "serialize.h"

using int_coord_t = unsigned long long int;

struct
crossings_t 
{

   static double start[3], end[3]; 

   unsigned dir = 1;
   bool init = 0; //1 when the ray is initilized
   double point[2]; //point x and z coords: the ones that prescribe the ray
 
   //std::vector<bool> flags; //in or out (in=1, out=0) ex: if flags[1]==1 -> [inters[1].first ; inters[2].first] is IN  
   std::vector<double> inters; //intersections
   std::vector<double> normals; //normals

};

struct
map_compare
{
  bool 
  operator ()(const std::array<double, 2> & arr1, const std::array<double, 2> & arr2) const;
};

struct
ray_cache_t 
{  
   
   //static double start[3], end[3]; //start and end point of the analytical surface
   static unsigned dir; //direction of the ray
   
   std::map<std::array<double, 2>, crossings_t, map_compare> xdir_rays; //y-coord and z-coord, correspondent ray
   std::map<std::array<double, 2>, crossings_t, map_compare> ydir_rays; //x-coord and z-coord, correspondent ray
   std::map<std::array<double, 2>, crossings_t, map_compare> zdir_rays; //x-coord and y-coord, correspondent ray

   static int_coord_t count_cache;
   static int_coord_t count_new;
   
   NS::NanoShaper ns;
   
   int num_req_rays = 0;
   std::set<std::array<double, 2>, map_compare> xdir_rays_list; //list of req rays
   std::set<std::array<double, 2>, map_compare> ydir_rays_list; //list of req rays
   std::set<std::array<double, 2>, map_compare> zdir_rays_list; //list of req rays

   crossings_t &
   operator() (double x0, double x1, unsigned dir = 1);
   
   void 
   fill_cache ();
   
   void 
   init_analytical_surf (const std::vector<NS::Atom> & atoms, const NS::surface_type & surf_type,
			  const double & surf_param, const double & stern_layer, const unsigned & num_threads);
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
