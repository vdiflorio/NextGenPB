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

#include "raytracer.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <unordered_set>
#include <tmesh_3d.h>

namespace
{
int64_t vertex_index (int i, int j, int k, const NS::PBAlignedGridDescriptor& grid)
{
  return static_cast<int64_t> (i) +
         static_cast<int64_t> (grid.ns_n[0]) *
         (static_cast<int64_t> (j) + static_cast<int64_t> (grid.ns_n[1]) * k);
}

int64_t edge_count (unsigned axis, const NS::PBAlignedGridDescriptor& grid)
{
  if (axis == 0)
    return static_cast<int64_t> (std::max (0, grid.ns_n[0] - 1)) * grid.ns_n[1] * grid.ns_n[2];

  if (axis == 1)
    return static_cast<int64_t> (grid.ns_n[0]) * std::max (0, grid.ns_n[1] - 1) * grid.ns_n[2];

  return static_cast<int64_t> (grid.ns_n[0]) * grid.ns_n[1] * std::max (0, grid.ns_n[2] - 1);
}

int64_t edge_index (unsigned axis, int i, int j, int k, const NS::PBAlignedGridDescriptor& grid)
{
  if (axis == 0)
    return static_cast<int64_t> (i) +
           static_cast<int64_t> (grid.ns_n[0] - 1) *
           (static_cast<int64_t> (j) + static_cast<int64_t> (grid.ns_n[1]) * k);

  if (axis == 1)
    return static_cast<int64_t> (i) +
           static_cast<int64_t> (grid.ns_n[0]) *
           (static_cast<int64_t> (j) + static_cast<int64_t> (grid.ns_n[1] - 1) * k);

  return static_cast<int64_t> (i) +
         static_cast<int64_t> (grid.ns_n[0]) *
         (static_cast<int64_t> (j) + static_cast<int64_t> (grid.ns_n[1]) * k);
}

template <class T>
void broadcast_vector (std::vector<T>& data)
{
  uint64_t size = static_cast<uint64_t> (data.size ());
  MPI_Bcast (&size, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  data.resize (static_cast<size_t> (size));

  if (size != 0)
    MPI_Bcast (data.data (), static_cast<int> (size * sizeof (T)), MPI_BYTE, 0, MPI_COMM_WORLD);
}
}

crossings_t &
ray_cache_t::operator () (double x0, double x1, unsigned direct)
{

  int rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  static crossings_t cr_t;

  static std::array<double, 2> start_point;
  start_point = {x0, x1};

  auto it0 = rays[direct].find (start_point);

  if (it0 != rays[direct].end () && it0->second.init == 1) {
    return rays[direct][start_point];
  }

  else if (it0 != rays[direct].end () && it0->second.init == 0) {
    count_new++;
    count_new_dir[direct]++;

    if (rank == 0)
      compute_ns_inters (it0->second);
  } else {
    // create a new ray
    crossings_t tmp;
    tmp.point[0] = x0;
    tmp.point[1] = x1;
    tmp.dir = direct;

    if (rank == 0)
      compute_ns_inters (tmp);

    rays[direct][start_point] = std::move (tmp);
    cr_t = rays[direct][start_point];
    count_new++;
    count_new_dir[direct]++;
  }

  return rays[direct][start_point];
}

void
ray_cache_t::fill_cache ()
{
  if (aligned_mode)
    return;

  int rank, size;
  MPI_Comm mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  for (unsigned idir = 0; idir < 3; ++idir) {
    std::vector<std::array<double, 2>> rays_vector (rays_list[idir].begin (), rays_list[idir].end ()); //copy the set content inside a vector

    num_req_rays[idir] = rays_vector.size (); //numb of req rays from each proc
    MPI_Barrier (mpicomm);
    // std::cout << "Sending " << num_req_rays[idir] << " rays requested from rank " << rank << std::endl;

    std::vector<unsigned char> local_ser_rays_vec = serialize::write (rays_vector); //vector of char for the rays_vector
    std::vector<unsigned char> global_ser_rays_vec; //vector with all the rays from all proc

    int ser_rays_len = local_ser_rays_vec.size (); //length of each vect of char
    int global_ser_rays_len; //variable with total number of rays (sum all the proc rays)

    auto proc_ser_rays_len = std::make_unique<int[]> (size);
    auto proc_num_req_rays = std::make_unique<int[]> (size);
    auto displ_ser_rays_vec = std::make_unique<int[]> (size);
    auto cum_rays_vec = std::make_unique<int[]> (size);


    std::vector<unsigned char> global_ser_map; //vector of char for the map of rank 0
    std::vector<unsigned char> local_ser_map; //vector of char for the local map to send to each rank

    int global_ser_map_len;
    int local_ser_map_len;

    auto proc_ser_map_len = std::make_unique<int[]> (size);
    auto displ_ser_map = std::make_unique<int[]> (size);


    MPI_Reduce (&ser_rays_len, &global_ser_rays_len, 1, MPI_INT, MPI_SUM, 0, mpicomm); //fill global_ser_rays_len

    MPI_Gather (&ser_rays_len, 1, MPI_INT, proc_ser_rays_len.get (), 1, MPI_INT, 0, mpicomm); //fill proc_ser_rays_len
    MPI_Gather (&num_req_rays[idir], 1, MPI_INT, proc_num_req_rays.get (), 1, MPI_INT, 0, mpicomm); //fill proc_num_req_rays

    if (rank == 0) {
      for (int i = 1; i < size; i++) {
        displ_ser_rays_vec[i] = displ_ser_rays_vec[i-1] + proc_ser_rays_len[i-1];
        cum_rays_vec[i] = cum_rays_vec[i-1] + proc_num_req_rays[i];
      }

      global_ser_rays_vec.resize (global_ser_rays_len);
    }

    // Send the requested rays to rank 0
    MPI_Gatherv (&local_ser_rays_vec[0], ser_rays_len, MPI_CHAR, &global_ser_rays_vec[0],
                 proc_ser_rays_len.get (), displ_ser_rays_vec.get (), MPI_CHAR, 0, mpicomm);

    std::map<std::array<double, 2>, crossings_t, map_compare> local_req_rays_map;

    if (rank == 0) {
      // Transform the char rays in vector<array<double,2>>
      serialize::read (global_ser_rays_vec, rays_vector);
      //std::copy (rays_vector.begin (), rays_vector.end (),
      //std::inserter(rays_list[idir], rays_list[idir].begin ())); //copy the req rays in a set, to avoid repetitions in for loop

      // Add the rays to rank 0 map:
      //for (auto it = rays_list[idir].begin (); it != rays_list[idir].end (); it ++)
      //crossings_t & ct = (*this)((*it)[0], (*it)[1]);

      for (int i = 1; i < size; i++) {
        local_req_rays_map.clear ();

        std::transform (rays_vector.begin () + cum_rays_vec[i-1],
                        rays_vector.begin () + cum_rays_vec[i],
                        std::inserter (local_req_rays_map, end (local_req_rays_map)),
        [this,idir] (std::array<double,2> arr) {
          return (std::pair<std::array<double,2>, crossings_t> (arr, (*this) (arr[0], arr[1],idir)));
        }); //return ((this->rays).find(arr));

        local_ser_map = write_map (local_req_rays_map);
        proc_ser_map_len[i] = local_ser_map.size ();
        global_ser_map.insert (global_ser_map.end (), local_ser_map.begin (), local_ser_map.end ());
        displ_ser_map[i] = displ_ser_map[i-1] + proc_ser_map_len[i-1];
      }

      global_ser_map_len = global_ser_map.size ();

    }

    MPI_Scatter (proc_ser_map_len.get (), 1, MPI_INT, &local_ser_map_len, 1, MPI_INT, 0, mpicomm);

    local_ser_map.resize (local_ser_map_len);
    MPI_Scatterv (&global_ser_map[0], proc_ser_map_len.get (), displ_ser_map.get (), MPI_CHAR,
                  &local_ser_map[0], local_ser_map_len, MPI_CHAR, 0, mpicomm);

    if (rank != 0) {
      local_req_rays_map.clear ();
      read_map (local_ser_map, local_req_rays_map);
      (this->rays)[idir].insert (local_req_rays_map.begin (), local_req_rays_map.end ());
    }

    //MPI_Barrier (mpicomm);
    // std::cout << "Rays created in rank " << rank << std::endl;
  }
}

/*
void
ray_cache_t::init_analytical_surf (const std::vector<NS::Atom> & atoms, const NS::surface_type & surf_type,
                                   const double & surf_param, const double & stern_layer, const unsigned & num_threads,
                                   const std::string* configFile)
{
  // NS::NanoShaper ns0 (atoms, surf_type, surf_param, stern_layer, num_threads);
  // ns = ns0;
  ns.initConstructor (atoms, surf_type, surf_param, stern_layer, num_threads,configFile);

  ns.setConfig<double> ("Grid_scale", 2.0);
  ns.setConfig<double> ("Self_Intersections_Grid_Coefficient", 1.5);
  ns.buildAnalyticalSurface();
  std::cout << "\n" << std::endl;

  std::vector<unsigned> idx;
  std::vector<double> coords;

  // right-upper corner
  ns.getGridSize (idx);
  ns.getGridPointCoordinates (idx[0]-1, idx[1]-1, idx[2]-1, coords);
  std::copy (coords.cbegin(), coords.cend(), crossings_t::end);

  // left-bottom corner
  ns.getGridPointCoordinates (0, 0, 0, coords);
  std::copy (coords.cbegin(), coords.cend(), crossings_t::start);
}


void
ray_cache_t::init_analytical_surf_ns (const std::vector<NS::Atom> & atoms, const NS::surface_type & surf_type,
                                      const double & surf_param, const double & stern_layer, const unsigned & num_threads,
                                      double* l_cr, double* r_cr, double scale, const std::string* configFile)
{

  ns.initConstructor (atoms, surf_type, surf_param, stern_layer, num_threads,configFile);
  // set here a consistent grid scale
  ns.setConfig<double> ("Grid_scale", scale );
  ns.setConfig<bool> ("Accurate_Triangulation",true);
  ns.setConfig<double> ("Self_Intersections_Grid_Coefficient", 1.5);

  ns.setConfig<bool> ("Build_epsilon_maps",true);

  // build the grid in the new mode
  ns.setConfig<bool> ("PB_grid_mode",true);

  // impose here the min and max of the cube
  ns.setConfig<double> ("xmin",l_cr[0]);
  ns.setConfig<double> ("ymin",l_cr[1]);
  ns.setConfig<double> ("zmin",l_cr[2]);

  ns.setConfig<double> ("xmax",r_cr[0]);
  ns.setConfig<double> ("ymax",r_cr[1]);
  ns.setConfig<double> ("zmax",r_cr[2]);
  ns.buildAnalyticalSurface();

  // remember to set this to true in order to collect rays intersections data
  ns.setCollectGridRays (true);
  ns.colourGrid();
  // retrieve intersections data from the map
  rays = ns.getRaysMap();

  // crossings_t::start[0] = l_cr[0];
  // crossings_t::start[1] = l_cr[1];
  // crossings_t::start[2] = l_cr[2];
  // crossings_t::end[0] = r_cr[0];
  // crossings_t::end[1] = r_cr[1];
  // crossings_t::end[2] = r_cr[2];

  l_c[0] = l_cr[0];
  l_c[1] = l_cr[1];
  l_c[2] = l_cr[2];
  r_c[0] = r_cr[0];
  r_c[1] = r_cr[1];
  r_c[2] = r_cr[2];
}


void
ray_cache_t::compute_ns_inters (crossings_t & ct)
{

  if (ct.dir == 0) {
    double start_ray[3] = {l_c[ct.dir], ct.point[0], ct.point[1]};
    bool compute_normals = true;


    if (ct.point[0] < l_c[1] || ct.point[1] < l_c[2] || ct.point[0] > r_c[1] || ct.point[1] > r_c[2]) { //out of the molecule
      ct.init = 1;
      return;
    }

    std::cout << "Sending new ray in x direction for NS!" << std::endl;
    ns.setDirection (ct.dir);
    std::vector<std::pair<double,double*>> ints_norms; //intersections and normals
    ns.castAxisOrientedRay (start_ray, r_c[ct.dir], ints_norms, ct.dir, compute_normals);

    if (ints_norms.size() != 0) {
      for (unsigned i = 0; i < ints_norms.size(); i++) {
        ct.inters.push_back (ints_norms[i].first);

        ct.normals.push_back (* (ints_norms[i].second));
        ct.normals.push_back (* (ints_norms[i].second+1));
        ct.normals.push_back (* (ints_norms[i].second+2));
      }
    }

    ct.init = 1; //ray is now initialized
  }

  if (ct.dir == 1) {
    double start_ray[3] = {ct.point[0], l_c[ct.dir], ct.point[1]};
    bool compute_normals = true;


    if (ct.point[0] < l_c[0] || ct.point[1] < l_c[2] || ct.point[0] > r_c[0] || ct.point[1] > r_c[2]) { //out of the molecule
      ct.init = 1;
      return;
    }
    std::cout << "Sending new ray in y direction for NS!" << std::endl;
    ns.setDirection (ct.dir);
    std::vector<std::pair<double,double*>> ints_norms; //intersections and normals
    ns.castAxisOrientedRay (start_ray, r_c[ct.dir], ints_norms, ct.dir, compute_normals);

    if (ints_norms.size() != 0) {
      for (unsigned i = 0; i < ints_norms.size(); i++) {
        ct.inters.push_back (ints_norms[i].first);

        ct.normals.push_back (* (ints_norms[i].second));
        ct.normals.push_back (* (ints_norms[i].second+1));
        ct.normals.push_back (* (ints_norms[i].second+2));
      }
    }

    ct.init = 1; //ray is now initialized
  }

  if (ct.dir == 2) {
    double start_ray[3] = {ct.point[0], ct.point[1], l_c[ct.dir]};
    bool compute_normals = true;


    if (ct.point[0] < l_c[0] || ct.point[1] < l_c[1] || ct.point[0] > r_c[0] || ct.point[1] > r_c[1]) { //out of the molecule
      ct.init = 1;
      return;
    }
    std::cout << "Sending new ray in z direction for NS!" << std::endl;
    ns.setDirection (ct.dir);
    std::vector<std::pair<double,double*>> ints_norms; //intersections and normals
    ns.castAxisOrientedRay (start_ray, r_c[ct.dir], ints_norms, ct.dir, compute_normals);

    if (ints_norms.size() != 0) {
      for (unsigned i = 0; i < ints_norms.size(); i++) {
        ct.inters.push_back (ints_norms[i].first);

        ct.normals.push_back (* (ints_norms[i].second));
        ct.normals.push_back (* (ints_norms[i].second+1));
        ct.normals.push_back (* (ints_norms[i].second+2));
      }
    }

    ct.init = 1; //ray is now initialized
  }

}
*/


void
ray_cache_t::init_analytical_surf_ns (const std::vector<NS::Atom> & atoms, const NS::surface_type & surf_type,
                                      const double & surf_param, const double & stern_layer, const unsigned & num_threads,
                                      double* l_cr, double* r_cr, double scale, const std::string* configFile)
{
  ns = std::make_unique<NS::NanoShaper> (atoms, surf_type, surf_param, stern_layer, num_threads,configFile);
  // set here a consistent grid scale
  ns->setConfig<double> ("Grid_scale", scale);
  ns->setConfig<bool> ("Accurate_Triangulation",true);
  ns->setConfig<double> ("Self_Intersections_Grid_Coefficient", 1.5);

  ns->setConfig<bool> ("Build_epsilon_maps",true);

  // build the grid in the new mode
  ns->setConfig<bool> ("PB_grid_mode",true);

  // impose here the min and max of the cube
  ns->setConfig<double> ("xmin",l_cr[0]);
  ns->setConfig<double> ("ymin",l_cr[1]);
  ns->setConfig<double> ("zmin",l_cr[2]);

  ns->setConfig<double> ("xmax",r_cr[0]);
  ns->setConfig<double> ("ymax",r_cr[1]);
  ns->setConfig<double> ("zmax",r_cr[2]);
  ns->buildAnalyticalSurface ();

  // remember to set this to true in order to collect rays intersections data
  ns->setCollectGridRays (true);
  ns->colourGrid ();
  // retrieve intersections data from the map

  rays = ns->getRaysMap ();

  l_c[0] = l_cr[0];
  l_c[1] = l_cr[1];
  l_c[2] = l_cr[2];
  r_c[0] = r_cr[0];
  r_c[1] = r_cr[1];
  r_c[2] = r_cr[2];
}

void
ray_cache_t::init_aligned_surface_ns (const std::vector<NS::Atom> & atoms, const NS::surface_type & surf_type,
                                      const double & surf_param, const double & stern_layer, const unsigned & num_threads,
                                      double* l_cr, double* r_cr, double scale, const std::string* configFile,
                                      bool do_triangulate, bool cavity_fill, double cavity_vol)
{
  ns = std::make_unique<NS::NanoShaper> (atoms, surf_type, surf_param, stern_layer, num_threads, configFile);
  ns->setConfig<double> ("Grid_scale", scale);
  ns->setConfig<bool> ("Accurate_Triangulation", true);
  ns->setConfig<double> ("Self_Intersections_Grid_Coefficient", 1.5);
  ns->setConfig<bool> ("Build_epsilon_maps", true);
  ns->setConfig<bool> ("PB_grid_mode", true);
  ns->setConfig<bool> ("PB_aligned_grid_mode", true);

  // Internal-cavity detection/filling. When enabled, enclosed voids with volume
  // <= cavity_vol [A^3] are filled in verticesInsidenessMap, so the exported
  // cavity-filled vertex_colors mark them inside (eps_in, no mobile ions).
  ns->setConfig<bool> ("Cavity_Detection_Filling", cavity_fill);
  ns->setConfig<double> ("Conditional_Volume_Filling_Value", cavity_vol);

  // PB-aligned (Approccio B): shift the NS input box up by +hside so the
  // marching-cubes triangulation lattice (which sits at DS->x - hside) lands
  // exactly on the PB/FEM nodes l_cr + i*h. colourGrid exports vertex_colors and
  // edge_crossings keyed on that staggered lattice (ns_origin = DS->x - hside =
  // l_cr), so vertex_color()/aligned_edge_crossings() formulas are unchanged.
  const double hside = 0.5 / scale;
  ns->setConfig<double> ("xmin", l_cr[0] + hside);
  ns->setConfig<double> ("ymin", l_cr[1] + hside);
  ns->setConfig<double> ("zmin", l_cr[2] + hside);
  ns->setConfig<double> ("xmax", r_cr[0] + hside);
  ns->setConfig<double> ("ymax", r_cr[1] + hside);
  ns->setConfig<double> ("zmax", r_cr[2] + hside);

  // Smoothing displaces vertList off the analytical crossing positions, which
  // would break the OFF <-> edge_crossings match (both must share vertList).
  ns->setConfig<bool> ("Smooth_Mesh", false);

  if (do_triangulate) {
    ns->setConfig<bool> ("Vertex_Atom_Info", true);
    ns->setConfig<bool> ("Compute_Vertex_Normals", true);
  }

  if (!ns->buildAnalyticalSurface ()) {
    std::cerr << "NanoShaper failed to build the analytical surface for PB aligned mode." << std::endl;
    aligned_mode = false;
    return;
  }

  ns->setCollectGridRays (true);
  if (!ns->colourGrid ()) {
    std::cerr << "NanoShaper failed to colour the grid for PB aligned mode." << std::endl;
    aligned_mode = false;
    return;
  }

  aligned_surface = ns->getPBAlignedSurfaceData ();
  aligned_mode = aligned_surface.valid;
  if (!aligned_mode)
    std::cerr << "NanoShaper did not export PB aligned surface data." << std::endl;

  if (do_triangulate) {
    double surf_area = 0.0;
    if (ns->triangulate (&surf_area))
      std::cout << "\nNanoShaper SES area: " << surf_area << " A^2\n";
    else
      std::cerr << "NanoShaper triangulation failed.\n";
  }

  l_c[0] = l_cr[0];
  l_c[1] = l_cr[1];
  l_c[2] = l_cr[2];
  r_c[0] = r_cr[0];
  r_c[1] = r_cr[1];
  r_c[2] = r_cr[2];
}

void
ray_cache_t::broadcast_aligned_surface ()
{
  // Only called when refine_surface will run before scatter_aligned_surface.
  // Broadcasts aligned_mode flag and vertex_colors so all ranks can answer
  // vertex_color() queries during parallel mesh refinement.
  // The grid descriptor and edge data are handled in scatter_aligned_surface.
  int enabled = aligned_mode ? 1 : 0;
  MPI_Bcast (&enabled, 1, MPI_INT, 0, MPI_COMM_WORLD);
  aligned_mode = enabled != 0;

  if (!aligned_mode)
    return;

  broadcast_vector (aligned_surface.vertex_colors);
  aligned_surface.valid = true;
}

void
ray_cache_t::scatter_aligned_surface (tmesh_3d& tmsh)
{
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  // With a single rank the full aligned_surface arrays are already local.
  if (size == 1)
    return;

  // Broadcast aligned_mode flag and grid descriptor to all ranks.
  // If broadcast_aligned_surface was called earlier (refine path), the
  // re-broadcast of grid is harmless (< 100 bytes).  If it was not called
  // (no-refine path), this is the first time non-zero ranks receive the data.
  {
    int enabled = aligned_mode ? 1 : 0;
    MPI_Bcast (&enabled, 1, MPI_INT, 0, MPI_COMM_WORLD);
    aligned_mode = enabled != 0;
  }
  if (!aligned_mode)
    return;

  MPI_Bcast (&aligned_surface.grid, static_cast<int> (sizeof (NS::PBAlignedGridDescriptor)),
             MPI_BYTE, 0, MPI_COMM_WORLD);

  // Non-zero ranks no longer need the broadcast vertex_colors (used only for
  // refine_surface).  Free them now, before scatter allocates working arrays,
  // to minimise peak RSS.  Rank 0 keeps the full arrays to serve scatter requests.
  // In the no-refine path vertex_colors was never broadcast, so this is a no-op.
  if (rank != 0)
    std::vector<uint8_t> ().swap (aligned_surface.vertex_colors);

  const NS::PBAlignedGridDescriptor& grid = aligned_surface.grid;

  // Edge topology matching pb_class.h edge_axis / edge2nodes
  static constexpr std::array<int, 12> edge_ax = {0,1,0,1,0,1,0,1,2,2,2,2};
  static constexpr std::array<int, 24> e2n = {
    0,1, 1,3, 2,3, 0,2, 4,5, 5,7, 6,7, 4,6, 0,4, 1,5, 3,7, 2,6
  };

  // Collect vertex and edge indices needed by this rank's local mesh
  std::unordered_set<int64_t> vset;
  std::array<std::unordered_set<int64_t>, 3> esets;

  for (auto q = tmsh.begin_quadrant_sweep (); q != tmsh.end_quadrant_sweep (); ++q) {
    for (int ii = 0; ii < 8; ++ii) {
      int idx[3]; bool ok = true;
      for (int d = 0; d < 3; ++d) {
        idx[d] = static_cast<int> (std::llround ((q->p (d, ii) - grid.ns_origin[d]) / grid.h));
        if (idx[d] < 0 || idx[d] >= grid.ns_n[d]) { ok = false; break; }
      }
      if (ok)
        vset.insert (vertex_index (idx[0], idx[1], idx[2], grid));
    }

    for (int j = 0; j < 12; ++j) {
      unsigned dir = static_cast<unsigned> (edge_ax[j]);
      int i1 = e2n[2 * j], i2 = e2n[2 * j + 1];
      double x1 = q->p (dir, i1), x2 = q->p (dir, i2);
      if (x2 < x1) std::swap (x1, x2);

      std::array<int, 2> fax;
      if (dir == 0) fax = {1, 2};
      else if (dir == 1) fax = {0, 2};
      else fax = {0, 1};

      int fixed[2]; bool ok = true;
      for (int q2 = 0; q2 < 2; ++q2) {
        fixed[q2] = static_cast<int> (std::llround ((q->p (fax[q2], i1) - grid.ns_origin[fax[q2]]) / grid.h));
        if (fixed[q2] < 0 || fixed[q2] >= grid.ns_n[fax[q2]]) { ok = false; break; }
      }
      if (!ok) continue;

      int a0 = static_cast<int> (std::llround ((x1 - grid.ns_origin[dir]) / grid.h));
      int a1 = static_cast<int> (std::llround ((x2 - grid.ns_origin[dir]) / grid.h));
      a0 = std::max (0, std::min (a0, grid.ns_n[dir] - 1));
      a1 = std::max (a0, std::min (a1, grid.ns_n[dir] - 1));

      for (int a = a0; a < a1; ++a) {
        int ci, cj, ck;
        if (dir == 0) { ci = a; cj = fixed[0]; ck = fixed[1]; }
        else if (dir == 1) { ci = fixed[0]; cj = a; ck = fixed[1]; }
        else { ci = fixed[0]; cj = fixed[1]; ck = a; }
        int64_t eid = edge_index (dir, ci, cj, ck, grid);
        if (eid >= 0 && eid < edge_count (dir, grid))
          esets[dir].insert (eid);
      }
    }
  }

  std::vector<int64_t> my_vlist (vset.begin (), vset.end ());
  std::array<std::vector<int64_t>, 3> my_elist;
  for (int d = 0; d < 3; ++d)
    my_elist[d].assign (esets[d].begin (), esets[d].end ());

  // Rank 0 builds its own local maps directly from the full arrays
  if (rank == 0) {
    for (int64_t vid : my_vlist)
      local_vertex_colors[vid] = aligned_surface.vertex_colors[static_cast<size_t> (vid)];
    for (unsigned d = 0; d < 3; ++d) {
      for (int64_t eid : my_elist[d]) {
        uint32_t beg = aligned_surface.edge_offsets[d][static_cast<size_t> (eid)];
        uint32_t end2 = aligned_surface.edge_offsets[d][static_cast<size_t> (eid) + 1];
        auto& vec = local_edge_crossings[d][eid];
        vec.assign (aligned_surface.edge_crossings[d].begin () + beg,
                    aligned_surface.edge_crossings[d].begin () + end2);
      }
    }
  }

  if (size > 1) {
    // ---- Vertex colour exchange ----
    // Rank 0 participates with empty send (count=0); it already built local maps above.
    {
      int64_t my_send_count = (rank != 0) ? static_cast<int64_t> (my_vlist.size ()) : 0;

      std::vector<int64_t> all_counts;
      if (rank == 0) all_counts.resize (size);
      MPI_Gather (&my_send_count, 1, MPI_INT64_T,
                  rank == 0 ? all_counts.data () : nullptr,
                  1, MPI_INT64_T, 0, MPI_COMM_WORLD);

      std::vector<int> vcounts, vdispls;
      std::vector<int64_t> global_vlist;
      int total_vreq = 0;
      if (rank == 0) {
        vcounts.resize (size);
        vdispls.resize (size, 0);
        for (int r = 0; r < size; ++r) vcounts[r] = static_cast<int> (all_counts[r]);
        for (int r = 1; r < size; ++r) vdispls[r] = vdispls[r-1] + vcounts[r-1];
        total_vreq = vdispls[size-1] + vcounts[size-1];
        global_vlist.resize (total_vreq);
      }

      if (rank == 0)
        MPI_Gatherv (MPI_IN_PLACE, 0, MPI_INT64_T,
                     global_vlist.data (), vcounts.data (), vdispls.data (),
                     MPI_INT64_T, 0, MPI_COMM_WORLD);
      else
        MPI_Gatherv (my_vlist.data (), static_cast<int> (my_vlist.size ()), MPI_INT64_T,
                     nullptr, nullptr, nullptr,
                     MPI_INT64_T, 0, MPI_COMM_WORLD);

      // Rank 0 looks up colours
      std::vector<uint8_t> global_vcolors;
      if (rank == 0) {
        global_vcolors.resize (total_vreq);
        for (int i = 0; i < total_vreq; ++i) {
          int64_t vid = global_vlist[i];
          global_vcolors[i] = (vid >= 0 && static_cast<size_t> (vid) < aligned_surface.vertex_colors.size ())
                              ? aligned_surface.vertex_colors[static_cast<size_t> (vid)] : 0;
        }
      }

      std::vector<uint8_t> my_vcolors;
      if (rank != 0) my_vcolors.resize (my_vlist.size ());

      if (rank == 0)
        MPI_Scatterv (global_vcolors.data (), vcounts.data (), vdispls.data (),
                      MPI_BYTE, MPI_IN_PLACE, 0, MPI_BYTE, 0, MPI_COMM_WORLD);
      else
        MPI_Scatterv (nullptr, nullptr, nullptr, MPI_BYTE,
                      my_vcolors.data (), static_cast<int> (my_vcolors.size ()), MPI_BYTE,
                      0, MPI_COMM_WORLD);

      if (rank != 0)
        for (size_t i = 0; i < my_vlist.size (); ++i)
          local_vertex_colors[my_vlist[i]] = my_vcolors[i];
    }

    // ---- Edge crossings exchange (per direction) ----
    for (unsigned dir = 0; dir < 3; ++dir) {
      auto& my_edir = my_elist[dir];

      int64_t my_send_count = (rank != 0) ? static_cast<int64_t> (my_edir.size ()) : 0;

      std::vector<int64_t> all_counts;
      if (rank == 0) all_counts.resize (size);
      MPI_Gather (&my_send_count, 1, MPI_INT64_T,
                  rank == 0 ? all_counts.data () : nullptr,
                  1, MPI_INT64_T, 0, MPI_COMM_WORLD);

      std::vector<int> ecounts, edispls;
      std::vector<int64_t> global_elist;
      int total_ereq = 0;
      if (rank == 0) {
        ecounts.resize (size);
        edispls.resize (size, 0);
        for (int r = 0; r < size; ++r) ecounts[r] = static_cast<int> (all_counts[r]);
        for (int r = 1; r < size; ++r) edispls[r] = edispls[r-1] + ecounts[r-1];
        total_ereq = edispls[size-1] + ecounts[size-1];
        global_elist.resize (total_ereq);
      }

      if (rank == 0)
        MPI_Gatherv (MPI_IN_PLACE, 0, MPI_INT64_T,
                     global_elist.data (), ecounts.data (), edispls.data (),
                     MPI_INT64_T, 0, MPI_COMM_WORLD);
      else
        MPI_Gatherv (my_edir.data (), static_cast<int> (my_edir.size ()), MPI_INT64_T,
                     nullptr, nullptr, nullptr,
                     MPI_INT64_T, 0, MPI_COMM_WORLD);

      // Rank 0: compute per-edge crossing counts and pack flat crossings
      std::vector<uint32_t> global_ec_counts;
      std::vector<NS::PBEdgeCrossing> global_ec_flat;
      std::vector<int> rank_xing_counts, rank_xing_displs;  // in bytes

      if (rank == 0) {
        global_ec_counts.resize (total_ereq, 0);
        for (int i = 0; i < total_ereq; ++i) {
          int64_t eid = global_elist[i];
          if (eid >= 0 && eid < edge_count (dir, grid)) {
            uint32_t beg  = aligned_surface.edge_offsets[dir][static_cast<size_t> (eid)];
            uint32_t end2 = aligned_surface.edge_offsets[dir][static_cast<size_t> (eid) + 1];
            global_ec_counts[i] = end2 - beg;
            for (uint32_t p = beg; p < end2; ++p)
              global_ec_flat.push_back (aligned_surface.edge_crossings[dir][p]);
          }
        }

        rank_xing_counts.resize (size, 0);
        rank_xing_displs.resize (size, 0);
        int flat_pos = 0;
        for (int r = 0; r < size; ++r) {
          rank_xing_displs[r] = flat_pos;
          int n_structs = 0;
          for (int e = edispls[r]; e < edispls[r] + ecounts[r]; ++e)
            n_structs += static_cast<int> (global_ec_counts[e]);
          rank_xing_counts[r] = n_structs * static_cast<int> (sizeof (NS::PBEdgeCrossing));
          flat_pos += rank_xing_counts[r];
        }
      }

      // Scatter per-edge crossing counts to non-zero ranks
      std::vector<uint32_t> my_ec_counts;
      if (rank != 0) my_ec_counts.resize (my_edir.size ());

      if (rank == 0)
        MPI_Scatterv (global_ec_counts.data (), ecounts.data (), edispls.data (),
                      MPI_UINT32_T, MPI_IN_PLACE, 0, MPI_UINT32_T, 0, MPI_COMM_WORLD);
      else
        MPI_Scatterv (nullptr, nullptr, nullptr, MPI_UINT32_T,
                      my_ec_counts.data (), static_cast<int> (my_ec_counts.size ()), MPI_UINT32_T,
                      0, MPI_COMM_WORLD);

      // Scatter actual crossings (as raw bytes)
      int my_total_xings = 0;
      if (rank != 0)
        for (auto c : my_ec_counts) my_total_xings += static_cast<int> (c);
      int my_recv_bytes = my_total_xings * static_cast<int> (sizeof (NS::PBEdgeCrossing));

      std::vector<NS::PBEdgeCrossing> my_crossings (my_total_xings);

      if (rank == 0)
        MPI_Scatterv (static_cast<void*> (global_ec_flat.data ()),
                      rank_xing_counts.data (), rank_xing_displs.data (),
                      MPI_BYTE, MPI_IN_PLACE, 0, MPI_BYTE, 0, MPI_COMM_WORLD);
      else
        MPI_Scatterv (nullptr, nullptr, nullptr, MPI_BYTE,
                      static_cast<void*> (my_crossings.data ()), my_recv_bytes, MPI_BYTE,
                      0, MPI_COMM_WORLD);

      if (rank != 0) {
        size_t offset = 0;
        for (size_t i = 0; i < my_edir.size (); ++i) {
          auto& vec = local_edge_crossings[dir][my_edir[i]];
          vec.resize (my_ec_counts[i]);
          std::copy (my_crossings.begin () + static_cast<ptrdiff_t> (offset),
                     my_crossings.begin () + static_cast<ptrdiff_t> (offset + my_ec_counts[i]),
                     vec.begin ());
          offset += my_ec_counts[i];
        }
      }
    }
  }

  use_local_maps = true;

  // Free full arrays on all ranks: non-zero ranks already freed vertex_colors
  // above; rank 0 frees everything here now that local maps are populated.
  std::vector<uint8_t> ().swap (aligned_surface.vertex_colors);
  for (int d = 0; d < 3; ++d) {
    std::vector<uint32_t> ().swap (aligned_surface.edge_offsets[d]);
    std::vector<NS::PBEdgeCrossing> ().swap (aligned_surface.edge_crossings[d]);
  }
  aligned_surface.valid = false;
}

int
ray_cache_t::vertex_color (double x, double y, double z) const
{
  const double coords[3] = {x, y, z};
  int idx[3] = {0, 0, 0};

  // Pass 1: locate the nearest NS node on every axis. A node outside the NS
  // lattice in ANY axis is definitely outside the molecule (return 0) and must
  // never reach the alignment check below — far-field octree nodes routinely
  // have one in-range coordinate while the others are far out of range, and the
  // alignment abort must not fire on those.
  for (int d = 0; d < 3; ++d) {
    const double scaled = (coords[d] - aligned_surface.grid.ns_origin[d]) / aligned_surface.grid.h;
    idx[d] = static_cast<int> (std::llround (scaled));

    if (idx[d] < 0 || idx[d] >= aligned_surface.grid.ns_n[d])
      return 0;
  }

  // Pass 2: only nodes fully inside the NS box must coincide with the lattice.
  for (int d = 0; d < 3; ++d) {
    const double snapped = aligned_surface.grid.ns_origin[d] + idx[d] * aligned_surface.grid.h;
    const double mismatch = std::fabs (coords[d] - snapped);
    if (mismatch > aligned_tol * aligned_surface.grid.h) {
      std::cerr << "PB node is not aligned with NanoShaper lattice: "
                << "(" << x << ", " << y << ", " << z << ") axis " << d
                << " mismatch " << mismatch << " > tol " << (aligned_tol * aligned_surface.grid.h)
                << " (h " << aligned_surface.grid.h << ")" << std::endl;
      MPI_Abort (MPI_COMM_WORLD, 1);
      return -1;
    }
  }

  const int64_t id = vertex_index (idx[0], idx[1], idx[2], aligned_surface.grid);

  if (use_local_maps) {
    auto it = local_vertex_colors.find (id);
    return (it != local_vertex_colors.end () && it->second) ? 1 : 0;
  }

  return aligned_surface.vertex_colors[static_cast<size_t> (id)] ? 1 : 0;
}

void
ray_cache_t::aligned_edge_crossings (unsigned dir, double x1, double x2,
                                     const std::array<double, 2>& ray,
                                     std::vector<NS::PBEdgeCrossing>& crossings) const
{
  crossings.clear ();

  if (dir > 2 || !aligned_mode)
    return;

  int fixed[2] = {0, 0};
  std::array<int, 2> fixed_axes;

  if (dir == 0)
    fixed_axes = {1, 2};
  else if (dir == 1)
    fixed_axes = {0, 2};
  else
    fixed_axes = {0, 1};

  // Pass 1: locate both fixed axes. If either is outside the NS lattice the edge
  // cannot cross the surface (return) and must not trip the alignment abort for
  // the other axis.
  for (int q = 0; q < 2; ++q) {
    const int axis = fixed_axes[q];
    const double scaled = (ray[q] - aligned_surface.grid.ns_origin[axis]) / aligned_surface.grid.h;
    fixed[q] = static_cast<int> (std::llround (scaled));

    if (fixed[q] < 0 || fixed[q] >= aligned_surface.grid.ns_n[axis])
      return;
  }

  // Pass 2: require the edge's fixed coordinates to sit on the NS lattice.
  for (int q = 0; q < 2; ++q) {
    const int axis = fixed_axes[q];
    const double snapped = aligned_surface.grid.ns_origin[axis] + fixed[q] * aligned_surface.grid.h;
    const double mismatch = std::fabs (ray[q] - snapped);
    if (mismatch > aligned_tol * aligned_surface.grid.h) {
      std::cerr << "PB edge is not aligned with NanoShaper lattice: axis " << axis
                << " coord " << ray[q] << " mismatch " << mismatch
                << " > tol " << (aligned_tol * aligned_surface.grid.h)
                << " (h " << aligned_surface.grid.h << ")" << std::endl;
      MPI_Abort (MPI_COMM_WORLD, 1);
      return;
    }
  }

  if (x2 < x1)
    std::swap (x1, x2);

  int a0 = static_cast<int> (std::llround ((x1 - aligned_surface.grid.ns_origin[dir]) / aligned_surface.grid.h));
  int a1 = static_cast<int> (std::llround ((x2 - aligned_surface.grid.ns_origin[dir]) / aligned_surface.grid.h));
  a0 = std::max (0, std::min (a0, aligned_surface.grid.ns_n[dir] - 1));
  a1 = std::max (a0, std::min (a1, aligned_surface.grid.ns_n[dir] - 1));

  for (int a = a0; a < a1; ++a) {
    int i = 0, j = 0, k = 0;

    if (dir == 0) {
      i = a;
      j = fixed[0];
      k = fixed[1];
    } else if (dir == 1) {
      i = fixed[0];
      j = a;
      k = fixed[1];
    } else {
      i = fixed[0];
      j = fixed[1];
      k = a;
    }

    const int64_t edge_id = edge_index (dir, i, j, k, aligned_surface.grid);
    if (edge_id < 0 || edge_id >= edge_count (dir, aligned_surface.grid))
      continue;

    if (use_local_maps) {
      auto it = local_edge_crossings[dir].find (edge_id);
      if (it != local_edge_crossings[dir].end ())
        for (const auto& c : it->second)
          crossings.push_back (c);
    } else {
      const auto begin = aligned_surface.edge_offsets[dir][static_cast<size_t> (edge_id)];
      const auto end   = aligned_surface.edge_offsets[dir][static_cast<size_t> (edge_id) + 1];
      for (uint32_t p = begin; p < end; ++p)
        crossings.push_back (aligned_surface.edge_crossings[dir][p]);
    }
  }
}


void
ray_cache_t::compute_ns_inters (crossings_t & ct)
{

  if (ct.dir == 0) {
    double start_ray[3] = {l_c[ct.dir], ct.point[0], ct.point[1]};
    bool compute_normals = true;


    if (ct.point[0] < l_c[1] || ct.point[1] < l_c[2] || ct.point[0] > r_c[1] || ct.point[1] > r_c[2]) { //out of the molecule
      ct.init = 1;
      return;
    }

    std::cout << "Sending new ray in x direction for NS!" << std::endl;
    ns->setDirection (ct.dir);
    std::vector<std::pair<double,double*>> ints_norms; //intersections and normals
    ns->castAxisOrientedRay (start_ray, r_c[ct.dir], ints_norms, ct.dir, compute_normals);

    if (ints_norms.size () != 0) {
      for (unsigned i = 0; i < ints_norms.size (); i++) {
        ct.inters.push_back (ints_norms[i].first);

        ct.normals.push_back (* (ints_norms[i].second));
        ct.normals.push_back (* (ints_norms[i].second+1));
        ct.normals.push_back (* (ints_norms[i].second+2));
      }
    }

    ct.init = 1; //ray is now initialized
  }

  if (ct.dir == 1) {
    double start_ray[3] = {ct.point[0], l_c[ct.dir], ct.point[1]};
    bool compute_normals = true;


    if (ct.point[0] < l_c[0] || ct.point[1] < l_c[2] || ct.point[0] > r_c[0] || ct.point[1] > r_c[2]) { //out of the molecule
      ct.init = 1;
      return;
    }

    std::cout << "Sending new ray in y direction for NS!" << std::endl;
    ns->setDirection (ct.dir);
    std::vector<std::pair<double,double*>> ints_norms; //intersections and normals
    ns->castAxisOrientedRay (start_ray, r_c[ct.dir], ints_norms, ct.dir, compute_normals);

    if (ints_norms.size () != 0) {
      for (unsigned i = 0; i < ints_norms.size (); i++) {
        ct.inters.push_back (ints_norms[i].first);

        ct.normals.push_back (* (ints_norms[i].second));
        ct.normals.push_back (* (ints_norms[i].second+1));
        ct.normals.push_back (* (ints_norms[i].second+2));
      }
    }

    ct.init = 1; //ray is now initialized
  }

  if (ct.dir == 2) {
    double start_ray[3] = {ct.point[0], ct.point[1], l_c[ct.dir]};
    bool compute_normals = true;


    if (ct.point[0] < l_c[0] || ct.point[1] < l_c[1] || ct.point[0] > r_c[0] || ct.point[1] > r_c[1]) { //out of the molecule
      ct.init = 1;
      return;
    }

    std::cout << "Sending new ray in z direction for NS!" << std::endl;
    ns->setDirection (ct.dir);
    std::vector<std::pair<double,double*>> ints_norms; //intersections and normals
    ns->castAxisOrientedRay (start_ray, r_c[ct.dir], ints_norms, ct.dir, compute_normals);

    if (ints_norms.size () != 0) {
      for (unsigned i = 0; i < ints_norms.size (); i++) {
        ct.inters.push_back (ints_norms[i].first);

        ct.normals.push_back (* (ints_norms[i].second));
        ct.normals.push_back (* (ints_norms[i].second+1));
        ct.normals.push_back (* (ints_norms[i].second+2));
      }
    }

    ct.init = 1; //ray is now initialized
  }

}


//! Convert an crossings_t object
//!  to a vector of bytes that can
//!  be saved to a binary file or transmitted through
//!  a channel such as a socket or sent as an MPI message.
std::vector<unsigned char>
ray_cache_t::write_ct (const crossings_t& ct)
{

  size_t N = sizeof (unsigned) + sizeof (bool) + sizeof (double)*2;
  std::vector<unsigned char> res (N, 0);

  std::copy (&ct.dir, &ct.dir+1, reinterpret_cast<unsigned*> (& (res[0])));
  std::copy (&ct.init, &ct.init+1, reinterpret_cast<bool*> (& (res[sizeof (unsigned)])));
  std::copy (ct.point, ct.point+2, reinterpret_cast<double*> (& (res[sizeof (bool)+sizeof (unsigned)])));

  std::vector<unsigned char> tmp = serialize::write (ct.inters);
  res.insert (res.end (), tmp.begin (), tmp.end ());
  tmp = serialize::write (ct.normals);
  res.insert (res.end (), tmp.begin (), tmp.end ());
  return res;
}

void
ray_cache_t::read_ct (const std::vector<unsigned char>& data, crossings_t& ct)
{

  size_t size_fixed = sizeof (unsigned) + sizeof (bool) + sizeof (double)*2;
  size_t size_non_fixed = 4*sizeof (double); // + sizeof(bool);
  size_t numel_flags = ( (data.size () * sizeof (unsigned char)) - size_fixed) / size_non_fixed;

  std::copy (data.begin (), data.begin () + sizeof (unsigned),
             reinterpret_cast<unsigned char*> (& (ct.dir)));

  auto init = data.begin () + sizeof (unsigned);
  std::copy (init, init + sizeof (bool),
             reinterpret_cast<unsigned char*> (& (ct.init)));

  init = init + sizeof (bool);
  std::copy (init, init + sizeof (double)*2,
             reinterpret_cast<unsigned char*> (& (ct.point)));

  if (numel_flags != 0) {
    init = init + sizeof (double)*2;
    //ct.flags.resize (numel_flags);
    //std::copy (init, init + sizeof(bool)*numel_flags,
    //reinterpret_cast<unsigned char*> (&(ct.flags[0])));

    //init = init + sizeof(bool)*numel_flags;
    ct.inters.resize (numel_flags);
    std::copy (init, init + sizeof (double)*numel_flags,
               reinterpret_cast<unsigned char*> (& (ct.inters[0])));

    init = init + sizeof (double)*numel_flags;
    ct.normals.resize (3*numel_flags);
    std::copy (init, data.end (),
               reinterpret_cast<unsigned char*> (& (ct.normals[0])));
  }
}

std::vector<unsigned char>
ray_cache_t::write_map (const std::map<std::array<double, 2>, crossings_t, map_compare>& container)
{
  size_t numel = container.size (); //number of elements inside the container
  size_t size = numel * (sizeof (double)*2); //array of double

  std::vector<unsigned char> res (size, 0);

  auto destf = reinterpret_cast<std::array<double, 2>*> (& (res[0]));

  for (auto ii = container.begin (); ii != container.end (); ++ii, ++destf)
    *destf = ii->first;

  for (auto ii = container.begin (); ii != container.end (); ++ii) {
    std::vector<unsigned char> tmp = write_ct (ii->second);
    res.insert (res.end (), tmp.begin (), tmp.end ());
  }

  std::vector<size_t> num_flags (numel + 1);
  int pos = 1;
  num_flags[0] = numel;

  for (auto it = container.cbegin (); it != container.cend (); it++) {
    num_flags[pos] = it->second.inters.size ();
    pos++;
  }

  std::vector<unsigned char> tmp2 = serialize::write (num_flags);
  res.insert (res.begin (), tmp2.begin (), tmp2.end ()); //at the begin the info about map size

  return res;
}

void
ray_cache_t::read_map (const std::vector<unsigned char>& data,
                       std::map<std::array<double, 2>, crossings_t, map_compare>& container)
{
  const size_t *psize = (reinterpret_cast<const size_t*> (& (data[0])));
  size_t num_maps = *psize;

  std::vector<size_t> num_flags (num_maps);
  size_t size_num_flags = sizeof (size_t) * num_maps + sizeof (size_t); //da dove iniziano i dati

  psize++;

  for (size_t jj = 0; jj < num_maps; jj++, psize++)
    num_flags[jj] = (*psize);

  auto destf = reinterpret_cast<const std::array<double, 2>*> (psize);

  for (size_t ii = 0; ii < num_maps; ++ii, ++destf)
    container[*destf];

  size_t jj = 0;
  auto destdir = reinterpret_cast<const unsigned*> (destf);

  for (auto ii = container.begin (); ii != container.end (); ++ii) {
    ii->second.dir = *destdir;
    destdir++;
    auto destinit = reinterpret_cast<const bool*> (destdir);
    ii->second.init = *destinit;
    destinit++;
    auto destpoint = reinterpret_cast<const double*> (destinit);
    ii->second.point[0] = *destpoint;
    destpoint++;
    ii->second.point[1] = *destpoint;
    destpoint++;

    if (num_flags[jj] != 0) {
      auto end = destpoint + num_flags[jj];
      ii->second.inters.resize (num_flags[jj]);
      std::copy (destpoint, end, ii->second.inters.begin ());
      destpoint = end;
      end = destpoint + 3*num_flags[jj];
      ii->second.normals.resize (3*num_flags[jj]);
      std::copy (destpoint, end, ii->second.normals.begin ());
      destpoint = end;
    }

    destdir = reinterpret_cast<const unsigned*> (destpoint);
    jj++;
  }
}
