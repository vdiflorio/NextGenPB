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
 * @file raytracer.h
 * @brief Ray-casting cache for molecular surface intersection detection.
 *
 * Provides `ray_cache_t`, which stores and lazily computes axis-aligned ray
 * intersections with the NanoShaper molecular surface.  Results are cached in
 * a map and distributed across MPI ranks via fill_cache().
 */
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

using int_coord_t = unsigned long long int; ///< Type for encoded ray coordinates

/**
 * @brief Cache of ray–surface intersections for the NanoShaper molecular surface.
 *
 * Axis-aligned rays are identified by a pair (fixed_coord_0, fixed_coord_1) plus
 * a direction (0=x, 1=y, 2=z).  Each ray stores a sorted list of crossing
 * distances (crossings_t).  Rays are computed on demand by calling operator()
 * and then synchronised across MPI ranks by fill_cache().
 */
struct
  ray_cache_t {

  rays_t rays; ///< All cached rays in 3 directions: rays[dir][(a,b)] → crossings

  static int_coord_t count_cache; ///< Total rays served from cache (all ranks)
  static int_coord_t count_new;   ///< Total rays computed fresh (all ranks)
  int count_new_dir[3]   = {0, 0, 0}; ///< New rays per direction
  int count_cache_dir[3] = {0, 0, 0}; ///< Cached rays per direction
  int count = 0;                       ///< Total ray lookups this run

  std::unique_ptr<NS::NanoShaper> ns; ///< NanoShaper surface instance

  double l_c[3] = {0, 0, 0}; ///< Domain bounding box lower corner [Å]
  double r_c[3] = {0, 0, 0}; ///< Domain bounding box upper corner [Å]

  int num_req_rays[3] = {0, 0, 0}; ///< Required ray counts per direction
  std::array<std::set<std::array<double, 2>, map_compare>, 3> rays_list; ///< Requested ray coordinates

  /// @brief Query (or compute) the surface crossings for the ray at (x0, x1) in direction @p dir.
  /// @param x0  First fixed coordinate of the ray [Å].
  /// @param x1  Second fixed coordinate of the ray [Å].
  /// @param dir Ray direction: 0=x, 1=y, 2=z.
  /// @return Reference to the sorted crossings list for this ray.
  crossings_t &
  operator () (double x0, double x1, unsigned dir = 1);

  /// @brief Compute all requested rays on rank 0 and broadcast results to all MPI ranks.
  void
  fill_cache ();

  /// @brief Initialise the NanoShaper surface from atomic data (analytical surface mode).
  /// @param atoms       List of atoms.
  /// @param surf_type   Surface type (SES, SAS, …).
  /// @param surf_param  Smoothing parameter.
  /// @param stern_layer Stern layer thickness [Å].
  /// @param num_threads Number of threads for NanoShaper.
  /// @param configFile  Optional path to a NanoShaper configuration file.
  void
  init_analytical_surf (const std::vector<NS::Atom> & atoms, const NS::surface_type & surf_type,
                        const double & surf_param, const double & stern_layer, const unsigned & num_threads,const std::string* configFile=nullptr);

  /// @brief Initialise the NanoShaper surface with a custom bounding box and scale (NS mode).
  void
  init_analytical_surf_ns (const std::vector<NS::Atom> & atoms, const NS::surface_type & surf_type,
                           const double & surf_param, const double & stern_layer, const unsigned & num_threads,
                           double* l_cr, double* r_cr, double scale,const std::string* configFile=nullptr);

  /// @brief Compute surface intersections for all rays in @p ct using NanoShaper.
  /// @param ct Crossing data container to populate.
  void
  compute_ns_inters (crossings_t & ct);

  /// @brief Serialise a ray map to a byte array (for MPI communication or checkpointing).
  std::vector<unsigned char>
  write_map (const std::map<std::array<double, 2>, crossings_t, map_compare>& container);

  /// @brief Deserialise a ray map from a byte array.
  void
  read_map (const std::vector<unsigned char>& data,
            std::map<std::array<double, 2>, crossings_t, map_compare>& container);

  /// @brief Serialise a single crossings_t object to a byte array.
  std::vector<unsigned char>
  write_ct (const crossings_t& ct);

  /// @brief Deserialise a single crossings_t object from a byte array.
  void
  read_ct (const std::vector<unsigned char>& data, crossings_t& ct);
};

#endif //RAYTRACER_H
