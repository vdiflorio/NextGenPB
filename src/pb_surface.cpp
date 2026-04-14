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
#include <quad_operators_3d.h>
#include <p8est.h>
#include <cmath>

double
poisson_boltzmann::levelsetfun (double x, double y, double z)
{
  double dist = 0.0;

  for (const NS::Atom& i : atoms) {

    dist += std::exp (surf_param * ( (std::pow (x - i.pos[0], 2) +
                                      std::pow (y - i.pos[1], 2) +
                                      std::pow (z - i.pos[2], 2)) /
                                     (i.radius *i.radius) - 1.0));

    if (dist > 1.5)
      break;

  }

  return dist;
}


double
poisson_boltzmann::is_in_ns_surf (ray_cache_t &ray_cache, double x, double y, double z, int dir)
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);
  double x1 = x;
  double x2 = y;
  double x3 = z;

  if (dir == 0) {
    x1 = y;
    x2 = z;
    x3 = x;
  } else if (dir == 1) {
    x1 = x;
    x2 = z;
    x3 = y;
  }

  crossings_t &ct = ray_cache (x1, x2, dir);

  if (!ct.init && rank != 0) {
    std::array<double, 2> ray = {x1, x2};
    ray_cache.rays[dir].erase (ray);
    return -1.;
  }


  int i = 0;

  if (ct.inters.size () == 0 || x3 < ct.inters[i])
    return 0; //if there are no inters or y_the coord is before the first intersection, the point is outside.

  while (i < ct.inters.size () && x3 >= ct.inters[i]) //go on until the inters is passed
    i++;

  return (i % 2);
}


bool
poisson_boltzmann::is_in (const NS::Atom &i,
                          tmesh_3d::quadrant_iterator q)
{
  double tol = p4esttol * (rr[0] - ll[0]);

  if (mesh_shape == 2)
    tol = p4esttol * (r_c[0] - l_c[0]);

  bool retval = false;
  double l, r, t, b, f, bk;

  l = q->p (0, 0);
  r = q->p (0, 7);

  f = q->p (1, 0);
  bk = q->p (1, 7);

  b = q->p (2, 0);
  t = q->p (2, 7);


  // for (int ii = 1; ii < 8; ++ii) {
  // l = q->p (0, ii) < l ? q->p (0, ii) : l;
  // r = q->p (0, ii) > r ? q->p (0, ii) : r;
  // f = q->p (1, ii) < f ? q->p (1, ii) : f;
  // bk = q->p (1, ii) > bk ? q->p (1, ii) : bk;
  // b = q->p (2, ii) < b ? q->p (2, ii) : b;
  // t = q->p (2, ii) > t ? q->p (2, ii) : t;

  // }

  retval = (i.pos[0] > l - tol) && (i.pos[0] <= r - tol); //make sure that the charge is assigned only once
  retval = retval && (i.pos[1] > f - tol) && (i.pos[1] <= bk - tol);
  retval = retval && (i.pos[2] > b - tol) && (i.pos[2] <= t - tol);

  return retval;
}

bool
poisson_boltzmann::is_in_ref (const NS::Atom &i,
                              tmesh_3d::quadrant_iterator q)
{

  bool retval = false;
  double l, r, t, b, f, bk;

  l = q->p (0, 0);
  r = q->p (0, 7);

  f = q->p (1, 0);
  bk = q->p (1, 7);

  b = q->p (2, 0);
  t = q->p (2, 7);


  // for (int ii = 1; ii < 8; ++ii) {
  // l = q->p (0, ii) < l ? q->p (0, ii) : l;
  // r = q->p (0, ii) > r ? q->p (0, ii) : r;
  // f = q->p (1, ii) < f ? q->p (1, ii) : f;
  // bk = q->p (1, ii) > bk ? q->p (1, ii) : bk;
  // b = q->p (2, ii) < b ? q->p (2, ii) : b;
  // t = q->p (2, ii) > t ? q->p (2, ii) : t;

  // }

  retval = (i.pos[0] >= l) && (i.pos[0] <= r);
  retval = retval && (i.pos[1] >= f) && (i.pos[1] <= bk);
  retval = retval && (i.pos[2] >= b) && (i.pos[2] <= t);

  return retval;
}

void
poisson_boltzmann::refine_surface (ray_cache_t &ray_cache)
{
  int rank, size;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);

  int num_cycles = 2;

  if (size == 1 || surf_type_num == 2)
    num_cycles = 1;

  int coars_ref_cycles = (maxlevel - unilevel) > (unilevel - minlevel) ? (maxlevel - unilevel) : (unilevel - minlevel);

  for (int kk = 0; kk < coars_ref_cycles; ++kk) {
    // REFINEMENT
    {
      distributed_vector rcoeff (tmsh.num_owned_nodes ());

      for (int jj = 0; jj < num_cycles; jj++) {
        ray_cache.num_req_rays[0] = 0; //zero at each ref/coarsen cycle
        ray_cache.num_req_rays[1] = 0; //zero at each ref/coarsen cycle
        ray_cache.num_req_rays[2] = 0; //zero at each ref/coarsen cycle
        ray_cache.rays_list[0].clear ();
        ray_cache.rays_list[1].clear ();
        ray_cache.rays_list[2].clear ();

        if (rank == 0 && jj == 0)
          std::cout << "Refinement: " << kk << std::endl;

        for (auto quadrant = tmsh.begin_quadrant_sweep ();
             quadrant != tmsh.end_quadrant_sweep ();
             ++quadrant) {

          for (int ii = 0; ii < 8; ++ii) {

            if (! quadrant->is_hanging (ii)) {
              if (surf_type_num == 2)
                rcoeff[quadrant->gt (ii)] = levelsetfun (quadrant->p (0, ii),
                                            quadrant->p (1, ii),
                                            quadrant->p (2, ii));
              else {
                for (int idir = 0; idir < 3; ++idir) {
                  rcoeff[quadrant->gt (ii)] = is_in_ns_surf (ray_cache,
                                              quadrant->p (0, ii),
                                              quadrant->p (1, ii),
                                              quadrant->p (2, ii), idir);

                  if (rcoeff[quadrant->gt (ii)] < -0.5) {
                    ray_cache.num_req_rays[idir]++;

                    std::array<double, 2> ray;

                    std::vector<int> direzioni = {0, 1, 2};
                    direzioni.erase (direzioni.begin () + idir);

                    for (unsigned i = 0; i < direzioni.size (); ++i) {
                      ray[i] = quadrant->p (direzioni[i], ii);
                    }

                    ray_cache.rays_list[idir].insert (ray);
                  }
                }
              }
            } else {
              for (int idir = 0; idir < 3; ++idir) {
                double pippo = is_in_ns_surf (ray_cache,
                                              quadrant->p (0, ii),
                                              quadrant->p (1, ii),
                                              quadrant->p (2, ii), idir);

                if (pippo < -0.5) {
                  ray_cache.num_req_rays[idir]++;

                  std::array<double, 2> ray;

                  std::vector<int> direzioni = {0, 1, 2};
                  direzioni.erase (direzioni.begin () + idir);

                  for (unsigned i = 0; i < direzioni.size (); ++i) {
                    ray[i] = quadrant->p (direzioni[i], ii);
                  }

                  ray_cache.rays_list[idir].insert (ray);
                }
              }

              for (int jj = 0; jj < quadrant->num_parents (ii); ++jj) {
                rcoeff[quadrant->gparent (jj, ii)] += 0.;
              }
            }
          }
        }

        MPI_Barrier (mpicomm);
        ray_cache.fill_cache ();
      }

      auto refinement = [&rcoeff, this]
      (tmesh_3d::quadrant_iterator q) -> int {
        int currentlevel = static_cast<int> (q->the_quadrant->level);
        int retval = 0;
        double min = 2.0;
        double max = 0.0;
        double tmp = 0.0;

        if (currentlevel >= this->maxlevel)
          retval = 0;
        else {
          for (int ii = 0; ii < 8; ++ii) {

            if (! q->is_hanging (ii)) {
              tmp = rcoeff[q->gt (ii)];

              if (tmp > max) max = tmp;

              if (tmp < min) min = tmp;
            }

          }

          if (this->surf_type_num == 2) {
            if (max > 1 && min < 1)
              retval = this->maxlevel - currentlevel;
            else
              for (const NS::Atom& i : atoms)
                if (is_in_ref (i, q)) {
                  retval = this->maxlevel - currentlevel;
                  break;
                }
          } else {
            if (max > 0.5 && min < 0.5)
              retval = this->maxlevel - currentlevel;
            else
              for (const NS::Atom& i : atoms)
                if (is_in_ref (i, q)) {
                  retval = this->maxlevel - currentlevel;
                  break;
                }
          }
        }

        return (retval);
      };

      tmsh.set_refine_marker (refinement);
      tmsh.refine (0, 1);
    }

    // COARSENING
    {
      distributed_vector rcoeff (tmsh.num_owned_nodes ());

      for (int jj = 0; jj < num_cycles; jj++) {
        ray_cache.num_req_rays[0] = 0; //zero at each ref/coarsen cycle
        ray_cache.num_req_rays[1] = 0; //zero at each ref/coarsen cycle
        ray_cache.num_req_rays[2] = 0; //zero at each ref/coarsen cycle
        ray_cache.rays_list[0].clear ();
        ray_cache.rays_list[1].clear ();
        ray_cache.rays_list[2].clear ();

        if (rank == 0 && jj == 0)
          std::cout << "Coarsening: " << kk << std::endl;

        for (auto quadrant = tmsh.begin_quadrant_sweep ();
             quadrant != tmsh.end_quadrant_sweep ();
             ++quadrant) {

          for (int ii = 0; ii < 8; ++ii) {

            if (! quadrant->is_hanging (ii)) {
              if (surf_type_num == 2)
                rcoeff[quadrant->gt (ii)] = levelsetfun (quadrant->p (0, ii),
                                            quadrant->p (1, ii),
                                            quadrant->p (2, ii));
              else {
                for (int idir = 0; idir < 3; ++idir) {
                  rcoeff[quadrant->gt (ii)] = is_in_ns_surf (ray_cache,
                                              quadrant->p (0, ii),
                                              quadrant->p (1, ii),
                                              quadrant->p (2, ii), idir);

                  if (rcoeff[quadrant->gt (ii)] < -0.5) {
                    ray_cache.num_req_rays[idir]++;

                    std::array<double, 2> ray;

                    std::vector<int> direzioni = {0, 1, 2};
                    direzioni.erase (direzioni.begin () + idir);

                    for (unsigned i = 0; i < direzioni.size (); ++i) {
                      ray[i] = quadrant->p (direzioni[i], ii);
                    }

                    // ray_cache.rays[idir].erase(ray);
                    ray_cache.rays_list[idir].insert (ray);
                  }
                }
              }
            } else {
              for (int idir = 0; idir < 3; ++idir) {
                double pippo = is_in_ns_surf (ray_cache,
                                              quadrant->p (0, ii),
                                              quadrant->p (1, ii),
                                              quadrant->p (2, ii), idir);

                if (pippo < -0.5) {
                  ray_cache.num_req_rays[idir]++;

                  std::array<double, 2> ray;

                  std::vector<int> direzioni = {0, 1, 2};
                  direzioni.erase (direzioni.begin () + idir);

                  for (unsigned i = 0; i < direzioni.size (); ++i) {
                    ray[i] = quadrant->p (direzioni[i], ii);
                  }

                  ray_cache.rays_list[idir].insert (ray);
                }
              }

              for (int jj = 0; jj < quadrant->num_parents (ii); ++jj) {
                rcoeff[quadrant->gparent (jj, ii)] += 0.;
              }
            }
          }
        }

        MPI_Barrier (mpicomm);
        ray_cache.fill_cache ();
      }

      auto coarsening = [&rcoeff, this]
      (tmesh_3d::quadrant_iterator q) -> int {
        int currentlevel = static_cast<int> (q->the_quadrant->level);
        int retval = 0;
        double min = 2.0;
        double max = 0.0;
        double tmp = 0.0;

        if (currentlevel <= this->minlevel)
          retval = 0;
        else {
          for (int ii = 0; ii < 8; ++ii) {

            if (! q->is_hanging (ii)) {
              tmp = rcoeff[q->gt (ii)];

              if (tmp > max) max = tmp;

              if (tmp < min) min = tmp;
            }

          }

          if (this->surf_type_num == 2) {
            if (min > 1 || max < 1)
              retval = currentlevel - this->minlevel;
          } else {
            if (min > 0.5 || max < 0.5)
              retval = currentlevel - this->minlevel;
          }

          for (const NS::Atom& i : atoms)
            if (is_in_ref (i, q)) {
              retval = 0;
              break;
            }
        }

        return (retval);
      };

      tmsh.set_coarsen_marker (coarsening);
      tmsh.coarsen (0, 1);
    }
  }
}



double
poisson_boltzmann::is_in_ns_surf_stern (ray_cache_t &ray_cache, double x, double y, double z, int dir)
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);
  double x1 = x;
  double x2 = y;
  double x3 = z;

  if (dir == 0) {
    x1 = y;
    x2 = z;
    x3 = x;
  } else if (dir == 1) {
    x1 = x;
    x2 = z;
    x3 = y;
  }

  crossings_t &ct = ray_cache (x1, x2, dir);

  if (!ct.init && rank != 0) {
    std::array<double, 2> ray = {x1, x2};
    ray_cache.rays[dir].erase (ray);
    return -1.;
  }

  int i = 0;
  int sign = 1;

  if (ct.inters.size () == 0 || x3 < (ct.inters[i] - stern_layer))
    return 0; //if there are no inters or y_the coord is before the first intersection, the point is outside.

  while (i < ct.inters.size () && x3 > (ct.inters[i] - stern_layer * sign)) {
    //go on until the inters is passed
    i++;
    sign *= -1;
  }

  return (i % 2);
}


void
poisson_boltzmann::create_markers (ray_cache_t &ray_cache)
{

  int size, rank;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);

  double eps_in = 4.0 * pi * e_0 * e_in * kb * T * Angs / (e *e); //adim e_in
  double eps_out = 4.0 * pi * e_0 * e_out * kb * T * Angs / (e *e); //adim e_out

  /////////////////////////////////////////////////////////
  //reactions
  double C_0 = 1.0e3 * N_av * ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0 * C_0 * Angs * Angs * e * e / (e_0 *e_out *kb *T);

  this->marker.assign (this->tmsh.num_local_quadrants (), 0.0); //marker = 0 -> in

  if (stern_layer_surf == 1) {
    this->marker_k.assign (this->tmsh.num_local_quadrants (), 1.0); //marker = 1 -> out stern
    this->reaction.assign (tmsh.num_local_quadrants (), eps_out *k2);
  }

  epsilon_nodes = std::make_unique<distributed_vector> (tmsh.num_owned_nodes (), mpicomm);
  epsilon_nodes->get_owned_data ().assign (tmsh.num_owned_nodes (), eps_out);

  reaction_nodes = std::make_unique<distributed_vector> (tmsh.num_owned_nodes (), mpicomm);
  reaction_nodes->get_owned_data ().assign (tmsh.num_owned_nodes (), eps_out *k2);

  int num_cycles = 2;

  if (size == 1) {
    num_cycles = 1;
  }

  int local_num;

  for (int jj = 0; jj < num_cycles; ++jj) {
    ray_cache.num_req_rays[0] = 0; //zero at each ref/coarsen cycle
    ray_cache.num_req_rays[1] = 0; //zero at each ref/coarsen cycle
    ray_cache.num_req_rays[2] = 0; //zero at each ref/coarsen cycle
    ray_cache.rays_list[0].clear ();
    ray_cache.rays_list[1].clear ();
    ray_cache.rays_list[2].clear ();


    for (auto quadrant = this->tmsh.begin_quadrant_sweep ();
         quadrant != this->tmsh.end_quadrant_sweep ();
         ++quadrant) {
      int num_int_nodes = 0;
      int num_int_nodes_stern = 0;
      int num_hanging[3] = {0, 0, 0};
      double x, y, z;

      for (int ii = 0; ii < 8; ++ii) {
        local_num = quadrant->gt (ii);
        x = quadrant->p (0, ii);
        y = quadrant->p (1, ii);
        z = quadrant->p (2, ii);

        if (! quadrant->is_hanging (ii)) {
          if (this->is_in_ns_surf (ray_cache, x, y, z, 2) > 0.5) { //inside the molecule
            ++num_int_nodes;
            ++num_int_nodes_stern;
            (*epsilon_nodes)[local_num] = eps_in;
            (*reaction_nodes)[local_num] = 0.0;
          } else if (this->is_in_ns_surf (ray_cache, x, y, z, 2) < -0.5) {
            ray_cache.num_req_rays[2]++;

            // std::array<double, 2> ray {x,y};

            ray_cache.rays_list[2].insert ({x, y});
          }

          if (this->is_in_ns_surf (ray_cache, x, y, z, 0) < -0.5) {
            ray_cache.num_req_rays[0]++;

            // std::array<double, 2> ray {y,z};

            ray_cache.rays_list[0].insert ({y, z});
          }

          if (this->is_in_ns_surf (ray_cache, x, y, z, 1) < -0.5) {
            ray_cache.num_req_rays[1]++;

            std::array<double, 2> ray {x, z};

            ray_cache.rays_list[1].insert ({x, z});
          }

          if (stern_layer_surf == 1) {
            if (this->is_in_ns_surf_stern (ray_cache, x, y, z, 2) > 0.5) { //inside the stern layer
              ++num_int_nodes_stern;
            }
          }


        } else
          for (int idir = 0; idir < 3; ++idir) {
            ++num_hanging[idir];

            if (this->is_in_ns_surf (ray_cache, x, y, z, idir) < -0.5) {
              ray_cache.num_req_rays[idir]++;
              std::array<double, 2> ray;

              std::vector<int> direzioni {0, 1, 2};
              direzioni.erase (direzioni.begin () + idir);

              for (unsigned i = 0; i < direzioni.size (); ++i) {
                ray[i] = quadrant->p (direzioni[i], ii);
              }

              ray_cache.rays_list[idir].insert (ray);

            }
          }
      }

      if (jj != 0 || num_cycles == 1) {
        if (num_int_nodes == 0) { //if there's no node inside the molecule
          this->marker[quadrant->get_forest_quad_idx ()] = 1.0; //quadrant is out
        } else if (num_int_nodes < (8 - num_hanging[2])) { //if the non hanging nodes are not all inside
          this->marker[quadrant->get_forest_quad_idx ()] = 1.0 / 2.0; //"border"
          border_quad.push_back (quadrant->get_forest_quad_idx ());
        }

        //else: all the nodes are inside: the quadrant is inside and the marker value is 0
        if (stern_layer_surf == 1) {
          for (int idir = 0; idir < 3; ++idir)
            if (num_int_nodes_stern != 0) { //if there is at least on node inside the stern layer along idir-axis
              this->marker_k[quadrant->get_forest_quad_idx ()] = 0.0; //quadrant is in
              this->reaction[quadrant->get_forest_quad_idx ()] = 0.0; //quadrant is in
            }
        }

      }

    }

    MPI_Barrier (mpicomm);
    ray_cache.fill_cache ();
  }

  if (size > 1) {
    bim3a_solution_with_ghosts (tmsh, *epsilon_nodes, replace_op);

    if (stern_layer_surf == 0)
      bim3a_solution_with_ghosts (tmsh, (*reaction_nodes), replace_op);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

std::array<double, 12>
poisson_boltzmann::cube_fraction_intersection (tmesh_3d::quadrant_iterator &quadrant,
    const ray_cache_t &ray_cache)


// v6_________e7_________v7
// /|                  /|
// e8 / |                 / |
// /  |             e6 /  |
// /   | e12           /   | e11
// v4/____|_____e5_______/v5  |
// |    |              |    |
// |  v2|______e3______|____|v3
// e9  |   /               |   /
// |  /            e10 |  /
// | /  e4             | / e2
// |/                  |/
// v0/_________e1________/v1
{
  std::array<double, 12> fraction = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};

  int dir;
  int i1, i2;
  double x1, x2;
  std::array<double, 2> ray;

  if (marker[quadrant->get_forest_quad_idx ()] == 0.5) {
    for (int j : {
           0, 2, 4, 6
         })
      //x-axis edge
    {
      dir = 0;
      i1 = edge2nodes[2 * j ];
      i2 = edge2nodes[2 * j + 1];
      x1 = quadrant->p (dir, i1);
      x2 = quadrant->p (dir, i2);
      std::vector<int> direzioni {0, 1, 2};
      direzioni.erase (direzioni.begin () + dir);

      for (unsigned i = 0; i < direzioni.size (); ++i) {
        ray[i] = quadrant->p (direzioni[i], i1);
      }

      auto it0 = ray_cache.rays[dir].find (ray);
      auto inters = it0->second.inters;

      for (int ii = 0; ii < inters.size (); ii++) {
        if (inters[ii] >= x1 && inters[ii] <= x2) {
          fraction[j] = (inters[ii] - x1) / (x2 - x1);
        }
      }
    }

    for (int j : {
           1, 3, 5, 7
         })
      //y-axis edge
    {
      dir = 1;
      i1 = edge2nodes[2 * j ];
      i2 = edge2nodes[2 * j + 1];
      x1 = quadrant->p (dir, i1);
      x2 = quadrant->p (dir, i2);
      std::vector<int> direzioni {0, 1, 2};
      direzioni.erase (direzioni.begin () + dir);

      for (unsigned i = 0; i < direzioni.size (); ++i) {
        ray[i] = quadrant->p (direzioni[i], i1);
      }

      auto it0 = ray_cache.rays[dir].find (ray);
      auto inters = it0->second.inters;

      for (int ii = 0; ii < inters.size (); ii++) {
        if (inters[ii] >= x1 && inters[ii] <= x2) {
          fraction[j] = (inters[ii] - x1) / (x2 - x1);
        }
      }
    }

    for (int j : {
           8, 9, 10, 11
         })
      //z-axis edge
    {
      dir = 2;
      i1 = edge2nodes[2 * j ];
      i2 = edge2nodes[2 * j + 1];
      x1 = quadrant->p (dir, i1);
      x2 = quadrant->p (dir, i2);
      std::vector<int> direzioni {0, 1, 2};
      direzioni.erase (direzioni.begin () + dir);

      for (unsigned i = 0; i < direzioni.size (); ++i) {
        ray[i] = quadrant->p (direzioni[i], i1);
      }

      auto it0 = ray_cache.rays[dir].find (ray);
      auto inters = it0->second.inters;

      for (int ii = 0; ii < inters.size (); ii++) {
        if (inters[ii] >= x1 && inters[ii] <= x2) {
          fraction[j] = (inters[ii] - x1) / (x2 - x1);
        }
      }
    }
  }

  return fraction;
}

void
poisson_boltzmann::normal_intersection (tmesh_3d::quadrant_iterator &quadrant,
                                        const ray_cache_t &ray_cache,
                                        int edge, std::array<double, 3> &norm,
                                        double &frac)
{
  int dir = edge_axis[edge];
  int i1 = edge2nodes[2 * edge ];
  int i2 = edge2nodes[2 * edge + 1];
  double x1 = quadrant->p (dir, i1),
         x2 = quadrant->p (dir, i2);

  std::array<double, 2> ray;
  std::vector<int> direzioni = {0, 1, 2};
  direzioni.erase (direzioni.begin () + dir);

  for (unsigned i = 0; i < direzioni.size (); ++i) {
    ray[i] = quadrant->p (direzioni[i], i1);
  }

  auto it0 = ray_cache.rays[dir].find (ray);
  auto normali = it0->second.normals;
  auto inters = it0->second.inters;

  frac = 0.5;

  for (int ii = 0; ii < inters.size (); ii++) {
    if (inters[ii] >= x1 && inters[ii] <= x2) {
      norm[0] = normali[0 + 3 * ii];
      norm[1] = normali[1 + 3 * ii];
      norm[2] = normali[2 + 3 * ii];
      frac = (inters[ii] - x1) / (x2 - x1);
    }
  }
}

int
poisson_boltzmann::classifyCube (tmesh_3d::quadrant_iterator &quadrant,
                                 double isolevel)
{
  int cubeindex = 0;
  int index = 1;
  double tmp = 0;
  constexpr double EPSILON = 1e-10;

  for (int ii : {
         0, 1, 3, 2, 4, 5, 7, 6
       }) {
    if (! quadrant->is_hanging (ii)) {
      if ( (*epsilon_nodes)[quadrant->gt (ii)] < (isolevel - EPSILON)) cubeindex |= index;
    } else {
      for (int jj = 0; jj < quadrant->num_parents (ii); ++jj) {
        tmp += (*epsilon_nodes)[quadrant->gparent (jj, ii)] / quadrant->num_parents (ii);
      }

      if (tmp < (isolevel - EPSILON)) cubeindex |= index;
    }

    tmp = 0;
    index *= 2;
  }

  // Cube is entirely in/out of the surface
  if (edgeTable[cubeindex] == 0)
    return -1;

  return cubeindex;
}

int
poisson_boltzmann::classifyCube_fast (tmesh_3d::quadrant_iterator &quadrant,
                                      double isolevel)
{
  int cubeindex = 0;
  int index = 1;
  double tmp = 0;
  constexpr double EPSILON = 1e-10;

  for (int ii : {
         0, 1, 3, 2, 4, 5, 7, 6
       }) {
    if ( (*epsilon_nodes)[quadrant->gt (ii)] < (isolevel - EPSILON)) cubeindex |= index;

    index *= 2;
  }

  // Cube is entirely in/out of the surface
  if (edgeTable[cubeindex] == 0)
    return -1;

  return cubeindex;
}

std::tuple<std::array<double, 8>, std::array<double, 8>, std::vector<int>, std::vector<int> >
poisson_boltzmann::classifyCube_flux (tmesh_3d::quadrant_iterator &quadrant,
                                      std::array<double, 8> &tmp_phi,
                                      std::array<double, 8> &tmp_eps)
{
  std::vector<int> edges {};
  std::vector<int> flux {};

  for (int ii = 0; ii < 8; ++ii) {
    if (! quadrant->is_hanging (ii)) {
      tmp_eps[ii] = (*epsilon_nodes)[quadrant->gt (ii)];
      tmp_phi[ii] = (*phi)[quadrant->gt (ii)];
    } else {
      tmp_eps[ii] = 0.0;
      tmp_phi[ii] = 0.0;

      for (int jj = 0; jj < quadrant->num_parents (ii); ++jj) {
        tmp_eps[ii] += (*epsilon_nodes)[quadrant->gparent (jj, ii)] / quadrant->num_parents (ii);
        tmp_phi[ii] += (*phi)[quadrant->gparent (jj, ii)] / quadrant->num_parents (ii);
      }
    }
  }

  for (int ii = 0; ii < 12; ++ii) {
    if (tmp_eps[edge2nodes[2 * ii]] < tmp_eps[edge2nodes[2 * ii + 1]]) {
      flux.push_back (1);
      edges.push_back (ii);
    } else if (tmp_eps[edge2nodes[2 * ii]] > tmp_eps[edge2nodes[2 * ii + 1]]) {
      flux.push_back (-1);
      edges.push_back (ii);
    }
  }

  return make_tuple (tmp_phi, tmp_eps, edges, flux);
}

std::tuple<std::array<double, 8>, std::array<double, 8>, std::vector<int>, std::vector<int> >
poisson_boltzmann::classifyCube_flux_fast (tmesh_3d::quadrant_iterator &quadrant,
    std::array<double, 8> &tmp_phi,
    std::array<double, 8> &tmp_eps)
{
  std::vector<int> edges {};
  std::vector<int> flux {};


  for (int ii = 0; ii < 8; ++ii) {
    tmp_eps[ii] = (*epsilon_nodes)[quadrant->gt (ii)];
    tmp_phi[ii] = (*phi)[quadrant->gt (ii)];
  }

  for (int ii = 0; ii < 12; ++ii) {
    if (tmp_eps[edge2nodes[2 * ii]] < tmp_eps[edge2nodes[2 * ii + 1]]) {
      flux.push_back (1);
      edges.push_back (ii);
    } else if (tmp_eps[edge2nodes[2 * ii]] > tmp_eps[edge2nodes[2 * ii + 1]]) {
      flux.push_back (-1);
      edges.push_back (ii);
    }
  }

  return make_tuple (tmp_phi, tmp_eps, edges, flux);
}

double wha (double eps1, double eps2, double frac)
{
  return 1.0 / (frac / eps1 + (1 - frac) / eps2);
}

double flux_dir (double eps1, double eps2)
{
  return eps1 < eps2 ? 1 : -1;
}

double phi0 (double eps1, double eps2,
             double phi1, double phi2, double frac)
{
  return phi1 + frac * eps2 * (phi2 - phi1) / (eps2 *frac + eps1* (1 - frac));
}

double areaTriangle (const std::array<std::array<double, 3>, 3> &triangle)
{

  double area;
  std::array<double, 3> ab;
  std::array<double, 3> ac;

  for (int i = 0; i < 3; ++i) {
    ab[i] = triangle[0][i] - triangle[1][i];
    ac[i] = triangle[0][i] - triangle[2][i];
  }

  area = 0.5 * std::hypot (ab[1] * ac[2] - ab[2] * ac[1],
                           ab[0] * ac[2] - ab[2] * ac[0],
                           ab[1] * ac[0] - ab[0] * ac[1]);
  return area;
}

double SphercalAreaTriangle (const std::array<std::array<double, 3>, 3> &triangle)
{

  double area;
  double rs = 2;
  double a, b, c, aa, bb, cc, E;
  double nA, nB, nC;
  std::array<double, 3> A = triangle[0];
  std::array<double, 3> B = triangle[1];
  std::array<double, 3> C = triangle[2];

  nA = std::hypot (A[0], A[1], A[2]);
  nB = std::hypot (B[0], B[1], B[2]);
  nC = std::hypot (C[0], C[1], C[2]);

  a = std::inner_product (std::begin (A), std::end (A), std::begin (B), 0.0);
  b = std::inner_product (std::begin (A), std::end (A), std::begin (C), 0.0);
  c = std::inner_product (std::begin (B), std::end (B), std::begin (C), 0.0);

  a = a / (nA *nB);
  b = b / (nA *nC);
  c = c / (nC *nB);

  a = std::acos (a);
  b = std::acos (b);
  c = std::acos (c);

  if (a > pi / 2.0)
    a = pi - a;

  if (b > pi / 2.0)
    b = pi - b;

  if (c > pi / 2.0)
    c = pi - c;

  aa = std::acos ( (cos (a) - cos (b) * cos (c)) / (sin (b) * sin (c)));
  bb = std::acos ( (cos (b) - cos (c) * cos (a)) / (sin (a) * sin (c)));
  cc = std::acos ( (cos (c) - cos (b) * cos (a)) / (sin (b) * sin (a)));

  E = aa + bb + cc - pi;

  area = rs * rs * E;

  return area;
}

int
poisson_boltzmann::getTriangles (int cubeindex,
                                 std::array<std::array<int, 3>, 5> &triangles)
{
  int ntriang = 0;
  int i;
  triangles.fill ({}); //set matrix to zero

  for (i = 0; triTable[cubeindex][i] != -1; i += 3) {
    // save all the assigned indexes
    triangles[ntriang][0] = triTable[cubeindex][i ];
    triangles[ntriang][1] = triTable[cubeindex][i + 1];
    triangles[ntriang][2] = triTable[cubeindex][i + 2];
    ntriang++;
  }

  return ntriang;
}

bool
poisson_boltzmann::controlla_coordinate (int i, const p8est_quadrant_t *quadrant)
{

  double tol = p4esttol * (rr[0] - ll[0]);

  if (mesh_shape == 2)
    tol = p4esttol * (r_c[0] - l_c[0]);


  bool retval = false;

  p8est_quadrant_t node;
  int ii, jj;
  double vxyz [3 * 8] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                        };

  for (ii = 0; ii < 8; ++ii) {
    p8est_quadrant_corner_node (quadrant, ii, &node);
    p8est_qcoord_to_vertex (tmsh.p8est->connectivity, 0,
                            node.x, node.y, node.z, & (vxyz[3 * ii]));
  }

  double l, r, t, b, f, bk;

  l = vxyz[0];
  r = vxyz[3 * 7];

  f = vxyz[1];
  bk = vxyz[3 * 7 + 1];

  b = vxyz[2];
  t = vxyz[3 * 7 + 2];

  // retval = (atoms[i].pos[0] > l- tol) && (atoms[i].pos[0] <= r- tol); //make sure that the charge is assigned only once
  // retval = retval && (atoms[i].pos[1] > f- tol) && (atoms[i].pos[1] <= bk- tol);
  // retval = retval && (atoms[i].pos[2] > b- tol) && (atoms[i].pos[2] <= t- tol);
  retval = (pos_atoms[i][0] > l - tol) && (pos_atoms[i][0] <= r - tol); //make sure that the charge is assigned only once
  retval = retval && (pos_atoms[i][1] > f - tol) && (pos_atoms[i][1] <= bk - tol);
  retval = retval && (pos_atoms[i][2] > b - tol) && (pos_atoms[i][2] <= t - tol);

  return retval;
}


int
poisson_boltzmann::cerca_atomo (p8est_t *p4est,
                                p4est_topidx_t which_tree,
                                p8est_quadrant_t *quadrant,
                                p4est_locidx_t local_num,
                                void *point)
{
  int *pt = (int *) point;
  // std::cout << "\n pt: " <<*pt <<"\n"<< std::endl;
  bool tf = controlla_coordinate (*pt, quadrant);

  if (tf) {
    if (local_num >= 0) {
      auto quadrant = this->tmsh.begin_quadrant_sweep ();
      quadrant[local_num];
      tmesh_3d::quadrant_t qi = tmsh.current_quadrant;
      lookup_table.emplace (*pt, qi);

    }

    return 1;
  }

  return 0;

}

void
poisson_boltzmann::search_points()
{
  size_t count = charge_atoms.size();
  auto base = std::make_unique<int[]>(count);

  for (size_t ii = 0; ii < count; ++ii)
    base[ii] = ii;

  sc_array_t* points = sc_array_new_data(
                         base.get(), sizeof(int), count);

  // Set the global pointer for the callback
  pb_global_wrapper = this;

  // call the search function
  p8est_search_local(
    tmsh.p8est,
    0,
    NULL,
    cerca_atomo_wrapper,
    points
  );
}

