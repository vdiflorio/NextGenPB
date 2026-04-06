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
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdio.h>

// Forward declarations of file-scope helpers defined in pb_surface.cpp
double wha (double eps1, double eps2, double frac);
double flux_dir (double eps1, double eps2);
double phi0 (double eps1, double eps2, double phi1, double phi2, double frac);
double areaTriangle (const std::array<std::array<double,3>,3> &triangle);

void
poisson_boltzmann::export_tmesh (ray_cache_t & ray_cache)
{
  int size, rank;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);
  bim3a_solution_with_ghosts (tmsh, (*epsilon_nodes), replace_op);
  tmsh.octbin_export ("eps_map_0", (*epsilon_nodes));
}

void
poisson_boltzmann::export_potential_map (ray_cache_t & ray_cache)
{
  int size, rank;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);
  tmsh.octbin_export ("potential_map_0", (*phi));
}

void
poisson_boltzmann::export_marked_tmesh ()
{
  tmsh.octbin_export_quadrant (markerfilename.c_str (), marker);

  if (stern_layer_surf == 1) {
    tmsh.octbin_export_quadrant ("mark_stern_0", marker_k);
  }
}

void
poisson_boltzmann::export_p4est ()
{
  tmsh.save (p4estfilename.c_str ());
}



void
poisson_boltzmann::write_potential_on_atoms_fast ()
{
  int rank, size;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  std::vector<std::string> local_lines;

  double phi_on_atom;
  double phi_hang_nodes = 0.0;

  // Costruisci stringhe localmente
  for (auto it = lookup_table.begin(); it != lookup_table.end(); ++it) {
    phi_on_atom = 0.0;
    double volume = (it->second.p (0, 7) - it->second.p (0, 0)) *
                    (it->second.p (1, 7) - it->second.p (1, 0)) *
                    (it->second.p (2, 7) - it->second.p (2, 0));

    for (int ii = 0; ii < 8; ++ii) {
      double weigth = std::abs ((pos_atoms[it->first][0] - it->second.p (0, 7-ii)) *
                                (pos_atoms[it->first][1] - it->second.p (1, 7-ii)) *
                                (pos_atoms[it->first][2] - it->second.p (2, 7-ii))) / volume;

      if (!it->second.is_hanging (ii)) {
        phi_on_atom += (*phi)[it->second.gt (ii)] * weigth;
      } else {
        phi_hang_nodes = 0.0;

        for (int jj = 0; jj < it->second.num_parents (ii); ++jj) {
          phi_hang_nodes += (*phi)[it->second.gparent (jj, ii)] / it->second.num_parents (ii);
        }

        phi_on_atom += phi_hang_nodes * weigth;
      }
    }

    std::ostringstream oss;
    oss << std::setw (8) << index_atoms[it->first]
        << std::fixed << std::setprecision (3)
        << std::setw (8) << pos_atoms[it->first][0]
        << std::setw (8) << pos_atoms[it->first][1]
        << std::setw (8) << pos_atoms[it->first][2]
        << std::fixed << std::setprecision (4)
        << "  " << phi_on_atom << "\n";

    local_lines.push_back (oss.str());
  }

  // Serializzazione delle stringhe
  std::string local_data;

  for (const auto& line : local_lines)
    local_data += line;

  int local_size = local_data.size();
  std::vector<int> all_sizes (size);

  MPI_Gather (&local_size, 1, MPI_INT, all_sizes.data(), 1, MPI_INT, 0, mpicomm);

  std::vector<int> displs (size);
  std::string global_data;

  if (rank == 0) {
    int total_size = 0;

    for (int i = 0; i < size; ++i) {
      displs[i] = total_size;
      total_size += all_sizes[i];
    }

    global_data.resize (total_size);
  }

  MPI_Gatherv (local_data.data(), local_size, MPI_CHAR,
               rank == 0 ? &global_data[0] : nullptr,
               all_sizes.data(), displs.data(), MPI_CHAR,
               0, mpicomm);

  // Solo rank 0 scrive sul file
  if (rank == 0) {
    std::ofstream phi_atoms ("phi_on_atoms.txt");
    phi_atoms << global_data;
    phi_atoms.close();
  }
}


void
poisson_boltzmann::energy (ray_cache_t & ray_cache)
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);

  if (rank == 0)
    std::cout << "\n================ [ Electrostatic Energy ] =================\n";


  // ===========================
  // Costanti fisiche e scalari
  // ===========================
  const double inv_4pi = 1.0 / (4.0 * pi);
  const double eps0 = e_0; // Permittività del vuoto
  const double eps_in = 4.0 * pi * eps0 * e_in * kb * T * Angs / (e * e);
  const double eps_out = 4.0 * pi * eps0 * e_out * kb * T * Angs / (e * e);

  const double C0 = 1.0e3 * N_av * ionic_strength; // [mol/m^3]
  const double k2 = 2.0 * C0 * Angs * Angs * e * e / (eps0 * e_out * kb * T);
  const double k = std::sqrt (k2);

  const double den_in = 1.0 / eps_in;
  const double constant_pol = (1.0 / eps_out - 1.0 / eps_in) * inv_4pi;
  const double constant_react = (1.0 / eps_out) * inv_4pi;

  // ===========================
  // Variabili locali
  // ===========================
  // double energy_pol = 0.0;
  // double energy_react = 0.0;
  // double coul_energy = 0.0;
  double charge_pol = 0.0;

  // --- Filtra gli atomi caricati ---
  std::vector<double> charge_atoms_tmp;
  std::vector<std::array<double, 3>> pos_atoms_tmp;
  // charge_atoms_tmp.reserve(charge_atoms.size());
  // pos_atoms_tmp.reserve(pos_atoms.size());

  for (size_t ii = 0; ii < charge_atoms.size(); ++ii) {
    if (std::fabs (charge_atoms[ii]) > 1.e-5) {
      charge_atoms_tmp.push_back (charge_atoms[ii]);
      pos_atoms_tmp.push_back (pos_atoms[ii]);
    }
  }

  const size_t num_atoms = charge_atoms_tmp.size();
  std::vector<double>().swap (charge_atoms);
  std::vector<std::array<double, 3>>().swap (pos_atoms);

  // --- Energia Coulombiana ---
  if (calc_coulombic == 1) {
    for (size_t i = 0; i < num_atoms; ++i) {
      const double qi = charge_atoms_tmp[i];
      const std::array<double,3>& ri = pos_atoms_tmp[i];

      for (size_t j = i + 1; j < num_atoms; ++j) {
        const std::array<double,3>& rj = pos_atoms_tmp[j];
        const double dx = ri[0] - rj[0];
        const double dy = ri[1] - rj[1];
        const double dz = ri[2] - rj[2];
        const double r = std::sqrt (dx * dx + dy * dy + dz * dz);
        this->coul_energy += qi * charge_atoms_tmp[j] / r;
      }
    }

    this->coul_energy *= den_in;
  }

  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  double first_int = 0.0, second_int = 0.0;

  std::array<double,3> h{0}, area_h{0};
  std::array<double,3> V, N;
  std::array<double,8> tmp_eps, tmp_phi;
  std::vector<int> edg, fl_dir;


  auto quadrant = this->tmsh.begin_quadrant_sweep ();

  // polarization energy
  if (calc_energy==1 || (calc_energy==2 && k < 1.e-5)) {
    for (const int ii : border_quad) {
      quadrant[ii];

      for (int d = 0; d < 3; ++d)
        h[d] = quadrant->p (d, 7) - quadrant->p (d, 0);

      area_h = {h[1]*h[2]/h[0]*0.25, h[0]*h[2]/h[1]*0.25, h[0]*h[1]/h[2]*0.25};

      std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux (quadrant, tmp_phi, tmp_eps);

      for (int ip = 0; ip < edg.size (); ++ip) {
        const int edge = edg[ip];
        const int axis = edge_axis[edge];
        const int i1 = edge2nodes[2 * edge];
        const int i2 = edge2nodes[2 * edge + 1];

        double fract = 0.0;
        normal_intersection (quadrant, ray_cache, edge, N, fract);

        V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
        V[axis] += fract * h[axis];

        const double tmp_flux =
          - (tmp_phi[i2] - tmp_phi[i1]) * wha (tmp_eps[i1], tmp_eps[i2], fract)
          * fl_dir[ip] * area_h[axis];

        charge_pol += tmp_flux;

        for (size_t ia = 0; ia < num_atoms; ++ia) {
          const std::array<double,3> &ra = pos_atoms_tmp[ia];
          const double dx = ra[0] - V[0];
          const double dy = ra[1] - V[1];
          const double dz = ra[2] - V[2];
          const double r = std::sqrt (dx * dx + dy * dy + dz * dz);
          const double inv_r3 = 1.0 / (r * r * r);
          const double qflux = tmp_flux / r;

          first_int += charge_atoms_tmp[ia] * qflux;
        }
      }
    }

    energy_pol = 0.5*constant_pol*first_int;
  }


  //polarization energy + ionic energy
  if (calc_energy==2 && k > 1.e-5) {
    int cubeindex = -1;
    std::array<std::array<double,3>,3> vert_triangles, norms_vert;
    std::array<double,3> dist_vert, phi_sup;
    int ntriang = 0;

    for (const int ii : border_quad) {
      quadrant[ii];
      cubeindex = classifyCube (quadrant, eps_out);

      for (int d = 0; d < 3; ++d)
        h[d] = quadrant->p (d, 7) - quadrant->p (d, 0);

      area_h = {h[1]*h[2]/h[0]*0.25, h[0]*h[2]/h[1]*0.25, h[0]*h[1]/h[2]*0.25};

      std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux (quadrant, tmp_phi, tmp_eps);
      ntriang = getTriangles (cubeindex, triangles);

      // --- flussi
      for (int ip = 0; ip < edg.size (); ++ip) {
        const int edge = edg[ip];
        const int axis = edge_axis[edge];
        const int i1 = edge2nodes[2 * edge];
        const int i2 = edge2nodes[2 * edge + 1];

        double fract = 0.0;
        normal_intersection (quadrant, ray_cache, edge, N, fract);

        V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
        V[axis] += fract * h[axis];

        const double tmp_flux =
          - (tmp_phi[i2] - tmp_phi[i1]) * wha (tmp_eps[i1], tmp_eps[i2], fract)
          * fl_dir[ip] * area_h[axis];

        charge_pol += tmp_flux;

        for (size_t ia = 0; ia < num_atoms; ++ia) {
          const std::array<double,3> &ra = pos_atoms_tmp[ia];
          const double dx = ra[0] - V[0];
          const double dy = ra[1] - V[1];
          const double dz = ra[2] - V[2];
          const double r = std::sqrt (dx * dx + dy * dy + dz * dz);
          const double inv_r3 = 1.0 / (r * r * r);

          first_int += charge_atoms_tmp[ia] * tmp_flux / r;
        }
      }

      // --- triangoli (componente ionica)
      for (int itri = 0; itri < ntriang; ++itri) {
        for (int jj = 0; jj < 3; ++jj) {
          const int edge = triangles[itri][jj];
          const int axis = edge_axis[edge];
          const int i1 = edge2nodes[2 * edge];
          const int i2 = edge2nodes[2 * edge + 1];

          double fract = 0.0;
          normal_intersection (quadrant, ray_cache, edge, N, fract);

          V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
          V[axis] += fract * h[axis];

          vert_triangles[jj] = V;
          norms_vert[jj] = N;

          phi_sup[jj] = phi0 (tmp_eps[i1], tmp_eps[i2], tmp_phi[i1], tmp_phi[i2], fract);
        }

        const double area = areaTriangle (vert_triangles);

        for (size_t ia = 0; ia < num_atoms; ++ia) {
          const double qi = charge_atoms_tmp[ia];
          const std::array<double,3> &ra = pos_atoms_tmp[ia];

          for (int kk = 0; kk < 3; ++kk) {
            dist_vert = {vert_triangles[kk][0] - ra[0],
                         vert_triangles[kk][1] - ra[1],
                         vert_triangles[kk][2] - ra[2]
                        };
            const double r2 = dist_vert[0]*dist_vert[0] + dist_vert[1]*dist_vert[1] + dist_vert[2]*dist_vert[2];
            const double r = std::sqrt (r2);
            const double inv_r3 = 1.0 / (r2 * r);
            const double inv_r5 = inv_r3 / r2;
            const double dot = dist_vert[0]*norms_vert[kk][0] + dist_vert[1]*norms_vert[kk][1] + dist_vert[2]*norms_vert[kk][2];
            const double factor = phi_sup[kk] * inv_4pi * area / 3.0;

            second_int += qi * phi_sup[kk] * dot * inv_r3 * inv_4pi * area / 3.0;
          }
        }
      }
    }

    energy_pol = 0.5 * constant_pol * first_int;
    energy_react = 0.5 * (second_int - first_int * constant_react);
  }



  auto reduce_double = [&] (double &x) {
    MPI_Reduce (rank == 0 ? MPI_IN_PLACE : &x, &x, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
  };

  auto reduce_vec = [&] (std::vector<double> &v) {
    MPI_Reduce (rank == 0 ? MPI_IN_PLACE : v.data(),
                rank == 0 ? v.data() : nullptr,
                (int)v.size(), MPI_DOUBLE, MPI_SUM, 0, mpicomm);
  };

  reduce_double (charge_pol);
  reduce_double (energy_pol);
  reduce_double (energy_react);

  if (rank == 0) {
    constexpr int label_width = 50;
    constexpr int precision = 16;

    std::cout << std::left << std::setw (label_width) << "  Net charge [e]:"
              << std::setprecision (precision) << net_charge << "\n";

    std::cout << std::left << std::setw (label_width) << "  Flux charge [e]:"
              << std::setprecision (precision) << charge_pol / (4.0 * pi) << "\n";

    std::cout << std::left << std::setw (label_width) << "  Polarization energy [kT]:"
              << std::setprecision (precision) << energy_pol << "\n";

    if (calc_energy == 2) {
      std::cout << std::left << std::setw (label_width) << "  Direct ionic energy [kT]:"
                << std::setprecision (precision) << energy_react << "\n";
    }

    if (calc_coulombic == 1) {
      std::cout << std::left << std::setw (label_width) << "  Coulombic energy [kT]:"
                << std::setprecision (precision) << coul_energy << "\n";
    }

    std::cout << std::left << std::setw (label_width) << "  Sum of electrostatic energy contributions [kT]:"
              << std::setprecision (precision)
              << (energy_pol + energy_react + coul_energy) << "\n";

    std::cout << "===========================================================\n";
  }
}

void
poisson_boltzmann::energy_fast (ray_cache_t & ray_cache)
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);

  if (rank == 0)
    std::cout << "\n================ [ Electrostatic Energy ] =================\n";

  // ===========================
  // Costanti fisiche e scalari
  // ===========================
  const double inv_4pi = 1.0 / (4.0 * pi);
  const double eps0 = e_0; // Permittività del vuoto
  const double eps_in = 4.0 * pi * eps0 * e_in * kb * T * Angs / (e * e);
  const double eps_out = 4.0 * pi * eps0 * e_out * kb * T * Angs / (e * e);

  const double C0 = 1.0e3 * N_av * ionic_strength; // [mol/m^3]
  const double k2 = 2.0 * C0 * Angs * Angs * e * e / (eps0 * e_out * kb * T);
  const double k = std::sqrt (k2);

  const double den_in = 1.0 / eps_in;
  const double constant_pol = (1.0 / eps_out - 1.0 / eps_in) * inv_4pi;
  const double constant_react = (1.0 / eps_out) * inv_4pi;

  // ===========================
  // Variabili locali
  // ===========================
  // double energy_pol = 0.0;
  // double energy_react = 0.0;
  // double coul_energy = 0.0;
  double charge_pol = 0.0;

  // --- Filtra gli atomi caricati ---
  std::vector<double> charge_atoms_tmp;
  std::vector<std::array<double, 3>> pos_atoms_tmp;
  // charge_atoms_tmp.reserve(charge_atoms.size());
  // pos_atoms_tmp.reserve(pos_atoms.size());

  for (size_t ii = 0; ii < charge_atoms.size(); ++ii) {
    if (std::fabs (charge_atoms[ii]) > 1.e-5) {
      charge_atoms_tmp.push_back (charge_atoms[ii]);
      pos_atoms_tmp.push_back (pos_atoms[ii]);
    }
  }

  const size_t num_atoms = charge_atoms_tmp.size();
  std::vector<double>().swap (charge_atoms);
  std::vector<std::array<double, 3>>().swap (pos_atoms);

  // --- Energia Coulombiana ---
  if (calc_coulombic == 1) {
    for (size_t i = 0; i < num_atoms; ++i) {
      const double qi = charge_atoms_tmp[i];
      const std::array<double,3>& ri = pos_atoms_tmp[i];

      for (size_t j = i + 1; j < num_atoms; ++j) {
        const std::array<double,3>& rj = pos_atoms_tmp[j];
        const double dx = ri[0] - rj[0];
        const double dy = ri[1] - rj[1];
        const double dz = ri[2] - rj[2];
        const double r = std::sqrt (dx * dx + dy * dy + dz * dz);
        this->coul_energy += qi * charge_atoms_tmp[j] / r;
      }
    }

    this->coul_energy *= den_in;
  }

  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  double first_int = 0.0, second_int = 0.0;

  std::array<double,3> h{0}, area_h{0};
  std::array<double,3> V, N;
  std::array<double,8> tmp_eps, tmp_phi;
  std::vector<int> edg, fl_dir;

  auto quadrant = this->tmsh.begin_quadrant_sweep ();

  if (!border_quad.empty()) {
    quadrant[border_quad[0]];

    for (int d = 0; d < 3; ++d)
      h[d] = quadrant->p (d, 7) - quadrant->p (d, 0);

    area_h = {h[1]*h[2]/h[0]*0.25, h[0]*h[2]/h[1]*0.25, h[0]*h[1]/h[2]*0.25};
  }

  // flux and polarization energy calculation
  if (calc_energy==1 || (calc_energy == 2 && k < 1.e-5)) {
    for (const int ii : border_quad) {
      quadrant[ii];
      std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux_fast (quadrant, tmp_phi, tmp_eps);

      for (int ip = 0; ip < edg.size (); ++ip) {
        const int edge = edg[ip];
        const int axis = edge_axis[edge];
        const int i1 = edge2nodes[2 * edge];
        const int i2 = edge2nodes[2 * edge + 1];

        double fract = 0.0;
        normal_intersection (quadrant, ray_cache, edge, N, fract);

        V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
        V[axis] += fract * h[axis];

        const double tmp_flux =
          - (tmp_phi[i2] - tmp_phi[i1]) * wha (tmp_eps[i1], tmp_eps[i2], fract)
          * fl_dir[ip] * area_h[axis];

        charge_pol += tmp_flux;

        for (size_t ia = 0; ia < num_atoms; ++ia) {
          const std::array<double,3> &ra = pos_atoms_tmp[ia];
          const double dx = ra[0] - V[0];
          const double dy = ra[1] - V[1];
          const double dz = ra[2] - V[2];
          const double r = std::sqrt (dx * dx + dy * dy + dz * dz);
          const double inv_r3 = 1.0 / (r * r * r);
          const double qflux = tmp_flux / r;

          first_int += charge_atoms_tmp[ia] * qflux;
        }
      }
    }

    this->energy_pol = 0.5*constant_pol*first_int;
  }

  //polarization energy + ionic energy
  if (calc_energy==2 && k > 1.e-5) {
    int cubeindex = -1;
    std::array<std::array<double,3>,3> vert_triangles, norms_vert;
    std::array<double,3> dist_vert, phi_sup;
    int ntriang = 0;

    for (const int ii : border_quad) {
      quadrant[ii];
      cubeindex = classifyCube_fast (quadrant, eps_out);
      std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux_fast (quadrant, tmp_phi, tmp_eps);
      ntriang = getTriangles (cubeindex, triangles);

      // --- flussi
      for (int ip = 0; ip < edg.size (); ++ip) {
        const int edge = edg[ip];
        const int axis = edge_axis[edge];
        const int i1 = edge2nodes[2 * edge];
        const int i2 = edge2nodes[2 * edge + 1];

        double fract = 0.0;
        normal_intersection (quadrant, ray_cache, edge, N, fract);

        V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
        V[axis] += fract * h[axis];

        const double tmp_flux =
          - (tmp_phi[i2] - tmp_phi[i1]) * wha (tmp_eps[i1], tmp_eps[i2], fract)
          * fl_dir[ip] * area_h[axis];

        charge_pol += tmp_flux;

        for (size_t ia = 0; ia < num_atoms; ++ia) {
          const std::array<double,3> &ra = pos_atoms_tmp[ia];
          const double dx = ra[0] - V[0];
          const double dy = ra[1] - V[1];
          const double dz = ra[2] - V[2];
          const double r = std::sqrt (dx * dx + dy * dy + dz * dz);
          const double inv_r3 = 1.0 / (r * r * r);

          first_int += charge_atoms_tmp[ia] * tmp_flux / r;
        }
      }

      // --- triangoli (componente ionica)
      for (int itri = 0; itri < ntriang; ++itri) {
        for (int jj = 0; jj < 3; ++jj) {
          const int edge = triangles[itri][jj];
          const int axis = edge_axis[edge];
          const int i1 = edge2nodes[2 * edge];
          const int i2 = edge2nodes[2 * edge + 1];

          double fract = 0.0;
          normal_intersection (quadrant, ray_cache, edge, N, fract);

          V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
          V[axis] += fract * h[axis];

          vert_triangles[jj] = V;
          norms_vert[jj] = N;

          phi_sup[jj] = phi0 (tmp_eps[i1], tmp_eps[i2], tmp_phi[i1], tmp_phi[i2], fract);
        }

        const double area = areaTriangle (vert_triangles);

        for (size_t ia = 0; ia < num_atoms; ++ia) {
          const double qi = charge_atoms_tmp[ia];
          const std::array<double,3> &ra = pos_atoms_tmp[ia];

          for (int kk = 0; kk < 3; ++kk) {
            dist_vert = {vert_triangles[kk][0] - ra[0],
                         vert_triangles[kk][1] - ra[1],
                         vert_triangles[kk][2] - ra[2]
                        };
            const double r2 = dist_vert[0]*dist_vert[0] + dist_vert[1]*dist_vert[1] + dist_vert[2]*dist_vert[2];
            const double r = std::sqrt (r2);
            const double inv_r3 = 1.0 / (r2 * r);
            const double inv_r5 = inv_r3 / r2;
            const double dot = dist_vert[0]*norms_vert[kk][0] + dist_vert[1]*norms_vert[kk][1] + dist_vert[2]*norms_vert[kk][2];
            const double factor = phi_sup[kk] * inv_4pi * area / 3.0;

            second_int += qi * phi_sup[kk] * dot * inv_r3 * inv_4pi * area / 3.0;
          }
        }
      }
    }

    this->energy_pol = 0.5 * constant_pol * first_int;
    this->energy_react = 0.5 * (second_int - first_int * constant_react);
  }


  auto reduce_double = [&] (double &x) {
    MPI_Reduce (rank == 0 ? MPI_IN_PLACE : &x, &x, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
  };

  auto reduce_vec = [&] (std::vector<double> &v) {
    MPI_Reduce (rank == 0 ? MPI_IN_PLACE : v.data(),
                rank == 0 ? v.data() : nullptr,
                (int)v.size(), MPI_DOUBLE, MPI_SUM, 0, mpicomm);
  };

  reduce_double (charge_pol);
  reduce_double (energy_pol);
  reduce_double (energy_react);

  // Print the result
  if (rank == 0) {
    constexpr int label_width = 50;
    constexpr int precision = 16;

    std::cout << std::left << std::setw (label_width) << "  Net charge [e]:"
              << std::setprecision (precision) << net_charge << "\n";

    std::cout << std::left << std::setw (label_width) << "  Flux charge [e]:"
              << std::setprecision (precision) << charge_pol / (4.0 * pi) << "\n";

    // std::cout << std::left << std::setw(label_width)
    // << "    Error w.r.t. net charge [%]:"
    // << std::setprecision(6)
    // << ((charge_pol / (4.0 * pi) - net_charge) / net_charge * 100.0) << "\n";

    std::cout << std::left << std::setw (label_width) << "  Polarization energy [kT]:"
              << std::setprecision (precision) << energy_pol << "\n";

    if (calc_energy == 2) {
      std::cout << std::left << std::setw (label_width) << "  Direct ionic energy [kT]:"
                << std::setprecision (precision) << energy_react << "\n";
    }

    if (calc_coulombic == 1) {
      std::cout << std::left << std::setw (label_width) << "  Coulombic energy [kT]:"
                << std::setprecision (precision) << coul_energy << "\n";
    }

    std::cout << std::left << std::setw (label_width) << "  Sum of electrostatic energy contributions [kT]:"
              << std::setprecision (precision)
              << (energy_pol + energy_react + coul_energy) << "\n";

    std::cout << "===========================================================\n";
  }
}

void
poisson_boltzmann::write_potential_on_surface (ray_cache_t & ray_cache)
{
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  const double eps_in = 4.0 * pi * e_0 * e_in * kb * T * Angs / (e * e);
  const double eps_out = 4.0 * pi * e_0 * e_out * kb * T * Angs / (e * e);
  const double C_0 = 1.0e3 * N_av * ionic_strength;
  const double k2 = 2.0 * C_0 * Angs * Angs * e * e / (e_0 * e_out * kb * T);
  const double k = std::sqrt (k2);

  std::array<double, 3> V, N, h, area_h;
  std::array<double, 8> tmp_eps, tmp_phi;
  std::vector<int> edg, fl_dir;
  std::array<std::array<double, 3>, 3> vert_triangles, norms_vert;
  std::array<double, 3> phi_sup;

  int cubeindex = -1;
  int edge, i1 = 0, i2 = 0, ntriang = 0;
  double fract;
  double tmp_phi_1 = 0.0, tmp_phi_2 = 0.0;
  double tmp_eps_1 = 0.0, tmp_eps_2 = 0.0;

  auto quadrant = this->tmsh.begin_quadrant_sweep();

  // File setup
  // std::string filename_nodes = "phi_nodes_" + pqrfilename + ".txt";
  // std::string filename_surf  = "phi_surf_" + pqrfilename + ".txt";

  // std::ofstream phi_nodes_txt(filename_nodes);
  // std::ofstream phi_surf_txt(filename_surf);
  // FILE* phi_nod_delphi = std::fopen("filename_nodes_delphi.txt", "w");
  // FILE* phi_sup_delphi = std::fopen("filename_sup_delphi.txt", "w");
  std::vector<std::string> phi_nodes_local;
  std::vector<std::string> phi_surf_local;
  // std::vector<std::string> phi_nodes_delphi_local;
  // std::vector<std::string> phi_sup_delphi_local;

  for (const int ii : border_quad) {
    quadrant[ii];
    cubeindex = classifyCube (quadrant, eps_out);

    // Compute edge lengths and area scale factors
    h[0] = quadrant->p (0, 7) - quadrant->p (0, 0);
    h[1] = quadrant->p (1, 7) - quadrant->p (1, 0);
    h[2] = quadrant->p (2, 7) - quadrant->p (2, 0);
    area_h[0] = h[1] * h[2] / h[0] * 0.25;
    area_h[1] = h[0] * h[2] / h[1] * 0.25;
    area_h[2] = h[0] * h[1] / h[2] * 0.25;

    std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux (quadrant, tmp_phi, tmp_eps);
    ntriang = getTriangles (cubeindex, triangles);

    for (int t = 0; t < ntriang; ++t) {
      for (int j = 0; j < 3; ++j) {
        edge = triangles[t][j];
        i1 = edge2nodes[2 * edge];
        i2 = edge2nodes[2 * edge + 1];

        // Intersection and vertex
        V[0] = quadrant->p (0, i1);
        V[1] = quadrant->p (1, i1);
        V[2] = quadrant->p (2, i1);

        normal_intersection (quadrant, ray_cache, edge, N, fract);
        V[edge_axis[edge]] += fract * h[edge_axis[edge]];

        vert_triangles[j] = V;
        norms_vert[j] = N;

        // Interpolate phi/epsilon at i1
        if (!quadrant->is_hanging (i1)) {
          tmp_phi_1 = (*phi)[quadrant->gt (i1)];
          tmp_eps_1 = (*epsilon_nodes)[quadrant->gt (i1)];
        } else {
          tmp_phi_1 = tmp_eps_1 = 0.0;
          int np = quadrant->num_parents (i1);

          for (int k = 0; k < np; ++k) {
            tmp_phi_1 += (*phi)[quadrant->gparent (k, i1)] / np;
            tmp_eps_1 += (*epsilon_nodes)[quadrant->gparent (k, i1)] / np;
          }
        }

        // Interpolate phi/epsilon at i2
        if (!quadrant->is_hanging (i2)) {
          tmp_phi_2 = (*phi)[quadrant->gt (i2)];
          tmp_eps_2 = (*epsilon_nodes)[quadrant->gt (i2)];
        } else {
          tmp_phi_2 = tmp_eps_2 = 0.0;
          int np = quadrant->num_parents (i2);

          for (int k = 0; k < np; ++k) {
            tmp_phi_2 += (*phi)[quadrant->gparent (k, i2)] / np;
            tmp_eps_2 += (*epsilon_nodes)[quadrant->gparent (k, i2)] / np;
          }
        }

        // Interpolate phi on surface
        phi_sup[j] = phi0 (tmp_eps_1, tmp_eps_2, tmp_phi_1, tmp_phi_2, fract);

        std::ostringstream oss;
        oss << std::scientific << std::setprecision (5)
            << quadrant->p (0, i1) << " "
            << quadrant->p (1, i1) << " "
            << quadrant->p (2, i1) << " "
            << tmp_phi_1 << "\n"
            << quadrant->p (0, i2) << " "
            << quadrant->p (1, i2) << " "
            << quadrant->p (2, i2) << " "
            << tmp_phi_2;
        phi_nodes_local.push_back (oss.str());
        oss.str ("");
        oss.clear();
        oss << std::scientific << std::setprecision (5)
            << V[0] << " "
            << V[1] << " "
            << V[2] << " "
            << phi_sup[j];
        phi_surf_local.push_back (oss.str());
        oss.str ("");
        oss.clear();
        // // Write to ASCII
        // phi_nodes_txt << quadrant->p(0, i1) << "  " << quadrant->p(1, i1) << "  " << quadrant->p(2, i1) << "  " << tmp_phi_1 << "\n";
        // phi_nodes_txt << quadrant->p(0, i2) << "  " << quadrant->p(1, i2) << "  " << quadrant->p(2, i2) << "  " << tmp_phi_2 << "\n";

        // phi_surf_txt  << V[0] << "  " << V[1] << "  " << V[2] << "  " << phi_sup[j] << "\n";

        // // Write to Delphi PDB-like format
        // std::fprintf(phi_nod_delphi,
        // "\nATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%8.4f%8.4f",
        // 1, "X", "XXX", " ", 0,
        // quadrant->p(0, i1), quadrant->p(1, i1), quadrant->p(2, i1), tmp_phi_1, tmp_phi_2);

        // std::fprintf(phi_nod_delphi,
        // "\nATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%8.4f%8.4f",
        // 1, "X", "XXX", " ", 0,
        // quadrant->p(0, i2), quadrant->p(1, i2), quadrant->p(2, i2), tmp_phi_1, tmp_phi_2);

        // std::fprintf(phi_sup_delphi,
        // "\nATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%8.4f%8.4f",
        // 1, "X", "XXX", " ", 0,
        // V[0], V[1], V[2], phi_sup[j], 0.0);
      }
    }
  }

  auto gather_and_write = [&] (const std::string& filename,
  const std::vector<std::string>& local_lines) {
    if (rank == 0) {
      std::ofstream ofs (filename);

      for (const auto& line : local_lines)
        ofs << line << "\n";

      for (int r = 1; r < size; ++r) {
        int n_lines;
        MPI_Recv (&n_lines, 1, MPI_INT, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int i = 0; i < n_lines; ++i) {
          char buf[512];
          MPI_Recv (buf, 512, MPI_CHAR, r, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          ofs << buf << "\n";
        }
      }

      ofs.close();
    } else {
      int n_lines = static_cast<int> (local_lines.size());
      MPI_Send (&n_lines, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

      for (const auto& line : local_lines) {
        MPI_Send (line.c_str(), static_cast<int> (line.size()) + 1,
                  MPI_CHAR, 0, 1, MPI_COMM_WORLD);
      }
    }
  };

  gather_and_write ("phi_nodes.txt", phi_nodes_local);
  gather_and_write ("phi_surf.txt", phi_surf_local);
  // Close files
  // phi_nodes_txt.close();
  // phi_surf_txt.close();
  // std::fclose(phi_nod_delphi);
  // std::fclose(phi_sup_delphi);
}

double
poisson_boltzmann::analytic_solution (double x, double y, double z)
{
  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/ (e*e); //adim e_out
  double eps_in = 4.0*pi*e_0*e_in*kb*T*Angs/ (e*e); //adim e_in
  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/ (e_0*e_out*kb*T);
  double k = std::sqrt (k2);
  double rs = 2.0;
  double pippo_out = eps_out* (1.0+k*rs);
  double pippo_in = eps_in*rs;
  double dist= 0.0;
  double pot = 0.0;

  dist = std::hypot (x, y, z);

  if (dist <= rs)
    pot = 1.0/ (pippo_out*rs) + (rs-dist)/ (pippo_in*dist);
  else
    pot = exp (k* (rs-dist))/ (dist*pippo_out);

  return pot;
}



void
poisson_boltzmann::analitic_potential ()
{
  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/ (e*e); //adim e_out
  double eps_in = 4.0*pi*e_0*e_in*kb*T*Angs/ (e*e); //adim e_in
  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/ (e_0*e_out*kb*T);
  double k = std::sqrt (k2);
  double rs = 2.0;
  double pippo_out = eps_out* (1.0+k*rs);
  double pippo_in = eps_in*rs;
  double dist= 0.0;


  distributed_vector phi_an (tmsh.num_owned_nodes (), mpicomm);
  // distributed_vector field_an (tmsh.num_owned_nodes (), mpicomm);
  bim3a_solution_with_ghosts (tmsh, phi_an, replace_op);
  // bim3a_solution_with_ghosts (tmsh, field_an, replace_op);

  for (auto quadrant = this->tmsh.begin_quadrant_sweep ();
       quadrant != this->tmsh.end_quadrant_sweep ();
       ++quadrant) {
    for (const NS::Atom& i : atoms) {
      if (std::fabs (i.charge) > std::numeric_limits<double>::epsilon ()) {
        for (int ii = 0; ii < 8; ++ii) {
          if (! quadrant->is_hanging (ii)) {
            dist = std::hypot ( (i.pos[0] - quadrant->p (0,ii)),
                                (i.pos[1] - quadrant->p (1,ii)),
                                (i.pos[2] - quadrant->p (2,ii)));

            if (dist <= rs) {
              phi_an[quadrant->gt (ii)] = 1.0/ (pippo_out*rs) +
                                          i.charge* (rs-dist)/ (pippo_in*dist);
              // field_an[quadrant->gt (ii)] = (1.0 + (10.0-dist)/dist)*
              // i.charge/(pippo_in*dist);
            } else
              phi_an[quadrant->gt (ii)] = i.charge*exp (k* (rs-dist))/ (dist*pippo_out);

            // field_an[quadrant->gt (ii)] = phi_an[quadrant->gt (ii)]*(k +1.0/dist);
          }

          // else{
          // for (int jj = 0; jj < quadrant->num_parents (ii); ++jj){
          // phi_an[quadrant->gparent (jj, ii)] += 0.0;
          // field_an[quadrant->gparent (jj, ii)] += 0.0;
          // }
          // }
        }
      }
    }
  }

  phi_an.assemble (replace_op);
  tmsh.octbin_export ("phi_an_0", phi_an);
  // field_an.assemble (replace_op);
  // tmsh.octbin_export ("field_an_0", field_an);
}

void
poisson_boltzmann::pot_field_fast (ray_cache_t & ray_cache)
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);

  if (rank == 0)
    std::cout << "\n================ [ Calculating Potential & Field Components ] =================\n";

  // ===========================
  // Costanti fisiche e scalari
  // ===========================
  const double inv_4pi = 1.0 / (4.0 * pi);
  const double eps0 = e_0; // Permittività del vuoto
  const double eps_in = 4.0 * pi * eps0 * e_in * kb * T * Angs / (e * e);
  const double eps_out = 4.0 * pi * eps0 * e_out * kb * T * Angs / (e * e);

  const double C0 = 1.0e3 * N_av * ionic_strength; // [mol/m^3]
  const double k2 = 2.0 * C0 * Angs * Angs * e * e / (eps0 * e_out * kb * T);
  const double k = std::sqrt (k2);

  const double den_in = 1.0 / eps_in;
  const double constant_pol = (1.0 / eps_out - 1.0 / eps_in) * inv_4pi;
  const double constant_react = (1.0 / eps_out) * inv_4pi;

  // ===========================
  // Variabili locali
  // ===========================
  // double energy_pol = 0.0;
  // double energy_react = 0.0;
  // double coul_energy = 0.0;
  double charge_pol = 0.0;
  const size_t num_atoms = charge_atoms.size();

  std::vector<double> phi_c, phi_p, phi_i;
  std::vector<double> field_cx, field_px, field_ix;
  std::vector<double> field_cy, field_py, field_iy;
  std::vector<double> field_cz, field_pz, field_iz;

  phi_c.assign (num_atoms, 0.0);
  field_cx.assign (num_atoms, 0.0);
  field_cy.assign (num_atoms, 0.0);
  field_cz.assign (num_atoms, 0.0);

  for (size_t i = 0; i < num_atoms; ++i) {
    const std::array<double,3> &ri = pos_atoms[i];
    const double qi = charge_atoms[i];

    for (size_t j = i + 1; j < num_atoms; ++j) {
      const std::array<double,3> &rj = pos_atoms[j];
      const double qj = charge_atoms[j];

      const double dx = ri[0] - rj[0];
      const double dy = ri[1] - rj[1];
      const double dz = ri[2] - rj[2];
      const double r2 = dx * dx + dy * dy + dz * dz;
      const double r = std::sqrt (r2);
      const double inv_r3 = 1.0 / (r2 * r);

      // Coulomb energy
      this->coul_energy += (qi * qj) / r * den_in;
      // Coulomb potential
      phi_c[i] += qj / r * den_in;
      phi_c[j] += qi / r * den_in;


      // Coulomb field
      std::array<double, 3> eij = {
        dx * inv_r3 * qj * den_in,
        dy * inv_r3 * qj * den_in,
        dz * inv_r3 * qj * den_in
      };
      std::array<double, 3> eji = {
        -dx * inv_r3 * qi * den_in,
        -dy * inv_r3 * qi * den_in,
        -dz * inv_r3 * qi * den_in
      };

      field_cx[i] += eij[0];
      field_cx[j] += eji[0];
      field_cy[i] += eij[1];
      field_cy[j] += eji[1];
      field_cz[i] += eij[2];
      field_cz[j] += eji[2];
    }
  }

  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  double first_int = 0.0, second_int = 0.0, distance = 0.0;

  std::array<double,3> h{0}, area_h{0};
  std::array<double,3> V, N;
  std::array<double,8> tmp_eps, tmp_phi;
  std::vector<int> edg, fl_dir;

  auto quadrant = this->tmsh.begin_quadrant_sweep ();

  if (!border_quad.empty()) {
    quadrant[border_quad[0]];

    for (int d = 0; d < 3; ++d)
      h[d] = quadrant->p (d, 7) - quadrant->p (d, 0);

    area_h = {h[1]*h[2]/h[0]*0.25, h[0]*h[2]/h[1]*0.25, h[0]*h[1]/h[2]*0.25};
  }

  // flux and polarization energy calculation
  if ((calc_field_term==1 || (calc_field_term == 2 && k < 1.e-5)) && calc_potential_term < 2) {
    phi_p.assign (num_atoms, 0.0);
    field_px.assign (num_atoms, 0.0);
    field_py.assign (num_atoms, 0.0);
    field_pz.assign (num_atoms, 0.0);

    for (const int ii : border_quad) {
      quadrant[ii];
      std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux_fast (quadrant, tmp_phi, tmp_eps);

      for (int ip = 0; ip < edg.size (); ++ip) {
        const int edge = edg[ip];
        const int axis = edge_axis[edge];
        const int i1 = edge2nodes[2 * edge];
        const int i2 = edge2nodes[2 * edge + 1];
        double fract = 0.0;
        normal_intersection (quadrant, ray_cache, edge, N, fract);
        V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
        V[axis] += fract * h[axis];
        const double tmp_flux =
          - (tmp_phi[i2] - tmp_phi[i1]) * wha (tmp_eps[i1], tmp_eps[i2], fract)
          * fl_dir[ip] * area_h[axis];

        charge_pol += tmp_flux;

        for (size_t ia = 0; ia < num_atoms; ++ia) {
          const std::array<double,3> &ra = pos_atoms[ia];
          const double dx = ra[0] - V[0];
          const double dy = ra[1] - V[1];
          const double dz = ra[2] - V[2];
          const double r = std::sqrt (dx * dx + dy * dy + dz * dz);
          const double inv_r3 = 1.0 / (r * r * r);
          const double qflux = tmp_flux / r;

          first_int += charge_atoms[ia] * qflux;
          phi_p[ia] += qflux * constant_pol;
          field_px[ia] += dx * inv_r3 * tmp_flux * constant_pol;
          field_py[ia] += dy * inv_r3 * tmp_flux * constant_pol;
          field_pz[ia] += dz * inv_r3 * tmp_flux * constant_pol;
        }
      }
    }

    this->energy_pol = 0.5*constant_pol*first_int;
  } else if ((calc_potential_term == 1 || (calc_potential_term == 2 && k < 1.e-5)) && calc_field_term == 0) {
    phi_p.assign (num_atoms, 0.0);

    for (const int ii : border_quad) {
      quadrant[ii];
      std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux_fast (quadrant, tmp_phi, tmp_eps);

      for (int ip = 0; ip < edg.size (); ++ip) {
        const int edge = edg[ip];
        const int axis = edge_axis[edge];
        const int i1 = edge2nodes[2 * edge];
        const int i2 = edge2nodes[2 * edge + 1];
        double fract = 0.0;
        normal_intersection (quadrant, ray_cache, edge, N, fract);
        V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
        V[axis] += fract * h[axis];
        const double tmp_flux =
          - (tmp_phi[i2] - tmp_phi[i1]) * wha (tmp_eps[i1], tmp_eps[i2], fract)
          * fl_dir[ip] * area_h[axis];

        charge_pol += tmp_flux;

        for (size_t ia = 0; ia < num_atoms; ++ia) {
          const std::array<double,3> &ra = pos_atoms[ia];
          const double dx = ra[0] - V[0];
          const double dy = ra[1] - V[1];
          const double dz = ra[2] - V[2];
          const double r = std::sqrt (dx * dx + dy * dy + dz * dz);
          const double inv_r3 = 1.0 / (r * r * r);
          const double qflux = tmp_flux / r;

          first_int += charge_atoms[ia] * qflux;
          phi_p[ia] += qflux * constant_pol;

        }
      }
    }

    this->energy_pol = 0.5*constant_pol*first_int;
  }

  auto allocate_potential_fields = [&] (void) {
    phi_p.assign (num_atoms, 0.0);
    phi_i.assign (num_atoms, 0.0);
    field_px.assign (num_atoms, 0.0);
    field_py.assign (num_atoms, 0.0);
    field_pz.assign (num_atoms, 0.0);
    field_ix.assign (num_atoms, 0.0);
    field_iy.assign (num_atoms, 0.0);
    field_iz.assign (num_atoms, 0.0);
  };

  //polarization energy + ionic energy
  if (calc_field_term==2 && k > 1.e-5) {
    allocate_potential_fields();
    int cubeindex = -1;
    std::array<std::array<double,3>,3> vert_triangles, norms_vert;
    std::array<double,3> dist_vert, phi_sup;
    int ntriang = 0;

    for (const int ii : border_quad) {
      quadrant[ii];
      cubeindex = classifyCube_fast (quadrant, eps_out);
      std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux_fast (quadrant, tmp_phi, tmp_eps);
      ntriang = getTriangles (cubeindex, triangles);

      // --- flussi
      for (int ip = 0; ip < edg.size (); ++ip) {
        const int edge = edg[ip];
        const int axis = edge_axis[edge];
        const int i1 = edge2nodes[2 * edge];
        const int i2 = edge2nodes[2 * edge + 1];

        double fract = 0.0;
        normal_intersection (quadrant, ray_cache, edge, N, fract);

        V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
        V[axis] += fract * h[axis];

        const double tmp_flux =
          - (tmp_phi[i2] - tmp_phi[i1]) * wha (tmp_eps[i1], tmp_eps[i2], fract)
          * fl_dir[ip] * area_h[axis];
        charge_pol += tmp_flux;

        for (size_t ia = 0; ia < num_atoms; ++ia) {
          const std::array<double,3> &ra = pos_atoms[ia];
          const double dx = ra[0] - V[0];
          const double dy = ra[1] - V[1];
          const double dz = ra[2] - V[2];
          const double r = std::sqrt (dx * dx + dy * dy + dz * dz);
          const double inv_r3 = 1.0 / (r * r * r);

          first_int += charge_atoms[ia] * tmp_flux / r;
          phi_p[ia] += tmp_flux / r * constant_pol;
          field_px[ia] += dx * inv_r3 * tmp_flux * constant_pol;
          field_py[ia] += dy * inv_r3 * tmp_flux * constant_pol;
          field_pz[ia] += dz * inv_r3 * tmp_flux * constant_pol;
        }
      }

      // --- triangoli (componente ionica)
      for (int itri = 0; itri < ntriang; ++itri) {
        for (int jj = 0; jj < 3; ++jj) {
          const int edge = triangles[itri][jj];
          const int axis = edge_axis[edge];
          const int i1 = edge2nodes[2 * edge];
          const int i2 = edge2nodes[2 * edge + 1];

          double fract = 0.0;
          normal_intersection (quadrant, ray_cache, edge, N, fract);

          V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
          V[axis] += fract * h[axis];

          vert_triangles[jj] = V;
          norms_vert[jj] = N;

          phi_sup[jj] = phi0 (tmp_eps[i1], tmp_eps[i2], tmp_phi[i1], tmp_phi[i2], fract);
        }

        const double area = areaTriangle (vert_triangles);

        for (size_t ia = 0; ia < num_atoms; ++ia) {
          const double qi = charge_atoms[ia];
          const std::array<double,3> &ra = pos_atoms[ia];

          for (int kk = 0; kk < 3; ++kk) {
            dist_vert = {vert_triangles[kk][0] - ra[0],
                         vert_triangles[kk][1] - ra[1],
                         vert_triangles[kk][2] - ra[2]
                        };
            const double r2 = dist_vert[0]*dist_vert[0] + dist_vert[1]*dist_vert[1] + dist_vert[2]*dist_vert[2];
            const double r = std::sqrt (r2);
            const double inv_r3 = 1.0 / (r2 * r);
            const double inv_r5 = inv_r3 / r2;
            const double dot = dist_vert[0]*norms_vert[kk][0] + dist_vert[1]*norms_vert[kk][1] + dist_vert[2]*norms_vert[kk][2];
            const double factor = phi_sup[kk] * inv_4pi * area / 3.0;

            second_int += qi * phi_sup[kk] * dot * inv_r3 * inv_4pi * area / 3.0;
            phi_i[ia] += phi_sup[kk] * dot * inv_r3 * inv_4pi * area / 3.0;

            field_ix[ia] += factor * (-3 * dist_vert[0] * inv_r5 * dot + norms_vert[kk][0] * inv_r3);
            field_iy[ia] += factor * (-3 * dist_vert[1] * inv_r5 * dot + norms_vert[kk][1] * inv_r3);
            field_iz[ia] += factor * (-3 * dist_vert[2] * inv_r5 * dot + norms_vert[kk][2] * inv_r3);
          }
        }
      }
    }

    this->energy_pol = 0.5 * constant_pol * first_int;
    this->energy_react = 0.5 * (second_int - first_int * constant_react);
  } else if ((calc_potential_term==2 && k > 1.e-5) && calc_field_term == 0) {
    phi_p.assign (num_atoms, 0.0);
    phi_i.assign (num_atoms, 0.0);
    int cubeindex = -1;
    std::array<std::array<double,3>,3> vert_triangles, norms_vert;
    std::array<double,3> dist_vert, phi_sup;
    int ntriang = 0;

    for (const int ii : border_quad) {
      quadrant[ii];
      cubeindex = classifyCube_fast (quadrant, eps_out);
      std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux_fast (quadrant, tmp_phi, tmp_eps);
      ntriang = getTriangles (cubeindex, triangles);

      // --- flussi
      for (int ip = 0; ip < edg.size (); ++ip) {
        const int edge = edg[ip];
        const int axis = edge_axis[edge];
        const int i1 = edge2nodes[2 * edge];
        const int i2 = edge2nodes[2 * edge + 1];

        double fract = 0.0;
        normal_intersection (quadrant, ray_cache, edge, N, fract);

        V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
        V[axis] += fract * h[axis];

        const double tmp_flux =
          - (tmp_phi[i2] - tmp_phi[i1]) * wha (tmp_eps[i1], tmp_eps[i2], fract)
          * fl_dir[ip] * area_h[axis];
        charge_pol += tmp_flux;

        for (size_t ia = 0; ia < num_atoms; ++ia) {
          const std::array<double,3> &ra = pos_atoms[ia];
          const double dx = ra[0] - V[0];
          const double dy = ra[1] - V[1];
          const double dz = ra[2] - V[2];
          const double r = std::sqrt (dx * dx + dy * dy + dz * dz);
          const double inv_r3 = 1.0 / (r * r * r);
          first_int += charge_atoms[ia] * tmp_flux / r;
          phi_p[ia] += tmp_flux / r * constant_pol;
        }
      }

      // --- triangoli (componente ionica)
      for (int itri = 0; itri < ntriang; ++itri) {
        for (int jj = 0; jj < 3; ++jj) {
          const int edge = triangles[itri][jj];
          const int axis = edge_axis[edge];
          const int i1 = edge2nodes[2 * edge];
          const int i2 = edge2nodes[2 * edge + 1];

          double fract = 0.0;
          normal_intersection (quadrant, ray_cache, edge, N, fract);

          V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
          V[axis] += fract * h[axis];

          vert_triangles[jj] = V;
          norms_vert[jj] = N;

          phi_sup[jj] = phi0 (tmp_eps[i1], tmp_eps[i2], tmp_phi[i1], tmp_phi[i2], fract);
        }

        const double area = areaTriangle (vert_triangles);

        for (size_t ia = 0; ia < num_atoms; ++ia) {
          const double qi = charge_atoms[ia];
          const std::array<double,3> &ra = pos_atoms[ia];

          for (int kk = 0; kk < 3; ++kk) {
            dist_vert = {vert_triangles[kk][0] - ra[0],
                         vert_triangles[kk][1] - ra[1],
                         vert_triangles[kk][2] - ra[2]
                        };
            const double r2 = dist_vert[0]*dist_vert[0] + dist_vert[1]*dist_vert[1] + dist_vert[2]*dist_vert[2];
            const double r = std::sqrt (r2);
            const double inv_r3 = 1.0 / (r2 * r);
            const double inv_r5 = inv_r3 / r2;
            const double dot = dist_vert[0]*norms_vert[kk][0] + dist_vert[1]*norms_vert[kk][1] + dist_vert[2]*norms_vert[kk][2];
            const double factor = phi_sup[kk] * inv_4pi * area / 3.0;
            second_int += qi * phi_sup[kk] * dot * inv_r3 * inv_4pi * area / 3.0;
            phi_i[ia] += phi_sup[kk] * dot * inv_r3 * inv_4pi * area / 3.0;
          }
        }
      }
    }

    this->energy_pol = 0.5 * constant_pol * first_int;
    this->energy_react = 0.5 * (second_int - first_int * constant_react);
  } else if ((calc_potential_term==2 && k > 1.e-5) && calc_field_term == 1) {
    phi_p.assign (num_atoms, 0.0);
    phi_i.assign (num_atoms, 0.0);
    field_px.assign (num_atoms, 0.0);
    field_py.assign (num_atoms, 0.0);
    field_pz.assign (num_atoms, 0.0);
    int cubeindex = -1;
    std::array<std::array<double,3>,3> vert_triangles, norms_vert;
    std::array<double,3> dist_vert, phi_sup;
    int ntriang = 0;

    for (const int ii : border_quad) {
      quadrant[ii];
      cubeindex = classifyCube_fast (quadrant, eps_out);
      std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux_fast (quadrant, tmp_phi, tmp_eps);
      ntriang = getTriangles (cubeindex, triangles);

      // --- flussi
      for (int ip = 0; ip < edg.size (); ++ip) {
        const int edge = edg[ip];
        const int axis = edge_axis[edge];
        const int i1 = edge2nodes[2 * edge];
        const int i2 = edge2nodes[2 * edge + 1];

        double fract = 0.0;
        normal_intersection (quadrant, ray_cache, edge, N, fract);

        V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
        V[axis] += fract * h[axis];

        const double tmp_flux =
          - (tmp_phi[i2] - tmp_phi[i1]) * wha (tmp_eps[i1], tmp_eps[i2], fract)
          * fl_dir[ip] * area_h[axis];
        charge_pol += tmp_flux;

        for (size_t ia = 0; ia < num_atoms; ++ia) {
          const std::array<double,3> &ra = pos_atoms[ia];
          const double dx = ra[0] - V[0];
          const double dy = ra[1] - V[1];
          const double dz = ra[2] - V[2];
          const double r = std::sqrt (dx * dx + dy * dy + dz * dz);
          const double inv_r3 = 1.0 / (r * r * r);
          first_int += charge_atoms[ia] * tmp_flux / r;
          phi_p[ia] += tmp_flux / r * constant_pol;
          field_px[ia] += dx * inv_r3 * tmp_flux * constant_pol;
          field_py[ia] += dy * inv_r3 * tmp_flux * constant_pol;
          field_pz[ia] += dz * inv_r3 * tmp_flux * constant_pol;
        }
      }

      // --- triangoli (componente ionica)
      for (int itri = 0; itri < ntriang; ++itri) {
        for (int jj = 0; jj < 3; ++jj) {
          const int edge = triangles[itri][jj];
          const int axis = edge_axis[edge];
          const int i1 = edge2nodes[2 * edge];
          const int i2 = edge2nodes[2 * edge + 1];

          double fract = 0.0;
          normal_intersection (quadrant, ray_cache, edge, N, fract);

          V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
          V[axis] += fract * h[axis];

          vert_triangles[jj] = V;
          norms_vert[jj] = N;

          phi_sup[jj] = phi0 (tmp_eps[i1], tmp_eps[i2], tmp_phi[i1], tmp_phi[i2], fract);
        }

        const double area = areaTriangle (vert_triangles);

        for (size_t ia = 0; ia < num_atoms; ++ia) {
          const double qi = charge_atoms[ia];
          const std::array<double,3> &ra = pos_atoms[ia];

          for (int kk = 0; kk < 3; ++kk) {
            dist_vert = {vert_triangles[kk][0] - ra[0],
                         vert_triangles[kk][1] - ra[1],
                         vert_triangles[kk][2] - ra[2]
                        };
            const double r2 = dist_vert[0]*dist_vert[0] + dist_vert[1]*dist_vert[1] + dist_vert[2]*dist_vert[2];
            const double r = std::sqrt (r2);
            const double inv_r3 = 1.0 / (r2 * r);
            const double inv_r5 = inv_r3 / r2;
            const double dot = dist_vert[0]*norms_vert[kk][0] + dist_vert[1]*norms_vert[kk][1] + dist_vert[2]*norms_vert[kk][2];
            const double factor = phi_sup[kk] * inv_4pi * area / 3.0;
            second_int += qi * phi_sup[kk] * dot * inv_r3 * inv_4pi * area / 3.0;
            phi_i[ia] += phi_sup[kk] * dot * inv_r3 * inv_4pi * area / 3.0;
          }
        }
      }
    }

    this->energy_pol = 0.5 * constant_pol * first_int;
    this->energy_react = 0.5 * (second_int - first_int * constant_react);
  }


  auto reduce_double = [&] (double &x) {
    MPI_Reduce (rank == 0 ? MPI_IN_PLACE : &x, &x, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
  };

  auto reduce_vec = [&] (std::vector<double> &v) {
    MPI_Reduce (rank == 0 ? MPI_IN_PLACE : v.data(),
                rank == 0 ? v.data() : nullptr,
                (int)v.size(), MPI_DOUBLE, MPI_SUM, 0, mpicomm);
  };

  reduce_double (charge_pol);
  reduce_double (energy_pol);
  reduce_double (energy_react);
  reduce_vec (phi_p);
  reduce_vec (phi_i);
  reduce_vec (field_px);
  reduce_vec (field_py);
  reduce_vec (field_pz);
  reduce_vec (field_ix);
  reduce_vec (field_iy);
  reduce_vec (field_iz);

  // Print the result
  if (rank == 0) {
    constexpr int label_width = 50;
    constexpr int precision = 16;


    if (calc_potential_term == 1 || calc_field_term == 1 || calc_potential_term == 2 || calc_field_term == 2) {
      std::cout << std::left << std::setw (label_width) << "  Polarization energy [kT]:"
                << std::setprecision (precision) << energy_pol << "\n";
    }

    if (calc_potential_term == 2 || calc_field_term == 2) {
      std::cout << std::left << std::setw (label_width) << "  Direct ionic energy [kT]:"
                << std::setprecision (precision) << energy_react << "\n";
    }


    std::cout << std::left << std::setw (label_width) << "  Coulombic energy [kT]:"
              << std::setprecision (precision) << coul_energy << "\n";


    std::cout << std::left << std::setw (label_width) << "  Sum of electrostatic energy contributions [kT]:"
              << std::setprecision (precision)
              << (energy_pol + energy_react + coul_energy) << "\n";

    std::cout << "===========================================================\n";


    // ============================
    // Scrittura su file risultati
    // ============================
    std::ofstream fout ("pot_field.dat");
    std::cout << "Writing potentials and fields to pot_field.dat\n";
    fout << "# index    x    y    z    ";

    if (calc_field_term==2 && k > 1.e-5) {
      fout <<"phi_c    phi_p    phi_i    Ex_c    Ey_c    Ez_c   Ex_p    Ey_p    Ez_p   Ex_i    Ey_i    Ez_i\n";

      for (size_t i = 0; i < num_atoms; ++i) {
        phi_i[i] -= phi_p[i] / constant_pol * constant_react;
        field_ix[i] -= field_px[i] / constant_pol * constant_react;
        field_iy[i] -= field_py[i] / constant_pol * constant_react;
        field_iz[i] -= field_pz[i] / constant_pol * constant_react;
        fout << std::setw (5) << i + 1 << "  "
             << std::setw (8) << pos_atoms[i][0] << "  "
             << std::setw (8) << pos_atoms[i][1] << "  "
             << std::setw (8) << pos_atoms[i][2] << "  "
             << std::setw (8) << phi_c[i] << "  "
             << std::setw (8) << phi_p[i] << "  "
             << std::setw (8) << phi_i[i] << "  ";
        fout << std::setw (8) << field_cx[i] << "  ";
        fout << std::setw (8) << field_cy[i] << "  ";
        fout << std::setw (8) << field_cz[i] << "  ";
        fout << std::setw (8) << field_px[i] << "  ";
        fout << std::setw (8) << field_py[i] << "  ";
        fout << std::setw (8) << field_pz[i] << "  ";
        fout << std::setw (8) << field_ix[i] << "  ";
        fout << std::setw (8) << field_iy[i] << "  ";
        fout << std::setw (8) << field_iz[i] << "  ";
        fout << "\n";
      }
    } else if ((calc_potential_term==2 && k > 1.e-5) && calc_field_term == 0) {
      fout <<"phi_c    phi_p    phi_i    Ex_c    Ey_c    Ez_c\n";

      for (size_t i = 0; i < num_atoms; ++i) {
        phi_i[i] -= phi_p[i] / constant_pol * constant_react;

        fout << std::setw (5) << i + 1 << "  "
             << std::setw (8) << pos_atoms[i][0] << "  "
             << std::setw (8) << pos_atoms[i][1] << "  "
             << std::setw (8) << pos_atoms[i][2] << "  "
             << std::setw (8) << phi_c[i] << "  "
             << std::setw (8) << phi_p[i] << "  "
             << std::setw (8) << phi_i[i] << "  ";
        fout << std::setw (8) << field_cx[i] << "  ";
        fout << std::setw (8) << field_cy[i] << "  ";
        fout << std::setw (8) << field_cz[i] << "  ";
        fout << "\n";
      }
    } else if ((calc_potential_term==2 && k > 1.e-5) && calc_field_term == 1) {
      fout <<"phi_c    phi_p    phi_i    Ex_c    Ey_c    Ez_c   Ex_p    Ey_p    Ez_p\n";

      for (size_t i = 0; i < num_atoms; ++i) {
        phi_i[i] -= phi_p[i] / constant_pol * constant_react;

        fout << std::setw (5) << i + 1 << "  "
             << std::setw (8) << pos_atoms[i][0] << "  "
             << std::setw (8) << pos_atoms[i][1] << "  "
             << std::setw (8) << pos_atoms[i][2] << "  "
             << std::setw (8) << phi_c[i] << "  "
             << std::setw (8) << phi_p[i] << "  "
             << std::setw (8) << phi_i[i] << "  ";
        fout << std::setw (8) << field_cx[i] << "  ";
        fout << std::setw (8) << field_cy[i] << "  ";
        fout << std::setw (8) << field_cz[i] << "  ";
        fout << std::setw (8) << field_px[i] << "  ";
        fout << std::setw (8) << field_py[i] << "  ";
        fout << std::setw (8) << field_pz[i] << "  ";

        fout << "\n";
      }
    } else if ((calc_potential_term == 1 || (calc_potential_term == 2 && k < 1.e-5)) && calc_field_term == 0) {
      fout <<"phi_c    phi_p    Ex_c    Ey_c    Ez_c\n";

      for (size_t i = 0; i < num_atoms; ++i) {
        fout << std::setw (5) << i + 1 << "  "
             << std::setw (8) << pos_atoms[i][0] << "  "
             << std::setw (8) << pos_atoms[i][1] << "  "
             << std::setw (8) << pos_atoms[i][2] << "  "
             << std::setw (8) << phi_c[i] << "  "
             << std::setw (8) << phi_p[i] << "  ";
        fout << std::setw (8) << field_cx[i] << "  ";
        fout << std::setw (8) << field_cy[i] << "  ";
        fout << std::setw (8) << field_cz[i] << "  ";
        fout << "\n";
      }
    } else if ((calc_field_term==1 || (calc_field_term == 2 && k < 1.e-5)) && calc_potential_term < 2) {
      fout <<"phi_c    phi_p    Ex_c    Ey_c    Ez_c   Ex_p    Ey_p    Ez_p\n";

      for (size_t i = 0; i < num_atoms; ++i) {
        fout << std::setw (5) << i + 1 << "  "
             << std::setw (8) << pos_atoms[i][0] << "  "
             << std::setw (8) << pos_atoms[i][1] << "  "
             << std::setw (8) << pos_atoms[i][2] << "  "
             << std::setw (8) << phi_c[i] << "  "
             << std::setw (8) << phi_p[i] << "  ";
        fout << std::setw (8) << field_cx[i] << "  ";
        fout << std::setw (8) << field_cy[i] << "  ";
        fout << std::setw (8) << field_cz[i] << "  ";
        fout << std::setw (8) << field_px[i] << "  ";
        fout << std::setw (8) << field_py[i] << "  ";
        fout << std::setw (8) << field_pz[i] << "  ";

        fout << "\n";
      }
    }


    fout.close();
    std::cout << "Atom potentials and fields written to 'pot_field.dat'\n";
  }
}

void
poisson_boltzmann::pot_field (ray_cache_t & ray_cache)
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);

  if (rank == 0)
    std::cout << "\n================ [ Calculating Potential & Field Components ] =================\n";

  // ===========================
  // Costanti fisiche e scalari
  // ===========================
  const double inv_4pi = 1.0 / (4.0 * pi);
  const double eps0 = e_0; // Permittività del vuoto
  const double eps_in = 4.0 * pi * eps0 * e_in * kb * T * Angs / (e * e);
  const double eps_out = 4.0 * pi * eps0 * e_out * kb * T * Angs / (e * e);

  const double C0 = 1.0e3 * N_av * ionic_strength; // [mol/m^3]
  const double k2 = 2.0 * C0 * Angs * Angs * e * e / (eps0 * e_out * kb * T);
  const double k = std::sqrt (k2);

  const double den_in = 1.0 / eps_in;
  const double constant_pol = (1.0 / eps_out - 1.0 / eps_in) * inv_4pi;
  const double constant_react = (1.0 / eps_out) * inv_4pi;

  // ===========================
  // Variabili locali
  // ===========================
  // double energy_pol = 0.0;
  // double energy_react = 0.0;
  // double coul_energy = 0.0;
  double charge_pol = 0.0;
  const size_t num_atoms = charge_atoms.size();

  std::vector<double> phi_c, phi_p, phi_i;
  std::vector<double> field_cx, field_px, field_ix;
  std::vector<double> field_cy, field_py, field_iy;
  std::vector<double> field_cz, field_pz, field_iz;

  phi_c.assign (num_atoms, 0.0);
  field_cx.assign (num_atoms, 0.0);
  field_cy.assign (num_atoms, 0.0);
  field_cz.assign (num_atoms, 0.0);

  for (size_t i = 0; i < num_atoms; ++i) {
    const std::array<double,3> &ri = pos_atoms[i];
    const double qi = charge_atoms[i];

    for (size_t j = i + 1; j < num_atoms; ++j) {
      const std::array<double,3> &rj = pos_atoms[j];
      const double qj = charge_atoms[j];

      const double dx = ri[0] - rj[0];
      const double dy = ri[1] - rj[1];
      const double dz = ri[2] - rj[2];
      const double r2 = dx * dx + dy * dy + dz * dz;
      const double r = std::sqrt (r2);
      const double inv_r3 = 1.0 / (r2 * r);

      // Coulomb energy
      this->coul_energy += (qi * qj) / r * den_in;
      // Coulomb potential
      phi_c[i] += qj / r * den_in;
      phi_c[j] += qi / r * den_in;


      // Coulomb field
      std::array<double, 3> eij = {
        dx * inv_r3 * qj * den_in,
        dy * inv_r3 * qj * den_in,
        dz * inv_r3 * qj * den_in
      };
      std::array<double, 3> eji = {
        -dx * inv_r3 * qi * den_in,
        -dy * inv_r3 * qi * den_in,
        -dz * inv_r3 * qi * den_in
      };

      field_cx[i] += eij[0];
      field_cx[j] += eji[0];
      field_cy[i] += eij[1];
      field_cy[j] += eji[1];
      field_cz[i] += eij[2];
      field_cz[j] += eji[2];
    }
  }

  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  double first_int = 0.0, second_int = 0.0, distance = 0.0;

  std::array<double,3> h{0}, area_h{0};
  std::array<double,3> V, N;
  std::array<double,8> tmp_eps, tmp_phi;
  std::vector<int> edg, fl_dir;

  auto quadrant = this->tmsh.begin_quadrant_sweep ();

  // flux and polarization energy calculation
  if ((calc_field_term==1 || (calc_field_term == 2 && k < 1.e-5)) && calc_potential_term < 2) {
    phi_p.assign (num_atoms, 0.0);
    field_px.assign (num_atoms, 0.0);
    field_py.assign (num_atoms, 0.0);
    field_pz.assign (num_atoms, 0.0);

    for (const int ii : border_quad) {
      quadrant[ii];

      for (int d = 0; d < 3; ++d)
        h[d] = quadrant->p (d, 7) - quadrant->p (d, 0);

      area_h = {h[1]*h[2]/h[0]*0.25, h[0]*h[2]/h[1]*0.25, h[0]*h[1]/h[2]*0.25};
      std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux (quadrant, tmp_phi, tmp_eps);

      for (int ip = 0; ip < edg.size (); ++ip) {
        const int edge = edg[ip];
        const int axis = edge_axis[edge];
        const int i1 = edge2nodes[2 * edge];
        const int i2 = edge2nodes[2 * edge + 1];
        double fract = 0.0;
        normal_intersection (quadrant, ray_cache, edge, N, fract);
        V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
        V[axis] += fract * h[axis];
        const double tmp_flux =
          - (tmp_phi[i2] - tmp_phi[i1]) * wha (tmp_eps[i1], tmp_eps[i2], fract)
          * fl_dir[ip] * area_h[axis];

        charge_pol += tmp_flux;

        for (size_t ia = 0; ia < num_atoms; ++ia) {
          const std::array<double,3> &ra = pos_atoms[ia];
          const double dx = ra[0] - V[0];
          const double dy = ra[1] - V[1];
          const double dz = ra[2] - V[2];
          const double r = std::sqrt (dx * dx + dy * dy + dz * dz);
          const double inv_r3 = 1.0 / (r * r * r);
          const double qflux = tmp_flux / r;

          first_int += charge_atoms[ia] * qflux;
          phi_p[ia] += qflux * constant_pol;
          field_px[ia] += dx * inv_r3 * tmp_flux * constant_pol;
          field_py[ia] += dy * inv_r3 * tmp_flux * constant_pol;
          field_pz[ia] += dz * inv_r3 * tmp_flux * constant_pol;
        }
      }
    }

    this->energy_pol = 0.5*constant_pol*first_int;
  } else if ((calc_potential_term == 1 || (calc_potential_term == 2 && k < 1.e-5)) && calc_field_term == 0) {
    phi_p.assign (num_atoms, 0.0);

    for (const int ii : border_quad) {
      quadrant[ii];

      for (int d = 0; d < 3; ++d)
        h[d] = quadrant->p (d, 7) - quadrant->p (d, 0);

      area_h = {h[1]*h[2]/h[0]*0.25, h[0]*h[2]/h[1]*0.25, h[0]*h[1]/h[2]*0.25};
      std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux (quadrant, tmp_phi, tmp_eps);

      for (int ip = 0; ip < edg.size (); ++ip) {
        const int edge = edg[ip];
        const int axis = edge_axis[edge];
        const int i1 = edge2nodes[2 * edge];
        const int i2 = edge2nodes[2 * edge + 1];
        double fract = 0.0;
        normal_intersection (quadrant, ray_cache, edge, N, fract);
        V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
        V[axis] += fract * h[axis];
        const double tmp_flux =
          - (tmp_phi[i2] - tmp_phi[i1]) * wha (tmp_eps[i1], tmp_eps[i2], fract)
          * fl_dir[ip] * area_h[axis];

        charge_pol += tmp_flux;

        for (size_t ia = 0; ia < num_atoms; ++ia) {
          const std::array<double,3> &ra = pos_atoms[ia];
          const double dx = ra[0] - V[0];
          const double dy = ra[1] - V[1];
          const double dz = ra[2] - V[2];
          const double r = std::sqrt (dx * dx + dy * dy + dz * dz);
          const double inv_r3 = 1.0 / (r * r * r);
          const double qflux = tmp_flux / r;

          first_int += charge_atoms[ia] * qflux;
          phi_p[ia] += qflux * constant_pol;

        }
      }
    }

    this->energy_pol = 0.5*constant_pol*first_int;
  }

  auto allocate_potential_fields = [&] (void) {
    phi_p.assign (num_atoms, 0.0);
    phi_i.assign (num_atoms, 0.0);
    field_px.assign (num_atoms, 0.0);
    field_py.assign (num_atoms, 0.0);
    field_pz.assign (num_atoms, 0.0);
    field_ix.assign (num_atoms, 0.0);
    field_iy.assign (num_atoms, 0.0);
    field_iz.assign (num_atoms, 0.0);
  };

  //polarization energy + ionic energy
  if (calc_field_term==2 && k > 1.e-5) {
    allocate_potential_fields();
    int cubeindex = -1;
    std::array<std::array<double,3>,3> vert_triangles, norms_vert;
    std::array<double,3> dist_vert, phi_sup;
    int ntriang = 0;

    for (const int ii : border_quad) {
      quadrant[ii];

      for (int d = 0; d < 3; ++d)
        h[d] = quadrant->p (d, 7) - quadrant->p (d, 0);

      area_h = {h[1]*h[2]/h[0]*0.25, h[0]*h[2]/h[1]*0.25, h[0]*h[1]/h[2]*0.25};
      cubeindex = classifyCube (quadrant, eps_out);
      std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux (quadrant, tmp_phi, tmp_eps);
      ntriang = getTriangles (cubeindex, triangles);

      // --- flussi
      for (int ip = 0; ip < edg.size (); ++ip) {
        const int edge = edg[ip];
        const int axis = edge_axis[edge];
        const int i1 = edge2nodes[2 * edge];
        const int i2 = edge2nodes[2 * edge + 1];

        double fract = 0.0;
        normal_intersection (quadrant, ray_cache, edge, N, fract);

        V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
        V[axis] += fract * h[axis];

        const double tmp_flux =
          - (tmp_phi[i2] - tmp_phi[i1]) * wha (tmp_eps[i1], tmp_eps[i2], fract)
          * fl_dir[ip] * area_h[axis];
        charge_pol += tmp_flux;

        for (size_t ia = 0; ia < num_atoms; ++ia) {
          const std::array<double,3> &ra = pos_atoms[ia];
          const double dx = ra[0] - V[0];
          const double dy = ra[1] - V[1];
          const double dz = ra[2] - V[2];
          const double r = std::sqrt (dx * dx + dy * dy + dz * dz);
          const double inv_r3 = 1.0 / (r * r * r);

          first_int += charge_atoms[ia] * tmp_flux / r;
          phi_p[ia] += tmp_flux / r * constant_pol;
          field_px[ia] += dx * inv_r3 * tmp_flux * constant_pol;
          field_py[ia] += dy * inv_r3 * tmp_flux * constant_pol;
          field_pz[ia] += dz * inv_r3 * tmp_flux * constant_pol;
        }
      }

      // --- triangoli (componente ionica)
      for (int itri = 0; itri < ntriang; ++itri) {
        for (int jj = 0; jj < 3; ++jj) {
          const int edge = triangles[itri][jj];
          const int axis = edge_axis[edge];
          const int i1 = edge2nodes[2 * edge];
          const int i2 = edge2nodes[2 * edge + 1];

          double fract = 0.0;
          normal_intersection (quadrant, ray_cache, edge, N, fract);

          V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
          V[axis] += fract * h[axis];

          vert_triangles[jj] = V;
          norms_vert[jj] = N;

          phi_sup[jj] = phi0 (tmp_eps[i1], tmp_eps[i2], tmp_phi[i1], tmp_phi[i2], fract);
        }

        const double area = areaTriangle (vert_triangles);

        for (size_t ia = 0; ia < num_atoms; ++ia) {
          const double qi = charge_atoms[ia];
          const std::array<double,3> &ra = pos_atoms[ia];

          for (int kk = 0; kk < 3; ++kk) {
            dist_vert = {vert_triangles[kk][0] - ra[0],
                         vert_triangles[kk][1] - ra[1],
                         vert_triangles[kk][2] - ra[2]
                        };
            const double r2 = dist_vert[0]*dist_vert[0] + dist_vert[1]*dist_vert[1] + dist_vert[2]*dist_vert[2];
            const double r = std::sqrt (r2);
            const double inv_r3 = 1.0 / (r2 * r);
            const double inv_r5 = inv_r3 / r2;
            const double dot = dist_vert[0]*norms_vert[kk][0] + dist_vert[1]*norms_vert[kk][1] + dist_vert[2]*norms_vert[kk][2];
            const double factor = phi_sup[kk] * inv_4pi * area / 3.0;

            second_int += qi * phi_sup[kk] * dot * inv_r3 * inv_4pi * area / 3.0;
            phi_i[ia] += phi_sup[kk] * dot * inv_r3 * inv_4pi * area / 3.0;

            field_ix[ia] += factor * (-3 * dist_vert[0] * inv_r5 * dot + norms_vert[kk][0] * inv_r3);
            field_iy[ia] += factor * (-3 * dist_vert[1] * inv_r5 * dot + norms_vert[kk][1] * inv_r3);
            field_iz[ia] += factor * (-3 * dist_vert[2] * inv_r5 * dot + norms_vert[kk][2] * inv_r3);
          }
        }
      }
    }

    this->energy_pol = 0.5 * constant_pol * first_int;
    this->energy_react = 0.5 * (second_int - first_int * constant_react);
  } else if ((calc_potential_term==2 && k > 1.e-5) && calc_field_term == 1) {
    phi_p.assign (num_atoms, 0.0);
    phi_i.assign (num_atoms, 0.0);
    field_px.assign (num_atoms, 0.0);
    field_py.assign (num_atoms, 0.0);
    field_pz.assign (num_atoms, 0.0);
    int cubeindex = -1;
    std::array<std::array<double,3>,3> vert_triangles, norms_vert;
    std::array<double,3> dist_vert, phi_sup;
    int ntriang = 0;

    for (const int ii : border_quad) {
      quadrant[ii];

      for (int d = 0; d < 3; ++d)
        h[d] = quadrant->p (d, 7) - quadrant->p (d, 0);

      area_h = {h[1]*h[2]/h[0]*0.25, h[0]*h[2]/h[1]*0.25, h[0]*h[1]/h[2]*0.25};
      cubeindex = classifyCube (quadrant, eps_out);
      std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux (quadrant, tmp_phi, tmp_eps);
      ntriang = getTriangles (cubeindex, triangles);

      // --- flussi
      for (int ip = 0; ip < edg.size (); ++ip) {
        const int edge = edg[ip];
        const int axis = edge_axis[edge];
        const int i1 = edge2nodes[2 * edge];
        const int i2 = edge2nodes[2 * edge + 1];

        double fract = 0.0;
        normal_intersection (quadrant, ray_cache, edge, N, fract);

        V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
        V[axis] += fract * h[axis];

        const double tmp_flux =
          - (tmp_phi[i2] - tmp_phi[i1]) * wha (tmp_eps[i1], tmp_eps[i2], fract)
          * fl_dir[ip] * area_h[axis];
        charge_pol += tmp_flux;

        for (size_t ia = 0; ia < num_atoms; ++ia) {
          const std::array<double,3> &ra = pos_atoms[ia];
          const double dx = ra[0] - V[0];
          const double dy = ra[1] - V[1];
          const double dz = ra[2] - V[2];
          const double r = std::sqrt (dx * dx + dy * dy + dz * dz);
          const double inv_r3 = 1.0 / (r * r * r);
          first_int += charge_atoms[ia] * tmp_flux / r;
          phi_p[ia] += tmp_flux / r * constant_pol;
          field_px[ia] += dx * inv_r3 * tmp_flux * constant_pol;
          field_py[ia] += dy * inv_r3 * tmp_flux * constant_pol;
          field_pz[ia] += dz * inv_r3 * tmp_flux * constant_pol;
        }
      }

      // --- triangoli (componente ionica)
      for (int itri = 0; itri < ntriang; ++itri) {
        for (int jj = 0; jj < 3; ++jj) {
          const int edge = triangles[itri][jj];
          const int axis = edge_axis[edge];
          const int i1 = edge2nodes[2 * edge];
          const int i2 = edge2nodes[2 * edge + 1];

          double fract = 0.0;
          normal_intersection (quadrant, ray_cache, edge, N, fract);

          V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
          V[axis] += fract * h[axis];

          vert_triangles[jj] = V;
          norms_vert[jj] = N;

          phi_sup[jj] = phi0 (tmp_eps[i1], tmp_eps[i2], tmp_phi[i1], tmp_phi[i2], fract);
        }

        const double area = areaTriangle (vert_triangles);

        for (size_t ia = 0; ia < num_atoms; ++ia) {
          const double qi = charge_atoms[ia];
          const std::array<double,3> &ra = pos_atoms[ia];

          for (int kk = 0; kk < 3; ++kk) {
            dist_vert = {vert_triangles[kk][0] - ra[0],
                         vert_triangles[kk][1] - ra[1],
                         vert_triangles[kk][2] - ra[2]
                        };
            const double r2 = dist_vert[0]*dist_vert[0] + dist_vert[1]*dist_vert[1] + dist_vert[2]*dist_vert[2];
            const double r = std::sqrt (r2);
            const double inv_r3 = 1.0 / (r2 * r);
            const double inv_r5 = inv_r3 / r2;
            const double dot = dist_vert[0]*norms_vert[kk][0] + dist_vert[1]*norms_vert[kk][1] + dist_vert[2]*norms_vert[kk][2];
            const double factor = phi_sup[kk] * inv_4pi * area / 3.0;
            second_int += qi * phi_sup[kk] * dot * inv_r3 * inv_4pi * area / 3.0;
            phi_i[ia] += phi_sup[kk] * dot * inv_r3 * inv_4pi * area / 3.0;
          }
        }
      }
    }

    this->energy_pol = 0.5 * constant_pol * first_int;
    this->energy_react = 0.5 * (second_int - first_int * constant_react);
  } else if ((calc_potential_term==2 && k > 1.e-5) && calc_field_term == 0) {
    phi_p.assign (num_atoms, 0.0);
    phi_i.assign (num_atoms, 0.0);
    int cubeindex = -1;
    std::array<std::array<double,3>,3> vert_triangles, norms_vert;
    std::array<double,3> dist_vert, phi_sup;
    int ntriang = 0;

    for (const int ii : border_quad) {
      quadrant[ii];

      for (int d = 0; d < 3; ++d)
        h[d] = quadrant->p (d, 7) - quadrant->p (d, 0);

      area_h = {h[1]*h[2]/h[0]*0.25, h[0]*h[2]/h[1]*0.25, h[0]*h[1]/h[2]*0.25};
      cubeindex = classifyCube (quadrant, eps_out);
      std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux (quadrant, tmp_phi, tmp_eps);
      ntriang = getTriangles (cubeindex, triangles);

      // --- flussi
      for (int ip = 0; ip < edg.size (); ++ip) {
        const int edge = edg[ip];
        const int axis = edge_axis[edge];
        const int i1 = edge2nodes[2 * edge];
        const int i2 = edge2nodes[2 * edge + 1];

        double fract = 0.0;
        normal_intersection (quadrant, ray_cache, edge, N, fract);

        V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
        V[axis] += fract * h[axis];

        const double tmp_flux =
          - (tmp_phi[i2] - tmp_phi[i1]) * wha (tmp_eps[i1], tmp_eps[i2], fract)
          * fl_dir[ip] * area_h[axis];
        charge_pol += tmp_flux;

        for (size_t ia = 0; ia < num_atoms; ++ia) {
          const std::array<double,3> &ra = pos_atoms[ia];
          const double dx = ra[0] - V[0];
          const double dy = ra[1] - V[1];
          const double dz = ra[2] - V[2];
          const double r = std::sqrt (dx * dx + dy * dy + dz * dz);
          const double inv_r3 = 1.0 / (r * r * r);
          first_int += charge_atoms[ia] * tmp_flux / r;
          phi_p[ia] += tmp_flux / r * constant_pol;
        }
      }

      // --- triangoli (componente ionica)
      for (int itri = 0; itri < ntriang; ++itri) {
        for (int jj = 0; jj < 3; ++jj) {
          const int edge = triangles[itri][jj];
          const int axis = edge_axis[edge];
          const int i1 = edge2nodes[2 * edge];
          const int i2 = edge2nodes[2 * edge + 1];

          double fract = 0.0;
          normal_intersection (quadrant, ray_cache, edge, N, fract);

          V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
          V[axis] += fract * h[axis];

          vert_triangles[jj] = V;
          norms_vert[jj] = N;

          phi_sup[jj] = phi0 (tmp_eps[i1], tmp_eps[i2], tmp_phi[i1], tmp_phi[i2], fract);
        }

        const double area = areaTriangle (vert_triangles);

        for (size_t ia = 0; ia < num_atoms; ++ia) {
          const double qi = charge_atoms[ia];
          const std::array<double,3> &ra = pos_atoms[ia];

          for (int kk = 0; kk < 3; ++kk) {
            dist_vert = {vert_triangles[kk][0] - ra[0],
                         vert_triangles[kk][1] - ra[1],
                         vert_triangles[kk][2] - ra[2]
                        };
            const double r2 = dist_vert[0]*dist_vert[0] + dist_vert[1]*dist_vert[1] + dist_vert[2]*dist_vert[2];
            const double r = std::sqrt (r2);
            const double inv_r3 = 1.0 / (r2 * r);
            const double inv_r5 = inv_r3 / r2;
            const double dot = dist_vert[0]*norms_vert[kk][0] + dist_vert[1]*norms_vert[kk][1] + dist_vert[2]*norms_vert[kk][2];
            const double factor = phi_sup[kk] * inv_4pi * area / 3.0;
            second_int += qi * phi_sup[kk] * dot * inv_r3 * inv_4pi * area / 3.0;
            phi_i[ia] += phi_sup[kk] * dot * inv_r3 * inv_4pi * area / 3.0;
          }
        }
      }
    }

    this->energy_pol = 0.5 * constant_pol * first_int;
    this->energy_react = 0.5 * (second_int - first_int * constant_react);
  }


  auto reduce_double = [&] (double &x) {
    MPI_Reduce (rank == 0 ? MPI_IN_PLACE : &x, &x, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
  };

  auto reduce_vec = [&] (std::vector<double> &v) {
    MPI_Reduce (rank == 0 ? MPI_IN_PLACE : v.data(),
                rank == 0 ? v.data() : nullptr,
                (int)v.size(), MPI_DOUBLE, MPI_SUM, 0, mpicomm);
  };

  reduce_double (charge_pol);
  reduce_double (energy_pol);
  reduce_double (energy_react);
  reduce_vec (phi_p);
  reduce_vec (phi_i);
  reduce_vec (field_px);
  reduce_vec (field_py);
  reduce_vec (field_pz);
  reduce_vec (field_ix);
  reduce_vec (field_iy);
  reduce_vec (field_iz);

  // Print the result
  if (rank == 0) {
    constexpr int label_width = 50;
    constexpr int precision = 16;


    if (calc_potential_term == 1 || calc_field_term == 1 || calc_potential_term == 2 || calc_field_term == 2) {
      std::cout << std::left << std::setw (label_width) << "  Polarization energy [kT]:"
                << std::setprecision (precision) << energy_pol << "\n";
    }

    if (calc_potential_term == 2 || calc_field_term == 2) {
      std::cout << std::left << std::setw (label_width) << "  Direct ionic energy [kT]:"
                << std::setprecision (precision) << energy_react << "\n";
    }


    std::cout << std::left << std::setw (label_width) << "  Coulombic energy [kT]:"
              << std::setprecision (precision) << coul_energy << "\n";


    std::cout << std::left << std::setw (label_width) << "  Sum of electrostatic energy contributions [kT]:"
              << std::setprecision (precision)
              << (energy_pol + energy_react + coul_energy) << "\n";

    std::cout << "===========================================================\n";


    // ============================
    // Scrittura su file risultati
    // ============================
    std::ofstream fout ("pot_field.dat");
    std::cout << "Writing potentials and fields to pot_field.dat\n";
    fout << "# index    x    y    z    ";

    if (calc_field_term==2 && k > 1.e-5) {
      fout <<"phi_c    phi_p    phi_i    Ex_c    Ey_c    Ez_c   Ex_p    Ey_p    Ez_p   Ex_i    Ey_i    Ez_i\n";

      for (size_t i = 0; i < num_atoms; ++i) {
        phi_i[i] -= phi_p[i] / constant_pol * constant_react;
        field_ix[i] -= field_px[i] / constant_pol * constant_react;
        field_iy[i] -= field_py[i] / constant_pol * constant_react;
        field_iz[i] -= field_pz[i] / constant_pol * constant_react;
        fout << std::setw (5) << i + 1 << "  "
             << std::setw (8) << pos_atoms[i][0] << "  "
             << std::setw (8) << pos_atoms[i][1] << "  "
             << std::setw (8) << pos_atoms[i][2] << "  "
             << std::setw (8) << phi_c[i] << "  "
             << std::setw (8) << phi_p[i] << "  "
             << std::setw (8) << phi_i[i] << "  ";
        fout << std::setw (8) << field_cx[i] << "  ";
        fout << std::setw (8) << field_cy[i] << "  ";
        fout << std::setw (8) << field_cz[i] << "  ";
        fout << std::setw (8) << field_px[i] << "  ";
        fout << std::setw (8) << field_py[i] << "  ";
        fout << std::setw (8) << field_pz[i] << "  ";
        fout << std::setw (8) << field_ix[i] << "  ";
        fout << std::setw (8) << field_iy[i] << "  ";
        fout << std::setw (8) << field_iz[i] << "  ";
        fout << "\n";
      }
    } else if ((calc_potential_term==2 && k > 1.e-5) && calc_field_term == 0) {
      fout <<"phi_c    phi_p    phi_i    Ex_c    Ey_c    Ez_c \n";

      for (size_t i = 0; i < num_atoms; ++i) {
        phi_i[i] -= phi_p[i] / constant_pol * constant_react;

        fout << std::setw (5) << i + 1 << "  "
             << std::setw (8) << pos_atoms[i][0] << "  "
             << std::setw (8) << pos_atoms[i][1] << "  "
             << std::setw (8) << pos_atoms[i][2] << "  "
             << std::setw (8) << phi_c[i] << "  "
             << std::setw (8) << phi_p[i] << "  "
             << std::setw (8) << phi_i[i] << "  ";
        fout << std::setw (8) << field_cx[i] << "  ";
        fout << std::setw (8) << field_cy[i] << "  ";
        fout << std::setw (8) << field_cz[i] << "  ";
        fout << "\n";
      }
    } else if ((calc_potential_term==2 && k > 1.e-5) && calc_field_term == 1) {
      fout <<"phi_c    phi_p    phi_i    Ex_c    Ey_c    Ez_c   Ex_p    Ey_p    Ez_p\n";

      for (size_t i = 0; i < num_atoms; ++i) {
        phi_i[i] -= phi_p[i] / constant_pol * constant_react;

        fout << std::setw (5) << i + 1 << "  "
             << std::setw (8) << pos_atoms[i][0] << "  "
             << std::setw (8) << pos_atoms[i][1] << "  "
             << std::setw (8) << pos_atoms[i][2] << "  "
             << std::setw (8) << phi_c[i] << "  "
             << std::setw (8) << phi_p[i] << "  "
             << std::setw (8) << phi_i[i] << "  ";
        fout << std::setw (8) << field_cx[i] << "  ";
        fout << std::setw (8) << field_cy[i] << "  ";
        fout << std::setw (8) << field_cz[i] << "  ";
        fout << std::setw (8) << field_px[i] << "  ";
        fout << std::setw (8) << field_py[i] << "  ";
        fout << std::setw (8) << field_pz[i] << "  ";

        fout << "\n";
      }
    } else if ((calc_potential_term == 1 || (calc_potential_term == 2 && k < 1.e-5)) && calc_field_term == 0) {
      fout <<"phi_c    phi_p    Ex_c    Ey_c    Ez_c\n";

      for (size_t i = 0; i < num_atoms; ++i) {
        fout << std::setw (5) << i + 1 << "  "
             << std::setw (8) << pos_atoms[i][0] << "  "
             << std::setw (8) << pos_atoms[i][1] << "  "
             << std::setw (8) << pos_atoms[i][2] << "  "
             << std::setw (8) << phi_c[i] << "  "
             << std::setw (8) << phi_p[i] << "  ";
        fout << std::setw (8) << field_cx[i] << "  ";
        fout << std::setw (8) << field_cy[i] << "  ";
        fout << std::setw (8) << field_cz[i] << "  ";
        fout << "\n";
      }
    } else if ((calc_field_term==1 || (calc_field_term == 2 && k < 1.e-5)) && calc_potential_term < 2) {
      fout <<"phi_c    phi_p    Ex_c    Ey_c    Ez_c   Ex_p    Ey_p    Ez_p\n";

      for (size_t i = 0; i < num_atoms; ++i) {
        fout << std::setw (5) << i + 1 << "  "
             << std::setw (8) << pos_atoms[i][0] << "  "
             << std::setw (8) << pos_atoms[i][1] << "  "
             << std::setw (8) << pos_atoms[i][2] << "  "
             << std::setw (8) << phi_c[i] << "  "
             << std::setw (8) << phi_p[i] << "  ";
        fout << std::setw (8) << field_cx[i] << "  ";
        fout << std::setw (8) << field_cy[i] << "  ";
        fout << std::setw (8) << field_cz[i] << "  ";
        fout << std::setw (8) << field_px[i] << "  ";
        fout << std::setw (8) << field_py[i] << "  ";
        fout << std::setw (8) << field_pz[i] << "  ";

        fout << "\n";
      }
    }


    fout.close();
    std::cout << "Atom potentials and fields written to 'pot_field.dat'\n";
  }
}
