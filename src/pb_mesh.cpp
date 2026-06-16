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
#include <map>
#include <tuple>
#include <p8est.h>
#include <random>
#include <cmath>
#include <numeric>
#include <limits>

void
poisson_boltzmann::create_mesh ()
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);

  double eps_out = 4.0 * pi * e_0 * e_out * kb * T * Angs / (e *e); //adim e_out
  double C_0 = 1.0e3 * N_av * ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0 * C_0 * Angs * Angs * e * e / (e_0 *e_out *kb *T);
  double k = std::sqrt (k2);

  int nx, ny, nz;
  double scale_tmp, scale_x, scale_y, scale_z;
  double lmax = 0;
  double l[3];
  scale_level = unilevel;
  double maxradius = *std::max_element (r_atoms.begin (), r_atoms.end ());

  auto comp_pos_x = [] (const std::array<double, 3> &a1, const std::array<double, 3> &a2) -> bool {
    return a1[0] < a2[0];
  };
  auto comp_pos_y = [] (const std::array<double, 3> &a1, const std::array<double, 3> &a2) -> bool {
    return a1[1] < a2[1];
  };
  auto comp_pos_z = [] (const std::array<double, 3> &a1, const std::array<double, 3> &a2) -> bool {
    return a1[2] < a2[2];
  };


  auto minmax_x = std::minmax_element (pos_atoms.begin (), pos_atoms.end (), comp_pos_x);
  auto minmax_y = std::minmax_element (pos_atoms.begin (), pos_atoms.end (), comp_pos_y);
  auto minmax_z = std::minmax_element (pos_atoms.begin (), pos_atoms.end (), comp_pos_z);



  net_charge = std::accumulate (charge_atoms.begin(), charge_atoms.end(), 0.0);
  int num_atoms = charge_atoms.size ();

  if (rank == 0) {
    std::cout << "\n========== [ System Information ] ==========\n";
    std::cout << "  Number of atoms    : " << num_atoms << '\n';
    std::cout << "  Size protein [Å]   : ";
    std::cout << "[" << (*minmax_x.second)[0] - (*minmax_x.first)[0] + 2 * maxradius << ", "
              << (*minmax_y.second)[1] - (*minmax_y.first)[1] + 2 * maxradius << ", "
              << (*minmax_z.second)[2] - (*minmax_z.first)[2] + 2 * maxradius << "]\n";

    if (std::fabs (net_charge - std::round (net_charge)) > 1.e-5)
      std::cerr << "  [WARNING] Net charge is not an integer: " << net_charge << '\n';

    std::cout << "  Solute epsilon     : " << e_in << '\n';
    std::cout << "  Solvent epsilon    : " << e_out << '\n';
    std::cout << "  Temperature        : " << T << " [K] \n";
    std::cout << "  Ionic strength     : " << ionic_strength << " [mol/L] \n";
    std::cout << "============================================\n\n";
  }

  if (mesh_shape != 2 && mesh_shape != MESH_SHAPE_MEM) {
    l_c[0] = (*minmax_x.first)[0] - maxradius - 2 * prb_radius;
    l_c[1] = (*minmax_y.first)[1] - maxradius - 2 * prb_radius;
    l_c[2] = (*minmax_z.first)[2] - maxradius - 2 * prb_radius;
    r_c[0] = (*minmax_x.second)[0] + maxradius + 2 * prb_radius;
    r_c[1] = (*minmax_y.second)[1] + maxradius + 2 * prb_radius;
    r_c[2] = (*minmax_z.second)[2] + maxradius + 2 * prb_radius;

    for (int kk = 0; kk < 3; ++kk) {
      l[kk] = (r_c[kk] - l_c[kk]);
      cc[kk] = (r_c[kk] + l_c[kk]) * 0.5;
      lmax = l[kk] > lmax ? l[kk] : lmax;
    }

    // For random displacement of the grid
    if (rand_center == 1) {
      std::random_device rd; // Will be used to obtain a seed for the random number engine
      std::mt19937 gen (rd ()); // Standard mersenne_twister_engine seeded with rd()
      std::uniform_real_distribution<> dis (-1. / scale * 0.5, 1. / scale * 0.5);
      double tmp_cc[3];

      for (int n = 0; n < 3; ++n) {
        tmp_cc[n] = cc[n];
        cc[n] = tmp_cc[n] + dis (gen);
      }

      if (rank == 0) {
        std::cout << "  Random grid displacement [Å]:\n";
        std::cout << "    x: " << tmp_cc[0] << " -> " << cc[0] << "  (delta " << std::abs (tmp_cc[0] - cc[0]) << ")\n";
        std::cout << "    y: " << tmp_cc[1] << " -> " << cc[1] << "  (delta " << std::abs (tmp_cc[1] - cc[1]) << ")\n";
        std::cout << "    z: " << tmp_cc[2] << " -> " << cc[2] << "  (delta " << std::abs (tmp_cc[2] - cc[2]) << ")\n";
      }
    }

    for (int kk = 0; kk < 3; ++kk) {
      ll[kk] = cc[kk] - lmax * 0.5;
      rr[kk] = cc[kk] + lmax * 0.5;
    }

    //stretched box with perfil1
    l_c[0] -= l[0] * 0.5 * (1.0 / perfil1 - 1);
    l_c[1] -= l[1] * 0.5 * (1.0 / perfil1 - 1);
    l_c[2] -= l[2] * 0.5 * (1.0 / perfil1 - 1);
    r_c[0] += l[0] * 0.5 * (1.0 / perfil1 - 1);
    r_c[1] += l[1] * 0.5 * (1.0 / perfil1 - 1);
    r_c[2] += l[2] * 0.5 * (1.0 / perfil1 - 1);

    //cubic box with perfil1
    ll[0] -= lmax * 0.5 * (1.0 / perfil1 - 1);
    ll[1] -= lmax * 0.5 * (1.0 / perfil1 - 1);
    ll[2] -= lmax * 0.5 * (1.0 / perfil1 - 1);
    rr[0] += lmax * 0.5 * (1.0 / perfil1 - 1);
    rr[1] += lmax * 0.5 * (1.0 / perfil1 - 1);
    rr[2] += lmax * 0.5 * (1.0 / perfil1 - 1);

    l_cr[0] = l_c[0];
    r_cr[0] = r_c[0];
    l_cr[1] = l_c[1];
    r_cr[1] = r_c[1];
    l_cr[2] = l_c[2];
    r_cr[2] = r_c[2];
  }

  if (mesh_shape == 0) {

    //cubic box with max perfil2
    double size = 1.0 / scale;
    scale_level = 0;

    for (int ii = 0; ii < 29; ++ii) {
      size *= 2;
      scale_level ++;

      if (lmax / size < perfil2)
        break;
    }

    minlevel = outlevel;

    ll[0] = cc[0] - size / 2;
    ll[1] = cc[1] - size / 2;
    ll[2] = cc[2] - size / 2;
    rr[0] = cc[0] + size / 2;
    rr[1] = cc[1] + size / 2;
    rr[2] = cc[2] + size / 2;


    nx = (int) ( (r_cr[0] - cc[0]) * scale + 0.5);
    ny = (int) ( (r_cr[1] - cc[1]) * scale + 0.5);
    nz = (int) ( (r_cr[2] - cc[2]) * scale + 0.5);

    //refined box
    l_cr[0] = cc[0] - nx * 1.0 / scale;
    l_cr[1] = cc[1] - ny * 1.0 / scale;
    l_cr[2] = cc[2] - nz * 1.0 / scale;
    r_cr[0] = cc[0] + nx * 1.0 / scale;
    r_cr[1] = cc[1] + ny * 1.0 / scale;
    r_cr[2] = cc[2] + nz * 1.0 / scale;

    double dist = size / 2 - lmax * 0.5;
    pot_bc = std::exp (-k *dist) / (dist *eps_out);

    if (rank == 0) {
      std::cout << "========== [ Domain Information ] ==========\n";
      std::cout << "  Scale:  " << scale << "\n";

      std::cout << "  Center of the System [Å]:";
      std::cout << "  [" << cc[0] << ", " << cc[1] << ", " << cc[2] << "]\n";

      std::cout << "  Perfil outer box:  " << perfil2 << "\n";
      std::cout << "  Complete Domain Box Size [Å]:\n";
      std::cout << "      x = [" << ll[0] << ", " << rr[0] << "]\n";
      std::cout << "      y = [" << ll[1] << ", " << rr[1] << "]\n";
      std::cout << "      z = [" << ll[2] << ", " << rr[2] << "]\n";

      std::cout << "  Perfil uniform grid:  " << perfil1 << "\n";
      std::cout << "  Uniform grid Size [Å]:\n";
      std::cout << "      x = [" << l_cr[0] << ", " << r_cr[0] << "]\n";
      std::cout << "      y = [" << l_cr[1] << ", " << r_cr[1] << "]\n";
      std::cout << "      z = [" << l_cr[2] << ", " << r_cr[2] << "]\n";

      std::cout << "  Number of Subdivisions in the Uniform grid:";
      std::cout << "  nx = " << nx * 2 << "  ny = " << ny * 2 << "  nz = " << nz * 2 << '\n';

      std::cout << "============================================\n";
    }

    simple_conn_num_vertices = 8;
    simple_conn_num_trees = 1;

    simple_conn_p = std::make_unique<double[]> (simple_conn_num_vertices * 3);
    simple_conn_t = std::make_unique<p4est_topidx_t[]> (simple_conn_num_vertices);

    auto tmp_p = {ll[0], ll[1], ll[2],
                  rr[0], ll[1], ll[2],
                  ll[0], rr[1], ll[2],
                  rr[0], rr[1], ll[2],
                  ll[0], ll[1], rr[2],
                  rr[0], ll[1], rr[2],
                  ll[0], rr[1], rr[2],
                  rr[0], rr[1], rr[2]
                 };
    auto tmp_t = {1, 2, 3, 4, 5, 6, 7, 8};

    std::copy (tmp_p.begin (), tmp_p.end (), simple_conn_p.get ());
    std::copy (tmp_t.begin (), tmp_t.end (), simple_conn_t.get ());

    for (int i = 0; i < 6; i++)
      bcells.push_back (std::make_pair (0, i));
  } else if (mesh_shape == 1) {
    l_cr[0] = ll[0];
    l_cr[1] = ll[1];
    l_cr[2] = ll[2];
    r_cr[0] = rr[0];
    r_cr[1] = rr[1];
    r_cr[2] = rr[2];
    scale = (1 << unilevel) / (rr[0] - ll[0]);
    minlevel = unilevel;

    if (refine_box == 1) {
      //cubic box with perfil2
      ll[0] = cc[0] - lmax * 0.5 * 1.0 / perfil2;
      ll[1] = cc[1] - lmax * 0.5 * 1.0 / perfil2;
      ll[2] = cc[2] - lmax * 0.5 * 1.0 / perfil2;
      rr[0] = cc[0] + lmax * 0.5 * 1.0 / perfil2;
      rr[1] = cc[1] + lmax * 0.5 * 1.0 / perfil2;
      rr[2] = cc[2] + lmax * 0.5 * 1.0 / perfil2;
      scale_tmp = (1 << unilevel) / (rr[0] - ll[0]);

      nx = (int) ( (r_cr[0] - cc[0]) * scale_tmp + 0.5);
      ny = (int) ( (r_cr[1] - cc[1]) * scale_tmp + 0.5);
      nz = (int) ( (r_cr[2] - cc[2]) * scale_tmp + 0.5);

      //refined stretched box
      l_cr[0] = cc[0] - nx * 1.0 / scale_tmp;
      l_cr[1] = cc[1] - ny * 1.0 / scale_tmp;
      l_cr[2] = cc[2] - nz * 1.0 / scale_tmp;
      r_cr[0] = cc[0] + nx * 1.0 / scale_tmp;
      r_cr[1] = cc[1] + ny * 1.0 / scale_tmp;
      r_cr[2] = cc[2] + nz * 1.0 / scale_tmp;
    }

    double dist = ((rr[0] - ll[0]) - lmax) * 0.5;
    pot_bc = std::exp (-k *dist) / (dist *eps_out);

    if (rank == 0) {
      std::cout << "========== [ Domain Information ] ==========\n";
      std::cout << "  Scale:  " << scale << "\n";

      std::cout << "  Center of the System [Å]:";
      std::cout << "  [" << cc[0] << ", " << cc[1] << ", " << cc[2] << "]\n";

      std::cout << "  Perfil box:  " << perfil1 << "\n";
      std::cout << "  Complete Domain Box Size [Å]:\n";
      std::cout << "      x = [" << ll[0] << ", " << rr[0] << "]\n";
      std::cout << "      y = [" << ll[1] << ", " << rr[1] << "]\n";
      std::cout << "      z = [" << ll[2] << ", " << rr[2] << "]\n";

      std::cout << "  Number of Subdivisions:";
      std::cout << "  nx = " << nx * 2 << "  ny = " << ny * 2 << "  nz = " << nz * 2 << '\n';

      std::cout << "============================================\n";
    }

    simple_conn_num_vertices = 8;
    simple_conn_num_trees = 1;

    simple_conn_p = std::make_unique<double[]> (simple_conn_num_vertices * 3);
    simple_conn_t = std::make_unique<p4est_topidx_t[]> (simple_conn_num_vertices);

    auto tmp_p = {ll[0], ll[1], ll[2],
                  rr[0], ll[1], ll[2],
                  ll[0], rr[1], ll[2],
                  rr[0], rr[1], ll[2],
                  ll[0], ll[1], rr[2],
                  rr[0], ll[1], rr[2],
                  ll[0], rr[1], rr[2],
                  rr[0], rr[1], rr[2]
                 };
    auto tmp_t = {1, 2, 3, 4, 5, 6, 7, 8};

    std::copy (tmp_p.begin (), tmp_p.end (), simple_conn_p.get ());
    std::copy (tmp_t.begin (), tmp_t.end (), simple_conn_t.get ());

    for (int i = 0; i < 6; i++)
      bcells.push_back (std::make_pair (0, i));
  } else if (mesh_shape == 2) {

    l_cr[0] = l_c[0];
    l_cr[1] = l_c[1];
    l_cr[2] = l_c[2];
    r_cr[0] = r_c[0];
    r_cr[1] = r_c[1];
    r_cr[2] = r_c[2];
    scale = (1 << unilevel) / (r_c[0] - l_c[0]);
    minlevel = unilevel;

    double dist = ((r_cr[0] - l_cr[0]) - lmax) * 0.5;
    pot_bc = std::exp (-k *dist) / (dist *eps_out);

    if (rank == 0) {
      std::cout << "========== [ Domain Information ] ==========\n";

      std::cout << "  Scale:  " << scale << "\n";

      std::cout << "  Center of the System [Å]:";
      std::cout << "  [" << (r_c[0] + l_c[0]) * 0.5 << ", "
                << (r_c[1] + l_c[1]) * 0.5 << ", "
                << (r_c[2] + l_c[2]) * 0.5 << "]\n";

      std::cout << "  Complete Domain Box Size [Å]:\n";
      std::cout << "      x = [" << l_cr[0] << ", " << r_cr[0] << "]\n";
      std::cout << "      y = [" << l_cr[1] << ", " << r_cr[1] << "]\n";
      std::cout << "      z = [" << l_cr[2] << ", " << r_cr[2] << "]\n";

      std::cout << "  Number of Subdivisions:";
      std::cout << "  nx = " << (1 << unilevel) << "  ny = " << (1 << unilevel) << "  nz = " << (1 << unilevel) << '\n';

      std::cout << "============================================\n";
    }

    if (refine_box == 1) {
      double nx_tmp, ny_tmp, nz_tmp;
      double cc_tmp[3];

      scale_x = (1 << unilevel) / (r_c[0] - l_c[0]);
      scale_y = (1 << unilevel) / (r_c[1] - l_c[1]);
      scale_z = (1 << unilevel) / (r_c[2] - l_c[2]);

      cc[0] = (r_c[0] + l_c[0]) * 0.5;
      cc[1] = (r_c[1] + l_c[1]) * 0.5;
      cc[2] = (r_c[2] + l_c[2]) * 0.5;

      cc_tmp[0] = (r_cr[0] + l_cr[0]) * 0.5;
      cc_tmp[1] = (r_cr[1] + l_cr[1]) * 0.5;
      cc_tmp[2] = (r_cr[2] + l_cr[2]) * 0.5;

      nx_tmp = (int) ( (cc[0] - cc_tmp[0]) * scale_x + 0.5);
      ny_tmp = (int) ( (cc[1] - cc_tmp[1]) * scale_y + 0.5);
      nz_tmp = (int) ( (cc[2] - cc_tmp[2]) * scale_z + 0.5);

      cc_tmp[0] = cc[0] + nx_tmp * 1.0 / scale_x;
      cc_tmp[1] = cc[1] + ny_tmp * 1.0 / scale_y;
      cc_tmp[2] = cc[2] + nz_tmp * 1.0 / scale_z;

      nx = (int) ( (r_cr[0] - cc_tmp[0]) * scale_x + 0.5);
      ny = (int) ( (r_cr[1] - cc_tmp[1]) * scale_y + 0.5);
      nz = (int) ( (r_cr[2] - cc_tmp[2]) * scale_z + 0.5);
      //refined stretched box
      l_cr[0] = cc_tmp[0] - nx * 1.0 / scale_x;
      l_cr[1] = cc_tmp[1] - ny * 1.0 / scale_y;
      l_cr[2] = cc_tmp[2] - nz * 1.0 / scale_z;
      r_cr[0] = cc_tmp[0] + nx * 1.0 / scale_x;
      r_cr[1] = cc_tmp[1] + ny * 1.0 / scale_y;
      r_cr[2] = cc_tmp[2] + nz * 1.0 / scale_z;

      if (rank == 0) {
        std::cout << "  Refined inner box [Å]:\n";
        std::cout << "      x = [" << l_cr[0] << ", " << r_cr[0] << "]\n";
        std::cout << "      y = [" << l_cr[1] << ", " << r_cr[1] << "]\n";
        std::cout << "      z = [" << l_cr[2] << ", " << r_cr[2] << "]\n";
      }
    }

    simple_conn_num_vertices = 8;
    simple_conn_num_trees = 1;

    simple_conn_p = std::make_unique<double[]> (simple_conn_num_vertices * 3);
    simple_conn_t = std::make_unique<p4est_topidx_t[]> (simple_conn_num_vertices);

    auto tmp_p = {l_c[0], l_c[1], l_c[2],
                  r_c[0], l_c[1], l_c[2],
                  l_c[0], r_c[1], l_c[2],
                  r_c[0], r_c[1], l_c[2],
                  l_c[0], l_c[1], r_c[2],
                  r_c[0], l_c[1], r_c[2],
                  l_c[0], r_c[1], r_c[2],
                  r_c[0], r_c[1], r_c[2]
                 };
    auto tmp_t = {1, 2, 3, 4, 5, 6, 7, 8};

    std::copy (tmp_p.begin (), tmp_p.end (), simple_conn_p.get ());
    std::copy (tmp_t.begin (), tmp_t.end (), simple_conn_t.get ());

    for (int i = 0; i < 6; i++)
      bcells.push_back (std::make_pair (0, i));

  } else if (mesh_shape == 3) {
    //cubic box with max perfil2
    double size = 1.0 / scale;
    double scale_min_box = 0.5;
    scale_level = 0;
    scale_level_min_box = 0;


    for (int ii = 0; ii < 28; ++ii) {
      size *= 2;
      scale_level ++;

      if (lmax / size < perfil2)
        break;
    }

    for (int ii = 0; ii < scale_level; ++ii) {
      scale_level_min_box ++;

      if ( (std::pow (2, scale_level_min_box) + 1) / size >= scale_min_box)
        break;
    }

    ll[0] = cc[0] - size / 2;
    ll[1] = cc[1] - size / 2;
    ll[2] = cc[2] - size / 2;
    rr[0] = cc[0] + size / 2;
    rr[1] = cc[1] + size / 2;
    rr[2] = cc[2] + size / 2;

    //refined box NS
    nx = (int) ( (r_cr[0] - cc[0]) * scale + 0.5);
    ny = (int) ( (r_cr[1] - cc[1]) * scale + 0.5);
    nz = (int) ( (r_cr[2] - cc[2]) * scale + 0.5);

    l_cr[0] = cc[0] - nx * 1.0 / scale;
    l_cr[1] = cc[1] - ny * 1.0 / scale;
    l_cr[2] = cc[2] - nz * 1.0 / scale;
    r_cr[0] = cc[0] + nx * 1.0 / scale;
    r_cr[1] = cc[1] + ny * 1.0 / scale;
    r_cr[2] = cc[2] + nz * 1.0 / scale;

    //refined box FOCUS

    double err_l;
    l_box[0] = cc_focusing[0] - std::round (n_grid / 2.) / scale;
    l_box[1] = cc_focusing[1] - std::round (n_grid / 2.) / scale;
    l_box[2] = cc_focusing[2] - std::round (n_grid / 2.) / scale;
    r_box[0] = cc_focusing[0] + std::round (n_grid / 2.) / scale;
    r_box[1] = cc_focusing[1] + std::round (n_grid / 2.) / scale;
    r_box[2] = cc_focusing[2] + std::round (n_grid / 2.) / scale;


    if ( (r_box[0] - l_box[0]) >= (r_cr[0] - l_cr[0])) {
      r_box[0] = r_cr[0];
      l_box[0] = l_cr[0];
    } else {
      if (l_box[0] <= l_cr[0]) {
        err_l = l_cr[0] - l_box[0];
        l_box[0] = l_box[0] + err_l;
        r_box[0] = r_box[0] + err_l;
      }

      if (r_box[0] >= r_cr[0]) {
        err_l = r_box[0] - r_cr[0];
        l_box[0] = l_box[0] - err_l;
        r_box[0] = r_box[0] - err_l;
      }
    }

    if ( (r_box[1] - l_box[1]) >= (r_cr[1] - l_cr[1])) {
      r_box[1] = r_cr[1];
      l_box[1] = l_cr[1];
    } else {
      if (l_box[1] <= l_cr[1]) {
        err_l = l_cr[1] - l_box[1];
        l_box[1] = l_box[1] + err_l;
        r_box[1] = r_box[1] + err_l;
      }

      if (r_box[1] >= r_cr[1]) {
        err_l = r_box[1] - r_cr[1];
        l_box[1] = l_box[1] - err_l;
        r_box[1] = r_box[1] - err_l;
      }
    }

    if ( (r_box[2] - l_box[2]) >= (r_cr[2] - l_cr[2])) {
      r_box[2] = r_cr[2];
      l_box[2] = l_cr[2];
    } else {
      if (l_box[2] <= l_cr[2]) {
        err_l = l_cr[2] - l_box[2];
        l_box[2] = l_box[2] + err_l;
        r_box[2] = r_box[2] + err_l;
      }

      if (r_box[2] >= r_cr[2]) {
        err_l = r_box[2] - r_cr[2];
        r_box[2] = r_box[2] - err_l;
      }
    }

    // compute outlevel
    unsigned int ratio_l_b;
    double len_box_foc = n_grid / scale;
    double len_box = rr[0] - ll[0];

    for (int ii = 1; ii <= 5; ++ii) {
      len_box /= 2.0;

      if (len_box <= len_box_foc) {
        ratio_l_b = ii;
        break;
      }
    }

    outlevel = ratio_l_b;
    minlevel = outlevel;

    if (rank == 0) {
      std::cout << "========== [ Domain Information ] ==========\n";


      std::cout << "  Center of the System [Å]:";
      std::cout << "  [" << cc[0] << ", " << cc[1] << ", " << cc[2] << "]\n";

      std::cout << "  Center of the focusing [Å]:";
      std::cout << "  [" << cc_focusing[0] << ", " << cc_focusing[1] << ", " << cc_focusing[2] << "]\n";

      std::cout << "  Scale in the box focusing:  " << scale << "\n";

      std::cout << "  Complete Domain Box Size [Å]:\n";
      std::cout << "      x = [" << ll[0] << ", " << rr[0] << "]\n";
      std::cout << "      y = [" << ll[1] << ", " << rr[1] << "]\n";
      std::cout << "      z = [" << ll[2] << ", " << rr[2] << "]\n";

      std::cout << "  Focusing Box Size [Å]:\n";
      std::cout << "      x = [" << l_box[0] << ", " << r_box[0] << "]\n";
      std::cout << "      y = [" << l_box[1] << ", " << r_box[1] << "]\n";
      std::cout << "      z = [" << l_box[2] << ", " << r_box[2] << "]\n";


      std::cout << "  Number of Subdivisions in the focusing Box:";
      std::cout << "  nx = " << n_grid << "  ny = " << n_grid << "  nz = " << n_grid << '\n';

      std::cout << "============================================\n";
    }


    simple_conn_num_vertices = 8;
    simple_conn_num_trees = 1;

    simple_conn_p = std::make_unique<double[]> (simple_conn_num_vertices * 3);
    simple_conn_t = std::make_unique<p4est_topidx_t[]> (simple_conn_num_vertices);

    auto tmp_p = {ll[0], ll[1], ll[2],
                  rr[0], ll[1], ll[2],
                  ll[0], rr[1], ll[2],
                  rr[0], rr[1], ll[2],
                  ll[0], ll[1], rr[2],
                  rr[0], ll[1], rr[2],
                  ll[0], rr[1], rr[2],
                  rr[0], rr[1], rr[2]
                 };
    auto tmp_t = {1, 2, 3, 4, 5, 6, 7, 8};

    std::copy (tmp_p.begin (), tmp_p.end (), simple_conn_p.get ());
    std::copy (tmp_t.begin (), tmp_t.end (), simple_conn_t.get ());

    for (int i = 0; i < 6; i++)
      bcells.push_back (std::make_pair (0, i));
  } else if (mesh_shape == 4) {

    //cubic box with max perfil2

    double size = 1.0 / scale_max;

    scale = scale_max;
    maxlevel = 0;

    for (int ii = 0; ii < 28; ++ii) {
      size *= 2;
      maxlevel ++;

      if (lmax / size < perfil2)
        break;
    }

    minlevel = (int) (maxlevel - std::sqrt (scale_max / scale_min));

    unilevel = (int) ( (maxlevel + minlevel) * 0.5);
    unilevel = unilevel + 1;
    outlevel = minlevel;
    scale_level = maxlevel;

    ll[0] = cc[0] - size / 2;
    ll[1] = cc[1] - size / 2;
    ll[2] = cc[2] - size / 2;
    rr[0] = cc[0] + size / 2;
    rr[1] = cc[1] + size / 2;
    rr[2] = cc[2] + size / 2;

    //refined box NS
    nx = (int) ( (r_cr[0] - cc[0]) * scale + 0.5);
    ny = (int) ( (r_cr[1] - cc[1]) * scale + 0.5);
    nz = (int) ( (r_cr[2] - cc[2]) * scale + 0.5);

    l_cr[0] = cc[0] - nx * 1.0 / scale;
    l_cr[1] = cc[1] - ny * 1.0 / scale;
    l_cr[2] = cc[2] - nz * 1.0 / scale;
    r_cr[0] = cc[0] + nx * 1.0 / scale;
    r_cr[1] = cc[1] + ny * 1.0 / scale;
    r_cr[2] = cc[2] + nz * 1.0 / scale;


    if (rank == 0) {
      std::cout << "========== [ Domain Information ] ==========\n";
      std::cout << "  Center of the System [Å]:";
      std::cout << "  [" << cc[0] << ", " << cc[1] << ", " << cc[2] << "]\n";
      std::cout << "  Complete Domain Box Size [Å]:\n";
      std::cout << "      x = [" << ll[0] << ", " << rr[0] << "]\n";
      std::cout << "      y = [" << ll[1] << ", " << rr[1] << "]\n";
      std::cout << "      z = [" << ll[2] << ", " << rr[2] << "]\n";
      std::cout << "  Levels: min = " << minlevel << "  max = " << maxlevel
                << "  uniform = " << unilevel << '\n';
      std::cout << "============================================\n";
    }

    simple_conn_num_vertices = 8;
    simple_conn_num_trees = 1;

    simple_conn_p = std::make_unique<double[]> (simple_conn_num_vertices * 3);
    simple_conn_t = std::make_unique<p4est_topidx_t[]> (simple_conn_num_vertices);

    auto tmp_p = {ll[0], ll[1], ll[2],
                  rr[0], ll[1], ll[2],
                  ll[0], rr[1], ll[2],
                  rr[0], rr[1], ll[2],
                  ll[0], ll[1], rr[2],
                  rr[0], ll[1], rr[2],
                  ll[0], rr[1], rr[2],
                  rr[0], rr[1], rr[2]
                 };
    auto tmp_t = {1, 2, 3, 4, 5, 6, 7, 8};

    std::copy (tmp_p.begin (), tmp_p.end (), simple_conn_p.get ());
    std::copy (tmp_t.begin (), tmp_t.end (), simple_conn_t.get ());

    for (int i = 0; i < 6; i++)
      bcells.push_back (std::make_pair (0, i));
  } else if (mesh_shape == 5) {
    l_cr[0] = cc[0] - len / 2;
    l_cr[1] = cc[1] - len / 2;
    l_cr[2] = cc[2] - len / 2;
    r_cr[0] = cc[0] + len / 2;
    r_cr[1] = cc[1] + len / 2;
    r_cr[2] = cc[2] + len / 2;


    if (unilevel == 0)
      scale = (num_trees[0]) / (len);
    else
      scale = (num_trees[0] * (1 << unilevel)) / (len);
    minlevel = unilevel;

    num_trees[1] = num_trees[0];
    num_trees[2] = num_trees[0];

    if (rank == 0) {
      std::cout << "========== [ Domain Information ] ==========\n";
      std::cout << "  Complete Domain Box Size [Å]:\n";
      std::cout << "      x = [" << l_cr[0] << ", " << r_cr[0] << "]\n";
      std::cout << "      y = [" << l_cr[1] << ", " << r_cr[1] << "]\n";
      std::cout << "      z = [" << l_cr[2] << ", " << r_cr[2] << "]\n";
      std::cout << "  Number of trees: " << num_trees[0] << '\n';
      std::cout << "============================================\n";
    }

    double bound_x = std::abs (r_cr[0] - l_cr[0]);
    double bound_y = std::abs (r_cr[1] - l_cr[1]);
    double bound_z = std::abs (r_cr[2] - l_cr[2]);
    double step[3] = { bound_x / num_trees[0],
                       bound_y / num_trees[1],
                       bound_z / num_trees[2]
                     };
    make_connectivity_3d (num_trees, step, simple_conn_p,
                          simple_conn_num_vertices, simple_conn_t,
                          simple_conn_num_trees, bcells);

    for (p4est_topidx_t i = 0; i < simple_conn_num_vertices; ++i) {
      p4est_topidx_t j = 0;
      simple_conn_p[3 * i + j++] += l_cr[0];
      simple_conn_p[3 * i + j++] += l_cr[1];
      simple_conn_p[3 * i + j] += l_cr[2];
    }
  } else if (mesh_shape == MESH_SHAPE_MEM) {

    // ── Membrane slab mesh ────────────────────────────────────────────────
    // Cubic box of side max(cell_length_x, cell_length_y).
    // The membrane fills the xy face; the slab in z covers both the membrane
    // and the protein (which may extend above/below the membrane).
    //
    // Level structure:
    // outlevel  = maxlevel - nlev_sol   (far solvent)
    // scale_level = maxlevel - nlev_mem  (membrane slab)
    // maxlevel                           (near molecular surface, loc_refinement)

    // --- lipid atom extents in all three directions ---
    auto minmax_lip_x = std::minmax_element (pos_lipid_atoms.begin (), pos_lipid_atoms.end (), comp_pos_x);
    auto minmax_lip_y = std::minmax_element (pos_lipid_atoms.begin (), pos_lipid_atoms.end (), comp_pos_y);
    auto minmax_lip_z = std::minmax_element (pos_lipid_atoms.begin (), pos_lipid_atoms.end (), comp_pos_z);

    double maxrad_lip = r_lipid_atoms.empty () ? 0.0
                        : *std::max_element (r_lipid_atoms.begin (), r_lipid_atoms.end ());
    double maxrad = std::max (maxradius, maxrad_lip);

    // --- z-slab bounds: covers both membrane and protein (protein may be taller) ---
    double z_prot_min = (*minmax_z.first)[2];
    double z_prot_max = (*minmax_z.second)[2];
    double z_lip_min = (*minmax_lip_z.first)[2];
    double z_lip_max = (*minmax_lip_z.second)[2];

    double z_slab_bot = std::min (z_prot_min, z_lip_min) - maxrad - 2.0 * prb_radius;
    double z_slab_top = std::max (z_prot_max, z_lip_max) + maxrad + 2.0 * prb_radius;
    double z_slab_cc = 0.5 * (z_slab_top + z_slab_bot);

    // --- xy extent of the membrane from lipid atoms + vdw radii ---
    // box_side is the minimum enclosing square: (max_pos + r_vdw) - (min_pos - r_vdw)
    double x_lip_min = (*minmax_lip_x.first)[0];
    double x_lip_max = (*minmax_lip_x.second)[0];
    double y_lip_min = (*minmax_lip_y.first)[1];
    double y_lip_max = (*minmax_lip_y.second)[1];


    double lip_lx = (x_lip_max + maxrad_lip) - (x_lip_min - maxrad_lip);
    double lip_ly = (y_lip_max + maxrad_lip) - (y_lip_min - maxrad_lip);

    // If lip_lx != lip_ly, use the smaller side so the square fits inside the
    // membrane patch on all ranks; atoms outside the square are trimmed after NS.
    double box_side = std::min (lip_lx, lip_ly);

    // --- box protein and membrane extents in all three directions ---
    l_prot[0] = (*minmax_x.first)[0] - maxradius - 2 * prb_radius;
    l_prot[1] = (*minmax_y.first)[1] - maxradius - 2 * prb_radius;
    l_prot[2] = (*minmax_z.first)[2] - maxradius - 2 * prb_radius;
    r_prot[0] = (*minmax_x.second)[0] + maxradius + 2 * prb_radius;
    r_prot[1] = (*minmax_y.second)[1] + maxradius + 2 * prb_radius;
    r_prot[2] = (*minmax_z.second)[2] + maxradius + 2 * prb_radius;

    l_mem[0] = x_lip_min - maxrad - 2 * prb_radius;
    l_mem[1] = y_lip_min - maxrad - 2 * prb_radius;
    l_mem[2] = z_lip_min - maxrad - 2 * prb_radius;
    r_mem[0] = x_lip_max + maxrad + 2 * prb_radius;
    r_mem[1] = y_lip_max + maxrad + 2 * prb_radius;
    r_mem[2] = z_lip_max + maxrad + 2 * prb_radius;


    // Center of the membrane patch in xy (overrides preamble's protein-only cc)
    cc[0] = 0.5 * (x_lip_max + x_lip_min);
    cc[1] = 0.5 * (y_lip_max + y_lip_min);
    cc[2] = z_slab_cc;

    // --- square outer box (cubic: same side in x, y, z) ---
    ll[0] = cc[0] - box_side / 2.0;
    rr[0] = cc[0] + box_side / 2.0;
    ll[1] = cc[1] - box_side / 2.0;
    rr[1] = cc[1] + box_side / 2.0;
    ll[2] = cc[2] - box_side / 2.0;
    rr[2] = cc[2] + box_side / 2.0;

    // --- levels ---
    maxlevel = (int) std::ceil (std::log2 (box_side *scale_max));
    scale = (double)(1 << maxlevel) / box_side;
    unilevel = maxlevel - nlev_sol;
    outlevel = maxlevel - nlev_sol;
    minlevel = outlevel;
    scale_level = maxlevel - nlev_mem;
    scale_level_min_box = maxlevel;

    // --- l_c / r_c: slab bounds for the refinement callback ---
    // xy spans the full box (condition always satisfied → only z filters)
    l_c[0] = ll[0];
    r_c[0] = rr[0];
    l_c[1] = ll[1];
    r_c[1] = rr[1];
    l_c[2] = z_slab_bot;
    r_c[2] = z_slab_top;

    // --- l_cr / r_cr: NS box, aligned to the grid (square in xy) ---
    int nx = (int) ((box_side / 2.0) * scale + 0.5);
    int ny = (int) ((box_side / 2.0) * scale + 0.5);
    int nz = (int) ((z_slab_top - z_slab_cc) * scale + 0.5);
    l_cr[0] = cc[0] - nx / scale;
    r_cr[0] = cc[0] + nx / scale;
    l_cr[1] = cc[1] - ny / scale;
    r_cr[1] = cc[1] + ny / scale;
    l_cr[2] = z_slab_cc - nz / scale;
    r_cr[2] = z_slab_cc + nz / scale;

    if (rank == 0) {
      std::cout << "========== [ Domain Information ] ==========\n";
      std::cout << "  Protein extents (with VdW):\n";
      std::cout << "    x: [" << l_prot[0] << ", " << r_prot[0] << "] Å\n";
      std::cout << "    y: [" << l_prot[1] << ", " << r_prot[1] << "] Å\n";
      std::cout << "    z: [" << l_prot[2] << ", " << r_prot[2] << "] Å\n";
      std::cout << "  Membrane patch extents (with VdW):\n";
      std::cout << "    x: [" << l_mem[0] << ", " << r_mem[0] << "] Å\n";
      std::cout << "    y: [" << l_mem[1] << ", " << r_mem[1] << "] Å\n";
      std::cout << "    z: [" << l_mem[2] << ", " << r_mem[2] << "] Å\n";
      std::cout << "  Lipid patch size: lx = " << lip_lx << " Å   ly = " << lip_ly << " Å\n";
      std::cout << "  Box side (min lx,ly)        : " << box_side << " Å\n";
      std::cout << "  Center (cx, cy, cz)         : "
                << cc[0] << "  " << cc[1] << "  " << cc[2] << " Å\n";
      std::cout << "  Outer box x                 : [" << ll[0] << ", " << rr[0] << "] Å\n";
      std::cout << "  Outer box y                 : [" << ll[1] << ", " << rr[1] << "] Å\n";
      std::cout << "  Outer box z                 : [" << ll[2] << ", " << rr[2] << "] Å\n";
      std::cout << "  Slab z                      : [" << z_slab_bot << ", " << z_slab_top
                << "] Å  (thickness " << z_slab_top - z_slab_bot << " Å)\n";
      std::cout << "  FEM domain (l_cr / r_cr) x  : [" << l_cr[0] << ", " << r_cr[0] << "] Å\n";
      std::cout << "  FEM domain (l_cr / r_cr) y  : [" << l_cr[1] << ", " << r_cr[1] << "] Å\n";
      std::cout << "  FEM domain (l_cr / r_cr) z  : [" << l_cr[2] << ", " << r_cr[2] << "] Å\n";
      std::cout << "  scale_max = " << scale_max
                << "   effective scale = " << scale << " cells/Å\n";
      std::cout << "  maxlevel = " << maxlevel
                << "   scale_level (slab) = " << scale_level
                << "   outlevel (solvent) = " << unilevel << '\n';

      if (unilevel >= scale_level || scale_level >= maxlevel)
        std::cerr << "  [ERROR] Level hierarchy violated: need outlevel < scale_level < maxlevel\n"
                  << "          Check nlev_mem and nlev_sol.\n";
      else
        std::cout << "  Level hierarchy OK: " << unilevel << " < " << scale_level
                  << " < " << maxlevel << '\n';

      std::cout << "============================================\n\n";
    }

    simple_conn_num_vertices = 8;
    simple_conn_num_trees = 1;

    simple_conn_p = std::make_unique<double[]> (simple_conn_num_vertices * 3);
    simple_conn_t = std::make_unique<p4est_topidx_t[]> (simple_conn_num_vertices);

    auto tmp_p = {ll[0], ll[1], ll[2],
                  rr[0], ll[1], ll[2],
                  ll[0], rr[1], ll[2],
                  rr[0], rr[1], ll[2],
                  ll[0], ll[1], rr[2],
                  rr[0], ll[1], rr[2],
                  ll[0], rr[1], rr[2],
                  rr[0], rr[1], rr[2]
                 };
    auto tmp_t = {1, 2, 3, 4, 5, 6, 7, 8};

    std::copy (tmp_p.begin (), tmp_p.end (), simple_conn_p.get ());
    std::copy (tmp_t.begin (), tmp_t.end (), simple_conn_t.get ());

    for (int i = 0; i < 6; i++)
      bcells.push_back (std::make_pair (0, i));

  } else {
    if (rank == 0) {
      std::cout << "x: " << ll[0] << ", " << rr[0] << std::endl;
      std::cout << "y: " << ll[1] << ", " << rr[1] << std::endl;
      std::cout << "z: " << ll[2] << ", " << rr[2] << "\n" << std::endl;
    }

    simple_conn_num_vertices = 8;
    simple_conn_num_trees = 1;

    simple_conn_p = std::make_unique<double[]> (simple_conn_num_vertices * 3);
    simple_conn_t = std::make_unique<p4est_topidx_t[]> (simple_conn_num_vertices);

    auto tmp_p = {ll[0], ll[1], ll[2],
                  rr[0], ll[1], ll[2],
                  ll[0], rr[1], ll[2],
                  rr[0], rr[1], ll[2],
                  ll[0], ll[1], rr[2],
                  rr[0], ll[1], rr[2],
                  ll[0], rr[1], rr[2],
                  rr[0], rr[1], rr[2]
                 };
    auto tmp_t = {1, 2, 3, 4, 5, 6, 7, 8};

    std::copy (tmp_p.begin (), tmp_p.end (), simple_conn_p.get ());
    std::copy (tmp_t.begin (), tmp_t.end (), simple_conn_t.get ());

    for (int i = 0; i < 6; i++)
      bcells.push_back (std::make_pair (0, i));
  }

  tmsh.read_connectivity (simple_conn_p.get (), simple_conn_num_vertices,
                          simple_conn_t.get (), simple_conn_num_trees);

}



void
poisson_boltzmann::init_tmesh ()
{
  for (auto i = 0; i < unilevel; ++i) {
    tmsh.set_refine_marker (uniform_refinement);
    tmsh.refine (0, 1);
  }
}


void
poisson_boltzmann::init_tmesh_with_refine_scale ()
{

  for (auto i = 0; i < outlevel; ++i) {
    tmsh.set_refine_marker (uniform_refinement);
    tmsh.refine (0, 1);
  }

  auto refinement = [this]
  (tmesh_3d::quadrant_iterator q) -> int {
    int currentlevel = static_cast<int> (q->the_quadrant->level);
    int retval = 0;
    double x1, y1, z1;

    if (currentlevel >= this->scale_level)
      retval = 0;
    else {
      for (int ii = 0; ii < 8; ++ii) {
        if (! q->is_hanging (ii)) {
          x1 = q -> p (0, ii);
          y1 = q -> p (1, ii);
          z1 = q -> p (2, ii);

          if ( (x1 >= this->l_c[0]) && (x1 <= this->r_c[0])
               && (y1 >= this->l_c[1]) && (y1 <= this->r_c[1])
               && (z1 >= this->l_c[2]) && (z1 <= this->r_c[2])) {
            retval = 1;
            break;
          }
        }
      }
    }

    return (retval);
  };

  for (auto i = 0; i < scale_level; ++i) {
    tmsh.set_refine_marker (refinement);
    tmsh.refine (0, 1);
  }

}


void
poisson_boltzmann::init_tmesh_with_refine_box_scale ()
{

  for (auto i = 0; i < outlevel; ++i) {
    tmsh.set_refine_marker (uniform_refinement);
    tmsh.refine (0, 1);
  }

  auto refinement = [this]
  (tmesh_3d::quadrant_iterator q) -> int {
    int currentlevel = static_cast<int> (q->the_quadrant->level);
    int retval = 0;
    double x1, y1, z1;

    if (currentlevel >= this->scale_level)
      retval = 0;
    else {
      for (int ii = 0; ii < 8; ++ii) {
        if (! q->is_hanging (ii)) {
          x1 = q -> p (0, ii);
          y1 = q -> p (1, ii);
          z1 = q -> p (2, ii);

          if ( (x1 >= this->l_box[0]) && (x1 <= this->r_box[0])
               && (y1 >= this->l_box[1]) && (y1 <= this->r_box[1])
               && (z1 >= this->l_box[2]) && (z1 <= this->r_box[2])) {
            retval = 1;
            break;
          }
        }
      }
    }

    return (retval);
  };

  for (auto i = 0; i < scale_level; ++i) {
    tmsh.set_refine_marker (refinement);
    tmsh.refine (0, 1);
  }
}

void
poisson_boltzmann::init_tmesh_mem_two_box ()
{

  for (auto i = 0; i < outlevel; ++i) {
    tmsh.set_refine_marker (uniform_refinement);
    tmsh.refine (0, 1);
  }


  auto refinement = [this]
  (tmesh_3d::quadrant_iterator q) -> int {
    int currentlevel = static_cast<int> (q->the_quadrant->level);
    int retval = 0;
    double x1, y1, z1;

    if (currentlevel >= this->maxlevel)
      retval = 0;
    else {
      for (int ii = 0; ii < 8; ++ii) {
        if (! q->is_hanging (ii)) {
          x1 = q -> p (0, ii);
          y1 = q -> p (1, ii);
          z1 = q -> p (2, ii);

          if ( (x1 >= this->l_prot[0]) && (x1 <= this->r_prot[0])
               && (y1 >= this->l_prot[1]) && (y1 <= this->r_prot[1])
               && (z1 >= this->l_prot[2]) && (z1 <= this->r_prot[2])) {
            retval = 1;
            break;
          } else if ((x1 >= this->l_mem[0]) && (x1 <= this->r_mem[0])
                     && (y1 >= this->l_mem[1]) && (y1 <= this->r_mem[1])
                     && (z1 >= this->l_mem[2]) && (z1 <= this->r_mem[2]) && (currentlevel <= this->scale_level)) {
            retval = 1;
            break;
          }
        }
      }
    }

    return (retval);
  };

  for (auto i = 0; i < maxlevel; ++i) {
    tmsh.set_refine_marker (refinement);
    tmsh.refine (0, 1);
  }
}

void
poisson_boltzmann::build_pbc_node_map ()
{
  if (!periodic_x && !periodic_y)
    return;

  int size, rank;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);

  double tol = 1e-10 * std::max ({r_cr[0] - l_cr[0], r_cr[1] - l_cr[1], r_cr[2] - l_cr[2]});

  // The true Dirichlet z-faces are the mesh z-extent (the Outer box ll/rr for the
  // membrane slab mesh, l_cr/r_cr otherwise), NOT the refined sub-box l_cr/r_cr:
  // for membrane runs these planes differ, so filtering on l_cr[2]/r_cr[2] would
  // leave the Dirichlet z± edge nodes inside the PBC maps and corrupt the strong
  // elimination (their huge-penalty rows get merged). Derive the bounds from the
  // mesh itself so both mesh shapes are handled correctly.
  double z_dir_min =  std::numeric_limits<double>::max ();
  double z_dir_max = -std::numeric_limits<double>::max ();
  for (auto q = tmsh.begin_quadrant_sweep (); q != tmsh.end_quadrant_sweep (); ++q)
    for (int ii = 0; ii < 8; ++ii) {
      if (q->is_hanging (ii)) continue;
      double z = q->p (2, ii);
      if (z < z_dir_min) z_dir_min = z;
      if (z > z_dir_max) z_dir_max = z;
    }
  if (size > 1) {
    MPI_Allreduce (MPI_IN_PLACE, &z_dir_min, 1, MPI_DOUBLE, MPI_MIN, mpicomm);
    MPI_Allreduce (MPI_IN_PLACE, &z_dir_max, 1, MPI_DOUBLE, MPI_MAX, mpicomm);
  }

  // Collect (t0, t1, global_node_id) for all non-Dirichlet nodes on a given face.
  // Packed as 3 doubles per entry. Rank-local deduplication via a seen-set.
  auto collect_face_nodes = [&] (int face_idx, int d0, int d1) {
    std::vector<double> buf;
    std::set<size_t> seen;
    for (auto q = tmsh.begin_quadrant_sweep (); q != tmsh.end_quadrant_sweep (); ++q) {
      auto fn = q->ef (face_idx);
      if (fn.empty ()) continue;
      for (int fi : fn) {
        double z = q->p (2, fi);
        if (std::abs (z - z_dir_min) < tol || std::abs (z - z_dir_max) < tol)
          continue;
        auto gid = static_cast<size_t> (q->gt (fi));
        if (!seen.insert (gid).second)
          continue;
        buf.push_back (q->p (d0, fi));
        buf.push_back (q->p (d1, fi));
        buf.push_back (static_cast<double> (gid));
      }
    }
    return buf;
  };

  // MPI_Allgatherv: all ranks send their local buffer, all ranks receive the merged result.
  auto allgatherv = [&] (const std::vector<double>& local) {
    int local_n = static_cast<int> (local.size ());
    std::vector<int> counts (size, 0), displs (size, 0);
    MPI_Allgather (&local_n, 1, MPI_INT, counts.data (), 1, MPI_INT, mpicomm);
    for (int r = 1; r < size; ++r)
      displs[r] = displs[r - 1] + counts[r - 1];
    int total_n = displs[size - 1] + counts[size - 1];
    std::vector<double> global (static_cast<size_t> (total_n));
    MPI_Allgatherv (local.data (), local_n, MPI_DOUBLE,
                    global.data (), counts.data (), displs.data (), MPI_DOUBLE, mpicomm);
    return global;
  };

  struct PBCPairDef {
    int face_left, face_right, d0, d1;
    const char* label;
    std::map<size_t, size_t>* map_out;
  };

  std::vector<PBCPairDef> pairs;
  if (periodic_x) pairs.push_back ({0, 1, 1, 2, "x±", &pbc_x_right_to_left});
  if (periodic_y) pairs.push_back ({2, 3, 0, 2, "y±", &pbc_y_right_to_left});

  if (rank == 0)
    std::cout << "  [PBC map] building node map...\n";

  for (auto& p : pairs) {
    auto left_all  = allgatherv (collect_face_nodes (p.face_left,  p.d0, p.d1));
    auto right_all = allgatherv (collect_face_nodes (p.face_right, p.d0, p.d1));

    // Build (t0, t1) → left_gid map.
    std::map<std::pair<double, double>, size_t> left_map;
    for (size_t i = 0; i + 3 <= left_all.size (); i += 3)
      left_map[{left_all[i], left_all[i + 1]}] = static_cast<size_t> (left_all[i + 2]);

    int matched = 0, unmatched = 0;
    p.map_out->clear ();
    for (size_t i = 0; i + 3 <= right_all.size (); i += 3) {
      auto key = std::make_pair (right_all[i], right_all[i + 1]);
      auto it  = left_map.find (key);
      if (it != left_map.end ()) {
        (*p.map_out)[static_cast<size_t> (right_all[i + 2])] = it->second;
        ++matched;
      } else {
        ++unmatched;
      }
    }

    if (rank == 0) {
      std::cout << "  [PBC map] " << p.label << ": "
                << matched << " matched";
      if (unmatched > 0)
        std::cout << ", " << unmatched << " UNMATCHED (mesh not conforming?)";
      std::cout << '\n';
    }
  }
}

void
poisson_boltzmann::check_pbc_face_conformity ()
{
  if (!periodic_x && !periodic_y)
    return;

  int size, rank;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);

  // Collect face quadrants on this rank as packed doubles:
  // { c0_min, c0_max, c1_min, c1_max, level } per quadrant (5 values).
  // d0, d1: tangential directions (bc.md face-node pattern is identical for all faces).
  auto collect_local = [&] (int face_idx, int d0, int d1) {
    std::vector<double> buf;
    for (auto q = tmsh.begin_quadrant_sweep (); q != tmsh.end_quadrant_sweep (); ++q) {
      auto fn = q->ef (face_idx);
      if (fn.empty ()) continue;
      buf.push_back (q->p (d0, fn[0]));  // c0_min
      buf.push_back (q->p (d0, fn[1]));  // c0_max
      buf.push_back (q->p (d1, fn[0]));  // c1_min
      buf.push_back (q->p (d1, fn[2]));  // c1_max
      buf.push_back (static_cast<double> (q->the_quadrant->level));
    }
    return buf;
  };

  // Gather packed buffers from all ranks to rank 0.
  auto gather_to_rank0 = [&] (const std::vector<double>& local) {
    int local_n = static_cast<int> (local.size ());
    std::vector<int> counts (size, 0), displs (size, 0);
    MPI_Gather (&local_n, 1, MPI_INT, counts.data (), 1, MPI_INT, 0, mpicomm);
    std::vector<double> global_buf;
    if (rank == 0) {
      for (int r = 1; r < size; ++r)
        displs[r] = displs[r - 1] + counts[r - 1];
      global_buf.resize (static_cast<size_t> (displs[size - 1] + counts[size - 1]));
    }
    MPI_Gatherv (local.data (), local_n, MPI_DOUBLE,
                 global_buf.data (), counts.data (), displs.data (),
                 MPI_DOUBLE, 0, mpicomm);
    return global_buf;
  };

  // Compare left vs right face. On rank 0 builds a bbox→level map for the right face
  // and checks each left-face quadrant against it.
  auto check_pair = [&] (const char* label, int fl, int fr, int d0, int d1) {
    auto left_buf  = gather_to_rank0 (collect_local (fl, d0, d1));
    auto right_buf = gather_to_rank0 (collect_local (fr, d0, d1));
    if (rank != 0)
      return;

    using Key = std::tuple<double, double, double, double>;
    std::map<Key, int> right_map;
    for (size_t i = 0; i + 5 <= right_buf.size (); i += 5)
      right_map[{right_buf[i], right_buf[i+1], right_buf[i+2], right_buf[i+3]}]
        = static_cast<int> (right_buf[i+4]);

    int n_match = 0, n_mismatch = 0;
    bool printed = false;
    for (size_t i = 0; i + 5 <= left_buf.size (); i += 5) {
      Key k {left_buf[i], left_buf[i+1], left_buf[i+2], left_buf[i+3]};
      int lev_l = static_cast<int> (left_buf[i+4]);
      auto it   = right_map.find (k);
      if (it != right_map.end () && it->second == lev_l) {
        ++n_match;
      } else {
        ++n_mismatch;
        if (!printed) {
          int lev_r = (it != right_map.end ()) ? it->second : -1;
          std::cout << "    example: d0=["
                    << left_buf[i] << "," << left_buf[i+1] << "] d1=["
                    << left_buf[i+2] << "," << left_buf[i+3] << "]"
                    << " level_left=" << lev_l
                    << " level_right=" << lev_r << "\n";
          printed = true;
        }
      }
    }
    int total = n_match + n_mismatch;
    std::cout << "  [PBC check] " << label << ": "
              << n_match << "/" << total;
    if (n_mismatch == 0)
      std::cout << " (CONFORMING)\n";
    else
      std::cout << " — " << n_mismatch << " mismatches (NON-CONFORMING)\n";
  };

  if (rank == 0)
    std::cout << "  [PBC check] face conformity:\n";
  if (periodic_x) check_pair ("x±", 0, 1, 1, 2);
  if (periodic_y) check_pair ("y±", 2, 3, 0, 2);
}

void
poisson_boltzmann::ensure_pbc_face_conformity ()
{
  if (!periodic_x && !periodic_y)
    return;

  int size, rank;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);

  struct PBCPair { int fl, fr, d0, d1; };
  std::vector<PBCPair> pairs;
  if (periodic_x) pairs.push_back ({0, 1, 1, 2});
  if (periodic_y) pairs.push_back ({2, 3, 0, 2});

  // face_idx → tangential directions, used in the refinement lambda
  std::map<int, std::pair<int,int>> face_dirs;
  for (auto& p : pairs) {
    face_dirs[p.fl] = {p.d0, p.d1};
    face_dirs[p.fr] = {p.d0, p.d1};
  }

  const double tol = 1e-10;

  // Collect packed quads on this rank: 5 doubles each (c0min, c0max, c1min, c1max, level)
  auto collect_local = [&] (int fi, int d0, int d1) {
    std::vector<double> buf;
    for (auto q = tmsh.begin_quadrant_sweep (); q != tmsh.end_quadrant_sweep (); ++q) {
      auto fn = q->ef (fi);
      if (fn.empty ()) continue;
      buf.push_back (q->p (d0, fn[0]));
      buf.push_back (q->p (d0, fn[1]));
      buf.push_back (q->p (d1, fn[0]));
      buf.push_back (q->p (d1, fn[2]));
      buf.push_back ((double) q->the_quadrant->level);
    }
    return buf;
  };

  // MPI collective: gather local double buffers to rank 0
  auto gather_to_0 = [&] (const std::vector<double>& local) {
    int local_n = (int) local.size ();
    std::vector<int> counts (size, 0), displs (size, 0);
    MPI_Gather (&local_n, 1, MPI_INT, counts.data (), 1, MPI_INT, 0, mpicomm);
    std::vector<double> global;
    if (rank == 0) {
      for (int r = 1; r < size; ++r) displs[r] = displs[r-1] + counts[r-1];
      global.resize ((size_t)(displs[size-1] + counts[size-1]));
    }
    MPI_Gatherv (local.data (), local_n, MPI_DOUBLE,
                 global.data (), counts.data (), displs.data (), MPI_DOUBLE, 0, mpicomm);
    return global;
  };

  // True if inner quad (by index in inner_buf) is strictly smaller and contained in outer quad.
  // "Strictly smaller" = inner d0-width < outer d0-width.
  auto strictly_inside = [&] (size_t ii, const std::vector<double>& ib,
                               size_t oj, const std::vector<double>& ob) {
    double i_w = ib[ii*5+1] - ib[ii*5];   // inner width in d0
    double o_w = ob[oj*5+1] - ob[oj*5];   // outer width in d0
    if (i_w >= o_w - tol) return false;    // inner not strictly smaller
    return ib[ii*5]   >= ob[oj*5]   - tol &&
           ib[ii*5+1] <= ob[oj*5+1] + tol &&
           ib[ii*5+2] >= ob[oj*5+2] - tol &&
           ib[ii*5+3] <= ob[oj*5+3] + tol;
  };

  int pass = 0;
  while (true) {
    // Find which face quads need refinement. All ranks join the MPI collectives;
    // only rank 0 builds the result.
    // A quad Q on face A is marked when any quad on the opposite face is strictly
    // contained in Q.bbox (opposite side finer → Q must subdivide to match).
    std::vector<double> to_ref;   // {face_idx, c0min, c0max, c1min, c1max} × n on rank 0

    for (auto& p : pairs) {
      auto left_buf  = gather_to_0 (collect_local (p.fl, p.d0, p.d1));
      auto right_buf = gather_to_0 (collect_local (p.fr, p.d0, p.d1));
      if (rank != 0) continue;

      size_t nl = left_buf.size () / 5;
      size_t nr = right_buf.size () / 5;

      for (size_t i = 0; i < nl; ++i) {
        for (size_t j = 0; j < nr; ++j) {
          if (strictly_inside (j, right_buf, i, left_buf)) {
            // right quad j is inside left quad i → left is coarser → refine left
            to_ref.push_back ((double) p.fl);
            to_ref.insert (to_ref.end (), {left_buf[i*5], left_buf[i*5+1],
                                           left_buf[i*5+2], left_buf[i*5+3]});
            break;
          }
        }
      }

      for (size_t j = 0; j < nr; ++j) {
        for (size_t i = 0; i < nl; ++i) {
          if (strictly_inside (i, left_buf, j, right_buf)) {
            // left quad i is inside right quad j → right is coarser → refine right
            to_ref.push_back ((double) p.fr);
            to_ref.insert (to_ref.end (), {right_buf[j*5], right_buf[j*5+1],
                                           right_buf[j*5+2], right_buf[j*5+3]});
            break;
          }
        }
      }
    }

    int n = (rank == 0) ? (int) (to_ref.size () / 5) : 0;
    MPI_Bcast (&n, 1, MPI_INT, 0, mpicomm);
    if (n == 0) break;

    to_ref.resize ((size_t)(n * 5), 0.0);
    MPI_Bcast (to_ref.data (), n * 5, MPI_DOUBLE, 0, mpicomm);

    // Build lookup for the refinement lambda (all ranks)
    using Key5 = std::tuple<int,double,double,double,double>;
    std::set<Key5> refine_set;
    for (int k = 0; k < n; ++k)
      refine_set.insert ({(int) to_ref[k*5], to_ref[k*5+1], to_ref[k*5+2],
                                              to_ref[k*5+3], to_ref[k*5+4]});

    auto balancer = [&] (tmesh_3d::quadrant_iterator q) -> int {
      for (auto& [fi, dirs] : face_dirs) {
        auto fn = q->ef (fi);
        if (fn.empty ()) continue;
        auto [d0, d1] = dirs;
        Key5 k {fi, q->p (d0, fn[0]), q->p (d0, fn[1]),
                    q->p (d1, fn[0]), q->p (d1, fn[2])};
        if (refine_set.count (k)) return 1;
      }
      return 0;
    };

    tmsh.set_refine_marker (balancer);
    tmsh.refine (0, 1);
    ++pass;

    if (rank == 0)
      std::cout << "  [PBC balance] pass " << pass
                << ": refined " << n << " face quad(s)\n";
  }

  if (rank == 0) {
    if (pass == 0)
      std::cout << "  [PBC balance] faces already conforming\n";
    else
      std::cout << "  [PBC balance] conformity reached in " << pass << " pass(es)\n";
  }
}

