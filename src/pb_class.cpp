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
#include "pb_membrane.h"
#include "GetPot"

#include <bim_distributed_vector.h>
#include <quad_operators_3d.h>
#include <mumps_class.h>
#include <lis_class.h>
// #include <lis_class_distributed.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdio.h>

#include <p8est.h>
#include <random>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <unordered_map>
#include <set>
#include <tuple>
#include <limits>


void
poisson_boltzmann::create_mesh ()
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);

  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/ (e*e); //adim e_out
  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/ (e_0*e_out*kb*T);
  double k = std::sqrt (k2);

  int nx, ny, nz;
  double scale_tmp, scale_x, scale_y, scale_z;
  double lmax = 0;
  double max_len = 0;
  double l[3];
  scale_level = unilevel;
  double maxradius = *std::max_element (r_atoms.begin (), r_atoms.end ());

  auto comp_pos_x = [] (const std::array<double, 3>& a1, const std::array<double, 3>& a2) -> bool {
    return a1[0] < a2[0];
  };
  auto comp_pos_y = [] (const std::array<double, 3>& a1, const std::array<double, 3>& a2) -> bool {
    return a1[1] < a2[1];
  };
  auto comp_pos_z = [] (const std::array<double, 3>& a1, const std::array<double, 3>& a2) -> bool {
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
    std::cout << "[" << (*minmax_x.second)[0] - (*minmax_x.first)[0] + 2*maxradius << ", "
              << (*minmax_y.second)[1] - (*minmax_y.first)[1] + 2*maxradius << ", "
              << (*minmax_z.second)[2] - (*minmax_z.first)[2] + 2*maxradius << "]\n";

    if (std::fabs (net_charge - std::round (net_charge)) > 1.e-5)
      std::cerr << "  [WARNING] Net charge is not an integer: " << net_charge << '\n';

    std::cout << "  Solute epsilon     : " << e_in << '\n';
    std::cout << "  Solvent epsilon    : " << e_out << '\n';
    std::cout << "  Temperature        : " << T << " [K] \n";
    std::cout << "  Ionic strength     : " << ionic_strength << " [mol/L] \n";
    std::cout << "============================================\n\n";
  }

  // MESH_SHAPE_MEM sets its own box (build_membrane_slab_mesh): the membrane
  // box is centred on the lipid patch, not on the protein.
  if (mesh_shape != 2 && mesh_shape != MESH_SHAPE_MEM) {
    l_c[0] = (*minmax_x.first)[0] - maxradius - 2*prb_radius;
    l_c[1] = (*minmax_y.first)[1] - maxradius - 2*prb_radius;
    l_c[2] = (*minmax_z.first)[2] - maxradius - 2*prb_radius;
    r_c[0] = (*minmax_x.second)[0] + maxradius + 2*prb_radius;
    r_c[1] = (*minmax_y.second)[1] + maxradius + 2*prb_radius;
    r_c[2] = (*minmax_z.second)[2] + maxradius + 2*prb_radius;

    for (int kk = 0; kk < 3; ++kk) {
      l[kk] = (r_c[kk] - l_c[kk]);
      cc[kk] = (r_c[kk] + l_c[kk])*0.5;
      lmax = l[kk] > lmax ? l[kk] : lmax;
    }

    // For random displacement of the grid
    if (rand_center == 1) {
      std::random_device rd; // Will be used to obtain a seed for the random number engine
      std::mt19937 gen (rd ()); // Standard mersenne_twister_engine seeded with rd()
      std::uniform_real_distribution<> dis (-1./scale*0.5, 1./scale*0.5);
      double tmp_cc[3];

      for (int n = 0; n < 3; ++n) {
        tmp_cc[n] = cc[n];
        cc[n] = tmp_cc[n] + dis (gen);
      }

      std::cout << tmp_cc[0] << "  " << cc[0] << "  " << std::abs (tmp_cc[0] -cc[0]) << std::endl;
      std::cout << tmp_cc[1] << "  " << cc[1] << "  " << std::abs (tmp_cc[1] -cc[1]) << std::endl;
      std::cout << tmp_cc[2] << "  " << cc[2] << "  " << std::abs (tmp_cc[2] -cc[2]) << std::endl;
    }

    for (int kk = 0; kk < 3; ++kk) {
      ll[kk] = cc[kk] - lmax*0.5;
      rr[kk] = cc[kk] + lmax*0.5;
    }

    //stretched box with perfil1
    l_c[0] -= l[0]*0.5* (1.0/perfil1 - 1);
    l_c[1] -= l[1]*0.5* (1.0/perfil1 - 1);
    l_c[2] -= l[2]*0.5* (1.0/perfil1 - 1);
    r_c[0] += l[0]*0.5* (1.0/perfil1 - 1);
    r_c[1] += l[1]*0.5* (1.0/perfil1 - 1);
    r_c[2] += l[2]*0.5* (1.0/perfil1 - 1);

    //cubic box with perfil1
    ll[0] -= lmax*0.5* (1.0/perfil1 - 1);
    ll[1] -= lmax*0.5* (1.0/perfil1 - 1);
    ll[2] -= lmax*0.5* (1.0/perfil1 - 1);
    rr[0] += lmax*0.5* (1.0/perfil1 - 1);
    rr[1] += lmax*0.5* (1.0/perfil1 - 1);
    rr[2] += lmax*0.5* (1.0/perfil1 - 1);

    l_cr[0] = l_c[0];
    r_cr[0] = r_c[0];
    l_cr[1] = l_c[1];
    r_cr[1] = r_c[1];
    l_cr[2] = l_c[2];
    r_cr[2] = r_c[2];
  }

  if (mesh_shape == 0) {

    //cubic box with max perfil2
    double size = 1.0/scale;
    scale_level = 0;

    for (int ii = 0; ii<29; ++ii) {
      size *= 2;
      scale_level ++;

      if (lmax/size < perfil2)
        break;
    }

    ll[0] = cc[0] - size/2;
    ll[1] = cc[1] - size/2;
    ll[2] = cc[2] - size/2;
    rr[0] = cc[0] + size/2;
    rr[1] = cc[1] + size/2;
    rr[2] = cc[2] + size/2;


    nx = (int) ( (r_cr[0] - cc[0])*scale + 0.5);
    ny = (int) ( (r_cr[1] - cc[1])*scale + 0.5);
    nz = (int) ( (r_cr[2] - cc[2])*scale + 0.5);

    //refined box
    l_cr[0] = cc[0] - nx*1.0/scale;
    l_cr[1] = cc[1] - ny*1.0/scale;
    l_cr[2] = cc[2] - nz*1.0/scale;
    r_cr[0] = cc[0] + nx*1.0/scale;
    r_cr[1] = cc[1] + ny*1.0/scale;
    r_cr[2] = cc[2] + nz*1.0/scale;

    double dist = size/2 - lmax*0.5;
    pot_bc = std::exp (-k*dist)/ (dist*eps_out);

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

    simple_conn_p = std::make_unique<double[]> (simple_conn_num_vertices*3);
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

    for (int i = 0; i<6; i++)
      bcells.push_back (std::make_pair (0, i));
  } else if (mesh_shape == 1) {
    l_cr[0] = ll[0];
    l_cr[1] = ll[1];
    l_cr[2] = ll[2];
    r_cr[0] = rr[0];
    r_cr[1] = rr[1];
    r_cr[2] = rr[2];
    scale = (1<<unilevel)/ (rr[0]-ll[0]);

    if (refine_box == 1) {
      //cubic box with perfil2
      ll[0] = cc[0] - lmax*0.5*1.0/perfil2;
      ll[1] = cc[1] - lmax*0.5*1.0/perfil2;
      ll[2] = cc[2] - lmax*0.5*1.0/perfil2;
      rr[0] = cc[0] + lmax*0.5*1.0/perfil2;
      rr[1] = cc[1] + lmax*0.5*1.0/perfil2;
      rr[2] = cc[2] + lmax*0.5*1.0/perfil2;
      scale_tmp = (1<<unilevel)/ (rr[0]-ll[0]);

      nx = (int) ( (r_cr[0] - cc[0])*scale_tmp + 0.5);
      ny = (int) ( (r_cr[1] - cc[1])*scale_tmp + 0.5);
      nz = (int) ( (r_cr[2] - cc[2])*scale_tmp + 0.5);

      //refined stretched box
      l_cr[0] = cc[0] - nx*1.0/scale_tmp;
      l_cr[1] = cc[1] - ny*1.0/scale_tmp;
      l_cr[2] = cc[2] - nz*1.0/scale_tmp;
      r_cr[0] = cc[0] + nx*1.0/scale_tmp;
      r_cr[1] = cc[1] + ny*1.0/scale_tmp;
      r_cr[2] = cc[2] + nz*1.0/scale_tmp;
    }

    double dist = ((rr[0]-ll[0]) - lmax)*0.5;
    pot_bc = std::exp (-k*dist)/ (dist*eps_out);

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

    simple_conn_p = std::make_unique<double[]> (simple_conn_num_vertices*3);
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

    for (int i = 0; i<6; i++)
      bcells.push_back (std::make_pair (0, i));
  } else if (mesh_shape == 2) {

    l_cr[0] = l_c[0];
    l_cr[1] = l_c[1];
    l_cr[2] = l_c[2];
    r_cr[0] = r_c[0];
    r_cr[1] = r_c[1];
    r_cr[2] = r_c[2];
    scale = (1<<unilevel)/ (r_c[0]-l_c[0]);

    double dist = ((r_cr[0]-l_cr[0]) - lmax)*0.5;
    pot_bc = std::exp (-k*dist)/ (dist*eps_out);

    if (rank == 0) {
      std::cout << "========== [ Domain Information ] ==========\n";

      std::cout << "  Scale:  " << scale << "\n";

      std::cout << "  Center of the System [Å]:";
      std::cout << "  [" << (r_c[0] + l_c[0])*0.5 << ", "
                << (r_c[1] + l_c[1])*0.5 << ", "
                << (r_c[2] + l_c[2])*0.5 << "]\n";

      std::cout << "  Complete Domain Box Size [Å]:\n";
      std::cout << "      x = [" << l_cr[0] << ", " << r_cr[0] << "]\n";
      std::cout << "      y = [" << l_cr[1] << ", " << r_cr[1] << "]\n";
      std::cout << "      z = [" << l_cr[2] << ", " << r_cr[2] << "]\n";

      std::cout << "  Number of Subdivisions:";
      std::cout << "  nx = " << (1<<unilevel) << "  ny = " << (1<<unilevel) << "  nz = " << (1<<unilevel) << '\n';

      std::cout << "============================================\n";
    }

    if (refine_box == 1) {
      double nx_tmp, ny_tmp, nz_tmp;
      double cc_tmp[3];

      scale_x = (1<<unilevel)/ (r_c[0]-l_c[0]);
      scale_y = (1<<unilevel)/ (r_c[1]-l_c[1]);
      scale_z = (1<<unilevel)/ (r_c[2]-l_c[2]);

      cc[0] = (r_c[0]+l_c[0])*0.5;
      cc[1] = (r_c[1]+l_c[1])*0.5;
      cc[2] = (r_c[2]+l_c[2])*0.5;

      cc_tmp[0] = (r_cr[0]+l_cr[0])*0.5;
      cc_tmp[1] = (r_cr[1]+l_cr[1])*0.5;
      cc_tmp[2] = (r_cr[2]+l_cr[2])*0.5;

      nx_tmp = (int) ( (cc[0] - cc_tmp[0])*scale_x + 0.5);
      ny_tmp = (int) ( (cc[1] - cc_tmp[1])*scale_y + 0.5);
      nz_tmp = (int) ( (cc[2] - cc_tmp[2])*scale_z + 0.5);

      cc_tmp[0] = cc[0] + nx_tmp*1.0/scale_x;
      cc_tmp[1] = cc[1] + ny_tmp*1.0/scale_y;
      cc_tmp[2] = cc[2] + nz_tmp*1.0/scale_z;

      nx = (int) ( (r_cr[0] - cc_tmp[0])*scale_x + 0.5);
      ny = (int) ( (r_cr[1] - cc_tmp[1])*scale_y + 0.5);
      nz = (int) ( (r_cr[2] - cc_tmp[2])*scale_z + 0.5);
      //refined stretched box
      l_cr[0] = cc_tmp[0] - nx*1.0/scale_x;
      l_cr[1] = cc_tmp[1] - ny*1.0/scale_y;
      l_cr[2] = cc_tmp[2] - nz*1.0/scale_z;
      r_cr[0] = cc_tmp[0] + nx*1.0/scale_x;
      r_cr[1] = cc_tmp[1] + ny*1.0/scale_y;
      r_cr[2] = cc_tmp[2] + nz*1.0/scale_z;

      std::cout << "xb: " << l_cr[0] << ", " << r_cr[0] << std::endl;
      std::cout << "yb: " << l_cr[1] << ", " << r_cr[1] << std::endl;
      std::cout << "zb: " << l_cr[2] << ", " << r_cr[2] << "\n" << std::endl;
    }

    simple_conn_num_vertices = 8;
    simple_conn_num_trees = 1;

    simple_conn_p = std::make_unique<double[]> (simple_conn_num_vertices*3);
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

    for (int i = 0; i<6; i++)
      bcells.push_back (std::make_pair (0, i));

  } else if (mesh_shape == 3) {
    //cubic box with max perfil2
    double size = 1.0/scale;
    double scale_min_box = 0.5;
    scale_level = 0;
    scale_level_min_box = 0;


    for (int ii = 0; ii<28; ++ii) {
      size *= 2;
      scale_level ++;

      if (lmax/size < perfil2)
        break;
    }

    for (int ii = 0; ii<scale_level; ++ii) {
      scale_level_min_box ++;

      if ( (std::pow (2,scale_level_min_box)+1)/size >= scale_min_box)
        break;
    }

    ll[0] = cc[0] - size/2;
    ll[1] = cc[1] - size/2;
    ll[2] = cc[2] - size/2;
    rr[0] = cc[0] + size/2;
    rr[1] = cc[1] + size/2;
    rr[2] = cc[2] + size/2;

    //refined box NS
    nx = (int) ( (r_cr[0] - cc[0])*scale + 0.5);
    ny = (int) ( (r_cr[1] - cc[1])*scale + 0.5);
    nz = (int) ( (r_cr[2] - cc[2])*scale + 0.5);

    l_cr[0] = cc[0] - nx*1.0/scale;
    l_cr[1] = cc[1] - ny*1.0/scale;
    l_cr[2] = cc[2] - nz*1.0/scale;
    r_cr[0] = cc[0] + nx*1.0/scale;
    r_cr[1] = cc[1] + ny*1.0/scale;
    r_cr[2] = cc[2] + nz*1.0/scale;

    //refined box FOCUS

    double err_l;
    l_box[0] = cc_focusing[0] - std::round (n_grid/2.)/scale;
    l_box[1] = cc_focusing[1] - std::round (n_grid/2.)/scale;
    l_box[2] = cc_focusing[2] - std::round (n_grid/2.)/scale;
    r_box[0] = cc_focusing[0] + std::round (n_grid/2.)/scale;
    r_box[1] = cc_focusing[1] + std::round (n_grid/2.)/scale;
    r_box[2] = cc_focusing[2] + std::round (n_grid/2.)/scale;


    if ( (r_box[0]- l_box[0])>= (r_cr[0] - l_cr[0])) {
      r_box[0] = r_cr[0];
      l_box[0] = l_cr[0];
    } else {
      if (l_box[0] <= l_cr[0]) {
        err_l = l_cr[0]-l_box[0];
        l_box[0] = l_box[0] + err_l;
        r_box[0] = r_box[0] + err_l;
      }

      if (r_box[0] >= r_cr[0]) {
        err_l = r_box[0]-r_cr[0];
        l_box[0] = l_box[0] - err_l;
        r_box[0] = r_box[0] - err_l;
      }
    }

    if ( (r_box[1]- l_box[1])>= (r_cr[1] - l_cr[1])) {
      r_box[1] = r_cr[1];
      l_box[1] = l_cr[1];
    } else {
      if (l_box[1] <= l_cr[1]) {
        err_l = l_cr[1]-l_box[1];
        l_box[1] = l_box[1] + err_l;
        r_box[1] = r_box[1] + err_l;
      }

      if (r_box[1] >= r_cr[1]) {
        err_l = r_box[1]-r_cr[1];
        l_box[1] = l_box[1] - err_l;
        r_box[1] = r_box[1] - err_l;
      }
    }

    if ( (r_box[2]- l_box[2])>= (r_cr[2] - l_cr[2])) {
      r_box[2] = r_cr[2];
      l_box[2] = l_cr[2];
    } else {
      if (l_box[2] <= l_cr[2]) {
        err_l = l_cr[2]-l_box[2];
        l_box[2] = l_box[2] + err_l;
        r_box[2] = r_box[2] + err_l;
      }

      if (r_box[2] >= r_cr[2]) {
        err_l = r_box[2]-r_cr[2];
        r_box[2] = r_box[2] - err_l;
      }
    }

    //calcolo outlevel
    unsigned int ratio_l_b;
    double len_box_foc = n_grid/scale;
    double len_box = rr[0]-ll[0];

    for (int ii = 1; ii<=5; ++ii) {
      len_box /= 2.0;

      if (len_box <= len_box_foc) {
        ratio_l_b =ii;
        break;
      }
    }

    outlevel = ratio_l_b;
    // outlevel = 1;
    ////////////////////////////////////////


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

    simple_conn_p = std::make_unique<double[]> (simple_conn_num_vertices*3);
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

    for (int i = 0; i<6; i++)
      bcells.push_back (std::make_pair (0, i));
  } else if (mesh_shape == 4) {

    //cubic box with max perfil2

    double size = 1.0/scale_max;

    scale = scale_max;
    maxlevel = 0;

    for (int ii = 0; ii<28; ++ii) {
      std::cout << size << "  " << maxlevel << std::endl;
      size *= 2;
      maxlevel ++;

      if (lmax/size < perfil2)
        break;
    }

    minlevel = (int) (maxlevel - std::sqrt (scale_max/scale_min));

    unilevel = (int) ( (maxlevel + minlevel)*0.5);
    unilevel = unilevel + 1;
    outlevel = minlevel;
    scale_level = maxlevel;

    ll[0] = cc[0] - size/2;
    ll[1] = cc[1] - size/2;
    ll[2] = cc[2] - size/2;
    rr[0] = cc[0] + size/2;
    rr[1] = cc[1] + size/2;
    rr[2] = cc[2] + size/2;

    //refined box NS
    nx = (int) ( (r_cr[0] - cc[0])*scale + 0.5);
    ny = (int) ( (r_cr[1] - cc[1])*scale + 0.5);
    nz = (int) ( (r_cr[2] - cc[2])*scale + 0.5);

    l_cr[0] = cc[0] - nx*1.0/scale;
    l_cr[1] = cc[1] - ny*1.0/scale;
    l_cr[2] = cc[2] - nz*1.0/scale;
    r_cr[0] = cc[0] + nx*1.0/scale;
    r_cr[1] = cc[1] + ny*1.0/scale;
    r_cr[2] = cc[2] + nz*1.0/scale;


    if (rank == 0) {
      std::cout << "cx: " << cc[0]
                << " , cy: " << cc[1]
                << " , cz: " << cc[2] << std::endl;

      std::cout << "x: " << ll[0] << ", " << rr[0] << std::endl;
      std::cout << "y: " << ll[1] << ", " << rr[1] << std::endl;
      std::cout << "z: " << ll[2] << ", " << rr[2] << std::endl;

      std::cout << minlevel <<" "<<maxlevel << "  " << unilevel <<std::endl;
    }

    simple_conn_num_vertices = 8;
    simple_conn_num_trees = 1;

    simple_conn_p = std::make_unique<double[]> (simple_conn_num_vertices*3);
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

    for (int i = 0; i<6; i++)
      bcells.push_back (std::make_pair (0, i));
  } else if (mesh_shape == 5) {
    l_cr[0] = cc[0] - len/2;
    l_cr[1] = cc[1] - len/2;
    l_cr[2] = cc[2] - len/2;
    r_cr[0] = cc[0] + len/2;
    r_cr[1] = cc[1] + len/2;
    r_cr[2] = cc[2] + len/2;


    if (unilevel == 0)
      scale = (num_trees[0])/ (len);
    else
      scale = (num_trees[0]* (1<<unilevel))/ (len);

    //////////////////////////////
    num_trees[1] = num_trees[0];
    num_trees[2] = num_trees[0];

    /////////////////////////////
    if (rank == 0) {
      std::cout << "x: " << l_cr[0] << ", " << r_cr[0] << std::endl;
      std::cout << "y: " << l_cr[1] << ", " << r_cr[1] << std::endl;
      std::cout << "z: " << l_cr[2] << ", " << r_cr[2] << "\n" << std::endl;
      std::cout << "Number of trees: " << num_trees[0] << std::endl;
    }

    double bound_x = std::abs (r_cr[0]-l_cr[0]);
    double bound_y = std::abs (r_cr[1]-l_cr[1]);
    double bound_z = std::abs (r_cr[2]-l_cr[2]);
    double step[3] = { bound_x/num_trees[0],
                       bound_y/num_trees[1],
                       bound_z/num_trees[2]
                     };
    double bound[3] = {bound_x, bound_y, bound_z};
    make_connectivity_3d (num_trees, step, simple_conn_p,
                          simple_conn_num_vertices, simple_conn_t,
                          simple_conn_num_trees, bcells);

    for (p4est_topidx_t i =0; i < simple_conn_num_vertices; ++i) {
      p4est_topidx_t j = 0;
      simple_conn_p[3*i + j++] += l_cr[0];
      simple_conn_p[3*i + j++] += l_cr[1];
      simple_conn_p[3*i + j] += l_cr[2];
    }
  } else if (mesh_shape == MESH_SHAPE_MEM) {
    build_membrane_slab_mesh (*this);
  } else {
    if (rank == 0) {
      std::cout << "x: " << ll[0] << ", " << rr[0] << std::endl;
      std::cout << "y: " << ll[1] << ", " << rr[1] << std::endl;
      std::cout << "z: " << ll[2] << ", " << rr[2] << "\n" << std::endl;
    }

    simple_conn_num_vertices = 8;
    simple_conn_num_trees = 1;

    simple_conn_p = std::make_unique<double[]> (simple_conn_num_vertices*3);
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

    for (int i = 0; i<6; i++)
      bcells.push_back (std::make_pair (0, i));
  }

  tmsh.read_connectivity (simple_conn_p.get (), simple_conn_num_vertices,
                          simple_conn_t.get (), simple_conn_num_trees);

}


// ============================================================================
// [PBC] Periodic boundary conditions on x+-/y+- (strong enforcement).
// Ported from branch `membrane` (src/pb_mesh.cpp), inline here per the
// membrane-clean layout decision (no module split). Mortar scaffold NOT
// ported: only the coordinate-based node map / conformity checks needed for
// strong enforcement.
// ============================================================================

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

    // right_all/left_all can contain the same physical node more than once:
    // a node shared between two face quadrants owned by different ranks gets
    // reported by each of them (the per-rank "seen" dedup above only catches
    // repeats within a single rank's own sweep). That inflates the raw
    // iteration count below, but is harmless here: every repeat of the same
    // right-face node resolves to the same left partner, so the overwrite in
    // map_out is idempotent. What actually matters is p.map_out->size (),
    // the number of distinct right-face nodes mapped.
    int unmatched = 0;
    p.map_out->clear ();
    for (size_t i = 0; i + 3 <= right_all.size (); i += 3) {
      auto key = std::make_pair (right_all[i], right_all[i + 1]);
      auto it  = left_map.find (key);
      if (it != left_map.end ())
        (*p.map_out)[static_cast<size_t> (right_all[i + 2])] = it->second;
      else
        ++unmatched;
    }

    if (rank == 0) {
      std::cout << "  [PBC map] " << p.label << ": "
                << p.map_out->size () << " matched";
      if (unmatched > 0)
        std::cout << ", " << unmatched << " UNMATCHED (mesh not conforming?)";
      std::cout << '\n';
    }
  }
}

// ============================================================================
// [PBC mortar] Weak periodic enforcement via Lagrange multipliers (K_static).
// Ported from branch `membrane` (src/pb_membrane.cpp assemble_mortar_block),
// stripped of the LIS-only scaffold (mortar rows written into A, diagonal
// regularization): the augmented saddle-point system is solved with MUMPS only,
// so the C^T block lives entirely in the dense mortar_C accumulator that
// mumps_compute_electric_potential converts to COO. Formulation: bc.md
// § "Derivazione del sistema aumentato (K_static)" onward.
// ============================================================================
void
poisson_boltzmann::assemble_mortar_system ()
{
  int size, rank;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);

  const int    nm       = (1 << minlevel);
  ndofm                 = (nm + 1) * (nm + 1);
  const int    n_pairs  = (periodic_x ? 1 : 0) + (periodic_y ? 1 : 0);
  const size_t N_global = static_cast<size_t> (tmsh.num_global_nodes ());

  const double Lx = r_cr[0] - l_cr[0];
  const double Ly = r_cr[1] - l_cr[1];
  const double Lz = r_cr[2] - l_cr[2];

  // Each pair carries its tangential directions (d0, d1), physical extents, and
  // the offsets of its DOFs: a global column base in A (N_global + cumulative)
  // and a row base in the dense C^T accumulator (cumulative).
  struct MortarFacePair {
    int    face_left, face_right, d0, d1;
    double d0_min, L0, d1_min, L1;
    size_t mortar_col_offset;    // global col in A: N_global + cumulative
    size_t mortar_dense_offset;  // row in mortar_C: cumulative
  };

  std::vector<MortarFacePair> pairs;
  size_t cumulative = 0;
  if (periodic_x) {
    pairs.push_back ({0, 1, 1, 2, l_cr[1], Ly, l_cr[2], Lz,
                      N_global + cumulative, cumulative});
    cumulative += static_cast<size_t> (ndofm);
  }
  if (periodic_y) {
    pairs.push_back ({2, 3, 0, 2, l_cr[0], Lx, l_cr[2], Lz,
                      N_global + cumulative, cumulative});
    cumulative += static_cast<size_t> (ndofm);
  }

  mortar_C.assign (static_cast<size_t> (n_pairs) * ndofm * N_global, 0.0);

  // Accumulate the C / C^T contributions of one boundary quadrant on one face.
  // sign = +1 on the left face, -1 on the right face. The integral of a FEM hat
  // against a mortar hat factorizes as a tensor product of two 1D mixed masses.
  auto assemble_block = [&] (tmesh_3d::quadrant_iterator q, int face_idx,
                             const MortarFacePair &pair, double sign) {
    const auto &fnodes = q->ef (face_idx);
    if (fnodes.empty ()) return;

    const int d0 = pair.d0, d1 = pair.d1;

    // Step 1 — tangential coords of the 4 face nodes.
    //   fnodes: [0]=(d0min,d1min) [1]=(d0max,d1min) [2]=(d0min,d1max) [3]=(d0max,d1max)
    double c0[2], c1[2];
    c0[0] = q->p (d0, fnodes[0]);
    c0[1] = q->p (d0, fnodes[1]);
    c1[0] = q->p (d1, fnodes[0]);
    c1[1] = q->p (d1, fnodes[2]);

    // Step 2 — Dirichlet filter: drop nodes on the z± planes (d1 == 2).
    const double z_lo = pair.d1_min;
    const double z_hi = pair.d1_min + pair.L1;
    const double tol  = pair.L1 * 1e-10;
    bool keep[4];
    for (int fi = 0; fi < 4; ++fi) {
      double zi = q->p (2, fnodes[fi]);
      keep[fi] = (zi > z_lo + tol) && (zi < z_hi - tol);
    }

    // Step 3 — mortar cell containing this quadrant (clamped to the last cell).
    const double hm0 = pair.L0 / nm;
    const double hm1 = pair.L1 / nm;
    const int j0 = std::min (static_cast<int> ((c0[0] - pair.d0_min) / hm0), nm - 1);
    const int j1 = std::min (static_cast<int> ((c1[0] - pair.d1_min) / hm1), nm - 1);

    // Step 4 — 1D mixed hat-function masses (2x2 each), exact.
    const double orig0 = pair.d0_min + j0 * hm0;
    const double orig1 = pair.d1_min + j1 * hm1;

    double psi_0[2][2], psi_1[2][2];
    psi_0[0][0] = (hm0 - (c0[0] - orig0)) / hm0;
    psi_0[0][1] = (hm0 - (c0[1] - orig0)) / hm0;
    psi_0[1][0] =        (c0[0] - orig0)  / hm0;
    psi_0[1][1] =        (c0[1] - orig0)  / hm0;
    psi_1[0][0] = (hm1 - (c1[0] - orig1)) / hm1;
    psi_1[0][1] = (hm1 - (c1[1] - orig1)) / hm1;
    psi_1[1][0] =        (c1[0] - orig1)  / hm1;
    psi_1[1][1] =        (c1[1] - orig1)  / hm1;

    const double h0 = c0[1] - c0[0], h1 = c1[1] - c1[0];
    double M_0[2][2], M_1[2][2];
    // int_a^b f*g dt = (b-a)/6 * (2 f_a g_a + f_a g_b + f_b g_a + 2 f_b g_b)
    for (int jl = 0; jl < 2; ++jl) {
      M_0[0][jl] = h0 / 6.0 * (2.0 * psi_0[jl][0] +       psi_0[jl][1]);
      M_0[1][jl] = h0 / 6.0 * (      psi_0[jl][0] + 2.0 * psi_0[jl][1]);
      M_1[0][jl] = h1 / 6.0 * (2.0 * psi_1[jl][0] +       psi_1[jl][1]);
      M_1[1][jl] = h1 / 6.0 * (      psi_1[jl][0] + 2.0 * psi_1[jl][1]);
    }

    // Step 5 — accumulate C (into A) and C^T (into mortar_C) by tensor product.
    static const int ii0[4] = {0, 1, 0, 1};
    static const int ii1[4] = {0, 0, 1, 1};
    for (int fi = 0; fi < 4; ++fi) {
      if (!keep[fi]) continue;
      const size_t phys = static_cast<size_t> (q->gt (fnodes[fi]));
      for (int dj0 = 0; dj0 < 2; ++dj0)
        for (int dj1 = 0; dj1 < 2; ++dj1) {
          const size_t k    = static_cast<size_t> ((j1 + dj1) * (nm + 1) + (j0 + dj0));
          const size_t gcol = pair.mortar_col_offset  + k;
          const size_t drow = pair.mortar_dense_offset + k;
          const double val  = sign * M_0[ii0[fi]][dj0] * M_1[ii1[fi]][dj1];
          (*A)[phys][gcol]                  += val;  // C block  (physical row, mortar col)
          mortar_C[drow * N_global + phys]  += val;  // C^T block (dense, for MUMPS rows)
        }
    }
  };

  for (auto &pair : pairs) {
    for (auto q = tmsh.begin_quadrant_sweep (); q != tmsh.end_quadrant_sweep (); ++q) {
      if (!q->ef (pair.face_left).empty ())
        assemble_block (q, pair.face_left,  pair, +1.0);
      if (!q->ef (pair.face_right).empty ())
        assemble_block (q, pair.face_right, pair, -1.0);
    }
  }

  // Sum the dense C^T contributions from every rank onto rank 0, which owns the
  // mortar rows of the augmented COO built at solve time.
  if (size > 1)
    MPI_Reduce (rank == 0 ? MPI_IN_PLACE : mortar_C.data (), mortar_C.data (),
                static_cast<int> (mortar_C.size ()), MPI_DOUBLE, MPI_SUM, 0, mpicomm);

  // Mortar DOFs on the z± edges of the mortar grid (j1 == 0 or nm) sit on the
  // Dirichlet boundary: no periodic constraint there. Zero their C^T rows now;
  // the identity that pins F = 0 for them is added with the COO rows at solve
  // time (see mumps_compute_electric_potential).
  int n_bd = 0;
  if (rank == 0)
    for (auto &pair : pairs)
      for (int j0 = 0; j0 <= nm; ++j0)
        for (int j1_bd : {0, nm}) {
          const size_t k    = static_cast<size_t> (j1_bd) * (nm + 1) + j0;
          const size_t drow = pair.mortar_dense_offset + k;
          for (size_t c = 0; c < N_global; ++c)
            mortar_C[drow * N_global + c] = 0.0;
          ++n_bd;
        }

  if (rank == 0)
    std::cout << "  [Mortar] pairs=" << n_pairs << "  ndofm=" << ndofm
              << " (nm=" << nm << ")  z-boundary DOFs zeroed=" << n_bd << '\n';
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
    for (size_t i = 0; i + 5 <= left_buf.size (); i += 5) {
      Key key {left_buf[i], left_buf[i+1], left_buf[i+2], left_buf[i+3]};
      auto it = right_map.find (key);
      if (it != right_map.end () && it->second == static_cast<int> (left_buf[i+4]))
        ++n_match;
      else
        ++n_mismatch;
    }

    std::cout << "  [PBC conformity] " << label << ": "
              << n_match << " matching, " << n_mismatch << " mismatching quadrant(s)\n";
  };

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

// ----------------------------------------------------------------------------
// [PBC] Strong enforcement.
//
// The identification happens during assembly: canonical[] is handed to the
// bimpp operators as the row and column ordering, so a right-face node's
// contributions land directly on its left partner, and the ones aimed at rows
// owned by another rank are shipped by A->assemble() like any other shared
// node. The rows left empty by the fold are closed with an identity and their
// value is restored afterwards by pbc_apply_periodic_copy().
//
// Building canonical[] needs no communication: it only depends on
// pbc_x/y_right_to_left, which build_pbc_node_map already made identical on
// every rank via MPI_Allgatherv.
// ----------------------------------------------------------------------------

void
poisson_boltzmann::pbc_build_canonical (std::vector<int> &canonical)
{
  const size_t N_phys = static_cast<size_t> (tmsh.num_global_nodes ());

  canonical.resize (N_phys);
  for (size_t i = 0; i < N_phys; ++i)
    canonical[i] = static_cast<int> (i);

  // y-pairs first, x-pairs second: a vertical-edge node (on both a periodic
  // x-face and a periodic y-face, e.g. x=Lx ∧ y=Ly) is a key in both maps.
  // The y-pass alone would leave it pointing at another right-face node
  // (x=Lx ∧ y=0); the x-pass then resolves *that* node's own canonical
  // target (already fixed by the y-pass), completing the chain down to the
  // single x=0 ∧ y=0 corner.
  for (auto const &kv : pbc_y_right_to_left)
    canonical[kv.first] = static_cast<int> (kv.second);
  for (auto const &kv : pbc_x_right_to_left)
    canonical[kv.first] = canonical[kv.second];

  int rank;
  MPI_Comm_rank (mpicomm, &rank);
  if (rank == 0) {
    size_t folded = 0;
    for (size_t i = 0; i < N_phys; ++i)
      if (canonical[i] != static_cast<int> (i))
        ++folded;
    std::cout << "  [PBC strong] " << folded
              << " right-face nodes folded onto their left partner\n";
  }
}

void
poisson_boltzmann::pbc_apply_periodic_copy ()
{
  if (!pbc_use_ordering || !phi)
    return;

  int size, rank;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);

  auto &data = phi->get_owned_data ();
  const int phi_is = phi->get_range_start ();
  const int phi_ie = phi_is + static_cast<int> (data.size ());

  // Reading a canonical node while writing a folded one is safe in any order:
  // canonical[l] == l for every target l, so no target is ever a destination.
  if (size == 1) {
    for (int r = phi_is; r < phi_ie; ++r)
      if (pbc_canonical[r] != r)
        data[static_cast<size_t> (r - phi_is)] =
          data[static_cast<size_t> (pbc_canonical[r] - phi_is)];
    return;
  }

  // The periodic partner sits on the opposite face of the box, so it is not a
  // ghost of r and bim3a_solution_with_ghosts cannot supply it. Only the
  // canonical left nodes are needed -- a few thousand values -- so one
  // Allreduce over that set is enough. pbc_canonical is identical on every
  // rank, hence so is the target list, and every global node is owned by
  // exactly one rank, so summing reconstructs the values exactly.
  std::vector<int> targets;
  for (size_t i = 0; i < pbc_canonical.size (); ++i)
    if (pbc_canonical[i] != static_cast<int> (i))
      targets.push_back (pbc_canonical[i]);
  std::sort (targets.begin (), targets.end ());
  targets.erase (std::unique (targets.begin (), targets.end ()), targets.end ());

  std::vector<double> buf (targets.size (), 0.0);
  for (size_t k = 0; k < targets.size (); ++k)
    if (targets[k] >= phi_is && targets[k] < phi_ie)
      buf[k] = data[static_cast<size_t> (targets[k] - phi_is)];
  MPI_Allreduce (MPI_IN_PLACE, buf.data (), static_cast<int> (buf.size ()),
                 MPI_DOUBLE, MPI_SUM, mpicomm);

  for (int r = phi_is; r < phi_ie; ++r) {
    if (pbc_canonical[r] == r)
      continue;
    auto it = std::lower_bound (targets.begin (), targets.end (), pbc_canonical[r]);
    data[static_cast<size_t> (r - phi_is)] =
      buf[static_cast<size_t> (it - targets.begin ())];
  }
}


double
poisson_boltzmann::levelsetfun (double x, double y, double z)
{
  double dist = 0.0;

  for (const NS::Atom& i : atoms) {

    dist += std::exp (surf_param * ( (std::pow (x - i.pos[0], 2) +
                                      std::pow (y - i.pos[1], 2) +
                                      std::pow (z - i.pos[2], 2)) /
                                     (i.radius*i.radius) - 1.0));

    if (dist > 1.5)
      break;

  }

  return dist;
}


double
poisson_boltzmann::is_in_ns_surf (ray_cache_t & ray_cache, double x, double y, double z, int dir)
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

  crossings_t & ct = ray_cache (x1, x2, dir);

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


int
poisson_boltzmann::parse_options (int argc, char **argv)
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);

  GetPot g (argc, argv);

  // =============================
  // 1. Leggi file dei parametri
  // =============================
  if (!g.search ("--prmfile") && !g.search ("--potfile")) {
    if (rank == 0)
      std::cout << "Warning: No parameters file selected, using the default one."
                << "\nTo select one use --prmfile or --potfile option followed by the desired file.\n";
  }

  // Cerca il file, dando precedenza a --prmfile se entrambi sono presenti
  if (g.search ("--prmfile")) {
    optionsfilename = g.next ("../../data/options.prm");
  } else if (g.search ("--potfile")) {
    optionsfilename = g.next ("../../data/options.prm");
  }

  if (rank == 0)
    std::cout << "Selected parameters file: " << optionsfilename << std::endl;

  std::ifstream optionsfile (optionsfilename);

  if (!optionsfile) {
    if (rank == 0)
      std::cerr << "Cannot find the options file" << std::endl;

    return 1;
  }

  GetPot g2 (optionsfilename.c_str());

  // =============================
  // 2. Leggi parametri input/
  // =============================
  const std::string input_section = "input/";
  std::string filename_from_file;

  filetype = g2 ((input_section + "filetype").c_str(), "pqr");
  pqrfilename = g2 ((input_section + "filename").c_str(), "../../data/1CCM.pqr");
  radiusfilename = g2 ((input_section + "radius_file").c_str(), "../../data/radius.siz");
  chargefilename = g2 ((input_section + "charge_file").c_str(), "../../data/charge.crg");
  write_pqr = g2 ((input_section + "write_pqr").c_str(), 0);
  name_pqr = g2 ((input_section + "name_pqr").c_str(), "output.pqr");



  // =============================
  // 3. Override da riga di comando: --pqrfile
  // =============================
  if (g.search ("--pqrfile")) {
    pqrfilename = g.next ("");
    filetype = "pqr"; // Forza il filetype a pqr se viene specificato un file
  }

  if (rank == 0)
    std::cout << "Selected molecule file:   " << pqrfilename << std::endl;

  std::ifstream pqrfile (pqrfilename);

  if (!pqrfile) {
    if (rank == 0) {
      std::cerr << "Cannot find the pqr file" << std::endl;
      return 1;
    }
  }

  const std::string mesh_options = "mesh/";
  maxlevel = g2 ( (mesh_options + "maxlevel").c_str (), 9);
  minlevel = g2 ( (mesh_options + "minlevel").c_str (), 3);
  unilevel = g2 ( (mesh_options + "unilevel").c_str (), 5);
  outlevel = g2 ( (mesh_options + "outlevel").c_str (), 1);
  loc_refinement = g2 ( (mesh_options + "loc_refinement").c_str (), 0);
  mesh_shape = g2 ( (mesh_options + "mesh_shape").c_str (), 1);
  refine_box = g2 ( (mesh_options + "refine_box").c_str (), 0);
  rand_center = g2 ( (mesh_options + "rand_center").c_str (), 0);
  aligned = g2 ( (mesh_options + "aligned").c_str (), 0);

  if (mesh_shape < 2) {
    perfil1 = g2 ( (mesh_options + "perfil1").c_str (), 0.8);
    perfil2 = g2 ( (mesh_options + "perfil2").c_str (), 0.2);
    scale = g2 ( (mesh_options + "scale").c_str (), 2.0);
  }

  if (mesh_shape == 2) {
    l_c[0] = g2 ( (mesh_options + "x1").c_str (), -128.0);
    r_c[0] = g2 ( (mesh_options + "x2").c_str (), 128.0);
    l_c[1] = g2 ( (mesh_options + "y1").c_str (), -128.0);
    r_c[1] = g2 ( (mesh_options + "y2").c_str (), 128.0);
    l_c[2] = g2 ( (mesh_options + "z1").c_str (), -128.0);
    r_c[2] = g2 ( (mesh_options + "z2").c_str (), 128.0);
    l_cr[0] = g2 ( (mesh_options + "refine_x1").c_str (), -64.0);
    r_cr[0] = g2 ( (mesh_options + "refine_x2").c_str (), 64.0);
    l_cr[1] = g2 ( (mesh_options + "refine_y1").c_str (), -64.0);
    r_cr[1] = g2 ( (mesh_options + "refine_y2").c_str (), 64.0);
    l_cr[2] = g2 ( (mesh_options + "refine_z1").c_str (), -64.0);
    r_cr[2] = g2 ( (mesh_options + "refine_z2").c_str (), 64.0);
  }

  if (mesh_shape == 3) {
    perfil1 = g2 ( (mesh_options + "perfil1").c_str (), 0.8);
    perfil2 = g2 ( (mesh_options + "perfil2").c_str (), 0.2);
    scale = g2 ( (mesh_options + "scale").c_str (), 2.0);
    cc_focusing[0] = g2 ( (mesh_options + "cx_foc").c_str (), 0.0);
    cc_focusing[1] = g2 ( (mesh_options + "cy_foc").c_str (), 0.0);
    cc_focusing[2] = g2 ( (mesh_options + "cz_foc").c_str (), 0.0);
    n_grid = g2 ( (mesh_options + "n_grid").c_str (), 8);
  }

  if (mesh_shape == 4) {
    perfil1 = g2 ( (mesh_options + "perfil1").c_str (), 0.8);
    perfil2 = g2 ( (mesh_options + "perfil2").c_str (), 0.2);
    scale_min = g2 ( (mesh_options + "scale_min").c_str (), 0.5);
    scale_max = g2 ( (mesh_options + "scale_max").c_str (), 2.0);
  }

  if (mesh_shape == 5) {
    num_trees [0] = g2 ( (mesh_options + "num_trees_x").c_str (), 10);
    num_trees [1] = g2 ( (mesh_options + "num_trees_y").c_str (), 10);
    num_trees [2] = g2 ( (mesh_options + "num_trees_z").c_str (), 10);
    len = g2 ( (mesh_options + "lato").c_str (), 50.0);
    perfil1 = g2 ( (mesh_options + "perfil1").c_str (), 0.8);
  }

  periodic_x = g2 ( (mesh_options + "periodic_x").c_str (), 0);
  periodic_y = g2 ( (mesh_options + "periodic_y").c_str (), 0);

  const std::string model_options = "model/";
  linearized = g2 ( (model_options + "linearized").c_str (), 1);
  bc = g2 ( (model_options + "bc_type").c_str (), 1);
  ionic_strength = g2 ( (model_options + "ionic_strength").c_str (), 0.145);
  e_in = g2 ( (model_options + "molecular_dielectric_constant").c_str (), 2.);
  e_out = g2 ( (model_options + "solvent_dielectric_constant").c_str (), 80.);
  T = g2 ( (model_options + "T").c_str (), 298.15);
  calc_energy = g2 ( (model_options + "calc_energy").c_str (), 2);
  calc_coulombic = g2 ( (model_options + "calc_coulombic").c_str (), 0);
  calc_potential_term = g2 ( (model_options + "calc_potential_terms").c_str (), 0);
  calc_field_term = g2 ( (model_options + "calc_field_terms").c_str (), 0);
  atoms_write = g2 ( (model_options + "atoms_write").c_str (), 0);
  surf_write = g2 ( (model_options + "surf_write").c_str (), 0);
  surf_write = g2 ( (model_options + "surf_write").c_str (), 0);
  map_type = g2 ( (model_options + "map_type").c_str (), "vtu");
  potential_map = g2 ( (model_options + "potential_map").c_str (), 0);
  eps_map = g2 ( (model_options + "eps_map").c_str (), 0);
  dataset_write = g2 ( (model_options + "dataset_write").c_str (), 0);
  const std::string surf_options = "surface/";
  surf_type_num = g2 ( (surf_options + "surface_type").c_str (), 0);

  if (surf_type_num == 1) surf_type = NS::skin;
  else if (surf_type_num == 0) surf_type = NS::ses;
  // else if (surf_type_num == 2) surf_type = NS::blobby;
  else surf_type = NS::ses;

  surf_param = g2 ( (surf_options + "surface_parameter").c_str (), 0.45);
  stern_layer_surf = g2 ( (surf_options + "stern_layer_surf").c_str (), 0);
  stern_layer = g2 ( (surf_options + "stern_layer_thickness").c_str (), 2.);
  num_threads = g2 ( (surf_options + "number_of_threads").c_str (), 1);

  const std::string alg_options = "algorithm/";
  linear_solver_name = g2 ( (alg_options + "linear_solver").c_str (), "lis");
  linear_solver_options = g2 ( (alg_options + "solver_options").c_str (), "-p ssor -ssor_omega 0.51 -i cgs -tol 1.e-6 -print 2 -conv_cond 2 -tol_w 0");
  pbc_mode = g2 ( (alg_options + "pbc_mode").c_str (), "strong");

  const std::string out_options = "output/";
  p4estfilename = g2 ( (out_options + "p4estfilename").c_str (), "poisson_boltzmann_p4est");
  markerfilename = g2 ( (out_options + "markerfilename").c_str (), "poisson_boltzmann_marker_0");
  surffilename = g2 ( (out_options + "surffilename").c_str (), "poisson_boltzmann_surface_0");

  // --- membrane ---
  const std::string mem_options = "membrane/";
  membrane_enabled = g2 ( (mem_options + "enabled").c_str (), 0);

  if (membrane_enabled) {
    const std::string mode_str = g2 ( (mem_options + "mode").c_str (), std::string ("ns"));

    if (mode_str == "ns")
      membrane_mode = MEM_MODE_NS;
    else if (mode_str == "implicit")
      membrane_mode = MEM_MODE_IMPLICIT;
    else if (mode_str == "implicit_charged")
      membrane_mode = MEM_MODE_IMPLICIT_CHARGED;
    else {
      if (rank == 0)
        std::cerr << "[ERROR] Unknown membrane mode '" << mode_str
                  << "'. Valid values: ns, implicit, implicit_charged.\n";

      return 1;
    }

    if (membrane_uses_lipids ()) {
      lipid_file = g2 ( (mem_options + "lipid_file").c_str (), std::string ("lipids.pqr"));
      lipid_filetype = g2 ( (mem_options + "lipid_filetype").c_str (), std::string ("pqr"));
    } else {
      // MEM_MODE_IMPLICIT reads no lipid file, so the slab has to be given.
      // Probe presence with two different defaults rather than treating 0 as
      // "unset": 0 is a legitimate membrane_center_z.
      auto is_set = [&g2] (const std::string &key) {
        return g2 (key.c_str (), 1.0) == g2 (key.c_str (), 2.0);
      };

      const std::string required[] = {mem_options + "cell_length_x",
                                      mem_options + "cell_length_y",
                                      mem_options + "membrane_thickness",
                                      mem_options + "membrane_center_z"
                                     };

      for (const std::string &key : required)
        if (! is_set (key)) {
          if (rank == 0)
            std::cerr << "[ERROR] membrane mode = implicit requires '" << key
                      << "': there are no lipids to derive the slab from.\n";

          return 1;
        }

      cell_length_x = g2 ( (mem_options + "cell_length_x").c_str (), 0.0);
      cell_length_y = g2 ( (mem_options + "cell_length_y").c_str (), 0.0);
      membrane_thickness = g2 ( (mem_options + "membrane_thickness").c_str (), 0.0);
      membrane_center_z = g2 ( (mem_options + "membrane_center_z").c_str (), 0.0);

      if (cell_length_x <= 0.0 || cell_length_y <= 0.0 || membrane_thickness <= 0.0) {
        if (rank == 0)
          std::cerr << "[ERROR] cell_length_x, cell_length_y and membrane_thickness "
                    << "must be > 0 (got " << cell_length_x << ", " << cell_length_y
                    << ", " << membrane_thickness << ").\n";

        return 1;
      }
    }
    // The membrane dielectric defaults to the solute one: with e_mem == e_in the
    // implicit modes reproduce the ns model (same low-dielectric region, only a
    // smoother boundary), which makes them comparable.
    e_mem = g2 ( (mem_options + "membrane_dielectric").c_str (), e_in);

    // In ns mode the lipids are merged into the solute and NanoShaper wraps
    // protein and lipids in a single surface, so create_markers cannot tell them
    // apart and assigns e_in to both. Honouring e_mem here needs a second,
    // lipid-only NS surface (see piano_cleanup.md, phase 2.5/M4). Warn rather
    // than silently ignoring the value.
    if (membrane_mode == MEM_MODE_NS && std::fabs (e_mem - e_in) > 1.e-12 && rank == 0)
      std::cerr << "[WARNING] membrane_dielectric = " << e_mem << " is ignored when "
                << "membrane mode = ns: the membrane dielectric is the solute one ("
                << e_in << ").\n          Use mode = implicit or implicit_charged "
                << "to control it.\n";

    // energy() locates the dielectric interface through the NanoShaper surface:
    // border_quad comes from is_in_ns_surf, classifyCube assumes a two-valued
    // epsilon field, and normal_intersection queries the ray cache. In the
    // implicit modes the membrane/solvent interface is a box, so it is in none of
    // those: the flux across it would be missed and the transmembrane border
    // quadrants mis-classified. The potential itself is unaffected (the assembly
    // only reads epsilon_nodes). Refuse rather than print a plausible wrong
    // number -- see piano_cleanup.md, phase 2.5/M5.
    if (membrane_mode != MEM_MODE_NS && calc_energy != 0) {
      if (rank == 0)
        std::cerr << "[ERROR] calc_energy = " << calc_energy << " is not supported with "
                  << "membrane mode = " << mode_str << ".\n"
                  << "        energy() assumes the dielectric interface is the NanoShaper "
                  << "surface, which is not true\n        for an implicit membrane box. "
                  << "Set calc_energy = 0.\n";

      return 1;
    }

    stern_membrane = g2 ( (mem_options + "stern_membrane").c_str (), 0);
    stern_membrane_d = g2 ( (mem_options + "stern_membrane_d").c_str (), 0.0);

    // Transmembrane applied potential, in kT/e (the internal unit of pot_bc,
    // no conversion). Default 0 = no applied potential: the whole DeltaV path
    // (phi_bar source and the z+ boundary value) is then skipped.
    applied_potential = g2 ( (mem_options + "applied_potential").c_str (), 0.0);

    // Membrane mode always uses MESH_SHAPE_MEM (slab mesh). The user cannot
    // override this: the slab geometry and level structure are determined by
    // the lipid/protein atom extents at mesh-build time.
    mesh_shape = MESH_SHAPE_MEM;
    scale_max = g2 ( (mesh_options + "scale_max").c_str (), 2.0);
    nlev_mem = g2 ( (mesh_options + "nlev_mem").c_str (), 2);
    nlev_sol = g2 ( (mesh_options + "nlev_sol").c_str (), 4);
    nlev_prot = g2 ( (mesh_options + "nlev_prot").c_str (), 1);

    // A membrane slab spans the whole xy face by construction, so the natural
    // default is periodic in x and y (same key as above: if the user set it
    // explicitly in [mesh], that value wins -- this only changes the default).
    periodic_x = g2 ( (mesh_options + "periodic_x").c_str (), 1);
    periodic_y = g2 ( (mesh_options + "periodic_y").c_str (), 1);
  }

  return 0;
}

// ====================================
// Autovettore dominante per matrice 3x3 simmetrica
// Usa il metodo di Jacobi (iterativo, robusto)
// ====================================
void compute_dominant_eigenvector (double cov[3][3], double axis[3])
{
  // Inizializza axis = (1,0,0)
  axis[0] = 1.0;
  axis[1] = 0.0;
  axis[2] = 0.0;

  // Potenza iterativa per il vettore principale
  for (int iter = 0; iter < 20; ++iter) {
    double x = cov[0][0]*axis[0] + cov[0][1]*axis[1] + cov[0][2]*axis[2];
    double y = cov[1][0]*axis[0] + cov[1][1]*axis[1] + cov[1][2]*axis[2];
    double z = cov[2][0]*axis[0] + cov[2][1]*axis[1] + cov[2][2]*axis[2];

    double norm = std::sqrt (x*x + y*y + z*z);

    if (norm < 1e-12) break;

    axis[0] = x / norm;
    axis[1] = y / norm;
    axis[2] = z / norm;
  }
}

void align_atoms_to_Z (std::vector<NS::Atom> &atoms)
{
  if (atoms.empty()) return;

  // 1. Centro geometrico (NON pesato)
  double center[3] = {0.0, 0.0, 0.0};

  for (const auto &a : atoms) {
    center[0] += a.pos[0];
    center[1] += a.pos[1];
    center[2] += a.pos[2];
  }

  center[0] /= atoms.size();
  center[1] /= atoms.size();
  center[2] /= atoms.size();

  // Traslazione
  for (auto &a : atoms) {
    a.pos[0] -= center[0];
    a.pos[1] -= center[1];
    a.pos[2] -= center[2];
  }

  // 2. Matrice di covarianza simmetrica
  double cov[3][3] = {{0.0}};

  for (const auto &a : atoms) {
    cov[0][0] += a.pos[0] * a.pos[0];
    cov[0][1] += a.pos[0] * a.pos[1];
    cov[0][2] += a.pos[0] * a.pos[2];
    cov[1][1] += a.pos[1] * a.pos[1];
    cov[1][2] += a.pos[1] * a.pos[2];
    cov[2][2] += a.pos[2] * a.pos[2];
  }

  cov[1][0] = cov[0][1];
  cov[2][0] = cov[0][2];
  cov[2][1] = cov[1][2];

  // 3. Autovettore principale (metodo potenza)
  double axis[3];
  compute_dominant_eigenvector (cov, axis);

  // Normalizza
  double norm = std::sqrt (axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);

  if (norm > 1e-12) {
    axis[0] /= norm;
    axis[1] /= norm;
    axis[2] /= norm;
  }

  // 4. Rotazione: porta axis su (0,0,1) usando rotazione di Rodrigues
  double z_axis[3] = {0.0, 0.0, 1.0};
  double v[3] = {
    axis[1]*z_axis[2] - axis[2]*z_axis[1],
    axis[2]*z_axis[0] - axis[0]*z_axis[2],
    axis[0]*z_axis[1] - axis[1]*z_axis[0]
  };
  double s = std::sqrt (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  double c = axis[0]*z_axis[0] + axis[1]*z_axis[1] + axis[2]*z_axis[2];

  double R[3][3];

  if (s < 1e-8) {
    R[0][0] = R[1][1] = R[2][2] = 1.0;
    R[0][1] = R[0][2] = R[1][0] = R[1][2] = R[2][0] = R[2][1] = 0.0;
  } else {
    double vx = v[0]/s, vy = v[1]/s, vz = v[2]/s;
    double k = 1.0 - c;
    R[0][0] = c + vx*vx*k;
    R[0][1] = vx*vy*k - vz*s;
    R[0][2] = vx*vz*k + vy*s;
    R[1][0] = vy*vx*k + vz*s;
    R[1][1] = c + vy*vy*k;
    R[1][2] = vy*vz*k - vx*s;
    R[2][0] = vz*vx*k - vy*s;
    R[2][1] = vz*vy*k + vx*s;
    R[2][2] = c + vz*vz*k;
  }

  // 5. Applica rotazione
  for (auto &a : atoms) {
    double x = a.pos[0], y = a.pos[1], z = a.pos[2];
    a.pos[0] = R[0][0]*x + R[0][1]*y + R[0][2]*z;
    a.pos[1] = R[1][0]*x + R[1][1]*y + R[1][2]*z;
    a.pos[2] = R[2][0]*x + R[2][1]*y + R[2][2]*z;
  }
}

void
poisson_boltzmann::read_atoms_from_pqr (std::basic_istream<char> &inputfile)
{
  static NS::Atom a;
  atoms.clear ();

  while (inputfile >> a)
    atoms.push_back (a);

  if (aligned == 1) {
    align_atoms_to_Z (atoms);
  }

  if (atoms.size() < 4) {
    auto comp = [] (const NS::Atom &a1, const NS::Atom &a2) -> bool {
      return a1.radius < a2.radius;
    };

    auto max_iter = std::max_element (atoms.begin(), atoms.end(), comp);
    std::array<double, 3> max_pos = {
      max_iter->pos[0],
      max_iter->pos[1],
      max_iter->pos[2]
    };

    const double epsilon = 0.001;

    std::array<std::array<int, 3>, 6> directions = {{
        {{-1, 0, 0}},
        {{+1, 0, 0}},
        {{0, -1, 0}},
        {{0, +1, 0}},
        {{0, 0, -1}},
        {{0, 0, +1}}
      }
    };

    for (int ii = 0; ii < 6; ++ii) {
      std::array<double, 3> new_pos = {
        max_pos[0] + epsilon * directions[ii][0],
        max_pos[1] + epsilon * directions[ii][1],
        max_pos[2] + epsilon * directions[ii][2]
      };
      NS::Atom dummy;
      dummy.pos[0] = new_pos[0];
      dummy.pos[1] = new_pos[1];
      dummy.pos[2] = new_pos[2];
      dummy.charge = 0.0;
      dummy.radius = 0.0;
      atoms.push_back (dummy);
    }
  }
}

void
poisson_boltzmann::read_atoms_from_class ()
{
  static std::array<double,3> pos;

  if (atoms_write == 1) {
    int atom_number = 1;

    for (const NS::Atom& i : atoms) {
      pos[0] = i.pos[0];
      pos[1] = i.pos[1];
      pos[2] = i.pos[2];
      pos_atoms.push_back (pos);
      index_atoms.push_back (atom_number);
      charge_atoms.push_back (i.charge);
      r_atoms.push_back (i.radius);
      atom_number++;
    }
  } else {
    for (const NS::Atom& i : atoms) {
      pos[0] = i.pos[0];
      pos[1] = i.pos[1];
      pos[2] = i.pos[2];
      pos_atoms.push_back (pos);
      charge_atoms.push_back (i.charge);
      r_atoms.push_back (i.radius);
    }
  }
}

void
poisson_boltzmann::broadcast_vectors ()
{
  int size_vec = charge_atoms.size ();
  MPI_Bcast (&size_vec, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Ogni processo alloca il vettore
  pos_atoms.resize (size_vec);
  charge_atoms.resize (size_vec);
  r_atoms.resize (size_vec);
  index_atoms.resize (size_vec);

  // Effettuare il broadcast del vettore
  MPI_Bcast (pos_atoms.data (), size_vec * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (charge_atoms.data (), size_vec, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (r_atoms.data (), size_vec, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (index_atoms.data (), size_vec, MPI_INT, 0, MPI_COMM_WORLD);
}

void
poisson_boltzmann::write_atoms_to_pqr (std::basic_ostream<char> &outputfile)
{
  int Atom_number = 1;

  outputfile << std::setw (10) << std::left << "fieldname" << std::setw (12)
             << std::left <<"Atom_number" << std::setw (12) << std::left << "Atom_name" << std::setw (16) << std::left
             << "Residue_name" << std::setw (16) << std::left << "Residue_number" << std::setw (10) << std::left << "X"
             << std::setw (10) << std::left << "Y" << std::setw (10) << std::left
             << "Z" << std::setw (10) << std::left << "Charge" << std::setw (10) << std::left << "Radius" << std::endl;

  for (auto & ii : atoms)
    outputfile << std::setw (10) << std::left << "ATOM" << std::setw (12)
               << std::left << Atom_number++ << std::setw (12) << std::left << ii.ai.name << std::setw (16) << std::left
               << ii.ai.resName << std::setw (16) << std::left << ii.ai.resNum << std::setw (10) << std::left << ii.pos[0]
               << std::setw (10) << std::left << ii.pos[1] << std::setw (10) << std::left
               << ii.pos[2] << std::setw (10) << std::left << ii.charge << std::setw (10) << std::left << ii.radius << std::endl;

}



std::basic_istream<char>&
operator>> (std::basic_istream<char>& inputfile, NS::Atom &a)
{
  int Atom_number;
  std::string Field_name;

  inputfile >> Field_name
            >> Atom_number
            >> a.ai.name
            >> a.ai.resName;

  std::string token;
  inputfile >> token;

  // Verifica se è un numero (resNum) oppure una stringa (chain)
  bool is_number = !token.empty() &&
                   (std::isdigit (token[0]) || token[0] == '-' || token[0] == '+');

  if (is_number) {
    // Era resNum - metti solo il token indietro nello stream
    for (auto it = token.rbegin(); it != token.rend(); ++it) {
      inputfile.putback (*it);
    }

    a.ai.chain.clear(); // Nessuna catena specificata
  } else {
    // Era chain - leggi il resNum successivo
    a.ai.chain = token;
  }

  // Ora leggi i dati numerici
  inputfile >> a.ai.resNum >> a.pos[0] >> a.pos[1] >> a.pos[2]
            >> a.charge >> a.radius;

  if (a.radius < 1.e-5)
    a.radius = 1.0;

  return inputfile;
}

std::basic_istream<char>&
operator>> (std::basic_istream<char>& inputfile, std::array<float,5> &a)
{
  int Atom_number;
  std::string Field_name;
  std::string name;
  std::string resName;
  int resNum;

  inputfile >> Field_name
            >> Atom_number
            >> name
            >> resName
            >> resNum
            >> a[0] // x_pos
            >> a[1] // y_pos
            >> a[2] // z_pos
            >> a[3] // charge
            >> a[4]; // radius

  if (a[4] < 1.e-5)
    a[4] = 1.0;

  return inputfile;
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

          // else if ((x1 >= this->l_c[0]) && (x1 <= this->r_c[0])
          // && (y1 >= this->l_c[1]) && (y1 <= this->r_c[1])
          // && (z1 >= this->l_c[2]) && (z1 <= this->r_c[2])
          // &&currentlevel < this->scale_level_min_box) {
          // retval = 1;
          // break;
          // }
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

bool
poisson_boltzmann::is_in (const NS::Atom& i,
                          tmesh_3d::quadrant_iterator q)
{
  double tol = p4esttol * (rr[0]-ll[0]);

  if (mesh_shape == 2)
    tol = p4esttol * (r_c[0]-l_c[0]);

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
poisson_boltzmann::is_in_ref (const NS::Atom& i,
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
poisson_boltzmann::refine_surface (ray_cache_t & ray_cache)
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

                    std::vector<int> direzioni = {0,1,2};
                    direzioni.erase (direzioni.begin ()+idir);

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

                  std::vector<int> direzioni = {0,1,2};
                  direzioni.erase (direzioni.begin ()+idir);

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

      auto refinement = [&rcoeff,this]
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

                    std::vector<int> direzioni = {0,1,2};
                    direzioni.erase (direzioni.begin ()+idir);

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

                  std::vector<int> direzioni = {0,1,2};
                  direzioni.erase (direzioni.begin ()+idir);

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

      auto coarsening = [&rcoeff,this]
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
poisson_boltzmann::is_in_ns_surf_stern (ray_cache_t & ray_cache, double x, double y, double z, int dir)
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

  crossings_t & ct = ray_cache (x1, x2, dir);

  if (!ct.init && rank != 0) {
    std::array<double, 2> ray = {x1, x2};
    ray_cache.rays[dir].erase (ray);
    return -1.;
  }

  int i = 0;
  int sign = 1;

  if (ct.inters.size () == 0 || x3 < (ct.inters[i]- stern_layer))
    return 0; //if there are no inters or y_the coord is before the first intersection, the point is outside.

  while (i < ct.inters.size () && x3 > (ct.inters[i] - stern_layer*sign)) {
    //go on until the inters is passed
    i++;
    sign *= -1;
  }

  return (i % 2);
}


void
poisson_boltzmann::create_markers (ray_cache_t & ray_cache)
{

  int size, rank;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);

  double eps_in = 4.0*pi*e_0*e_in*kb*T*Angs/ (e*e); //adim e_in
  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/ (e*e); //adim e_out
  double eps_mem = 4.0*pi*e_0*e_mem*kb*T*Angs/ (e*e); //adim e_mem

  // Implicit membrane: the slab is a box, not a molecular surface. It spans the
  // whole xy face, so only the z bounds filter. Nodes outside the molecular
  // surface but inside the slab get the membrane dielectric and no ionic
  // screening (no ions in the hydrophobic core).
  const bool implicit_membrane = membrane_enabled && membrane_mode != MEM_MODE_NS;

  /////////////////////////////////////////////////////////
  //reactions
  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/ (e_0*e_out*kb*T);

  this->marker.assign (this->tmsh.num_local_quadrants (), 0.0); //marker = 0 -> in

  if (stern_layer_surf == 1) {
    this->marker_k.assign (this->tmsh.num_local_quadrants (), 1.0); //marker = 1 -> out stern
    this->reaction.assign (tmsh.num_local_quadrants (), eps_out*k2);
  }

  epsilon_nodes = std::make_unique<distributed_vector> (tmsh.num_owned_nodes (),mpicomm);
  epsilon_nodes->get_owned_data ().assign (tmsh.num_owned_nodes (), eps_out);

  reaction_nodes = std::make_unique<distributed_vector> (tmsh.num_owned_nodes (),mpicomm);
  reaction_nodes->get_owned_data ().assign (tmsh.num_owned_nodes (), eps_out*k2);

  // Applied potential: regional bulk reference phi_bar, 0 in the lower bath and
  // applied_potential in the upper one, split at the lipid-slab midplane. Inside
  // the membrane and the solute c = 0, so the value phi_bar takes there never
  // reaches the system. z_mem_bot/z_mem_top are the physical (van der Waals)
  // slab, not the l_mem/r_mem refinement boxes; the midplane is the same either
  // way since their padding is symmetric. Consumed as the RHS source +c*phi_bar
  // in assemple_system_matrix. Left null when there is no applied potential, and
  // every downstream use tests for that, so V = 0 is a strict no-op.
  const bool use_phi_bar =
    membrane_enabled && std::fabs (applied_potential) > 1.e-12;
  const double z_mid = 0.5 * (z_mem_bot + z_mem_top);

  if (use_phi_bar) {
    phi_bar = std::make_unique<distributed_vector> (tmsh.num_owned_nodes (), mpicomm);
    phi_bar->get_owned_data ().assign (tmsh.num_owned_nodes (), 0.0);

    // Worth printing: the slab bounds are echoed only in the implicit modes, so
    // in ns mode this is the one place the midplane splitting the two baths --
    // which is where the whole applied-potential setup hinges -- is visible.
    if (rank == 0)
      std::cout << "  Applied potential: phi_bar = 0 below z = " << z_mid
                << " Å, " << applied_potential << " kT/e above\n";
  }

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
      double x,y,z;

      for (int ii = 0; ii < 8; ++ii) {
        local_num =quadrant->gt (ii);
        x = quadrant->p (0, ii);
        y = quadrant->p (1, ii);
        z = quadrant->p (2, ii);

        if (! quadrant->is_hanging (ii)) {
          if (use_phi_bar)
            (*phi_bar)[local_num] = (z > z_mid) ? applied_potential : 0.0;

          if (this->is_in_ns_surf (ray_cache, x, y, z, 2) > 0.5) { //inside the molecule
            ++num_int_nodes;
            ++num_int_nodes_stern;
            (*epsilon_nodes)[local_num] = eps_in;
            (*reaction_nodes)[local_num] = 0.0;
          } else if (this->is_in_ns_surf (ray_cache, x, y, z, 2) < -0.5) {
            ray_cache.num_req_rays[2]++;

            // std::array<double, 2> ray {x,y};

            ray_cache.rays_list[2].insert ({x,y});
          } else if (implicit_membrane && z >= z_mem_bot && z <= z_mem_top) {
            //outside the molecule, inside the membrane slab
            (*epsilon_nodes)[local_num] = eps_mem;
            (*reaction_nodes)[local_num] = 0.0;
          }

          if (this->is_in_ns_surf (ray_cache, x, y, z, 0) < -0.5) {
            ray_cache.num_req_rays[0]++;

            // std::array<double, 2> ray {y,z};

            ray_cache.rays_list[0].insert ({y,z});
          }

          if (this->is_in_ns_surf (ray_cache, x, y, z, 1) < -0.5) {
            ray_cache.num_req_rays[1]++;

            std::array<double, 2> ray {x,z};

            ray_cache.rays_list[1].insert ({x,z});
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

              std::vector<int> direzioni {0,1,2};
              direzioni.erase (direzioni.begin ()+idir);

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
          this->marker[quadrant->get_forest_quad_idx ()] = 1.0/2.0; //"border"
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

  if (size >1) {
    bim3a_solution_with_ghosts (tmsh, *epsilon_nodes, replace_op);

    if (stern_layer_surf == 0)
      bim3a_solution_with_ghosts (tmsh, (*reaction_nodes), replace_op);

    if (use_phi_bar)
      bim3a_solution_with_ghosts (tmsh, *phi_bar, replace_op);
  }
}

void
poisson_boltzmann::create_density_map (ray_cache_t & ray_cache)
{
  int size, rank;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);

  // ------------------------------------------------------------
  // Early exit: if rho_fixed is already initialized, skip the whole procedure.
  // The density map must be created only once.
  // ------------------------------------------------------------
  if (this->rho_fixed) {
    if (rank == 0)
      std::cout << "[INFO] Density map already initialized. Skipping.\n";

    return;
  }

  // ------------------------------------------------------------
  // Allocate and initialize nodal density vector (rho)
  // ------------------------------------------------------------
  this->rho_fixed = std::make_unique<distributed_vector> (tmsh.num_owned_nodes(), mpicomm);
  this->rho_fixed->get_owned_data().assign (tmsh.num_owned_nodes(), 0.0);

  // ------------------------------------------------------------
  // Vector of ones at nodes
  // ------------------------------------------------------------
  this->ones = std::make_unique<distributed_vector> (tmsh.num_owned_nodes(), mpicomm);
  this->ones->get_owned_data().assign (tmsh.num_owned_nodes(), 1.0);

  // ------------------------------------------------------------
  // Per-cell constant field used for computing patch volumes
  // ------------------------------------------------------------
  this->const_ones.assign (tmsh.num_local_quadrants(), 1.0);

  // ------------------------------------------------------------
  // vol_patch will store the nodal patch volumes obtained via BIM
  // ------------------------------------------------------------
  std::unique_ptr<distributed_vector> vol_patch =
    std::make_unique<distributed_vector> (tmsh.num_owned_nodes(), mpicomm);

  // ------------------------------------------------------------
  // Update ghost values if running in parallel
  // ------------------------------------------------------------
  if (size > 1)
    bim3a_solution_with_ghosts (tmsh, *ones, replace_op);

  // ------------------------------------------------------------
  // Compute patch volumes by integrating the constant field = 1
  // ------------------------------------------------------------
  bim3a_rhs (tmsh, const_ones, *ones, *vol_patch);

  if (size > 1)
    vol_patch->assemble();

  // ------------------------------------------------------------
  // Locate the atom positions inside the adaptive mesh
  // ------------------------------------------------------------
  search_points();

  // ------------------------------------------------------------
  // Distribute the atomic charges to the surrounding nodes
  // using a linear volumetric approximation.
  // ------------------------------------------------------------
  for (auto it = lookup_table.begin(); it != lookup_table.end(); ++it) {

    // Compute cell volume (Cartesian and axis-aligned)
    double volume =
      (it->second.p (0, 7) - it->second.p (0, 0)) *
      (it->second.p (1, 7) - it->second.p (1, 0)) *
      (it->second.p (2, 7) - it->second.p (2, 0));

    for (int ii = 0; ii < 8; ++ii) {

      // Linear interpolation weight based on opposite corner distances
      double weight = std::abs (
                        (pos_atoms[it->first][0] - it->second.p (0, 7 - ii)) *
                        (pos_atoms[it->first][1] - it->second.p (1, 7 - ii)) *
                        (pos_atoms[it->first][2] - it->second.p (2, 7 - ii))) / volume;

      // Regular node (not hanging)
      if (!it->second.is_hanging (ii)) {
        (*rho_fixed)[it->second.gt (ii)] +=
          charge_atoms[it->first] * 4.0 * pi * weight /
          (*vol_patch)[it->second.gt (ii)];
      }
      // Hanging node → distribute to parents
      else {
        for (int jj = 0; jj < it->second.num_parents (ii); ++jj) {
          double denom =
            it->second.num_parents (ii) *
            (*vol_patch)[it->second.gparent (jj, ii)];
          (*rho_fixed)[it->second.gparent (jj, ii)] +=
            charge_atoms[it->first] * 4.0 * pi * weight / denom;
        }
      }
    }
  }

  // Release temporary vector
  vol_patch.reset();

  // ------------------------------------------------------------
  // Sync ghost nodes of rho_fixed in parallel runs
  // ------------------------------------------------------------
  if (size > 1)
    bim3a_solution_with_ghosts (tmsh, *rho_fixed);
}

void
poisson_boltzmann::assemple_system_matrix (ray_cache_t & ray_cache)
{
  int size, rank;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);

  // ------------------------------------------------------------
  // 1) CHECK / BUILD REQUIRED DATA STRUCTURES
  // ------------------------------------------------------------

  // If the density map has not been created yet,
  // the system matrix cannot be assembled. Build it first.
  if (!rho_fixed) {
    create_density_map (ray_cache);
  }

  // Function computing fractional volume intersections (for cut-cells)
  auto func_frac = [&] (tmesh_3d::quadrant_iterator& quadrant) {
    return cube_fraction_intersection (quadrant, ray_cache);
  };

  // PBC: diagnostic conformity check + right->left node map, before assembling
  // the matrix. No-op if neither periodic_x nor periodic_y.
  check_pbc_face_conformity ();
  build_pbc_node_map ();

  // Without a periodic pair this stays false and every bimpp call below gets
  // default_ord, i.e. the exact call of the non-periodic code path.
  // Mortar mode also uses default_ord: the periodic identification there is
  // weak (Lagrange multipliers assembled by assemble_mortar_system), not the
  // strong node folding driven by the ordering.
  pbc_use_ordering = (periodic_x || periodic_y) && (pbc_mode != "mortar");

  if (pbc_use_ordering)
    pbc_build_canonical (pbc_canonical);

  // Maps every periodic right-face node onto its left partner; identity on
  // every other node, the Dirichlet z± nodes included (build_pbc_node_map
  // deliberately keeps those out of the maps, see bc.md § "Spigoli
  // orizzontali"). Never called unless pbc_use_ordering.
  ordering periodic_ord = [this] (p4est_gloidx_t g) -> size_t {
    return static_cast<size_t> (pbc_canonical[g]);
  };
  const ordering &ord = pbc_use_ordering ? periodic_ord : default_ord;

  // Allocate sparse matrix A and RHS vector
  A = std::make_unique<distributed_sparse_matrix> (mpicomm);
  A->set_ranges (tmsh.num_owned_nodes());

  rhs = std::make_unique<distributed_vector> (tmsh.num_owned_nodes(), mpicomm);

  // ------------------------------------------------------------
  // 2) ASSEMBLE THE LINEAR SYSTEM (RHS + STIFFNESS MATRIX)
  // ------------------------------------------------------------

  // Assemble RHS using previously computed fixed charge density
  bim3a_rhs (tmsh, const_ones, *rho_fixed, *rhs, ord);

  // rho_fixed and const_ones are no longer needed after building RHS
  rho_fixed.reset();
  std::vector<double>().swap (const_ones);

  // Assemble Laplace operator (with fractional cell treatment)
  bim3a_laplacian_frac (tmsh, *epsilon_nodes, *A, func_frac, ord, ord);

  // Add reaction term (Stern layer present or fractional formulation)
  if (stern_layer_surf == 1) {
    // Standard Stern-layer reaction
    bim3a_reaction (tmsh, reaction, *ones, *A, ord, ord);
  } else {
    // Fractional reaction for intersected cells
    bim3a_reaction_frac (tmsh, (*reaction_nodes), *ones, *A, func_frac, ord, ord);
  }

  // Applied potential: RHS source +c*phi_bar. The membrane reaction term is
  // c (phi - phi_bar); c*phi is the matrix block assembled just above, while
  // c*phi_bar is a known source and belongs on the RHS. Built with the same
  // operator, the same fractional weights and the same ordering as that block,
  // so it is exactly M_c*phi_bar (the box-method reaction matrix is lumped).
  // Using ord is not cosmetic: with strong PBC it folds the right-face rows
  // onto their left partners, and default_ord here would drop the term on rows
  // the rest of the system no longer uses.
  //
  // Accumulated straight into *rhs rather than into a temporary: bim3a_rhs*
  // adds (+=) at *global* indices, so a rank also writes into nodes owned by
  // its neighbours, and those land in the vector's non-local map. Only
  // assemble() ships them to the owner, and get_owned_data() does not see them
  // -- summing a temporary's owned data into *rhs would silently drop every
  // cross-rank contribution (invisible on one rank, where nothing is
  // non-local). Writing into *rhs lets the single rhs->assemble() below reduce
  // this term together with the rho_fixed one assembled the same way.
  if (phi_bar) {
    if (stern_layer_surf == 1)
      bim3a_rhs (tmsh, reaction, *phi_bar, *rhs, ord);
    else
      bim3a_rhs_frac (tmsh, *reaction_nodes, *phi_bar, *rhs, func_frac, ord);

    phi_bar.reset ();
  }

  // Reaction-related vectors are no longer required
  reaction_nodes.reset();
  ones.reset();
  std::vector<double>().swap (reaction);
  std::vector<double>().swap (marker);

  // ------------------------------------------------------------
  // 2b) MORTAR ASSEMBLY (weak PBC, before Dirichlet)
  // ------------------------------------------------------------
  // bc.md § "Il lavoro da fare in NGPB": the mortar blocks are written into the
  // free (pre-Dirichlet) system, so the Dirichlet pass below zeroes the C-block
  // rows of any z± node just like it does for A_free.
  if (pbc_mode == "mortar" && (periodic_x || periodic_y))
    assemble_mortar_system ();

  // ------------------------------------------------------------
  // 3) APPLY BOUNDARY CONDITIONS
  // ------------------------------------------------------------
  // Build Dirichlet boundary list
  dirichlet_bcs3 bcs;

  if (bc == 1) { // Dirichlet BC, homogeneous unless an applied potential is set
    if (std::fabs (pot_bc) > 1.e-5 && rank == 0)
      std::cerr << "[WARNING] Boundary conditions may be inaccurate!!\n";

    // Applied potential: the two z faces become the electrodes, z- (face 4)
    // grounded and z+ (face 5) held at applied_potential. This is orthogonal to
    // pbc_mode: the z faces are Dirichlet under strong and mortar alike, only
    // x/y are periodic, so the face filter below is untouched. Without an
    // applied potential every face keeps the homogeneous value.
    const bool membrane_voltage =
      membrane_enabled && std::fabs (applied_potential) > 1.e-12;

    if (membrane_voltage && rank == 0)
      std::cout << "  Membrane voltage: z- = 0, z+ = " << applied_potential
                << " kT/e\n";

    for (auto const& ibc : bcells) {
      if (periodic_x && (ibc.second == 0 || ibc.second == 1)) continue;
      if (periodic_y && (ibc.second == 2 || ibc.second == 3)) continue;

      const double bc_val =
        (membrane_voltage && ibc.second == 5) ? applied_potential : 0.0;

      bcs.emplace_back (ibc.first, ibc.second,
      [bc_val] (double, double, double) {
        return bc_val;
      });
    }

    bim3a_dirichlet_bc (tmsh, bcs, *A, *rhs, ord);
  }

  if (bc == 2) { // Coulombic Dirichlet BC
    for (auto const& ibc : bcells) {
      if (periodic_x && (ibc.second == 0 || ibc.second == 1)) continue;
      if (periodic_y && (ibc.second == 2 || ibc.second == 3)) continue;
      bcs.emplace_back (ibc.first, ibc.second,
      [&] (double x, double y, double z) {
        return coulomb_boundary_conditions (x, y, z);
      });
    }

    bim3a_dirichlet_bc (tmsh, bcs, *A, *rhs, ord);
  }

  if (bc == 3) { // Analytic Dirichlet BC (sphere test case)
    for (auto const& ibc : bcells) {
      if (periodic_x && (ibc.second == 0 || ibc.second == 1)) continue;
      if (periodic_y && (ibc.second == 2 || ibc.second == 3)) continue;
      bcs.emplace_back (ibc.first, ibc.second,
      [&] (double x, double y, double z) {
        return analytic_solution (x, y, z);
      });
    }

    bim3a_dirichlet_bc (tmsh, bcs, *A, *rhs, ord);
  }

  // ------------------------------------------------------------
  // 3b) CLOSE THE ROWS OF THE FOLDED NODES
  // ------------------------------------------------------------
  // The periodic ordering folded every right-face node r onto its left partner
  // canonical[r], on rows *and* columns: row r and column r are both untouched,
  // so phi_r would be undetermined and A singular. Close those rows with a
  // plain identity, phi_r = 0, and recover the real value with
  // pbc_apply_periodic_copy() after the solve.
  //
  // NOT with the constraint row phi_r - phi_l = 0 (a -1 in column l): that
  // extra off-diagonal entry has no symmetric counterpart, since column r
  // stays empty. It turns A into a block-triangular, non-symmetric matrix and
  // CG -- which is what LIS runs by default, legitimately, because the system
  // is otherwise symmetric (measured: ||A - A^T||_max = 4.4e-16) -- then has no
  // convergence guarantee and silently returns garbage on phi_r. MUMPS
  // factorizes and does not care, which is why the two solvers disagreed.
  // Keeping the identity leaves A symmetric and the solve unaffected.
  if (pbc_use_ordering) {
    const int row_start = A->range_start ();
    const int row_end   = A->range_end ();
    for (int r = row_start; r < row_end; ++r) {
      if (pbc_canonical[r] == r)
        continue;
      (*A)[r][r] = 1.0;
      (*rhs)[r]  = 0.0;
    }
  }

  // ------------------------------------------------------------
  // 4) PARALLEL ASSEMBLY (MPI)
  // ------------------------------------------------------------

  if (size > 1) {
    A->assemble();
    rhs->assemble();
  }
}


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
poisson_boltzmann::mumps_compute_electric_potential (ray_cache_t & ray_cache)
{
  int rank, size;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);

  mumps mumps_solver;

  std::vector<double> vals;
  std::vector<int> irow, jcol;

  // Owned rows only (flag=false, the default): the matrix stays distributed
  // across ranks, with or without PBC -- the periodic identification is
  // already baked into it by the assembly ordering.
  (*A).aij (vals, irow, jcol, mumps_solver.get_index_base ());

  // ------------------------------------------------------------------------
  // Mortar (weak PBC): solve the augmented saddle-point system.
  //
  // The physical rows already carry the C block (mortar columns >= N_global)
  // from assemble_mortar_system. Here rank 0 appends the C^T rows from the
  // dense mortar_C accumulator plus an identity for the z-boundary multipliers,
  // gathers the RHS and pads it with zeros for the mortar DOFs, and MUMPS
  // solves the whole N_phys + n_pairs*ndofm system with a centralized RHS
  // (solution written back in place on rank 0, then scattered to phi).
  // MUMPS only: the augmented system is indefinite.
  // ------------------------------------------------------------------------
  if (pbc_mode == "mortar" && (periodic_x || periodic_y)) {
    const int    base     = mumps_solver.get_index_base ();
    const int    nm       = (1 << minlevel);
    const int    n_pairs  = (periodic_x ? 1 : 0) + (periodic_y ? 1 : 0);
    const size_t N_global = static_cast<size_t> (tmsh.num_global_nodes ());
    const int    n_system = static_cast<int> (N_global)
                          + n_pairs * ndofm;

    // Zone D -- rank 0 appends the mortar rows (C^T) and z-boundary identity.
    if (rank == 0) {
      const size_t n_mortar_rows = static_cast<size_t> (n_pairs) * ndofm;
      for (size_t drow = 0; drow < n_mortar_rows; ++drow) {
        const int global_row = static_cast<int> (N_global + drow) + base;
        for (size_t c = 0; c < N_global; ++c) {
          const double v = mortar_C[drow * N_global + c];
          if (v != 0.0) {
            irow.push_back (global_row);
            jcol.push_back (static_cast<int> (c) + base);
            vals.push_back (v);
          }
        }
        // z-boundary multipliers (j1 == 0 or nm) were zeroed in assembly: pin
        // them with an identity so the augmented matrix is non-singular.
        const int j1 = static_cast<int> (drow % static_cast<size_t> (ndofm)) / (nm + 1);
        if (j1 == 0 || j1 == nm) {
          irow.push_back (global_row);
          jcol.push_back (global_row);
          vals.push_back (1.0);
        }
      }
    }
    std::vector<double>().swap (mortar_C);

    // Gather the physical RHS onto rank 0, padded to n_system with zeros for
    // the mortar DOFs. rhs_counts/displs are reused for the solution scatter.
    std::vector<double> rhs_aug;
    std::vector<int>    rhs_counts (size, 0), rhs_displs (size, 0);
    const int local_n = static_cast<int> (tmsh.num_owned_nodes ());
    MPI_Gather (&local_n, 1, MPI_INT, rhs_counts.data (), 1, MPI_INT, 0, mpicomm);
    if (rank == 0) {
      for (int r = 1; r < size; ++r)
        rhs_displs[r] = rhs_displs[r - 1] + rhs_counts[r - 1];
      rhs_aug.assign (static_cast<size_t> (n_system), 0.0);
    }
    MPI_Gatherv (rhs->get_owned_data ().data (), local_n, MPI_DOUBLE,
                 rhs_aug.data (), rhs_counts.data (), rhs_displs.data (),
                 MPI_DOUBLE, 0, mpicomm);

    mumps_solver.set_lhs_distributed ();
    mumps_solver.set_distributed_lhs_structure (n_system, irow, jcol);
    mumps_solver.set_distributed_lhs_data (vals);
    if (rank == 0)
      mumps_solver.set_rhs (rhs_aug);   // centralized rhs, solution in place

    rhs.reset ();
    A.reset ();

    if (rank == 0)
      std::cout << "  [Mortar] augmented system size for MUMPS: " << n_system
                << "  (N_phys=" << N_global << " + " << n_pairs << "x" << ndofm << ")\n";

    std::cout << "mumps_solver.analyze () = "   << mumps_solver.analyze ()   << std::endl;
    std::cout << "mumps_solver.factorize () = " << mumps_solver.factorize () << std::endl;
    std::cout << "mumps_solver.solve () = "     << mumps_solver.solve ()     << std::endl;

    // Solution: rhs_aug[0 .. N_phys) on rank 0 is the potential u; scatter it
    // back to the owned ranges (mortar multipliers rhs_aug[N_phys ..) dropped).
    phi = std::make_unique<distributed_vector> (tmsh.num_owned_nodes ());
    MPI_Scatterv (rhs_aug.data (), rhs_counts.data (), rhs_displs.data (), MPI_DOUBLE,
                  phi->get_owned_data ().data (), local_n, MPI_DOUBLE, 0, mpicomm);

    if (size > 1)
      bim3a_solution_with_ghosts (tmsh, *phi, replace_op);

    mumps_solver.cleanup ();
    return;
  }

  const int n_system = static_cast<int> (tmsh.num_global_nodes ());

  mumps_solver.set_lhs_distributed ();
  mumps_solver.set_distributed_lhs_structure (n_system, irow, jcol);
  mumps_solver.set_distributed_lhs_data (vals);

  mumps_solver.set_rhs_distributed (*rhs);
  rhs.reset ();
  A.reset ();

  std::cout << "mumps_solver.analyze () = "
            << mumps_solver.analyze ()
            << std::endl;
  std::cout << "mumps_solver.factorize () = "
            << mumps_solver.factorize ()
            << std::endl;
  std::cout << "mumps_solver.solve () = "
            << mumps_solver.solve ()
            << std::endl;

  phi = std::make_unique<distributed_vector> (tmsh.num_owned_nodes ());
  (*phi) = mumps_solver.get_distributed_solution ();

  // The folded rows were solved as phi_r = 0; fill in their real value.
  pbc_apply_periodic_copy ();

  if (size > 1)
    bim3a_solution_with_ghosts (tmsh, *phi, replace_op);

  mumps_solver.cleanup ();
}



void
poisson_boltzmann::lis_compute_electric_potential (ray_cache_t & ray_cache)
{
  int rank, size;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);

  LIS_INT is, ie;

  std::vector<double> vals;
  std::vector<int> irow, jcol;

  (*A).csr (vals, jcol, irow);

  LIS_INT ln = static_cast<LIS_INT> (rhs->get_owned_data ().size ());
  LIS_VECTOR rhs_lis;
  lis_vector_create (mpicomm, &rhs_lis);
  lis_vector_set_size (rhs_lis, ln, 0);
  lis_vector_get_range (rhs_lis, &is, &ie);

  for (LIS_INT k = is; k < ie; ++k)
    lis_vector_set_value (LIS_INS_VALUE, k, rhs->get_owned_data ()[k - is], rhs_lis);

  rhs.reset ();

  LIS_VECTOR phi_lis;
  lis_vector_create (mpicomm, &phi_lis);
  lis_vector_set_size (phi_lis, ln, 0);
  lis_vector_get_range (phi_lis, &is, &ie);

  LIS_INT nnz = (*A).owned_nnz ();
  LIS_INT n   = static_cast<LIS_INT> (tmsh.num_owned_nodes ());

  A.reset ();

  LIS_MATRIX A_lis;
  lis_matrix_create (mpicomm, &A_lis);
  lis_matrix_set_size (A_lis, n, 0);
  lis_matrix_set_csr (nnz, &irow[0], &jcol[0], &vals[0], A_lis);
  lis_matrix_assemble (A_lis);

  LIS_SOLVER solver;
  lis_solver_create (&solver);
  std::string opts = linear_solver_options;
  lis_solver_set_option (&opts[0], solver);

  LIS_INT solve_err = lis_solve (A_lis, rhs_lis, phi_lis, solver);

  // Report what the solver actually did: the return code alone is not enough,
  // since a run that stops on maxiter still returns a plausible-looking phi.
  {
    LIS_INT  status = -1, iter = -1, nsol = -1, precon = -1;
    LIS_REAL resnorm = -1.0, rhs_nrm1 = -1.0;
    char     solvername[128] = "?", preconname[128] = "?";

    lis_solver_get_status (solver, &status);
    lis_solver_get_iter (solver, &iter);
    lis_solver_get_residualnorm (solver, &resnorm);

    // ||b||_1: with -conv_cond 2 (LIS_CONV_COND_NRM1_B) the stopping test is
    // ||r||_1 <= tol * ||b||_1, so the residual norm above is only meaningful
    // against this scale.
    lis_vector_nrm1 (rhs_lis, &rhs_nrm1);

    if (lis_solver_get_solver (solver, &nsol) == LIS_SUCCESS)
      lis_solver_get_solvername (nsol, solvername);

    if (lis_solver_get_precon (solver, &precon) == LIS_SUCCESS)
      lis_solver_get_preconname (precon, preconname);

    if (rank == 0) {
      std::cout << "\n========== [ LIS Solver Report ] ==========\n";
      std::cout << "  Options string     : \"" << linear_solver_options << "\"\n";
      std::cout << "  Method             : " << solvername << '\n';
      std::cout << "  Preconditioner     : " << preconname << '\n';
      std::cout << "  Iterations         : " << iter << '\n';
      std::cout << "  Residual norm      : " << std::scientific << resnorm << '\n';
      std::cout << "  ||b||_1            : " << std::scientific << rhs_nrm1 << '\n';
      std::cout << "  Status             : " << status << '\n';
      std::cout << "===========================================\n\n";

      if (solve_err != LIS_SUCCESS || status != LIS_SUCCESS)
        std::cerr << "WARNING: LIS did not converge (return code " << solve_err
                  << ", status " << status << ", " << iter << " iterations, residual "
                  << resnorm << "). The potential below is NOT trustworthy.\n";
    }
  }

  lis_solver_destroy (solver);
  lis_vector_destroy (rhs_lis);

  phi = std::make_unique<distributed_vector> (tmsh.num_owned_nodes (), mpicomm);
  lis_vector_get_values (phi_lis, is, ln, phi->get_owned_data ().data ());
  lis_vector_destroy (phi_lis);

  // The folded rows were solved as phi_r = 0; fill in their real value.
  pbc_apply_periodic_copy ();

  // The two faces of a periodic pair must carry the exact same solution --
  // that is the definition of strong enforcement.
  if (pbc_use_ordering) {
    static const char* fname[6] = {"x-", "x+", "y-", "y+", "z-", "z+"};
    const int phi_is = static_cast<int> (phi->get_range_start ());
    const int phi_ie = phi_is + static_cast<int> (phi->get_owned_data ().size ());
    if (rank == 0)
      std::cout << "  [PBC ordering] face phi ranges:\n";
    for (int face : {0, 1, 2, 3}) {
      if (face < 2 && !periodic_x) continue;
      if (face >= 2 && !periodic_y) continue;
      double lmin =  std::numeric_limits<double>::max ();
      double lmax = -std::numeric_limits<double>::max ();
      for (auto q = tmsh.begin_quadrant_sweep (); q != tmsh.end_quadrant_sweep (); ++q)
        for (auto local_idx : q->ef (face)) {
          int g = q->gt (local_idx);
          if (g >= phi_is && g < phi_ie) {
            double v = phi->get_owned_data ()[static_cast<size_t> (g - phi_is)];
            lmin = std::min (lmin, v);
            lmax = std::max (lmax, v);
          }
        }
      double gmin, gmax;
      MPI_Reduce (&lmin, &gmin, 1, MPI_DOUBLE, MPI_MIN, 0, mpicomm);
      MPI_Reduce (&lmax, &gmax, 1, MPI_DOUBLE, MPI_MAX, 0, mpicomm);
      if (rank == 0)
        std::cout << "    face " << fname[face] << "  phi=[" << gmin << ", " << gmax << "]\n";
    }
  }

  if (size > 1)
    bim3a_solution_with_ghosts (tmsh, *phi, replace_op);
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


/////////////////////////////////////////////////////////////////////////////////////////////////////

std::array<double,12>
poisson_boltzmann::cube_fraction_intersection (tmesh_3d::quadrant_iterator& quadrant,
    const ray_cache_t & ray_cache)


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
  std::array<double,12> fraction = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};

  int dir;
  int i1, i2;
  double x1, x2;
  std::array<double,2> ray;

  if (marker[quadrant->get_forest_quad_idx ()] == 0.5) {
    for (int j: {
           0,2,4,6
         })
      //x-axis edge
    {
      dir = 0;
      i1 = edge2nodes[2 * j ];
      i2 = edge2nodes[2 * j + 1];
      x1 = quadrant->p (dir, i1);
      x2 = quadrant->p (dir, i2);
      std::vector<int> direzioni {0,1,2};
      direzioni.erase (direzioni.begin ()+dir);

      for (unsigned i = 0; i < direzioni.size (); ++i) {
        ray[i] = quadrant->p (direzioni[i], i1);
      }

      auto it0 = ray_cache.rays[dir].find (ray);
      if (it0 == ray_cache.rays[dir].end ()) continue;
      auto inters = it0->second.inters;

      for (int ii =0; ii<inters.size (); ii++) {
        if (inters[ii]>= x1 && inters[ii] <=x2) {
          fraction[j] = (inters[ii] - x1)/ (x2 - x1);
        }
      }
    }

    for (int j: {
           1,3,5,7
         })
      //y-axis edge
    {
      dir = 1;
      i1 = edge2nodes[2 * j ];
      i2 = edge2nodes[2 * j + 1];
      x1 = quadrant->p (dir, i1);
      x2 = quadrant->p (dir, i2);
      std::vector<int> direzioni {0,1,2};
      direzioni.erase (direzioni.begin ()+dir);

      for (unsigned i = 0; i < direzioni.size (); ++i) {
        ray[i] = quadrant->p (direzioni[i], i1);
      }

      auto it0 = ray_cache.rays[dir].find (ray);
      if (it0 == ray_cache.rays[dir].end ()) continue;
      auto inters = it0->second.inters;

      for (int ii =0; ii<inters.size (); ii++) {
        if (inters[ii]>= x1 && inters[ii] <=x2) {
          fraction[j] = (inters[ii] - x1)/ (x2 - x1);
        }
      }
    }

    for (int j: {
           8,9,10,11
         })
      //z-axis edge
    {
      dir = 2;
      i1 = edge2nodes[2 * j ];
      i2 = edge2nodes[2 * j + 1];
      x1 = quadrant->p (dir, i1);
      x2 = quadrant->p (dir, i2);
      std::vector<int> direzioni {0,1,2};
      direzioni.erase (direzioni.begin ()+dir);

      for (unsigned i = 0; i < direzioni.size (); ++i) {
        ray[i] = quadrant->p (direzioni[i], i1);
      }

      auto it0 = ray_cache.rays[dir].find (ray);
      if (it0 == ray_cache.rays[dir].end ()) continue;
      auto inters = it0->second.inters;

      for (int ii =0; ii<inters.size (); ii++) {
        if (inters[ii]>= x1 && inters[ii] <=x2) {
          fraction[j] = (inters[ii] - x1)/ (x2 - x1);
        }
      }
    }
  }

  return fraction;
}

void
poisson_boltzmann::normal_intersection (tmesh_3d::quadrant_iterator& quadrant,
                                        const ray_cache_t & ray_cache,
                                        int edge, std::array<double,3> &norm,
                                        double &frac)
{
  int dir = edge_axis[edge];
  int i1 = edge2nodes[2*edge ];
  int i2 = edge2nodes[2*edge + 1];
  double x1 = quadrant->p (dir, i1),
         x2 = quadrant->p (dir, i2);

  std::array<double,2> ray;
  std::vector<int> direzioni = {0,1,2};
  direzioni.erase (direzioni.begin ()+dir);

  for (unsigned i = 0; i < direzioni.size (); ++i) {
    ray[i] = quadrant->p (direzioni[i], i1);
  }

  frac = 0.5;

  auto it0 = ray_cache.rays[dir].find (ray);

  // In ns mode the caller only visits border_quad, where the ray always exists.
  // The membrane full sweep (implicit modes) may reach an edge with no cached
  // ray; keep the default frac = 0.5 rather than dereferencing end().
  if (it0 == ray_cache.rays[dir].end ())
    return;

  auto normali = it0->second.normals;
  auto inters = it0->second.inters;

  for (int ii =0; ii<inters.size (); ii++) {
    if (inters[ii]>= x1 && inters[ii] <=x2) {
      norm[0] = normali[0 + 3*ii];
      norm[1] = normali[1 + 3*ii];
      norm[2] = normali[2 + 3*ii];
      frac = (inters[ii] - x1)/ (x2 - x1);
    }
  }
}

int
poisson_boltzmann::classifyCube (tmesh_3d::quadrant_iterator& quadrant,
                                 double isolevel)
{
  int cubeindex = 0;
  int index = 1;
  double tmp = 0;
  constexpr double EPSILON = 1e-10;

  for (int ii : {
         0,1,3,2,4,5,7,6
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
poisson_boltzmann::classifyCube_fast (tmesh_3d::quadrant_iterator& quadrant,
                                      double isolevel)
{
  int cubeindex = 0;
  int index = 1;
  double tmp = 0;
  constexpr double EPSILON = 1e-10;

  for (int ii : {
         0,1,3,2,4,5,7,6
       }) {
    if ( (*epsilon_nodes)[quadrant->gt (ii)] < (isolevel - EPSILON)) cubeindex |= index;

    index *= 2;
  }

  // Cube is entirely in/out of the surface
  if (edgeTable[cubeindex] == 0)
    return -1;

  return cubeindex;
}

std::tuple<std::array<double,8>, std::array<double,8>, std::vector<int>,std::vector<int> >
poisson_boltzmann::classifyCube_flux (tmesh_3d::quadrant_iterator& quadrant,
                                      std::array<double,8>& tmp_phi,
                                      std::array<double,8>& tmp_eps)
{
  std::vector<int> edges {};
  std::vector<int> flux {};

  for (int ii = 0; ii < 8; ++ii) {
    if (! quadrant->is_hanging (ii)) {
      tmp_eps[ii]= (*epsilon_nodes)[quadrant->gt (ii)];
      tmp_phi[ii]= (*phi)[quadrant->gt (ii)];
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
    if (tmp_eps[edge2nodes[2*ii]] < tmp_eps[edge2nodes[2*ii +1]]) {
      flux.push_back (1);
      edges.push_back (ii);
    } else if (tmp_eps[edge2nodes[2*ii]] > tmp_eps[edge2nodes[2*ii +1]]) {
      flux.push_back (-1);
      edges.push_back (ii);
    }
  }

  return make_tuple (tmp_phi, tmp_eps, edges, flux);
}

std::tuple<std::array<double,8>, std::array<double,8>, std::vector<int>,std::vector<int> >
poisson_boltzmann::classifyCube_flux_fast (tmesh_3d::quadrant_iterator& quadrant,
    std::array<double,8>& tmp_phi,
    std::array<double,8>& tmp_eps)
{
  std::vector<int> edges {};
  std::vector<int> flux {};


  for (int ii = 0; ii < 8; ++ii) {
    tmp_eps[ii]= (*epsilon_nodes)[quadrant->gt (ii)];
    tmp_phi[ii]= (*phi)[quadrant->gt (ii)];
  }

  for (int ii = 0; ii < 12; ++ii) {
    if (tmp_eps[edge2nodes[2*ii]] < tmp_eps[edge2nodes[2*ii +1]]) {
      flux.push_back (1);
      edges.push_back (ii);
    } else if (tmp_eps[edge2nodes[2*ii]] > tmp_eps[edge2nodes[2*ii +1]]) {
      flux.push_back (-1);
      edges.push_back (ii);
    }
  }

  return make_tuple (tmp_phi, tmp_eps, edges, flux);
}

double wha (double eps1, double eps2, double frac)
{
  return 1.0/ (frac/eps1 + (1-frac)/eps2);
}

double flux_dir (double eps1, double eps2)
{
  return eps1 < eps2 ? 1 : -1;
}

double phi0 (double eps1, double eps2,
             double phi1, double phi2, double frac)
{
  return phi1 + frac*eps2* (phi2-phi1)/ (eps2*frac + eps1* (1-frac));
}

double areaTriangle (const std::array<std::array<double,3>,3> &triangle)
{

  double area;
  std::array<double,3> ab;
  std::array<double,3> ac;

  for (int i = 0; i < 3; ++i) {
    ab[i] = triangle[0][i] - triangle[1][i];
    ac[i] = triangle[0][i] - triangle[2][i];
  }

  area = 0.5* std::hypot (ab[1]*ac[2] - ab[2]*ac[1],
                          ab[0]*ac[2] - ab[2]*ac[0],
                          ab[1]*ac[0] - ab[0]*ac[1]);
  return area;
}

double SphercalAreaTriangle (const std::array<std::array<double,3>,3> &triangle)
{

  double area;
  double rs =2;
  double a,b,c, aa, bb, cc, E;
  double nA,nB,nC;
  std::array<double,3> A = triangle[0];
  std::array<double,3> B = triangle[1];
  std::array<double,3> C = triangle[2];

  nA = std::hypot (A[0], A[1], A[2]);
  nB = std::hypot (B[0], B[1], B[2]);
  nC = std::hypot (C[0], C[1], C[2]);

  a = std::inner_product (std::begin (A), std::end (A), std::begin (B), 0.0);
  b = std::inner_product (std::begin (A), std::end (A), std::begin (C), 0.0);
  c = std::inner_product (std::begin (B), std::end (B), std::begin (C), 0.0);

  a = a/ (nA*nB);
  b = b/ (nA*nC);
  c = c/ (nC*nB);

  a = std::acos (a);
  b = std::acos (b);
  c = std::acos (c);

  if (a>pi/2.0)
    a = pi - a;

  if (b>pi/2.0)
    b = pi - b;

  if (c>pi/2.0)
    c = pi - c;

  aa = std::acos ( (cos (a)-cos (b)*cos (c))/ (sin (b)*sin (c)));
  bb = std::acos ( (cos (b)-cos (c)*cos (a))/ (sin (a)*sin (c)));
  cc = std::acos ( (cos (c)-cos (b)*cos (a))/ (sin (b)*sin (a)));

  E = aa + bb + cc - pi;

  area = rs*rs * E;

  return area;
}

int
poisson_boltzmann::getTriangles (int cubeindex,
                                 std::array<std::array<int,3>,5> &triangles)
{
  int ntriang = 0;
  int i;
  triangles.fill ({}); //set matrix to zero

  // classifyCube returns -1 for a cube that no isosurface crosses (all eight
  // corners on the same side): there is no triangle to build, and indexing
  // triTable with it walks off both tables and overflows the caller's array.
  // In ns mode this never happens on a border quadrant -- being outside the
  // molecular surface and being solvent are the same thing there, so a border
  // cube always has corners on both sides.  With an implicit membrane they are
  // not: a cube straddling the protein surface *inside* the slab has every
  // corner below eps_out (protein or membrane) and is entirely non-solvent.
  // Zero triangles is the right answer for it -- the solvent boundary really is
  // absent there -- and the flux loop over the epsilon jumps is unaffected.
  if (cubeindex < 0)
    return 0;

  for (i=0; triTable[cubeindex][i]!=-1; i+=3) {
    // save all the assigned indexes
    triangles[ntriang][0] = triTable[cubeindex][i ];
    triangles[ntriang][1] = triTable[cubeindex][i+1];
    triangles[ntriang][2] = triTable[cubeindex][i+2];
    ntriang++;
  }

  return ntriang;
}

// On-shell free energy of the membrane problem:
//
//   dG* = 1/2 sum_i q_i phi*(r_i) + 1/2 Vbar int_{S_int^up} D*_n,solv dS
//
// This is a *different* formula from the one energy()/energy_fast() evaluate,
// not an extra term added on top of it: it is the on-shell value of the LPB
// functional written with the regional bulk reference phi_bar, so the two are
// not comparable and the pre-membrane path is deliberately left untouched.
//
// S_int^up is the internal boundary of the upper bath, i.e. everything bounding
// the region where phi_bar = Vbar except the z+ electrode and the lateral
// periodic faces -- not just the flat top of the lipid slab.  In this geometry
// the protein sticks well above the membrane, so its surface is a large part of
// it.  Discretely the region is exactly {z > z_mid}, the same split that sets
// phi_bar in create_markers, which makes the flux term consistent by
// construction with the RHS source that produced phi.
//
// Only called when a transmembrane potential is applied, so at Vbar = 0 nothing
// here runs and the existing energies stay bit-identical.
void
poisson_boltzmann::energy_membrane (ray_cache_t & ray_cache)
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);

  if (rank == 0)
    std::cout << "\n============ [ Membrane On-Shell Energy ] =================\n";

  // ---------------------------------------------------------------------
  // Term 1:  1/2 sum_i q_i phi*(r_i)
  // ---------------------------------------------------------------------
  // Trilinear interpolation of phi at the charge positions, with the same
  // weights as write_potential_on_atoms_fast (complementary corner volumes over
  // the cell volume) and the same averaging of hanging nodes over their parents.
  // lookup_table is filled in create_rhs via search_points and assigns every
  // atom to exactly one rank (controlla_coordinate rejects the upper faces), so
  // summing the local contributions and reducing is not double counting.
  double energy_charge = 0.0;

  for (auto it = lookup_table.begin (); it != lookup_table.end (); ++it) {
    const double volume = (it->second.p (0, 7) - it->second.p (0, 0)) *
                          (it->second.p (1, 7) - it->second.p (1, 0)) *
                          (it->second.p (2, 7) - it->second.p (2, 0));

    double phi_on_atom = 0.0;

    for (int ii = 0; ii < 8; ++ii) {
      const double weight = std::abs ((pos_atoms[it->first][0] - it->second.p (0, 7 - ii)) *
                                      (pos_atoms[it->first][1] - it->second.p (1, 7 - ii)) *
                                      (pos_atoms[it->first][2] - it->second.p (2, 7 - ii))) / volume;

      if (! it->second.is_hanging (ii)) {
        phi_on_atom += (*phi)[it->second.gt (ii)] * weight;
      } else {
        double phi_hang_nodes = 0.0;

        for (int jj = 0; jj < it->second.num_parents (ii); ++jj)
          phi_hang_nodes += (*phi)[it->second.gparent (jj, ii)] / it->second.num_parents (ii);

        phi_on_atom += phi_hang_nodes * weight;
      }
    }

    energy_charge += charge_atoms[it->first] * phi_on_atom;
  }

  energy_charge *= 0.5;

  // ---------------------------------------------------------------------
  // Term 2:  1/2 Vbar int_{S_int^up} D*_n,solv dS
  // ---------------------------------------------------------------------
  // Same flux loop as energy(), minus the inner sweep over the atoms (which is
  // what makes energy() expensive): classifyCube_flux already returns exactly
  // the quantity needed.  Its fl_dir is +1 when epsilon grows along the edge
  // axis, so the normal it uses points from low to high epsilon, i.e. from the
  // membrane/solute into the solvent -- the n_hat of the notes -- and
  //
  //   tmp_flux = -(phi_2 - phi_1) * eps_harmonic * fl_dir * area_h[axis]
  //            = (D . n_hat) * (face area / 4)
  //
  // is D_n,solv on the dual-cell face element.  The epsilon field carries a
  // factor 4*pi (see eps_in/eps_out above), hence the same 1/(4*pi) that turns
  // charge_pol into "Flux charge [e]".
  //
  // The upper bath is the region where phi_bar = Vbar, i.e. z > z_mid with the
  // very same midplane used in create_markers, so the crossing points are
  // filtered on their z coordinate.  The whole-interface flux is accumulated
  // too, as a diagnostic to compare against the Flux charge of energy().
  const double z_mid = 0.5 * (z_mem_bot + z_mem_top);

  double flux_up = 0.0, flux_all = 0.0;

  std::array<double,3> h {0}, area_h {0};
  std::array<double,3> V, N;
  std::array<double,8> tmp_eps, tmp_phi;
  std::vector<int> edg, fl_dir;

  auto quadrant = this->tmsh.begin_quadrant_sweep ();

  // One quadrant's contribution to the two flux sums.  In ns mode every crossing
  // is on the NanoShaper surface, so the fraction always comes from the ray
  // cache.  In the implicit modes the dielectric interface is partly a box: the
  // slab top/bottom planes.  A dielectric jump can only sit on such a plane along
  // a z edge whose endpoints straddle it, and there the fraction is analytic
  // (linear in z, the plane is flat); every other jump is still the protein
  // NanoShaper surface and uses the ray cache as before.  Inside the protein
  // footprint both sides of z_mem_top/bot are e_in, so no jump is selected there
  // and the analytic branch only ever fires on genuine slab/solvent interfaces.
  auto accumulate_quadrant = [&] (bool implicit_box) {
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
      bool on_box_plane = false;

      // edge2nodes lists the lower corner first, so p(axis,i2) - p(axis,i1) = h.
      if (implicit_box && axis == 2) {
        const double z1 = quadrant->p (2, i1);

        for (const double zp : {z_mem_bot, z_mem_top})
          if (z1 < zp && zp < z1 + h[2]) {
            fract = (zp - z1) / h[2];
            on_box_plane = true;
          }
      }

      if (! on_box_plane)
        normal_intersection (quadrant, ray_cache, edge, N, fract);

      V = {quadrant->p (0, i1), quadrant->p (1, i1), quadrant->p (2, i1)};
      V[axis] += fract * h[axis];

      const double tmp_flux =
        - (tmp_phi[i2] - tmp_phi[i1]) * wha (tmp_eps[i1], tmp_eps[i2], fract)
        * fl_dir[ip] * area_h[axis];

      flux_all += tmp_flux;

      if (V[2] > z_mid)
        flux_up += tmp_flux;
    }
  };

  if (membrane_mode == MEM_MODE_NS) {
    // ns: the interface is the NanoShaper surface, tracked by border_quad.
    for (const int ii : border_quad) {
      quadrant[ii];
      accumulate_quadrant (false);
    }
  } else {
    // implicit / implicit_charged: the slab planes are not in border_quad (which
    // is built from the protein surface), so sweep the whole mesh and let
    // classifyCube_flux pick up every dielectric jump.
    for (quadrant = this->tmsh.begin_quadrant_sweep ();
         quadrant != this->tmsh.end_quadrant_sweep ();
         ++quadrant)
      accumulate_quadrant (true);
  }

  flux_up /= (4.0 * pi);
  flux_all /= (4.0 * pi);

  auto reduce_double = [&] (double &x) {
    MPI_Reduce (rank == 0 ? MPI_IN_PLACE : &x, &x, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
  };

  reduce_double (energy_charge);
  reduce_double (flux_up);
  reduce_double (flux_all);

  const double energy_flux = 0.5 * applied_potential * flux_up;

  if (rank == 0) {
    constexpr int label_width = 50;
    constexpr int precision = 16;

    std::cout << std::left << std::setw (label_width) << "  Applied potential [kT/e]:"
              << std::setprecision (precision) << applied_potential << "\n";

    std::cout << std::left << std::setw (label_width) << "  Upper-bath midplane z_mid [A]:"
              << std::setprecision (precision) << z_mid << "\n";

    std::cout << std::left << std::setw (label_width) << "  Flux through S_int^up [e]:"
              << std::setprecision (precision) << flux_up << "\n";

    std::cout << std::left << std::setw (label_width) << "  Flux through the whole interface [e]:"
              << std::setprecision (precision) << flux_all << "\n";

    std::cout << std::left << std::setw (label_width) << "  Term 1, 1/2 sum q_i phi(r_i) [kT]:"
              << std::setprecision (precision) << energy_charge << "\n";

    std::cout << std::left << std::setw (label_width) << "  Term 2, 1/2 Vbar int D_n dS [kT]:"
              << std::setprecision (precision) << energy_flux << "\n";

    std::cout << std::left << std::setw (label_width) << "  Membrane on-shell energy dG* [kT]:"
              << std::setprecision (precision) << (energy_charge + energy_flux) << "\n";

    std::cout << "===========================================================\n";
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
poisson_boltzmann::coulomb_boundary_conditions (double x, double y, double z)
{
  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/ (e*e); //adim e_out
  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/ (e_0*e_out*kb*T);

  double dist = 0.0;
  double pot = 0.0;
  double k = std::sqrt (k2);

  for (const NS::Atom& i : atoms) {
    dist = std::hypot ( (i.pos[0] - x), (i.pos[1] - y), (i.pos[2] - z));
    pot += i.charge*exp (-k*dist)/ (dist*eps_out);
  }

  return pot;
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

bool
poisson_boltzmann::controlla_coordinate (int i, const p8est_quadrant_t *quadrant)
{

  double tol = p4esttol * (rr[0]-ll[0]);

  if (mesh_shape == 2)
    tol = p4esttol * (r_c[0]-l_c[0]);


  bool retval = false;

  p8est_quadrant_t node;
  int ii, jj;
  double vxyz [3 * 8] = {0,0,0, 0,0,0, 0,0,0, 0,0,0,
                         0,0,0, 0,0,0, 0,0,0, 0,0,0
                        };

  for (ii = 0; ii < 8; ++ii) {
    p8est_quadrant_corner_node (quadrant, ii, &node);
    p8est_qcoord_to_vertex (tmsh.p8est->connectivity, 0,
                            node.x, node.y, node.z, & (vxyz[3 * ii]));
  }

  double l, r, t, b, f, bk;

  l = vxyz[0];
  r = vxyz[3*7];

  f = vxyz[1];
  bk = vxyz[3*7 +1];

  b = vxyz[2];
  t = vxyz[3*7 +2];

  // retval = (atoms[i].pos[0] > l- tol) && (atoms[i].pos[0] <= r- tol); //make sure that the charge is assigned only once
  // retval = retval && (atoms[i].pos[1] > f- tol) && (atoms[i].pos[1] <= bk- tol);
  // retval = retval && (atoms[i].pos[2] > b- tol) && (atoms[i].pos[2] <= t- tol);
  retval = (pos_atoms[i][0] > l- tol) && (pos_atoms[i][0] <= r- tol); //make sure that the charge is assigned only once
  retval = retval && (pos_atoms[i][1] > f- tol) && (pos_atoms[i][1] <= bk- tol);
  retval = retval && (pos_atoms[i][2] > b- tol) && (pos_atoms[i][2] <= t- tol);

  return retval;
}


int
poisson_boltzmann::cerca_atomo (p8est_t * p4est,
                                p4est_topidx_t which_tree,
                                p8est_quadrant_t * quadrant,
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
void
poisson_boltzmann::write_dataset (ray_cache_t & ray_cache)
{
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0)
    std::cout << "\n================ [ Computing vertex quantities ] =================\n";

  const double eps_in  = 4.0 * pi * e_0 * e_in  * kb * T * Angs / (e * e);
  const double eps_out = 4.0 * pi * e_0 * e_out * kb * T * Angs / (e * e);

  std::array<double, 3> N;
  std::array<double, 3> V;
  std::array<double, 3> h;
  std::array<double, 8> tmp_eps;
  std::array<double, 8> tmp_phi;
  std::vector<int> edg;
  std::vector<int> fl_dir;

  int cubeindex = -1;
  std::unordered_map<EdgeKey, VertexData, EdgeHash> edgeMap;

  auto quadrant = this->tmsh.begin_quadrant_sweep ();

  if (!border_quad.empty ()) {
    quadrant[border_quad[0]];
    h[0] = quadrant->p (0, 7) - quadrant->p (0, 0);
    h[1] = quadrant->p (1, 7) - quadrant->p (1, 0);
    h[2] = quadrant->p (2, 7) - quadrant->p (2, 0);
  }

  for (const int ii : border_quad) {
    quadrant[ii];
    cubeindex = classifyCube_fast (quadrant, eps_out);
    std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux_fast (quadrant, tmp_phi, tmp_eps);

    for (int ip = 0; ip < (int)edg.size (); ++ip) {
      const int edge = edg[ip];
      const int axis = edge_axis[edge];
      const int i1   = edge2nodes[2 * edge];
      const int i2   = edge2nodes[2 * edge + 1];

      EdgeKey triEdges{
        static_cast<int>(std::min (quadrant->gt (i1), quadrant->gt (i2))),
        static_cast<int>(std::max (quadrant->gt (i1), quadrant->gt (i2)))
      };

      auto & vertexData = edgeMap[triEdges];

      vertexData.axis = axis;

      double fract = 0.0;
      normal_intersection (quadrant, ray_cache, edge, N, fract);

      for (int jj = 0; jj < 3; ++jj) {
        int nu = (axis + jj) % 3;
        vertexData.N[jj] = N[nu];
      }

      vertexData.pos1[0] = quadrant->p (0, i1);
      vertexData.pos1[1] = quadrant->p (1, i1);
      vertexData.pos1[2] = quadrant->p (2, i1);

      vertexData.pos0[0] = vertexData.pos1[0];
      vertexData.pos0[1] = vertexData.pos1[1];
      vertexData.pos0[2] = vertexData.pos1[2];
      vertexData.pos0[edge_axis[edge]] += fract * h[edge_axis[edge]];

      vertexData.phi1 = tmp_phi[i1];
      vertexData.phi2 = tmp_phi[i2];
      vertexData.eps1 = tmp_eps[i1];
      vertexData.eps2 = tmp_eps[i2];

      vertexData.alpha = fract;
      vertexData.phi0  = phi0 (tmp_eps[i1], tmp_eps[i2], tmp_phi[i1], tmp_phi[i2], fract);

      const int i1_nn1   = edge2nodes_nn1[2 * edge];
      const int i1_nn2   = edge2nodes_nn1[2 * edge + 1];
      const int i2_nn1   = edge2nodes_nn2[2 * edge];
      const int i2_nn2   = edge2nodes_nn2[2 * edge + 1];
      const int ind1_nn1 = edge2index_nn1[2 * edge];
      const int ind1_nn2 = edge2index_nn1[2 * edge + 1];
      const int ind2_nn1 = edge2index_nn2[2 * edge];
      const int ind2_nn2 = edge2index_nn2[2 * edge + 1];

      vertexData.phi1_nn[ind1_nn1] = tmp_phi[i1_nn1];
      vertexData.phi1_nn[ind1_nn2] = tmp_phi[i1_nn2];
      vertexData.eps1_nn[ind1_nn1] = tmp_eps[i1_nn1];
      vertexData.eps1_nn[ind1_nn2] = tmp_eps[i1_nn2];
      vertexData.phi2_nn[ind2_nn1] = tmp_phi[i2_nn1];
      vertexData.phi2_nn[ind2_nn2] = tmp_phi[i2_nn2];
      vertexData.eps2_nn[ind2_nn1] = tmp_eps[i2_nn1];
      vertexData.eps2_nn[ind2_nn2] = tmp_eps[i2_nn2];
    }
  }

  {
    std::ostringstream fname;
    fname << "vertexdata_rank" << rank << ".csv";
    std::ofstream ofs (fname.str ());
    if (!ofs) {
      std::cerr << "Error: cannot open " << fname.str () << " for writing\n";
      return;
    }

    ofs << "phi1,eps1,alpha,N_nu,N_nu1,N_nu2,"
        << "phi2,eps2,"
        << "phi_perp_1_p_1,eps_perp_1_p_1,"
        << "phi_perp_1_m_1,eps_perp_1_m_1,"
        << "phi_perp_2_p_1,eps_perp_2_p_1,"
        << "phi_perp_2_m_1,eps_perp_2_m_1,"
        << "phi_perp_1_p_2,eps_perp_1_p_2,"
        << "phi_perp_1_m_2,eps_perp_1_m_2,"
        << "phi_perp_2_p_2,eps_perp_2_p_2,"
        << "phi_perp_2_m_2,eps_perp_2_m_2,"
        << "x0,y0,z0,phi0,x1,y1,z1,axis\n";

    ofs << std::scientific << std::setprecision (8);

    for (const auto & [k, vd] : edgeMap) {
      ofs << vd.phi1 << "," << vd.eps1 << "," << vd.alpha << ","
          << vd.N[0] << "," << vd.N[1] << "," << vd.N[2] << ","
          << vd.phi2 << "," << vd.eps2 << ",";

      ofs << vd.phi1_nn[0] << "," << vd.eps1_nn[0] << ","
          << vd.phi1_nn[1] << "," << vd.eps1_nn[1] << ","
          << vd.phi1_nn[2] << "," << vd.eps1_nn[2] << ","
          << vd.phi1_nn[3] << "," << vd.eps1_nn[3] << ",";

      ofs << vd.phi2_nn[0] << "," << vd.eps2_nn[0] << ","
          << vd.phi2_nn[1] << "," << vd.eps2_nn[1] << ","
          << vd.phi2_nn[2] << "," << vd.eps2_nn[2] << ","
          << vd.phi2_nn[3] << "," << vd.eps2_nn[3] << ",";

      ofs << vd.pos0[0] << "," << vd.pos0[1] << "," << vd.pos0[2] << "," << vd.phi0 << ",";
      ofs << vd.pos1[0] << "," << vd.pos1[1] << "," << vd.pos1[2] << "," << vd.axis << "\n";
    }
  }

  MPI_Barrier (MPI_COMM_WORLD);

  if (rank == 0) {
    std::ofstream final ("vertexdata.csv");
    if (!final) {
      std::cerr << "Error: cannot open vertexdata.csv for writing\n";
      return;
    }

    final << "phi1,eps1,alpha,N_nu,N_nu1,N_nu2,"
          << "phi2,eps2,"
          << "phi_perp_1_p_1,eps_perp_1_p_1,"
          << "phi_perp_1_m_1,eps_perp_1_m_1,"
          << "phi_perp_2_p_1,eps_perp_2_p_1,"
          << "phi_perp_2_m_1,eps_perp_2_m_1,"
          << "phi_perp_1_p_2,eps_perp_1_p_2,"
          << "phi_perp_1_m_2,eps_perp_1_m_2,"
          << "phi_perp_2_p_2,eps_perp_2_p_2,"
          << "phi_perp_2_m_2,eps_perp_2_m_2,"
          << "x0,y0,z0,phi0,x1,y1,z1,axis\n";

    final << std::scientific << std::setprecision (8);

    for (int r = 0; r < size; ++r) {
      std::ostringstream fname;
      fname << "vertexdata_rank" << r << ".csv";
      std::ifstream ifs (fname.str ());
      if (!ifs) continue;

      std::string line;
      bool first_line = true;
      while (std::getline (ifs, line)) {
        if (first_line) { first_line = false; continue; }
        final << line << "\n";
      }
      ifs.close ();
      std::filesystem::remove (fname.str ());
    }

    std::cout << "\nAll data merged into vertexdata.csv\n";
  }
}
