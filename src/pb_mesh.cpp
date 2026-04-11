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
#include <p8est.h>
#include <random>
#include <cmath>
#include <numeric>

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

  if (mesh_shape !=2) {
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

    // ── Membrane slab mesh ────────────────────────────────────────────────
    // Cubic box of side max(cell_length_x, cell_length_y).
    // The membrane fills the xy face; the slab in z covers both the membrane
    // and the protein (which may extend above/below the membrane).
    //
    // Level structure:
    //   outlevel  = maxlevel - nlev_sol   (far solvent)
    //   scale_level = maxlevel - nlev_mem  (membrane slab)
    //   maxlevel                           (near molecular surface, loc_refinement)

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
    double z_lip_min  = pos_lipid_atoms.empty () ? z_prot_min : (*minmax_lip_z.first)[2];
    double z_lip_max  = pos_lipid_atoms.empty () ? z_prot_max : (*minmax_lip_z.second)[2];

    double z_slab_bot = std::min (z_prot_min, z_lip_min) - maxrad - 2.0 * prb_radius;
    double z_slab_top = std::max (z_prot_max, z_lip_max) + maxrad + 2.0 * prb_radius;
    double z_slab_cc  = 0.5 * (z_slab_top + z_slab_bot);

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

    // Center of the membrane patch in xy (overrides preamble's protein-only cc)
    cc[0] = 0.5 * (x_lip_max + x_lip_min);
    cc[1] = 0.5 * (y_lip_max + y_lip_min);
    cc[2] = z_slab_cc;

    // --- square outer box (cubic: same side in x, y, z) ---
    ll[0] = cc[0] - box_side / 2.0;  rr[0] = cc[0] + box_side / 2.0;
    ll[1] = cc[1] - box_side / 2.0;  rr[1] = cc[1] + box_side / 2.0;
    ll[2] = cc[2] - box_side / 2.0;  rr[2] = cc[2] + box_side / 2.0;

    // --- levels ---
    maxlevel    = (int) std::ceil (std::log2 (box_side * scale_max));
    scale       = (double)(1 << maxlevel) / box_side;
    unilevel    = maxlevel - nlev_sol;
    outlevel    = maxlevel - nlev_sol;
    scale_level = maxlevel - nlev_mem;
    scale_level_min_box = maxlevel;

    // --- l_c / r_c: slab bounds for the refinement callback ---
    // xy spans the full box (condition always satisfied → only z filters)
    l_c[0] = ll[0];      r_c[0] = rr[0];
    l_c[1] = ll[1];      r_c[1] = rr[1];
    l_c[2] = z_slab_bot; r_c[2] = z_slab_top;

    // --- l_cr / r_cr: NS box, aligned to the grid (square in xy) ---
    int nx = (int) ((box_side / 2.0) * scale + 0.5);
    int ny = (int) ((box_side / 2.0) * scale + 0.5);
    int nz = (int) ((z_slab_top - z_slab_cc) * scale + 0.5);
    l_cr[0] = cc[0] - nx / scale;     r_cr[0] = cc[0] + nx / scale;
    l_cr[1] = cc[1] - ny / scale;     r_cr[1] = cc[1] + ny / scale;
    l_cr[2] = z_slab_cc - nz / scale; r_cr[2] = z_slab_cc + nz / scale;

    if (rank == 0) {
      std::cout << "cx: " << cc[0] << "  cy: " << cc[1] << "  cz: " << cc[2] << '\n';
      std::cout << "x: "  << ll[0] << " , " << rr[0] << '\n';
      std::cout << "y: "  << ll[1] << " , " << rr[1] << '\n';
      std::cout << "z: "  << ll[2] << " , " << rr[2] << '\n';
      std::cout << "z_slab: " << z_slab_bot << " , " << z_slab_top << '\n';
      std::cout << "maxlevel: "    << maxlevel
                << "  scale_level: " << scale_level
                << "  unilevel: "    << unilevel << '\n';
    }

    simple_conn_num_vertices = 8;
    simple_conn_num_trees    = 1;

    simple_conn_p = std::make_unique<double[]> (simple_conn_num_vertices * 3);
    simple_conn_t = std::make_unique<p4est_topidx_t[]> (simple_conn_num_vertices);

    auto tmp_p = {ll[0], ll[1], ll[2],
                  rr[0], ll[1], ll[2],
                  ll[0], rr[1], ll[2],
                  rr[0], rr[1], ll[2],
                  ll[0], ll[1], rr[2],
                  rr[0], ll[1], rr[2],
                  ll[0], rr[1], rr[2],
                  rr[0], rr[1], rr[2]};
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

