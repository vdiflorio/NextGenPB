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


void
poisson_boltzmann::create_mesh ()
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);

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


  
  double net_charge = std::accumulate(charge_atoms.begin(), charge_atoms.end(), 0.0);
  int num_atoms = charge_atoms.size ();

  if (rank == 0) {
    std::cout << "\n========== [ System Information ] ==========\n";
    std::cout << "  Number of atoms    : " << num_atoms << '\n';
    std::cout << "  Size protein [Å]   : ";
    std::cout << "[" << (*minmax_x.second)[0] - (*minmax_x.first)[0] + 2*maxradius << ", " 
                     << (*minmax_y.second)[0] - (*minmax_y.first)[0] + 2*maxradius << ", " 
                     << (*minmax_z.second)[0] - (*minmax_z.first)[0] + 2*maxradius << "]\n"; 
    std::cout << "  Net charge         : " << std::scientific << net_charge << std::defaultfloat << '\n';
    std::cout << "  Solvent epsilon    : " << e_in << '\n';
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

  if (!g.search ("--pqrfile")) {
    if (rank == 0) {
      std::cout << "Warning: No pqr file selected, using the default one." <<
                "\nTo select one use --pqrfile option followed by the desired one." << std::endl;
    }
  }

  pqrfilename = g.next ("../../data/1CCM.pqr");

  //Check that the pqr file exists
  if (rank == 0)
    std::cout << "Selected pqr file:        " << pqrfilename << std::endl;

  std::ifstream pqrfile (pqrfilename);

  if (!pqrfile) {
    if (rank == 0) {
      std::cerr << "Cannot find the pqr file" << std::endl;
      return 1;
    }
  }

  if (!g.search ("--potfile")) {
    if (rank == 0) {
      std::cout << "Warning: No parameters file selected, using the default one." <<
                "\nTo select one use --potfile option followed by the desired one." << std::endl;
    }
  }

  optionsfilename = g.next ("../../data/options.pot");

  //Check that the pot file exists
  if (rank == 0)
    std::cout << "Selected parameters file: " << optionsfilename << std::endl;

  std::ifstream optionsfile (optionsfilename);

  if (!optionsfile) {
    if (rank == 0) {
      std::cerr << "Cannot find the options file" << std::endl;
      return 1;
    }
  }

  //Read the options from the file
  GetPot g2 (optionsfilename.c_str ());

  const std::string mesh_options = "mesh/";
  maxlevel = g2 ( (mesh_options + "maxlevel").c_str (), 9);
  minlevel = g2 ( (mesh_options + "minlevel").c_str (), 3);
  unilevel = g2 ( (mesh_options + "unilevel").c_str (), 5);
  outlevel = g2 ( (mesh_options + "outlevel").c_str (), 1);
  loc_refinement = g2 ( (mesh_options + "loc_refinement").c_str (), 0);
  mesh_shape = g2 ( (mesh_options + "mesh_shape").c_str (), 1);
  refine_box = g2 ( (mesh_options + "refine_box").c_str (), 0);
  rand_center = g2 ( (mesh_options + "rand_center").c_str (), 0);

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

  const std::string model_options = "model/";
  linearized = g2 ( (model_options + "linearized").c_str (), 1);
  bc = g2 ( (model_options + "bc_type").c_str (), 1);
  ionic_strength = g2 ( (model_options + "ionic_strength").c_str (), 0.145);
  e_in = g2 ( (model_options + "molecular_dielectric_constant").c_str (), 2.);
  e_out = g2 ( (model_options + "solvent_dielectric_constant").c_str (), 80.);
  T = g2 ( (model_options + "T").c_str (), 298.15);
  calc_energy = g2 ( (model_options + "calc_energy").c_str (), 2);
  atoms_write = g2 ( (model_options + "atoms_write").c_str (), 0);
  surf_write = g2 ( (model_options + "surf_write").c_str (), 0);
  surf_write = g2 ( (model_options + "surf_write").c_str (), 0);
  map_type = g2 ( (model_options + "map_type").c_str (), "vtu");
  potential_map = g2 ( (model_options + "potential_map").c_str (), 0);
  eps_map = g2 ( (model_options + "eps_map").c_str (), 0);
  const std::string surf_options = "surface/";
  int surf_type_num = g2 ( (surf_options + "surface_type").c_str (), 1);

  if (surf_type_num == 1) surf_type = NS::skin;
  else if (surf_type_num == 0) surf_type = NS::ses;
  else if (surf_type_num == 2) surf_type = NS::blobby;
  else surf_type = NS::ses;

  surf_param = g2 ( (surf_options + "surface_parameter").c_str (), 0.45);
  stern_layer_surf = g2 ( (surf_options + "stern_layer_surf").c_str (), 0);
  stern_layer = g2 ( (surf_options + "stern_layer_thickness").c_str (), 2.);
  num_threads = g2 ( (surf_options + "number_of_threads").c_str (), 1);

  const std::string alg_options = "algorithm/";
  linear_solver_name = g2 ( (alg_options + "linear_solver").c_str (), "lis");
  linear_solver_options = g2 ( (alg_options + "solver_options").c_str (), "-i cg -p ilu [0] -tol 1.e-12");

  const std::string out_options = "output/";
  p4estfilename = g2 ( (out_options + "p4estfilename").c_str (), "poisson_boltzmann_p4est");
  markerfilename = g2 ( (out_options + "markerfilename").c_str (), "poisson_boltzmann_marker_0");
  surffilename = g2 ( (out_options + "surffilename").c_str (), "poisson_boltzmann_surface_0");

  return 0;
}

void
poisson_boltzmann::print_options ()
{
  std::cout << "\nChosen options: " << std::endl;

  if (mesh_shape == 0) {
    std::cout << "Mesh shape = cubic at max " << perfil2 << " perfil" << " and refined box at "<< perfil1 << " perfil" << std::endl;
    std::cout << "scale: " << scale << std::endl;
    std::cout << "with the following domain vertices: " << std::endl;
  } else if (mesh_shape == 1) {
    std::cout << "Mesh shape = cubic at "<< perfil1 << "perfil" << std::endl;
    std::cout << "minlevel = " << minlevel << "\nmaxlevel = " << maxlevel << std::endl;
    std::cout << "unilevel = " << unilevel << std::endl;
    std::cout << "with the following domain vertices: " << std::endl;
  } else if (mesh_shape == 2) {
    std::cout << "Manual setting of the mesh vertices:: " << std::endl;
    std::cout << "minlevel = " << minlevel << "\nmaxlevel = " << maxlevel << std::endl;
    std::cout << "unilevel = " << unilevel << std::endl;
    std::cout << "with the following domain vertices: " << std::endl;
  } else if (mesh_shape == 3) {
    std::cout << "Focusing centered at:" << "\n x = " << cc_focusing[0]
              << "\n y = " << cc_focusing[1] << "\n z = " << cc_focusing[2] << std::endl;
    std::cout << "scale: " << scale << std::endl;
  } else
    std::cout << "Manual setting of the mesh vertices:" << std::endl;


  std::cout << "Linearized model = " << linearized << std::endl;

  if (bc == 1)
    std::cout << "Dirichlet boundary conditions" << std::endl;
  else if (bc == 2)
    std::cout << "Coulombic boundary conditions" << std::endl;
  else
    std::cout << "Neumann boundary conditions" << std::endl;

  std::cout << "e_in = " << e_in << "\ne_out = " << e_out
            << "\nionic_strenght = " << ionic_strength << std::endl;

  std::cout << "Linear solver = " << linear_solver_name << std::endl;
  std::cout << "Linear solver options = " << linear_solver_options << std::endl;

  std::cout << "Chosen surface: " << surf_type << std::endl;
  std::cout << "Surface parameter: " << surf_param << std::endl;

  if (stern_layer_surf == 1)
    std::cout << "Stern layer thickness: " << stern_layer << std::endl;

  std::cout << "Number of threads for nanoshaper: " << num_threads << std::endl;
}

void
poisson_boltzmann::read_atoms_from_pqr (std::basic_istream<char> &inputfile)
{
  static NS::Atom a;
  atoms.clear ();

  while (inputfile >> a)
    atoms.push_back (a);
  
  if (atoms.size() < 4) {
    auto comp = [] (const NS::Atom &a1, const NS::Atom &a2) -> bool {
      return a1.radius < a2.radius;
    };

    auto max_iter = std::max_element(atoms.begin(), atoms.end(), comp);
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
    }};

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
      atoms.push_back(dummy);
    }
  }
}

void
poisson_boltzmann::read_atoms_from_class ()
{
  static std::array<double,3> pos;

  for (const NS::Atom& i : atoms) {
    pos[0] = i.pos[0];
    pos[1] = i.pos[1];
    pos[2] = i.pos[2];
    pos_atoms.push_back (pos);
    charge_atoms.push_back (i.charge);
    r_atoms.push_back (i.radius);
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

  // Effettuare il broadcast del vettore
  MPI_Bcast (pos_atoms.data (), size_vec * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (charge_atoms.data (), size_vec, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (r_atoms.data (), size_vec, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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
            >> a.ai.resName
            >> a.ai.resNum
            >> a.pos[0]
            >> a.pos[1]
            >> a.pos[2]
            >> a.charge
            >> a.radius;

  if (a.radius < 1.e-5)
    a.radius = 1.0;

  a.radius2 = a.radius*a.radius;
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

  if (size == 1 || surf_type == 2)
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
              if (surf_type == 2)
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

          if (this->surf_type == 2) {
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
              if (surf_type == 2)
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

          if (this->surf_type == 2) {
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

  /////////////////////////////////////////////////////////
  //reactions
  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/ (e_0*e_out*kb*T);

  this->marker.assign (this->tmsh.num_local_quadrants (), 0.0); //marker = 0 -> in

  if (stern_layer_surf == 1) {
    this->marker_k.assign (this->tmsh.num_local_quadrants (), 1.0); //marker = 1 -> out stern
    this->reaction.assign (tmsh.num_local_quadrants (), eps_out*k2);
  }

  // markn = std::make_unique<distributed_vector> (tmsh.num_owned_nodes ());
  // markn->get_owned_data ().assign (tmsh.num_owned_nodes (), 0.0); //markn = 0 -> out
  epsilon_nodes = std::make_unique<distributed_vector> (tmsh.num_owned_nodes (),mpicomm);
  epsilon_nodes->get_owned_data ().assign (tmsh.num_owned_nodes (), eps_out);

  reaction_nodes = std::make_unique<distributed_vector> (tmsh.num_owned_nodes (),mpicomm);
  reaction_nodes->get_owned_data ().assign (tmsh.num_owned_nodes (), eps_out*k2);

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
          if (this->is_in_ns_surf (ray_cache, x, y, z, 2) > 0.5) { //inside the molecule
            ++num_int_nodes;
            ++num_int_nodes_stern;
            (*epsilon_nodes)[local_num] = eps_in;
            (*reaction_nodes)[local_num] = 0.0;
          } else if (this->is_in_ns_surf (ray_cache, x, y, z, 2) < -0.5) {
            ray_cache.num_req_rays[2]++;

            // std::array<double, 2> ray {x,y};

            ray_cache.rays_list[2].insert ({x,y});
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

  // if (rank == 0)
  //   std::cout << "\nStarting MUMPS solution" << std::endl;

  // diffusion
  double eps_in = 4.0*pi*e_0*e_in*kb*T*Angs/ (e*e); //adim e_in
  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/ (e*e); //adim e_out

  /////////////////////////////////////////////////////////
  //reactions
  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/ (e_0*e_out*kb*T);


  if (size >1) {
    bim3a_solution_with_ghosts (tmsh, *epsilon_nodes, replace_op);

    if (stern_layer_surf == 0)
      bim3a_solution_with_ghosts (tmsh, (*reaction_nodes), replace_op);
  }


  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  rho_fixed = std::make_unique<distributed_vector> (tmsh.num_owned_nodes (), mpicomm);
  rho_fixed->get_owned_data ().assign (tmsh.num_owned_nodes (), 0.0);

  std::unique_ptr<distributed_vector> ones =
    std::make_unique<distributed_vector> (tmsh.num_owned_nodes (),mpicomm);
  ones->get_owned_data ().assign (tmsh.num_owned_nodes (), 1.0);

  std::vector<double> const_ones (tmsh.num_local_quadrants (), 1.0);

  std::unique_ptr<distributed_vector> vol_patch =
    std::make_unique<distributed_vector> (tmsh.num_owned_nodes (),mpicomm);

  if (size >1)
    bim3a_solution_with_ghosts (tmsh, *ones, replace_op);

  bim3a_rhs (tmsh, const_ones, *ones, *vol_patch);

  if (size >1)
    vol_patch->assemble ();

  MPI_Barrier (mpicomm);
  auto start_rho = std::chrono::steady_clock::now ();


  search_points ();

  for (auto it = lookup_table.begin (); it!=lookup_table.end (); ++it) {
    //linear approx:
    double volume = (it->second.p (0, 7) - it->second.p (0, 0)) *
                    (it->second.p (1, 7) - it->second.p (1, 0)) *
                    (it->second.p (2, 7) - it->second.p (2, 0)); //volume


    {
      for (int ii = 0; ii < 8; ++ii) {
        double weigth = std::abs ( (pos_atoms[it->first][0] - it->second.p (0, 7-ii))*
                                   (pos_atoms[it->first][1] - it->second.p (1, 7-ii))*
                                   (pos_atoms[it->first][2] - it->second.p (2, 7-ii))) / volume;

        if (! it->second.is_hanging (ii))
          (*rho_fixed)[it->second.gt (ii)] += charge_atoms[it->first]*4.0*pi*weigth / (*vol_patch)[it->second.gt (ii)];
        else
          for (int jj = 0; jj < it->second.num_parents (ii); ++jj) {
            double denom = it->second.num_parents (ii) * (*vol_patch)[it->second.gparent (jj, ii)];
            (*rho_fixed)[it->second.gparent (jj, ii)] += charge_atoms[it->first]*4.0*pi*weigth / denom;
          }

      }
    }
  }

  vol_patch.reset ();

  MPI_Barrier (mpicomm);
  auto end_rho = std::chrono::steady_clock::now ();

  if (rank==0) {
    std::cout << "\nTime to calculate rho:  "
              << std::chrono::duration_cast<std::chrono::milliseconds> (end_rho- start_rho).count ()
              << " ms"
              <<std::endl;
  }

  //////////////////////////////////////////////////////////////////
  auto func_frac = [&] (tmesh_3d::quadrant_iterator& quadrant) {
    return cube_fraction_intersection (quadrant,ray_cache);
  };

  std::unique_ptr<distributed_sparse_matrix> A =
    std::make_unique<distributed_sparse_matrix> (mpicomm);
  A->set_ranges (tmsh.num_owned_nodes ());
  std::unique_ptr<distributed_vector> rhs =
    std::make_unique<distributed_vector> (tmsh.num_owned_nodes (), mpicomm);

  bim3a_laplacian_frac (tmsh, *epsilon_nodes, *A, func_frac);


  if (stern_layer_surf == 1) { //se c'è lo stern leyer
    bim3a_reaction (tmsh, reaction, *ones, *A);
  } else {
    bim3a_reaction_frac (tmsh, (*reaction_nodes), *ones, *A, func_frac);
  }

  if (size >1) {
    bim3a_solution_with_ghosts (tmsh, *rho_fixed);
  }

  bim3a_rhs (tmsh, const_ones, *rho_fixed, *rhs);

  rho_fixed.reset ();
  reaction_nodes.reset ();
  ones.reset ();
  std::vector<double> ().swap (const_ones);
  std::vector<double> ().swap (reaction);
  std::vector<double> ().swap (marker);

  // Set boundary conditions.
  dirichlet_bcs3 bcs;

  if (bc == 1) { //hom Dir bc
    for (auto const & ibc : bcells) {
      auto cella = ibc.first;
      auto lato = ibc.second;
      bcs.push_back (std::make_tuple (cella, lato,
      [] (double x, double y, double z) {
        return 0;
      }));
    }

    bim3a_dirichlet_bc (tmsh, bcs, *A, *rhs);
  }

  if (bc == 2) { //coulombic Dir bc
    for (auto const & ibc : bcells) {
      auto cella = ibc.first;
      auto lato = ibc.second;
      bcs.push_back (std::make_tuple (cella, lato,
      [&] (double x, double y, double z) {
        return coulomb_boundary_conditions (x,y,z);
      }));
    }

    bim3a_dirichlet_bc (tmsh, bcs, *A, *rhs);
  }

  if (bc == 3) { //analytic Dir bc for a sphere of R = 2AA and q=1
    for (auto const & ibc : bcells) {
      auto cella = ibc.first;
      auto lato = ibc.second;
      bcs.push_back (std::make_tuple (cella, lato,
      [&] (double x, double y, double z) {
        return analytic_solution (x,y,z);
      }));
    }

    bim3a_dirichlet_bc (tmsh, bcs, *A, *rhs);
  }

  if (size > 1) {
    A->assemble ();
    rhs->assemble ();
  }

  mumps mumps_solver;

  std::vector<double> vals;
  std::vector<int> irow, jcol;

  (*A).aij (vals, irow, jcol, mumps_solver.get_index_base ());

  mumps_solver.set_lhs_distributed ();
  mumps_solver.set_distributed_lhs_structure (tmsh.num_global_nodes (), irow, jcol);
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
  if (size > 1)
    bim3a_solution_with_ghosts (tmsh, *phi, replace_op);

  ///////

  mumps_solver.cleanup ();
}



void
poisson_boltzmann::lis_compute_electric_potential (ray_cache_t & ray_cache)
{
  int rank, size;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);

  // if (rank == 0)
  //   std::cout << "\nStarting LIS solution" << std::endl;
    

  // diffusion
  double eps_in = 4.0*pi*e_0*e_in*kb*T*Angs/ (e*e); //adim e_in
  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/ (e*e); //adim e_out

  /////////////////////////////////////////////////////////
  //reactions
  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/ (e_0*e_out*kb*T);

  if (size >1) {
    bim3a_solution_with_ghosts (tmsh, *epsilon_nodes, replace_op);

    if (stern_layer_surf == 0)
      bim3a_solution_with_ghosts (tmsh, (*reaction_nodes), replace_op);
  }


  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  rho_fixed = std::make_unique<distributed_vector> (tmsh.num_owned_nodes (), mpicomm);
  rho_fixed->get_owned_data ().assign (tmsh.num_owned_nodes (), 0.0);

  std::unique_ptr<distributed_vector> ones =
    std::make_unique<distributed_vector> (tmsh.num_owned_nodes (),mpicomm);
  ones->get_owned_data ().assign (tmsh.num_owned_nodes (), 1.0);

  std::vector<double> const_ones (tmsh.num_local_quadrants (), 1.0);

  std::unique_ptr<distributed_vector> vol_patch =
    std::make_unique<distributed_vector> (tmsh.num_owned_nodes (),mpicomm);

  if (size >1)
    bim3a_solution_with_ghosts (tmsh, *ones, replace_op);

  bim3a_rhs (tmsh, const_ones, *ones, *vol_patch);

  if (size >1)
    vol_patch->assemble ();

  MPI_Barrier (mpicomm);
  auto start_rho = std::chrono::steady_clock::now ();


  search_points ();


  for (auto it = lookup_table.begin (); it!=lookup_table.end (); ++it) {
    //linear approx:
    double volume = (it->second.p (0, 7) - it->second.p (0, 0)) *
                    (it->second.p (1, 7) - it->second.p (1, 0)) *
                    (it->second.p (2, 7) - it->second.p (2, 0)); //volume


    {
      for (int ii = 0; ii < 8; ++ii) {
        double weigth = std::abs ( (pos_atoms[it->first][0] - it->second.p (0, 7-ii))*
                                   (pos_atoms[it->first][1] - it->second.p (1, 7-ii))*
                                   (pos_atoms[it->first][2] - it->second.p (2, 7-ii))) / volume;

        if (! it->second.is_hanging (ii))
          (*rho_fixed)[it->second.gt (ii)] += charge_atoms[it->first]*4.0*pi*weigth / (*vol_patch)[it->second.gt (ii)];
        else
          for (int jj = 0; jj < it->second.num_parents (ii); ++jj) {
            double denom = it->second.num_parents (ii) * (*vol_patch)[it->second.gparent (jj, ii)];
            (*rho_fixed)[it->second.gparent (jj, ii)] += charge_atoms[it->first]*4.0*pi*weigth / denom;
          }

      }
    }
  }

  vol_patch.reset ();

  MPI_Barrier (mpicomm);
  auto end_rho = std::chrono::steady_clock::now ();

  if (rank==0) {
    std::cout << "\nTime to calculate rho : "
              << std::chrono::duration_cast<std::chrono::milliseconds> (end_rho- start_rho).count ()
              << " ms"
              <<std::endl;
  }

  //////////////////////////////////////////////////////////////////
  auto func_frac = [&] (tmesh_3d::quadrant_iterator& quadrant) {
    return cube_fraction_intersection (quadrant,ray_cache);
  };

  // MPI_Barrier (mpicomm);
  // auto start_matrix = std::chrono::steady_clock::now ();
  std::unique_ptr<distributed_sparse_matrix> A =
    std::make_unique<distributed_sparse_matrix> (mpicomm);
  A->set_ranges (tmsh.num_owned_nodes ());
  std::unique_ptr<distributed_vector> rhs =
    std::make_unique<distributed_vector> (tmsh.num_owned_nodes (), mpicomm);

  bim3a_laplacian_frac (tmsh, *epsilon_nodes, *A, func_frac);


  if (stern_layer_surf == 1) { //se c'è lo stern leyer
    bim3a_reaction (tmsh, reaction, *ones, *A);
  } else {
    bim3a_reaction_frac (tmsh, (*reaction_nodes), *ones, *A, func_frac);
  }

  if (size >1) {
    bim3a_solution_with_ghosts (tmsh, *rho_fixed);
  }

  bim3a_rhs (tmsh, const_ones, *rho_fixed, *rhs);

  rho_fixed.reset ();
  reaction_nodes.reset ();
  ones.reset ();
  std::vector<double> ().swap (const_ones);
  std::vector<double> ().swap (reaction);
  std::vector<double> ().swap (marker);

  // Set boundary conditions.
  dirichlet_bcs3 bcs;

  if (bc == 1) { //hom Dir bc
    for (auto const & ibc : bcells) {
      auto cella = ibc.first;
      auto lato = ibc.second;
      bcs.push_back (std::make_tuple (cella, lato,
      [] (double x, double y, double z) {
        return 0;
      }));
    }

    bim3a_dirichlet_bc (tmsh, bcs, *A, *rhs);
  }

  if (bc == 2) { //coulombic Dir bc
    for (auto const & ibc : bcells) {
      auto cella = ibc.first;
      auto lato = ibc.second;
      bcs.push_back (std::make_tuple (cella, lato,
      [&] (double x, double y, double z) {
        return coulomb_boundary_conditions (x,y,z);
      }));
    }

    bim3a_dirichlet_bc (tmsh, bcs, *A, *rhs);
  }

  if (bc == 3) { //analytic Dir bc for a sphere of R = 2AA and q=1
    for (auto const & ibc : bcells) {
      auto cella = ibc.first;
      auto lato = ibc.second;
      bcs.push_back (std::make_tuple (cella, lato,
      [&] (double x, double y, double z) {
        return analytic_solution (x,y,z);
      }));
    }

    bim3a_dirichlet_bc (tmsh, bcs, *A, *rhs);
  }

  if (size > 1) {
    A->assemble ();
    rhs->assemble ();
  }

  // std::cout << "\nMatrix and rhs assembled!"<< std::endl;
  
  MPI_Barrier (mpicomm);
  // auto end_matrix = std::chrono::steady_clock::now ();
  // if (rank==0) {
  //   std::cout << "\nTime to calculate Matrix and rhs:  "
  //             << std::chrono::duration_cast<std::chrono::milliseconds> (end_matrix- start_matrix).count ()
  //             << " ms"
  //             <<std::endl<<std::endl;
  // }


  auto start_sol = std::chrono::steady_clock::now ();


  //CSR
  std::vector<double> vals;
  std::vector<int> irow, jcol;


  (*A).csr (vals, jcol, irow);

  // lis RHS
  LIS_INT i, is, ie, n_rhs, ln;
  LIS_VECTOR rhs_lis;
  //n_rhs = tmsh.num_global_nodes();
  ln = rhs->get_owned_data ().size ();

  lis_vector_create (mpicomm, &rhs_lis);
  lis_vector_set_size (rhs_lis, ln, 0);
  lis_vector_get_range (rhs_lis, &is, &ie);

  for (i=is; i<ie; i++)
    lis_vector_set_value (LIS_INS_VALUE, i, rhs->get_owned_data ()[i-is], rhs_lis);

  //cleaning of rhs
  rhs.reset ();

  //lis_vector_print(rhs_lis);
  // lis PHI
  LIS_VECTOR phi_lis;

  lis_vector_create (mpicomm, &phi_lis);
  lis_vector_set_size (phi_lis, ln, 0);
  //lis_vector_set_size(phi_lis, 0, n_rhs);
  lis_vector_get_range (phi_lis, &is, &ie);

  // lis MATRIX
  LIS_INT n, nnz; //n: matrix dim ; nnz: numb of non zero elems
  LIS_INT *index; //array of integer containing the col index of non zero elems
  LIS_INT *ptr; //array of integer with starting points of rows
  LIS_SCALAR *value; //array of double stores non-zero elements of matrix A along the row
  LIS_MATRIX A_lis; //array of integer containing the col index of non zero elems

  nnz = (*A).owned_nnz ();
  n = tmsh.num_owned_nodes ();

  A.reset ();

  lis_matrix_create (mpicomm, &A_lis);

  lis_matrix_set_size (A_lis, n, 0);
  ptr = &irow[0];
  index = &jcol[0];
  value = &vals[0];

  lis_matrix_set_csr (nnz, ptr, index, value, A_lis);


  lis_matrix_assemble (A_lis);

  //Solve linear system
  LIS_SOLVER solver;

  lis_solver_create (&solver);

  std::string opts = linear_solver_options;

  lis_solver_set_option (&opts[0], solver);

  lis_solve (A_lis, rhs_lis, phi_lis, solver);

  lis_solver_destroy (solver);
  lis_vector_destroy (rhs_lis);

  phi = std::make_unique<distributed_vector> (tmsh.num_owned_nodes (), mpicomm);

  lis_vector_get_values (phi_lis, is, ln, phi->get_owned_data ().data ());

  lis_vector_destroy (phi_lis);


  /*
  distributed_vector phi (tmsh.num_owned_nodes (), mpicomm);
  {
    lis_distributed solver;
    A.csr (vals, jcol, irow, solver.get_index_base ());
    solver.set_lhs_structure (tmsh.num_global_nodes (), irow, jcol);
    solver.set_lhs_data (vals);
    solver.set_rhs (rhs.get_owned_data ());
    solver.set_initial_guess (phi.get_owned_data ());
    solver.set_iterative_method ("Conjugate Gradient");
    solver.set_preconditioner ("ilu");
    solver.set_convergence_condition ("nrm1_b");
    solver.analyze ();
    solver.factorize ();
    solver.solve ();
    solver.cleanup ();
  }
  */
  MPI_Barrier (mpicomm);

  auto end_sol = std::chrono::steady_clock::now ();

  if (rank==0) {
    std::cout << "\nTime to solve linear problem:  "
              << std::chrono::duration_cast<std::chrono::milliseconds> (end_sol- start_sol).count ()
              << " ms"
              <<std::endl;
  }

  if (size > 1)
    bim3a_solution_with_ghosts (tmsh, *phi, replace_op);

}


void
poisson_boltzmann::write_potential_on_atoms_fast ()
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);

  std::ofstream phi_atoms;

  std::string filename = "phi_on_atoms_";
  std::string extension = ".txt";
  filename += std::to_string (rank);
  filename += extension;
  phi_atoms.open (filename.c_str ());

  double phi_on_atom;
  double phi_hang_nodes = 0.0;

  for (auto it = lookup_table.begin (); it!=lookup_table.end (); ++it) {
    phi_on_atom = 0.0;
    //linear approx:
    double volume = (it->second.p (0, 7) - it->second.p (0, 0)) *
                    (it->second.p (1, 7) - it->second.p (1, 0)) *
                    (it->second.p (2, 7) - it->second.p (2, 0)); //volume



    {
      for (int ii = 0; ii < 8; ++ii) {
        double weigth = std::abs ( (pos_atoms[it->first][0] - it->second.p (0, 7-ii))*
                                   (pos_atoms[it->first][1] - it->second.p (1, 7-ii))*
                                   (pos_atoms[it->first][2] - it->second.p (2, 7-ii))) / volume;

        if (! it->second.is_hanging (ii))
          phi_on_atom += (*phi)[it->second.gt (ii)]*weigth;
        else {
          phi_hang_nodes = 0.0;

          for (int jj = 0; jj < it->second.num_parents (ii); ++jj)
            phi_hang_nodes += (*phi)[it->second.gparent (jj, ii)]/it->second.num_parents (ii);


          phi_on_atom += phi_hang_nodes*weigth;
        }
      }
    }

    phi_atoms << std::fixed << std::setprecision (3)
              << std::setw (8) << pos_atoms[it->first][0]
              << std::setw (8) << pos_atoms[it->first][1]
              << std::setw (8) << pos_atoms[it->first][2]
              << std::fixed << std::setprecision (4)
              << "  " << phi_on_atom << std::endl;
  }

  /*for (auto it = lookup_table.begin(); it!=lookup_table.end(); ++it) {
    phi_on_atom = 0.0;
      //linear approx:
      double volume = (it->second.p (0, 7) - it->second.p (0, 0)) *
                      (it->second.p (1, 7) - it->second.p (1, 0)) *
                      (it->second.p (2, 7) - it->second.p (2, 0)); //volume



      {
        for (int ii = 0; ii < 8; ++ii) {
          double weigth = std::abs ((atoms[it->first].pos[0] - it->second.p (0, 7-ii))*
                                    (atoms[it->first].pos[1] - it->second.p (1, 7-ii))*
                                    (atoms[it->first].pos[2] - it->second.p (2, 7-ii))) / volume;

          if (! it->second.is_hanging (ii))
            phi_on_atom += (*phi)[it->second.gt (ii)]*weigth;
          else {
            phi_hang_nodes = 0.0;
            for (int jj = 0; jj < it->second.num_parents (ii); ++jj)
              phi_hang_nodes += (*phi)[it->second.gparent (jj, ii)]/it->second.num_parents (ii);


            phi_on_atom += phi_hang_nodes*weigth;
          }
        }
      }

      phi_atoms << std::fixed << std::setprecision(3)
                << std::setw(8) << atoms[it->first].pos[0]
                << std::setw(8) << atoms[it->first].pos[1]
                << std::setw(8) << atoms[it->first].pos[2]
                << std::fixed << std::setprecision(4)
                << "  " << phi_on_atom << std::endl;
  }*/

  phi_atoms.close ();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
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

  auto it0 = ray_cache.rays[dir].find (ray);
  auto normali = it0->second.normals;
  auto inters = it0->second.inters;

  frac = 0.5;

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

  for (int ii : {
         0,1,3,2,4,5,7,6
       }) {
    if (! quadrant->is_hanging (ii)) {
      if ( (*epsilon_nodes)[quadrant->gt (ii)] < isolevel) cubeindex |= index;
    } else {
      for (int jj = 0; jj < quadrant->num_parents (ii); ++jj) {
        tmp += (*epsilon_nodes)[quadrant->gparent (jj, ii)] / quadrant->num_parents (ii);
      }

      if (tmp < isolevel) cubeindex |= index;
    }

    tmp = 0;
    index *= 2;
  }

  // Cube is entirely in/out of the surface
  if (edgeTable[cubeindex] == 0)
    return -1;

  return cubeindex;
}

std::tuple<std::array<double,8>, std::array<double,8>, std::vector<double>,std::vector<int> >
poisson_boltzmann::classifyCube_flux (tmesh_3d::quadrant_iterator& quadrant)
{
  std::vector<double> edges {};
  std::vector<int> flux {};
  std::array<double,8> tmp_eps;
  std::array<double,8> tmp_phi;

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

  for (i=0; triTable[cubeindex][i]!=-1; i+=3) {
    // save all the assigned indexes
    triangles[ntriang][0] = triTable[cubeindex][i ];
    triangles[ntriang][1] = triTable[cubeindex][i+1];
    triangles[ntriang][2] = triTable[cubeindex][i+2];
    ntriang++;
  }

  return ntriang;
}

void
poisson_boltzmann::energy (ray_cache_t & ray_cache)
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);

  if (rank == 0)
    std::cout << "\nStarting energy calculation with surface integrals" << std::endl;


  double eps_in = 4.0*pi*e_0*e_in*kb*T*Angs/ (e*e); //adim e_in
  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/ (e*e); //adim e_out

  double net_charge = 0.0;

  for (const NS::Atom& i : atoms) {
    net_charge += i.charge;
  }

  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  double fract;
  std::array<double,3> V;
  std::array<double,3> N;
  std::array<double,3> dist_vert;
  std::array<double,3> h;
  std::array<double,3> area_h;

  std::array<double,8> tmp_eps;
  std::array<double,8> tmp_phi;
  std::vector<double> edg;
  std::vector<int> fl_dir;

  int cubeindex = -1;

  double charge_pol = 0.0;

  double constant_pol = 0.5* (1.0/eps_out - 1.0/eps_in)/ (4.0*pi);
  double constant_react = 1.0/ (8*pi*eps_out);
  double distance = 0.0;
  double product = 0.0;
  double first_int = 0.0;
  double second_int = 0.0;
  double tmp_flux;
  int i1 = 0, i2 = 0;
  double tmp_phi_1 = 0.0, tmp_phi_2 = 0.0,
         tmp_eps_1 = 0.0, tmp_eps_2 = 0.0;




  // flux and polarization energy calculation

  for (auto quadrant = this->tmsh.begin_quadrant_sweep ();
       quadrant != this->tmsh.end_quadrant_sweep ();
       ++quadrant) {

    if (marker[quadrant->get_forest_quad_idx ()] == 0.5) {
      h[0] = quadrant->p (0, 7) - quadrant->p (0, 0);
      h[1] = quadrant->p (1, 7) - quadrant->p (1, 0);
      h[2] = quadrant->p (2, 7) - quadrant->p (2, 0);
      area_h[0] = h[1]*h[2]/h[0] * 0.25;
      area_h[1] = h[0]*h[2]/h[1] * 0.25;
      area_h[2] = h[0]*h[1]/h[2] * 0.25;

      std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux (quadrant);

      for (int ip = 0; ip < edg.size (); ++ip) {

        tmp_flux = 0.0;

        i1 = edge2nodes[2 * edg[ip] ];
        i2 = edge2nodes[2 * edg[ip] + 1];

        normal_intersection (quadrant, ray_cache, edg[ip], N,fract);
        V[0] = quadrant->p (0, i1);
        V[1] = quadrant->p (1, i1);
        V[2] = quadrant->p (2, i1);
        V[edge_axis[edg[ip]]] += fract*h[edge_axis[edg[ip]]];


        tmp_flux = - (tmp_phi[i2] - tmp_phi[i1]) * wha (tmp_eps[i1],tmp_eps[i2], fract)*
                   fl_dir[ip] * area_h[edge_axis[edg[ip]]];
        charge_pol += tmp_flux;


        for (const NS::Atom& i : atoms) {
          distance = std::hypot (i.pos[0]-V[0], i.pos[1]-V[1], i.pos[2]-V[2]);
          first_int += i.charge*tmp_flux/distance;
        }
      }
    }
  }

  double energy_pol = constant_pol*first_int;

  //direct reaction energy
  double energy_react = 0.0;
  int ntriang = 0;
  int edge;
  std::array<std::array<double,3>,3> vert_triangles;
  std::array<std::array<double,3>,3> norms_vert;
  std::array<double,3> phi_sup;

  double area = 0.0;

  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/ (e_0*e_out*kb*T);
  double k = std::sqrt (k2);

  if (calc_energy>=2 && k >1.e-5) {
    // // Open the write file
    /*
    std::ofstream phi_nodes_txt;
    std::ofstream phi_surf_txt;
    FILE* phi_nod_delphi;
    FILE* phi_sup_delphi;
    std::string filename_nodes = "phi_nodes_";
    std::string filename_nodes_delphi = "phi_nodes_delphi_";
    std::string filename_surf = "phi_surf_";
    std::string filename_sup_delphi = "phi_sup_delphi_";
    std::string extension = ".txt";
    filename_nodes += std::to_string (bc);
    filename_nodes += "_";
    filename_surf += std::to_string (bc);
    filename_surf += "_";
    filename_nodes_delphi += std::to_string (bc);
    filename_nodes_delphi += "_";
    filename_sup_delphi += std::to_string (bc);
    filename_sup_delphi += "_";
    filename_nodes += pqrfilename;
    filename_surf += pqrfilename;
    filename_nodes += extension;
    filename_surf += extension;
    filename_nodes_delphi += pqrfilename;
    filename_sup_delphi += pqrfilename;
    filename_nodes_delphi += extension;
    filename_sup_delphi += extension;

    phi_nodes_txt.open (filename_nodes.c_str ());
    phi_surf_txt.open (filename_surf.c_str ());

    phi_sup_delphi = std::fopen ("filename_sup_delphi.txt", "w");
    phi_nod_delphi = std::fopen ("filename_nodes_delphi.txt", "w");
    /* */
    /////////////////////////////////////////////////

    for (auto quadrant = this->tmsh.begin_quadrant_sweep ();
         quadrant != this->tmsh.end_quadrant_sweep ();
         ++quadrant) {
      cubeindex = classifyCube (quadrant, eps_out);

      if (cubeindex > -1) {
        h[0] = quadrant->p (0, 7) - quadrant->p (0, 0);
        h[1] = quadrant->p (1, 7) - quadrant->p (1, 0);
        h[2] = quadrant->p (2, 7) - quadrant->p (2, 0);

        ntriang = getTriangles (cubeindex, triangles);

        for (int ii = 0; ii < ntriang; ++ii) {
          for (int jj = 0; jj < 3; ++jj) {
            tmp_eps_1 = 0.0;
            tmp_eps_2 = 0.0;
            tmp_phi_1 = 0.0;
            tmp_phi_2 = 0.0;
            edge = triangles[ii][jj];
            i1 = edge2nodes[2 * edge ];
            i2 = edge2nodes[2 * edge + 1];
            V[0] = quadrant->p (0, i1);
            V[1] = quadrant->p (1, i1);
            V[2] = quadrant->p (2, i1);

            normal_intersection (quadrant, ray_cache, edge, N,fract);
            V[edge_axis[edge]] += fract*h[edge_axis[edge]];

            vert_triangles[jj] = V;
            norms_vert[jj] = N;

            if (! quadrant->is_hanging (i1)) {
              tmp_phi_1 = (*phi)[quadrant->gt (i1)];
              tmp_eps_1 = (*epsilon_nodes)[quadrant->gt (i1)];
            } else {
              for (int jj = 0; jj < quadrant->num_parents (i1); ++jj) {
                tmp_phi_1 += (*phi)[quadrant->gparent (jj, i1)] / quadrant->num_parents (i1);
                tmp_eps_1 += (*epsilon_nodes)[quadrant->gparent (jj, i1)] / quadrant->num_parents (i1);
              }
            }

            if (! quadrant->is_hanging (i2)) {
              tmp_phi_2 = (*phi)[quadrant->gt (i2)];
              tmp_eps_2 = (*epsilon_nodes)[quadrant->gt (i2)];
            } else {
              for (int jj = 0; jj < quadrant->num_parents (i2); ++jj) {
                tmp_phi_2 += (*phi)[quadrant->gparent (jj, i2)] / quadrant->num_parents (i2);
                tmp_eps_2 += (*epsilon_nodes)[quadrant->gparent (jj, i2)] / quadrant->num_parents (i2);
              }
            }

            phi_sup[jj]= phi0 (tmp_eps_1, tmp_eps_2, tmp_phi_1, tmp_phi_2, fract);

            // // writing potential on surf and nodes

            /*
            phi_nodes_txt << quadrant->p (0, i1) << "  "
                          << quadrant->p (1, i1) << "  "
                          << quadrant->p (2, i1) << "  "
                          << tmp_phi_1 << std::endl;
            phi_nodes_txt << quadrant->p (0, i2) << "  "
                          << quadrant->p (1, i2) << "  "
                          << quadrant->p (2, i2) << "  "
                          << tmp_phi_2 << std::endl;

            phi_surf_txt << V[0] << "  " << V[1] << "  " << V[2] << "  " << phi_sup[jj] << std::endl;

            std::fprintf (phi_nod_delphi,"\nATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%8.4f%8.4f",1,"X","XXX"," ",0,
                          quadrant->p (0, i1),quadrant->p (1, i1),quadrant->p (2, i1),tmp_phi_1,tmp_phi_2);
            std::fprintf (phi_nod_delphi,"\nATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%8.4f%8.4f",1,"X","XXX"," ",0,
                          quadrant->p (0, i2),quadrant->p (1, i2),quadrant->p (2, i2),tmp_phi_1,tmp_phi_2);
            std::fprintf (phi_sup_delphi,"\nATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%8.4f%8.4f",1,"X","XXX"," ",0,
                          V[0],V[1],V[2],phi_sup[jj],0.0);
            /* */
            /////////////////////////////////////////////////
          }

          area = areaTriangle (vert_triangles);
          // area = SphercalAreaTriangle (vert_triangles);

          for (const NS::Atom& i : atoms) {
            for (int kk = 0; kk < 3; ++kk) {
              dist_vert[0] = vert_triangles[kk][0]-i.pos[0];
              dist_vert[1] = vert_triangles[kk][1]-i.pos[1];
              dist_vert[2] = vert_triangles[kk][2]-i.pos[2];
              distance = std::hypot (dist_vert[0], dist_vert[1], dist_vert[2]);
              product = dist_vert[0]*norms_vert[kk][0] +
                        dist_vert[1]*norms_vert[kk][1] +
                        dist_vert[2]*norms_vert[kk][2];
              second_int += i.charge*phi_sup[kk]*product/ (4.0*pi*distance*distance*distance)*area/3;
            }
          }
        }
      }

    }

    /*
    phi_nodes_txt.close ();
    phi_surf_txt.close ();
    fclose (phi_nod_delphi);
    fclose (phi_sup_delphi);
    /* */

    energy_react = 0.5*second_int - first_int*constant_react;
  }

  //coulombic energy
  double coul_energy = 0.0;

  if (calc_energy==3) {
    double den_in = 1.0/ (eps_in);
    int i_atom = 0;
    int j_atom = 0;

    if (rank == 0) {
      for (const NS::Atom& i : atoms) {
        if (std::fabs (i.charge) > 0.0) {
          for (const NS::Atom& j : atoms) {
            if (std::fabs (j.charge) > 0.0) {
              if (j_atom > i_atom) {
                distance = std::hypot ( (i.pos[0] - j.pos[0]),
                                        (i.pos[1] - j.pos[1]),
                                        (i.pos[2] - j.pos[2]));
                coul_energy += i.charge*j.charge/distance;
              }
            }

            j_atom ++;
          }
        }

        i_atom ++;
        j_atom = 0;
      }

      coul_energy = coul_energy*den_in;
    }
  }



  if (rank == 0) {
    MPI_Reduce (MPI_IN_PLACE, &charge_pol, 1, MPI_DOUBLE, MPI_SUM, 0,
                mpicomm);
    MPI_Reduce (MPI_IN_PLACE, &energy_pol, 1, MPI_DOUBLE, MPI_SUM, 0,
                mpicomm);
    MPI_Reduce (MPI_IN_PLACE, &energy_react, 1, MPI_DOUBLE, MPI_SUM, 0,
                mpicomm);
  } else {
    MPI_Reduce (&charge_pol, &charge_pol, 1, MPI_DOUBLE, MPI_SUM, 0,
                mpicomm);
    MPI_Reduce (&energy_pol, &energy_pol, 1, MPI_DOUBLE, MPI_SUM, 0,
                mpicomm);
    MPI_Reduce (&energy_react, &energy_react, 1, MPI_DOUBLE, MPI_SUM, 0,
                mpicomm);
  }

  // Print the result
  if (rank == 0) {
    std::cout << std::endl;
    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "Net charge: "
              << std::setprecision (16)<<net_charge
              << std::endl;
    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << std::endl;


    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "Polarization charge: "
              << std::setprecision (16)<<charge_pol/ (4.0*pi)
              << "  Errore %:" << (charge_pol/ (4.0*pi) - net_charge)/net_charge*100
              << std::endl;
    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << std::endl;

    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "Polarization energy: "
              << std::setprecision (16)<<energy_pol
              << std::endl;
    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << std::endl;


    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "Direct ionic energy: "
              << std::setprecision (16)<<energy_react
              << std::endl;
    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << std::endl;

    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "Coulumbic energy: "
              << std::setprecision (16)<<coul_energy
              << std::endl;
    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << std::endl;

    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "Total energy: "
              << std::setprecision (16)<<energy_pol + energy_react + coul_energy
              << std::endl;
    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;
  }
}

void
poisson_boltzmann::energy_fast (ray_cache_t & ray_cache)
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);

  if (rank == 0)
    std::cout << "\n================ [ Electrostatic Energy ] =================\n";


  double eps_in = 4.0*pi*e_0*e_in*kb*T*Angs/ (e*e); //adim e_in
  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/ (e*e); //adim e_out

  double net_charge = 0.0;

  // Store charged atoms
  std::vector<double> charge_atoms_tmp;
  std::vector<std::array<double,3>> pos_atoms_tmp;

  // Store charged atoms
  for (int ii = 0; ii < charge_atoms.size (); ++ii) {
    if (std::fabs (charge_atoms[ii]) > 1.e-5) {
      net_charge += charge_atoms[ii];
      charge_atoms_tmp.push_back (std::move (charge_atoms[ii]));
      pos_atoms_tmp.push_back (std::move (pos_atoms[ii]));
    }
  }

  std::vector<double> ().swap (charge_atoms);
  std::vector<std::array<double,3>> ().swap (pos_atoms);


  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  double fract;
  std::array<double,3> V;
  std::array<double,3> N;
  std::array<double,3> dist_vert;
  std::array<double,3> h;
  std::array<double,3> area_h;

  std::array<double,8> tmp_eps;
  std::array<double,8> tmp_phi;
  std::vector<double> edg;
  std::vector<int> fl_dir;

  int cubeindex = -1;

  double charge_pol = 0.0;

  double constant_pol = 0.5* (1.0/eps_out - 1.0/eps_in)/ (4.0*pi);
  double constant_react = 1.0/ (8*pi*eps_out);
  double distance = 0.0;
  double product = 0.0;
  double first_int = 0.0;
  double second_int = 0.0;
  double tmp_flux;
  int i1 = 0, i2 = 0;
  double tmp_phi_1 = 0.0, tmp_phi_2 = 0.0,
         tmp_eps_1 = 0.0, tmp_eps_2 = 0.0;


  int ntriang = 0;
  int edge;
  std::array<std::array<double,3>,3> vert_triangles;
  std::array<std::array<double,3>,3> norms_vert;
  std::array<double,3> phi_sup;

  double area = 0.0;

  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/ (e_0*e_out*kb*T);
  double k = std::sqrt (k2);
  double energy_pol = 0.0;
  double energy_react = 0.0;
  double coul_energy = 0.0;

  auto quadrant = this->tmsh.begin_quadrant_sweep ();

  // flux and polarization energy calculation
  if (calc_energy==1) {
    for (const int ii : border_quad) {
      quadrant[ii];
      h[0] = quadrant->p (0, 7) - quadrant->p (0, 0);
      h[1] = quadrant->p (1, 7) - quadrant->p (1, 0);
      h[2] = quadrant->p (2, 7) - quadrant->p (2, 0);
      area_h[0] = h[1]*h[2]/h[0] * 0.25;
      area_h[1] = h[0]*h[2]/h[1] * 0.25;
      area_h[2] = h[0]*h[1]/h[2] * 0.25;

      std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux (quadrant);

      for (int ip = 0; ip < edg.size (); ++ip) {
        tmp_flux = 0.0;
        i1 = edge2nodes[2 * edg[ip] ];
        i2 = edge2nodes[2 * edg[ip] + 1];
        normal_intersection (quadrant, ray_cache, edg[ip], N,fract);
        V[0] = quadrant->p (0, i1);
        V[1] = quadrant->p (1, i1);
        V[2] = quadrant->p (2, i1);
        V[edge_axis[edg[ip]]] += fract*h[edge_axis[edg[ip]]];
        tmp_flux = - (tmp_phi[i2] - tmp_phi[i1]) * wha (tmp_eps[i1],tmp_eps[i2], fract)*
                   fl_dir[ip] * area_h[edge_axis[edg[ip]]];
        charge_pol += tmp_flux;

        for (int ii = 0; ii < charge_atoms_tmp.size (); ++ii) {
          distance = std::hypot (pos_atoms_tmp[ii][0]-V[0], pos_atoms_tmp[ii][1]-V[1], pos_atoms_tmp[ii][2]-V[2]);
          first_int += charge_atoms_tmp[ii]*tmp_flux/distance;
        }

      }
    }

    energy_pol = constant_pol*first_int;
  }

  //direct reaction energy
  if (calc_energy==2) {
    for (const int ii : border_quad) {
      quadrant[ii];
      cubeindex = classifyCube (quadrant, eps_out);

      h[0] = quadrant->p (0, 7) - quadrant->p (0, 0);
      h[1] = quadrant->p (1, 7) - quadrant->p (1, 0);
      h[2] = quadrant->p (2, 7) - quadrant->p (2, 0);
      area_h[0] = h[1]*h[2]/h[0] * 0.25;
      area_h[1] = h[0]*h[2]/h[1] * 0.25;
      area_h[2] = h[0]*h[1]/h[2] * 0.25;

      std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux (quadrant);
      ntriang = getTriangles (cubeindex, triangles);

      for (int ip = 0; ip < edg.size (); ++ip) {

        tmp_flux = 0.0;

        i1 = edge2nodes[2 * edg[ip] ];
        i2 = edge2nodes[2 * edg[ip] + 1];

        normal_intersection (quadrant, ray_cache, edg[ip], N,fract);
        V[0] = quadrant->p (0, i1);
        V[1] = quadrant->p (1, i1);
        V[2] = quadrant->p (2, i1);
        V[edge_axis[edg[ip]]] += fract*h[edge_axis[edg[ip]]];
        tmp_flux = - (tmp_phi[i2] - tmp_phi[i1]) * wha (tmp_eps[i1],tmp_eps[i2], fract)*
                   fl_dir[ip] * area_h[edge_axis[edg[ip]]];
        charge_pol += tmp_flux;

        for (int ii = 0; ii < charge_atoms_tmp.size (); ++ii) {
          distance = std::hypot (pos_atoms_tmp[ii][0]-V[0], pos_atoms_tmp[ii][1]-V[1], pos_atoms_tmp[ii][2]-V[2]);
          first_int += charge_atoms_tmp[ii]*tmp_flux/distance;
        }
      }

      if (k > 1.e-5)
        for (int ii = 0; ii < ntriang; ++ii) {
          for (int jj = 0; jj < 3; ++jj) {
            tmp_eps_1 = 0.0;
            tmp_eps_2 = 0.0;
            tmp_phi_1 = 0.0;
            tmp_phi_2 = 0.0;
            edge = triangles[ii][jj];
            i1 = edge2nodes[2 * edge ];
            i2 = edge2nodes[2 * edge + 1];
            V[0] = quadrant->p (0, i1);
            V[1] = quadrant->p (1, i1);
            V[2] = quadrant->p (2, i1);

            normal_intersection (quadrant, ray_cache, edge, N,fract);
            V[edge_axis[edge]] += fract*h[edge_axis[edge]];

            vert_triangles[jj] = V;
            norms_vert[jj] = N;

            if (! quadrant->is_hanging (i1)) {
              tmp_phi_1 = (*phi)[quadrant->gt (i1)];
              tmp_eps_1 = (*epsilon_nodes)[quadrant->gt (i1)];
            } else {
              for (int jj = 0; jj < quadrant->num_parents (i1); ++jj) {
                tmp_phi_1 += (*phi)[quadrant->gparent (jj, i1)] / quadrant->num_parents (i1);
                tmp_eps_1 += (*epsilon_nodes)[quadrant->gparent (jj, i1)] / quadrant->num_parents (i1);
              }
            }

            if (! quadrant->is_hanging (i2)) {
              tmp_phi_2 = (*phi)[quadrant->gt (i2)];
              tmp_eps_2 = (*epsilon_nodes)[quadrant->gt (i2)];
            } else {
              for (int jj = 0; jj < quadrant->num_parents (i2); ++jj) {
                tmp_phi_2 += (*phi)[quadrant->gparent (jj, i2)] / quadrant->num_parents (i2);
                tmp_eps_2 += (*epsilon_nodes)[quadrant->gparent (jj, i2)] / quadrant->num_parents (i2);
              }
            }

            phi_sup[jj]= phi0 (tmp_eps_1, tmp_eps_2, tmp_phi_1, tmp_phi_2, fract);

          }

          area = areaTriangle (vert_triangles);
          // area = SphercalAreaTriangle (vert_triangles);

          for (int ii = 0; ii < charge_atoms_tmp.size (); ++ii) {
            for (int kk = 0; kk < 3; ++kk) {
              dist_vert[0] = vert_triangles[kk][0]-pos_atoms_tmp[ii][0];
              dist_vert[1] = vert_triangles[kk][1]-pos_atoms_tmp[ii][1];
              dist_vert[2] = vert_triangles[kk][2]-pos_atoms_tmp[ii][2];
              distance = std::hypot (dist_vert[0], dist_vert[1], dist_vert[2]);
              product = dist_vert[0]*norms_vert[kk][0] +
                        dist_vert[1]*norms_vert[kk][1] +
                        dist_vert[2]*norms_vert[kk][2];
              second_int += charge_atoms_tmp[ii]*phi_sup[kk]*product/ (4.0*pi*distance*distance*distance)*area/3;
            }
          }
        }
    }

    energy_pol = constant_pol*first_int;
    energy_react = 0.5*second_int - first_int*constant_react;
  }

  //coulombic energy
  if (calc_energy==3) {
    for (const int ii : border_quad) {
      quadrant[ii];
      cubeindex = classifyCube (quadrant, eps_out);

      h[0] = quadrant->p (0, 7) - quadrant->p (0, 0);
      h[1] = quadrant->p (1, 7) - quadrant->p (1, 0);
      h[2] = quadrant->p (2, 7) - quadrant->p (2, 0);
      area_h[0] = h[1]*h[2]/h[0] * 0.25;
      area_h[1] = h[0]*h[2]/h[1] * 0.25;
      area_h[2] = h[0]*h[1]/h[2] * 0.25;
      std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux (quadrant);
      ntriang = getTriangles (cubeindex, triangles);

      for (int ip = 0; ip < edg.size (); ++ip) {
        tmp_flux = 0.0;
        i1 = edge2nodes[2 * edg[ip] ];
        i2 = edge2nodes[2 * edg[ip] + 1];
        normal_intersection (quadrant, ray_cache, edg[ip], N,fract);
        V[0] = quadrant->p (0, i1);
        V[1] = quadrant->p (1, i1);
        V[2] = quadrant->p (2, i1);
        V[edge_axis[edg[ip]]] += fract*h[edge_axis[edg[ip]]];
        tmp_flux = - (tmp_phi[i2] - tmp_phi[i1]) * wha (tmp_eps[i1],tmp_eps[i2], fract)*
                   fl_dir[ip] * area_h[edge_axis[edg[ip]]];
        charge_pol += tmp_flux;

        for (int ii = 0; ii < charge_atoms_tmp.size (); ++ii) {
          distance = std::hypot (pos_atoms_tmp[ii][0]-V[0], pos_atoms_tmp[ii][1]-V[1], pos_atoms_tmp[ii][2]-V[2]);
          first_int += charge_atoms_tmp[ii]*tmp_flux/distance;
        }
      }

      if (k > 1.e-5)
        for (int ii = 0; ii < ntriang; ++ii) {
          for (int jj = 0; jj < 3; ++jj) {
            tmp_eps_1 = 0.0;
            tmp_eps_2 = 0.0;
            tmp_phi_1 = 0.0;
            tmp_phi_2 = 0.0;
            edge = triangles[ii][jj];
            i1 = edge2nodes[2 * edge ];
            i2 = edge2nodes[2 * edge + 1];
            V[0] = quadrant->p (0, i1);
            V[1] = quadrant->p (1, i1);
            V[2] = quadrant->p (2, i1);

            normal_intersection (quadrant, ray_cache, edge, N,fract);
            V[edge_axis[edge]] += fract*h[edge_axis[edge]];

            vert_triangles[jj] = V;
            norms_vert[jj] = N;

            if (! quadrant->is_hanging (i1)) {
              tmp_phi_1 = (*phi)[quadrant->gt (i1)];
              tmp_eps_1 = (*epsilon_nodes)[quadrant->gt (i1)];
            } else {
              for (int jj = 0; jj < quadrant->num_parents (i1); ++jj) {
                tmp_phi_1 += (*phi)[quadrant->gparent (jj, i1)] / quadrant->num_parents (i1);
                tmp_eps_1 += (*epsilon_nodes)[quadrant->gparent (jj, i1)] / quadrant->num_parents (i1);
              }
            }

            if (! quadrant->is_hanging (i2)) {
              tmp_phi_2 = (*phi)[quadrant->gt (i2)];
              tmp_eps_2 = (*epsilon_nodes)[quadrant->gt (i2)];
            } else {
              for (int jj = 0; jj < quadrant->num_parents (i2); ++jj) {
                tmp_phi_2 += (*phi)[quadrant->gparent (jj, i2)] / quadrant->num_parents (i2);
                tmp_eps_2 += (*epsilon_nodes)[quadrant->gparent (jj, i2)] / quadrant->num_parents (i2);
              }
            }

            phi_sup[jj]= phi0 (tmp_eps_1, tmp_eps_2, tmp_phi_1, tmp_phi_2, fract);
          }

          area = areaTriangle (vert_triangles);
          // area = SphercalAreaTriangle (vert_triangles);

          for (int ii = 0; ii < charge_atoms_tmp.size (); ++ii) {
            for (int kk = 0; kk < 3; ++kk) {
              dist_vert[0] = vert_triangles[kk][0]-pos_atoms_tmp[ii][0];
              dist_vert[1] = vert_triangles[kk][1]-pos_atoms_tmp[ii][1];
              dist_vert[2] = vert_triangles[kk][2]-pos_atoms_tmp[ii][2];
              distance = std::hypot (dist_vert[0], dist_vert[1], dist_vert[2]);
              product = dist_vert[0]*norms_vert[kk][0] +
                        dist_vert[1]*norms_vert[kk][1] +
                        dist_vert[2]*norms_vert[kk][2];
              second_int += charge_atoms_tmp[ii]*phi_sup[kk]*product/ (4.0*pi*distance*distance*distance)*area/3;
            }
          }
        }

    }

    double den_in = 1.0 / eps_in;

    // Calculate energy
    for (size_t i = 0; i < charge_atoms_tmp.size (); ++i) {
      for (size_t j = i + 1; j < charge_atoms_tmp.size (); ++j) {
        distance = std::hypot ( (pos_atoms_tmp[i][0] - pos_atoms_tmp[j][0]),
                                (pos_atoms_tmp[i][1] - pos_atoms_tmp[j][1]),
                                (pos_atoms_tmp[i][2] - pos_atoms_tmp[j][2]));
        // if (distance > 0) {  // Check to avoid division by zero
        coul_energy += (charge_atoms_tmp[i] * charge_atoms_tmp[j]) / distance;
        // }
      }
    }

    coul_energy *= den_in;
    energy_pol = constant_pol*first_int;
    energy_react = 0.5*second_int - first_int*constant_react;
  }

  if (rank == 0) {
    MPI_Reduce (MPI_IN_PLACE, &charge_pol, 1, MPI_DOUBLE, MPI_SUM, 0,
                mpicomm);
    MPI_Reduce (MPI_IN_PLACE, &energy_pol, 1, MPI_DOUBLE, MPI_SUM, 0,
                mpicomm);
    MPI_Reduce (MPI_IN_PLACE, &energy_react, 1, MPI_DOUBLE, MPI_SUM, 0,
                mpicomm);
  } else {
    MPI_Reduce (&charge_pol, &charge_pol, 1, MPI_DOUBLE, MPI_SUM, 0,
                mpicomm);
    MPI_Reduce (&energy_pol, &energy_pol, 1, MPI_DOUBLE, MPI_SUM, 0,
                mpicomm);
    MPI_Reduce (&energy_react, &energy_react, 1, MPI_DOUBLE, MPI_SUM, 0,
                mpicomm);
  }

  if (rank == 0) {
    constexpr int label_width = 50;
    constexpr int precision = 16;

    std::cout << std::left << std::setw(label_width) << "  Net charge [e]:"
              << std::setprecision(precision) << net_charge << "\n";

    std::cout << std::left << std::setw(label_width) << "  Flux charge [e]:"
              << std::setprecision(precision) << charge_pol / (4.0 * pi) << "\n";

    std::cout << std::left << std::setw(label_width)
              << "    Error w.r.t. net charge [%]:"
              << std::setprecision(6)
              << ((charge_pol / (4.0 * pi) - net_charge) / net_charge * 100.0) << "\n";

    std::cout << std::left << std::setw(label_width) << "  Polarization energy [kT]:"
              << std::setprecision(precision) << energy_pol << "\n";

    std::cout << std::left << std::setw(label_width) << "  Direct ionic energy [kT]:"
              << std::setprecision(precision) << energy_react << "\n";

    std::cout << std::left << std::setw(label_width) << "  Coulombic energy [kT]:"
              << std::setprecision(precision) << coul_energy << "\n";

    std::cout << std::left << std::setw(label_width) << "  Sum of electrostatic energy contributions [kT]:"
              << std::setprecision(precision)
              << (energy_pol + energy_react + coul_energy) << "\n";

    std::cout << "===========================================================\n";
  }
}


void
poisson_boltzmann::write_potential_on_surface (ray_cache_t & ray_cache)
{
  double eps_in = 4.0*pi*e_0*e_in*kb*T*Angs/ (e*e); //adim e_in
  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/ (e*e); //adim e_out

  double fract;
  std::array<double,3> V;
  std::array<double,3> N;
  std::array<double,3> dist_vert;
  std::array<double,3> h;
  std::array<double,3> area_h;

  std::array<double,8> tmp_eps;
  std::array<double,8> tmp_phi;
  std::vector<double> edg;
  std::vector<int> fl_dir;

  int cubeindex = -1;

  int i1 = 0, i2 = 0;
  double tmp_phi_1 = 0.0, tmp_phi_2 = 0.0,
         tmp_eps_1 = 0.0, tmp_eps_2 = 0.0;

  int ntriang = 0;
  int edge;
  std::array<std::array<double,3>,3> vert_triangles;
  std::array<std::array<double,3>,3> norms_vert;
  std::array<double,3> phi_sup;

  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/ (e_0*e_out*kb*T);
  double k = std::sqrt (k2);

  auto quadrant = this->tmsh.begin_quadrant_sweep ();

  std::ofstream phi_nodes_txt;
  std::ofstream phi_surf_txt;
  FILE* phi_nod_delphi;
  FILE* phi_sup_delphi;
  std::string filename_nodes = "phi_nodes_";
  std::string filename_nodes_delphi = "phi_nodes_delphi_";
  std::string filename_surf = "phi_surf_";
  std::string filename_sup_delphi = "phi_sup_delphi_";
  std::string extension = ".txt";
  filename_nodes += pqrfilename;
  filename_surf += pqrfilename;
  filename_nodes += extension;
  filename_surf += extension;

  phi_nodes_txt.open (filename_nodes.c_str ());
  phi_surf_txt.open (filename_surf.c_str ());

  phi_sup_delphi = std::fopen ("filename_sup_delphi.txt", "w");
  phi_nod_delphi = std::fopen ("filename_nodes_delphi.txt", "w");

  for (const int ii : border_quad) {
    quadrant[ii];
    cubeindex = classifyCube (quadrant, eps_out);

    h[0] = quadrant->p (0, 7) - quadrant->p (0, 0);
    h[1] = quadrant->p (1, 7) - quadrant->p (1, 0);
    h[2] = quadrant->p (2, 7) - quadrant->p (2, 0);
    area_h[0] = h[1]*h[2]/h[0] * 0.25;
    area_h[1] = h[0]*h[2]/h[1] * 0.25;
    area_h[2] = h[0]*h[1]/h[2] * 0.25;

    std::tie (tmp_phi, tmp_eps, edg, fl_dir) = classifyCube_flux (quadrant);
    ntriang = getTriangles (cubeindex, triangles);

    for (int ii = 0; ii < ntriang; ++ii) {
      for (int jj = 0; jj < 3; ++jj) {
        tmp_eps_1 = 0.0;
        tmp_eps_2 = 0.0;
        tmp_phi_1 = 0.0;
        tmp_phi_2 = 0.0;
        edge = triangles[ii][jj];
        i1 = edge2nodes[2 * edge ];
        i2 = edge2nodes[2 * edge + 1];
        V[0] = quadrant->p (0, i1);
        V[1] = quadrant->p (1, i1);
        V[2] = quadrant->p (2, i1);

        normal_intersection (quadrant, ray_cache, edge, N,fract);
        V[edge_axis[edge]] += fract*h[edge_axis[edge]];

        vert_triangles[jj] = V;
        norms_vert[jj] = N;

        if (! quadrant->is_hanging (i1)) {
          tmp_phi_1 = (*phi)[quadrant->gt (i1)];
          tmp_eps_1 = (*epsilon_nodes)[quadrant->gt (i1)];
        } else {
          for (int jj = 0; jj < quadrant->num_parents (i1); ++jj) {
            tmp_phi_1 += (*phi)[quadrant->gparent (jj, i1)] / quadrant->num_parents (i1);
            tmp_eps_1 += (*epsilon_nodes)[quadrant->gparent (jj, i1)] / quadrant->num_parents (i1);
          }
        }

        if (! quadrant->is_hanging (i2)) {
          tmp_phi_2 = (*phi)[quadrant->gt (i2)];
          tmp_eps_2 = (*epsilon_nodes)[quadrant->gt (i2)];
        } else {
          for (int jj = 0; jj < quadrant->num_parents (i2); ++jj) {
            tmp_phi_2 += (*phi)[quadrant->gparent (jj, i2)] / quadrant->num_parents (i2);
            tmp_eps_2 += (*epsilon_nodes)[quadrant->gparent (jj, i2)] / quadrant->num_parents (i2);
          }
        }

        phi_sup[jj]= phi0 (tmp_eps_1, tmp_eps_2, tmp_phi_1, tmp_phi_2, fract);

        phi_nodes_txt << quadrant->p (0, i1) << "  "
                      << quadrant->p (1, i1) << "  "
                      << quadrant->p (2, i1) << "  "
                      << tmp_phi_1 << std::endl;
        phi_nodes_txt << quadrant->p (0, i2) << "  "
                      << quadrant->p (1, i2) << "  "
                      << quadrant->p (2, i2) << "  "
                      << tmp_phi_2 << std::endl;

        phi_surf_txt << V[0] << "  " << V[1] << "  " << V[2] << "  " << phi_sup[jj] << std::endl;

        std::fprintf (phi_nod_delphi,"\nATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%8.4f%8.4f",1,"X","XXX"," ",0,
                      quadrant->p (0, i1),quadrant->p (1, i1),quadrant->p (2, i1),tmp_phi_1,tmp_phi_2);
        std::fprintf (phi_nod_delphi,"\nATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%8.4f%8.4f",1,"X","XXX"," ",0,
                      quadrant->p (0, i2),quadrant->p (1, i2),quadrant->p (2, i2),tmp_phi_1,tmp_phi_2);
        std::fprintf (phi_sup_delphi,"\nATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f%8.4f%8.4f",1,"X","XXX"," ",0,
                      V[0],V[1],V[2],phi_sup[jj],0.0);


      }
    }
  }

  phi_nodes_txt.close ();
  phi_surf_txt.close ();
  fclose (phi_nod_delphi);
  fclose (phi_sup_delphi);
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
poisson_boltzmann::search_points ()
{
  // size_t count = atoms.size();
  size_t count = charge_atoms.size ();
  auto base =std::make_unique<int[]> (count);

  for (size_t ii = 0; ii < count; ++ii) {
    base[ii] = ii;
  }

  sc_array_t *points = sc_array_new_data (base.get (), sizeof (int), count);


  p8est_search_local (tmsh.p8est, 0, NULL, cerca_atomo_wrapper, points);

}

