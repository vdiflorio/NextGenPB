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
#include <cstdio>
#include <fstream>
#include <stdio.h>

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
  mesh_shape = g2 ( (mesh_options + "mesh_shape").c_str (), 0);
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
  const std::string surf_options = "surface/";
  surf_type_num = g2 ( (surf_options + "surface_type").c_str (), 0);

  if (surf_type_num == 1) surf_type = NS::skin;
  else if (surf_type_num == 0) surf_type = NS::ses;
  // else if (surf_type_num == 2) surf_type = NS::blobby;
  else surf_type = NS::ses;

  surf_param = g2 ( (surf_options + "surface_parameter").c_str (), 1.4);
  stern_layer_surf = g2 ( (surf_options + "stern_layer_surf").c_str (), 0);
  stern_layer = g2 ( (surf_options + "stern_layer_thickness").c_str (), 2.);
  num_threads = g2 ( (surf_options + "number_of_threads").c_str (), 1);

  const std::string alg_options = "algorithm/";
  linear_solver_name = g2 ( (alg_options + "linear_solver").c_str (), "lis");
  linear_solver_options = g2 ( (alg_options + "solver_options").c_str (), "-p ssor -ssor_omega 0.51 -i cgs -tol 1.e-6 -print 2 -conv_cond 2 -tol_w 0");

  const std::string out_options = "output/";
  p4estfilename = g2 ( (out_options + "p4estfilename").c_str (), "poisson_boltzmann_p4est");
  markerfilename = g2 ( (out_options + "markerfilename").c_str (), "poisson_boltzmann_marker_0");
  surffilename = g2 ( (out_options + "surffilename").c_str (), "poisson_boltzmann_surface_0");

  // --- membrane ---
  const std::string mem_options = "membrane/";
  membrane_enabled = g2 ( (mem_options + "enabled").c_str (), 0);
  if (membrane_enabled) {
    lipid_file     = g2 ( (mem_options + "lipid_file").c_str (),    std::string ("lipids.pqr"));
    lipid_filetype = g2 ( (mem_options + "lipid_filetype").c_str (), std::string ("pqr"));
    periodic_x     = g2 ( (mem_options + "periodic_x").c_str (), 0);
    periodic_y     = g2 ( (mem_options + "periodic_y").c_str (), 0);
    cell_length_x  = g2 ( (mem_options + "cell_length_x").c_str (), 0.0);
    cell_length_y  = g2 ( (mem_options + "cell_length_y").c_str (), 0.0);
    e_mem          = g2 ( (mem_options + "membrane_dielectric").c_str (), 2.0);
    kappa_in       = g2 ( (mem_options + "kappa_in").c_str (), 0.0);
    kappa_out      = g2 ( (mem_options + "kappa_out").c_str (), ionic_strength);
    stern_membrane   = g2 ( (mem_options + "stern_membrane").c_str (), 0);
    stern_membrane_d = g2 ( (mem_options + "stern_membrane_d").c_str (), 0.0);

    // When the membrane is enabled and the user has not explicitly chosen a
    // mesh shape (default 0), automatically switch to mesh_shape 4 (adaptive
    // scaling), which is better suited for protein-in-membrane systems.
    // perfil1/perfil2 are already read above (mesh_shape 0 falls in the
    // mesh_shape < 2 block); only scale_min/scale_max need to be read here.
    if (mesh_shape == 0) {
      mesh_shape = 4;
      scale_min = g2 ( (mesh_options + "scale_min").c_str (), 0.5);
      scale_max = g2 ( (mesh_options + "scale_max").c_str (), 2.0);
    }
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

