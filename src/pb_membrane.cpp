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

#include "pb_membrane.h"
#include "readpdb.h"

#include <fstream>
#include <iostream>
#include <mpi.h>

// ─── Lipid I/O ───────────────────────────────────────────────────────────────

void
read_lipids (poisson_boltzmann& pb)
{
  pb.lipid_atoms.clear ();
  pb.pos_lipid_atoms.clear ();
  pb.charge_lipid_atoms.clear ();
  pb.r_lipid_atoms.clear ();

  if (pb.lipid_filetype == "pqr") {
    std::ifstream f (pb.lipid_file);
    if (!f)
      throw std::runtime_error ("read_lipids: cannot open " + pb.lipid_file);

    NS::Atom a;
    while (f >> a)
      pb.lipid_atoms.push_back (a);

  } else if (pb.lipid_filetype == "pdb") {
    read_pdb (pb.lipid_file, pb.lipid_atoms);
    // TODO: load radii and charges from parameter files if needed

  } else {
    throw std::runtime_error ("read_lipids: unknown filetype '" + pb.lipid_filetype + "'");
  }

  for (const NS::Atom& a : pb.lipid_atoms) {
    pb.pos_lipid_atoms.push_back ({a.pos[0], a.pos[1], a.pos[2]});
    pb.charge_lipid_atoms.push_back (a.charge);
    pb.r_lipid_atoms.push_back (a.radius);
  }

  int rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  if (rank == 0)
    std::cout << "  Lipid atoms read     : " << pb.lipid_atoms.size () << '\n';
}

void
broadcast_lipid_vectors (poisson_boltzmann& pb)
{
  int n = static_cast<int> (pb.charge_lipid_atoms.size ());
  MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  pb.pos_lipid_atoms.resize (n);
  pb.charge_lipid_atoms.resize (n);
  pb.r_lipid_atoms.resize (n);

  MPI_Bcast (pb.pos_lipid_atoms.data (),    n * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (pb.charge_lipid_atoms.data (), n,     MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (pb.r_lipid_atoms.data (),      n,     MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Rebuild lipid_atoms on non-zero ranks from the broadcast flat arrays
  pb.lipid_atoms.resize (n);
  for (int i = 0; i < n; ++i) {
    pb.lipid_atoms[i].pos[0] = pb.pos_lipid_atoms[i][0];
    pb.lipid_atoms[i].pos[1] = pb.pos_lipid_atoms[i][1];
    pb.lipid_atoms[i].pos[2] = pb.pos_lipid_atoms[i][2];
    pb.lipid_atoms[i].charge = pb.charge_lipid_atoms[i];
    pb.lipid_atoms[i].radius = pb.r_lipid_atoms[i];
  }
}

// ─── NanoShaper supercell (PBC workaround) ───────────────────────────────────

std::vector<NS::Atom>
build_ns_supercell (const poisson_boltzmann& pb)
{
  // Start with the protein atoms
  std::vector<NS::Atom> supercell = pb.atoms;

  // [PLACEHOLDER] Replicate lipid atoms across periodic images in xy.
  // The replication strategy (number of images, domain extension passed to NS)
  // is TBD.  For now, append the lipid atoms of the central cell only.
  for (const NS::Atom& a : pb.lipid_atoms)
    supercell.push_back (a);

  return supercell;
}

// ─── Assembly placeholder ────────────────────────────────────────────────────

void
apply_mixing_mass_bc (poisson_boltzmann& /*pb*/, ray_cache_t& /*ray_cache*/)
{
  // [PLACEHOLDER] Mixing mass / periodic BC contribution to the FEM matrix.
  // To be implemented once the method is defined.
}
