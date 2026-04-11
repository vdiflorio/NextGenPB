/*
 *  Copyright (C) 2021-2026 Vincenzo Di Florio
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

#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <utility>
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
    read_pdb     (pb.lipid_file,      pb.lipid_atoms);
    load_radii   (pb.radiusfilename,  pb.lipid_atoms);
    load_charges (pb.chargefilename,  pb.lipid_atoms);

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

// ─── NanoShaper supercell ────────────────────────────────────────────────────

std::vector<NS::Atom>
build_ns_supercell (const poisson_boltzmann& pb)
{
  // Start with the protein atoms
  std::vector<NS::Atom> supercell = pb.atoms;

  // Append lipid atoms of the central cell only.
  // Periodicity in xy is handled by applying periodic BCs on the potential,
  // not by replicating the geometry into NanoShaper.
  for (const NS::Atom& a : pb.lipid_atoms)
    supercell.push_back (a);

  return supercell;
}

// ─── Post-NanoShaper lipid trimming ──────────────────────────────────────────

void
trim_lipid_atoms (poisson_boltzmann& pb)
{
  // Remove lipid atoms that lie outside the computational domain in xy.
  // These atoms were included in the NanoShaper supercell so that the
  // membrane surface closes at the periodic boundary, but they must not
  // contribute to the FEM assembly.
  const double x_lo = pb.l_cr[0], x_hi = pb.r_cr[0];
  const double y_lo = pb.l_cr[1], y_hi = pb.r_cr[1];

  std::vector<NS::Atom> kept;
  kept.reserve (pb.lipid_atoms.size ());
  for (const NS::Atom& a : pb.lipid_atoms)
    if (a.pos[0] >= x_lo && a.pos[0] <= x_hi &&
        a.pos[1] >= y_lo && a.pos[1] <= y_hi)
      kept.push_back (a);

  int rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  if (rank == 0)
    std::cout << "  Lipid atoms after xy trim: " << kept.size ()
              << " (removed " << pb.lipid_atoms.size () - kept.size () << ")\n";

  pb.lipid_atoms = std::move (kept);

  // Rebuild flat vectors to stay in sync with lipid_atoms
  pb.pos_lipid_atoms.clear ();
  pb.charge_lipid_atoms.clear ();
  pb.r_lipid_atoms.clear ();
  pb.pos_lipid_atoms.reserve (pb.lipid_atoms.size ());
  pb.charge_lipid_atoms.reserve (pb.lipid_atoms.size ());
  pb.r_lipid_atoms.reserve (pb.lipid_atoms.size ());
  for (const NS::Atom& a : pb.lipid_atoms) {
    pb.pos_lipid_atoms.push_back    ({a.pos[0], a.pos[1], a.pos[2]});
    pb.charge_lipid_atoms.push_back (a.charge);
    pb.r_lipid_atoms.push_back      (a.radius);
  }
}

void
zero_boundary_residue_charges (poisson_boltzmann& pb)
{
  // Zero the charge of every atom belonging to a membrane residue that has
  // at least one atom within (max_radius + probe_radius) of the xy domain
  // boundary.  This prevents singularities near the periodic faces and
  // ensures postprocessing consistency.
  if (pb.lipid_atoms.empty ()) return;

  const double x_lo = pb.l_cr[0], x_hi = pb.r_cr[0];
  const double y_lo = pb.l_cr[1], y_hi = pb.r_cr[1];

  double max_r = 0.0;
  for (const NS::Atom& a : pb.lipid_atoms)
    max_r = std::max (max_r, static_cast<double> (a.radius));
  const double cutoff = max_r + pb.surf_param;

  // Collect residue IDs (chain string + resNum) that straddle the boundary
  using ResID = std::pair<std::string, int>;
  std::set<ResID> boundary_res;
  for (const NS::Atom& a : pb.lipid_atoms)
    if (a.pos[0] < x_lo + cutoff || a.pos[0] > x_hi - cutoff ||
        a.pos[1] < y_lo + cutoff || a.pos[1] > y_hi - cutoff)
      boundary_res.insert ({a.ai.chain, a.ai.resNum});

  // Zero charges of all atoms in those residues
  std::size_t n_zeroed = 0;
  for (std::size_t i = 0; i < pb.lipid_atoms.size (); ++i) {
    NS::Atom& a = pb.lipid_atoms[i];
    if (boundary_res.count ({a.ai.chain, a.ai.resNum})) {
      a.charge                    = 0.0;
      pb.charge_lipid_atoms[i]    = 0.0;
      ++n_zeroed;
    }
  }

  int rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  if (rank == 0)
    std::cout << "  Boundary membrane residues zeroed: " << boundary_res.size ()
              << " residues (" << n_zeroed << " atoms)\n";
}

// ─── Assembly placeholder ────────────────────────────────────────────────────

void
apply_mixing_mass_bc (poisson_boltzmann& /*pb*/, ray_cache_t& /*ray_cache*/)
{
  // [PLACEHOLDER] Mixing mass / periodic BC contribution to the FEM matrix.
  // To be implemented once the method is defined.
}
