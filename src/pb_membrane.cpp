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
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <utility>
#include <mpi.h>

// ─── Lipid I/O ───────────────────────────────────────────────────────────────

void
read_lipids (poisson_boltzmann &pb)
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
    load_radii (pb.radiusfilename, pb.lipid_atoms);
    load_charges (pb.chargefilename, pb.lipid_atoms);

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
broadcast_lipid_vectors (poisson_boltzmann &pb)
{
  int n = static_cast<int> (pb.charge_lipid_atoms.size ());
  MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  pb.pos_lipid_atoms.resize (n);
  pb.charge_lipid_atoms.resize (n);
  pb.r_lipid_atoms.resize (n);

  MPI_Bcast (pb.pos_lipid_atoms.data (), n * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (pb.charge_lipid_atoms.data (), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (pb.r_lipid_atoms.data (), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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
build_ns_supercell (const poisson_boltzmann &pb)
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
trim_lipid_atoms (poisson_boltzmann &pb)
{
  // Remove lipid atoms that lie outside the computational domain in xy.
  // These atoms were included in the NanoShaper supercell so that the
  // membrane surface closes at the periodic boundary, but they must not
  // contribute to the FEM assembly.
  const double x_lo = pb.l_cr[0], x_hi = pb.r_cr[0];
  const double y_lo = pb.l_cr[1], y_hi = pb.r_cr[1];

  std::vector<NS::Atom> kept;
  kept.reserve (pb.lipid_atoms.size ());
  std::size_t out_x = 0, out_y = 0, out_xy = 0;

  for (const NS::Atom& a : pb.lipid_atoms) {
    bool in_x = (a.pos[0] >= x_lo && a.pos[0] <= x_hi);
    bool in_y = (a.pos[1] >= y_lo && a.pos[1] <= y_hi);

    if (in_x && in_y) {
      kept.push_back (a);
    } else {
      if (!in_x && !in_y) ++out_xy;
      else if (!in_x) ++out_x;
      else ++out_y;
    }
  }

  int rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    std::cout << "\n=== [ Lipid atom xy-trim debug ] ===\n";
    std::cout << "  FEM domain  x : [" << x_lo << ", " << x_hi << "] Å\n";
    std::cout << "  FEM domain  y : [" << y_lo << ", " << y_hi << "] Å\n";
    std::cout << "  Lipid atoms before trim : " << pb.lipid_atoms.size () << '\n';
    std::cout << "  Removed (outside x only): " << out_x << '\n';
    std::cout << "  Removed (outside y only): " << out_y << '\n';
    std::cout << "  Removed (outside x & y) : " << out_xy << '\n';
    std::cout << "  Lipid atoms after trim  : " << kept.size () << '\n';
    std::cout << "====================================\n\n";
  }

  pb.lipid_atoms = std::move (kept);

  // Rebuild flat vectors to stay in sync with lipid_atoms
  pb.pos_lipid_atoms.clear ();
  pb.charge_lipid_atoms.clear ();
  pb.r_lipid_atoms.clear ();
  pb.pos_lipid_atoms.reserve (pb.lipid_atoms.size ());
  pb.charge_lipid_atoms.reserve (pb.lipid_atoms.size ());
  pb.r_lipid_atoms.reserve (pb.lipid_atoms.size ());

  for (const NS::Atom& a : pb.lipid_atoms) {
    pb.pos_lipid_atoms.push_back ({a.pos[0], a.pos[1], a.pos[2]});
    pb.charge_lipid_atoms.push_back (a.charge);
    pb.r_lipid_atoms.push_back (a.radius);
  }
}

void
zero_boundary_residue_charges (poisson_boltzmann &pb)
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

  // Collect residue IDs (chain string + resNum) that straddle the boundary.
  // Also store (resName, chain, resNum) for the diagnostic print.
  using ResID = std::pair<std::string, int>;
  std::set<ResID> boundary_res;
  // Map ResID → resName for the diagnostic listing
  std::map<ResID, std::string> res_name_map;

  for (const NS::Atom& a : pb.lipid_atoms) {
    if (a.pos[0] < x_lo + cutoff || a.pos[0] > x_hi - cutoff ||
        a.pos[1] < y_lo + cutoff || a.pos[1] > y_hi - cutoff) {
      ResID rid {a.ai.chain, a.ai.resNum};
      boundary_res.insert (rid);
      res_name_map.emplace (rid, a.ai.resName);
    }
  }

  // Zero charges of all atoms in those residues
  std::size_t n_zeroed = 0;

  for (std::size_t i = 0; i < pb.lipid_atoms.size (); ++i) {
    NS::Atom &a = pb.lipid_atoms[i];

    if (boundary_res.count ({a.ai.chain, a.ai.resNum})) {
      a.charge = 0.0;
      pb.charge_lipid_atoms[i] = 0.0;
      ++n_zeroed;
    }
  }

  int rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    std::cout << "\n=== [ Boundary residue charge-zeroing debug ] ===\n";
    std::cout << "  Cutoff (max_r + probe_r): " << cutoff << " Å\n";
    std::cout << "  FEM domain  x : [" << x_lo << ", " << x_hi << "] Å\n";
    std::cout << "  FEM domain  y : [" << y_lo << ", " << y_hi << "] Å\n";
    std::cout << "  Zeroed residues (" << boundary_res.size () << "):\n";

    for (const auto& rid : boundary_res) {
      const std::string &rname = res_name_map.at (rid);
      std::cout << "    " << rname
                << "  chain=" << (rid.first.empty () ? "-" : rid.first)
                << "  resNum=" << rid.second << '\n';
    }

    std::cout << "  Total atoms zeroed: " << n_zeroed << '\n';
    std::cout << "=================================================\n\n";
  }
}

// ─── Merge lipid charges into the solute charge set ──────────────────────────

void
merge_lipid_charges_into_solute (poisson_boltzmann &pb)
{
  // Inject the lipid charges into the main solute arrays (pos_atoms /
  // charge_atoms) so that they act as field sources both for the potential
  // (density map / RHS) and for the energy / flux postprocessing.
  //
  // Must be called AFTER trim_lipid_atoms() and zero_boundary_residue_charges():
  // out-of-domain lipids are already removed and boundary residues already
  // zeroed, so only the surviving charged lipid atoms are appended.  Safe on
  // all MPI ranks: the lipid arrays are identical across ranks (broadcast +
  // deterministic trim/zeroing), and pos_atoms/charge_atoms are already present
  // on every rank at this point.
  if (pb.lipid_atoms.empty ()) return;

  std::size_t n_added = 0;

  for (std::size_t i = 0; i < pb.charge_lipid_atoms.size (); ++i) {
    if (std::fabs (pb.charge_lipid_atoms[i]) > 1.e-5) {
      pb.pos_atoms.push_back (pb.pos_lipid_atoms[i]);
      pb.charge_atoms.push_back (pb.charge_lipid_atoms[i]);

      // index_atoms is only populated when atoms_write == 1, and is indexed by
      // pos_atoms position when writing phi on atoms: it has to grow with
      // pos_atoms or the lipid entries read past its end.
      if (pb.atoms_write == 1)
        pb.index_atoms.push_back (static_cast<int> (pb.index_atoms.size ()) + 1);

      ++n_added;
    }
  }

  // net_charge was computed in create_mesh() from the protein only; refresh it
  // so postprocessing compares the Gauss flux against the full source charge.
  pb.net_charge = std::accumulate (pb.charge_atoms.begin (),
                                   pb.charge_atoms.end (), 0.0);

  int rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    std::cout << "\n=== [ Lipid charge merge ] ===\n";
    std::cout << "  Charged lipid atoms added to source : " << n_added << '\n';
    std::cout << "  Net charge (protein + lipids) [e]   : " << pb.net_charge << '\n';
    std::cout << "==============================\n\n";
  }
}

// ─── Membrane slab mesh (MESH_SHAPE_MEM) ─────────────────────────────────────

void
build_membrane_slab_mesh (poisson_boltzmann &pb)
{
  int rank;
  MPI_Comm_rank (pb.mpicomm, &rank);

  // ── Membrane slab mesh ────────────────────────────────────────────────
  // Cubic box whose side is the smaller of the two lipid patch extents in xy.
  // The membrane fills the xy face; the slab in z covers both the membrane
  // and the protein (which may extend above/below the membrane).
  //
  // Level structure:
  // outlevel  = maxlevel - nlev_sol   (far solvent)
  // scale_level = maxlevel - nlev_mem  (membrane slab)
  // maxlevel                           (near molecular surface, loc_refinement)

  auto comp_pos_x = [] (const std::array<double, 3>& a1, const std::array<double, 3>& a2) -> bool {
    return a1[0] < a2[0];
  };
  auto comp_pos_y = [] (const std::array<double, 3>& a1, const std::array<double, 3>& a2) -> bool {
    return a1[1] < a2[1];
  };
  auto comp_pos_z = [] (const std::array<double, 3>& a1, const std::array<double, 3>& a2) -> bool {
    return a1[2] < a2[2];
  };

  // --- protein extents: recomputed here rather than reusing create_mesh()'s
  //     preamble locals, which are not visible from this module ---
  auto minmax_x = std::minmax_element (pb.pos_atoms.begin (), pb.pos_atoms.end (), comp_pos_x);
  auto minmax_y = std::minmax_element (pb.pos_atoms.begin (), pb.pos_atoms.end (), comp_pos_y);
  auto minmax_z = std::minmax_element (pb.pos_atoms.begin (), pb.pos_atoms.end (), comp_pos_z);
  double maxradius = *std::max_element (pb.r_atoms.begin (), pb.r_atoms.end ());

  // --- lipid atom extents in all three directions ---
  auto minmax_lip_x = std::minmax_element (pb.pos_lipid_atoms.begin (), pb.pos_lipid_atoms.end (), comp_pos_x);
  auto minmax_lip_y = std::minmax_element (pb.pos_lipid_atoms.begin (), pb.pos_lipid_atoms.end (), comp_pos_y);
  auto minmax_lip_z = std::minmax_element (pb.pos_lipid_atoms.begin (), pb.pos_lipid_atoms.end (), comp_pos_z);

  double maxrad_lip = pb.r_lipid_atoms.empty () ? 0.0
                      : *std::max_element (pb.r_lipid_atoms.begin (), pb.r_lipid_atoms.end ());
  double maxrad = std::max (maxradius, maxrad_lip);

  // --- z-slab bounds: covers both membrane and protein (protein may be taller) ---
  double z_prot_min = (*minmax_z.first)[2];
  double z_prot_max = (*minmax_z.second)[2];
  double z_lip_min = (*minmax_lip_z.first)[2];
  double z_lip_max = (*minmax_lip_z.second)[2];

  double z_slab_bot = std::min (z_prot_min, z_lip_min) - maxrad - 2.0 * pb.prb_radius;
  double z_slab_top = std::max (z_prot_max, z_lip_max) + maxrad + 2.0 * pb.prb_radius;
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
  pb.l_prot[0] = (*minmax_x.first)[0] - maxradius - 2 * pb.prb_radius;
  pb.l_prot[1] = (*minmax_y.first)[1] - maxradius - 2 * pb.prb_radius;
  pb.l_prot[2] = (*minmax_z.first)[2] - maxradius - 2 * pb.prb_radius;
  pb.r_prot[0] = (*minmax_x.second)[0] + maxradius + 2 * pb.prb_radius;
  pb.r_prot[1] = (*minmax_y.second)[1] + maxradius + 2 * pb.prb_radius;
  pb.r_prot[2] = (*minmax_z.second)[2] + maxradius + 2 * pb.prb_radius;

  pb.l_mem[0] = x_lip_min - maxrad - 2 * pb.prb_radius;
  pb.l_mem[1] = y_lip_min - maxrad - 2 * pb.prb_radius;
  pb.l_mem[2] = z_lip_min - maxrad - 2 * pb.prb_radius;
  pb.r_mem[0] = x_lip_max + maxrad + 2 * pb.prb_radius;
  pb.r_mem[1] = y_lip_max + maxrad + 2 * pb.prb_radius;
  pb.r_mem[2] = z_lip_max + maxrad + 2 * pb.prb_radius;

  // Dielectric slab for the implicit modes: van der Waals envelope of the lipids,
  // which is what the NanoShaper surface approximates in ns mode. l_mem/r_mem are
  // not usable here: they add maxrad + 2*prb_radius per side because they steer
  // mesh refinement, and would give a slab ~4*prb_radius thicker than the ns one.
  pb.z_mem_bot = z_lip_min - maxrad_lip;
  pb.z_mem_top = z_lip_max + maxrad_lip;

  // Center of the membrane patch in xy (create_mesh()'s protein-only box block
  // is skipped for MESH_SHAPE_MEM, so cc is set here for the first time)
  pb.cc[0] = 0.5 * (x_lip_max + x_lip_min);
  pb.cc[1] = 0.5 * (y_lip_max + y_lip_min);
  pb.cc[2] = z_slab_cc;

  // --- square outer box (cubic: same side in x, y, z) ---
  pb.ll[0] = pb.cc[0] - box_side / 2.0;
  pb.rr[0] = pb.cc[0] + box_side / 2.0;
  pb.ll[1] = pb.cc[1] - box_side / 2.0;
  pb.rr[1] = pb.cc[1] + box_side / 2.0;
  pb.ll[2] = pb.cc[2] - box_side / 2.0;
  pb.rr[2] = pb.cc[2] + box_side / 2.0;

  // --- levels ---
  pb.maxlevel = (int) std::ceil (std::log2 (box_side * pb.scale_max));
  pb.scale = (double) (1 << pb.maxlevel) / box_side;
  pb.unilevel = pb.maxlevel - pb.nlev_sol;
  pb.outlevel = pb.maxlevel - pb.nlev_sol;
  pb.minlevel = pb.outlevel;
  pb.scale_level = pb.maxlevel - pb.nlev_mem;
  pb.scale_level_min_box = pb.maxlevel;

  // --- l_c / r_c: slab bounds for the refinement callback ---
  // xy spans the full box (condition always satisfied → only z filters)
  pb.l_c[0] = pb.ll[0];
  pb.r_c[0] = pb.rr[0];
  pb.l_c[1] = pb.ll[1];
  pb.r_c[1] = pb.rr[1];
  pb.l_c[2] = z_slab_bot;
  pb.r_c[2] = z_slab_top;

  // --- l_cr / r_cr: NS box, aligned to the grid (square in xy) ---
  int nx = (int) ((box_side / 2.0) * pb.scale + 0.5);
  int ny = (int) ((box_side / 2.0) * pb.scale + 0.5);
  int nz = (int) ((z_slab_top - z_slab_cc) * pb.scale + 0.5);
  pb.l_cr[0] = pb.cc[0] - nx / pb.scale;
  pb.r_cr[0] = pb.cc[0] + nx / pb.scale;
  pb.l_cr[1] = pb.cc[1] - ny / pb.scale;
  pb.r_cr[1] = pb.cc[1] + ny / pb.scale;
  pb.l_cr[2] = z_slab_cc - nz / pb.scale;
  pb.r_cr[2] = z_slab_cc + nz / pb.scale;

  if (rank == 0) {
    std::cout << "========== [ Domain Information ] ==========\n";
    std::cout << "  Protein extents (with VdW):\n";
    std::cout << "    x: [" << pb.l_prot[0] << ", " << pb.r_prot[0] << "] Å\n";
    std::cout << "    y: [" << pb.l_prot[1] << ", " << pb.r_prot[1] << "] Å\n";
    std::cout << "    z: [" << pb.l_prot[2] << ", " << pb.r_prot[2] << "] Å\n";
    std::cout << "  Membrane patch extents (with VdW):\n";
    std::cout << "    x: [" << pb.l_mem[0] << ", " << pb.r_mem[0] << "] Å\n";
    std::cout << "    y: [" << pb.l_mem[1] << ", " << pb.r_mem[1] << "] Å\n";
    std::cout << "    z: [" << pb.l_mem[2] << ", " << pb.r_mem[2] << "] Å\n";
    std::cout << "  Lipid patch size: lx = " << lip_lx << " Å   ly = " << lip_ly << " Å\n";
    std::cout << "  Box side (min lx,ly)        : " << box_side << " Å\n";
    std::cout << "  Center (cx, cy, cz)         : "
              << pb.cc[0] << "  " << pb.cc[1] << "  " << pb.cc[2] << " Å\n";
    std::cout << "  Outer box x                 : [" << pb.ll[0] << ", " << pb.rr[0] << "] Å\n";
    std::cout << "  Outer box y                 : [" << pb.ll[1] << ", " << pb.rr[1] << "] Å\n";
    std::cout << "  Outer box z                 : [" << pb.ll[2] << ", " << pb.rr[2] << "] Å\n";
    std::cout << "  Slab z                      : [" << z_slab_bot << ", " << z_slab_top
              << "] Å  (thickness " << z_slab_top - z_slab_bot << " Å)\n";
    std::cout << "  FEM domain (l_cr / r_cr) x  : [" << pb.l_cr[0] << ", " << pb.r_cr[0] << "] Å\n";
    std::cout << "  FEM domain (l_cr / r_cr) y  : [" << pb.l_cr[1] << ", " << pb.r_cr[1] << "] Å\n";
    std::cout << "  FEM domain (l_cr / r_cr) z  : [" << pb.l_cr[2] << ", " << pb.r_cr[2] << "] Å\n";
    std::cout << "  scale_max = " << pb.scale_max
              << "   effective scale = " << pb.scale << " cells/Å\n";
    std::cout << "  maxlevel = " << pb.maxlevel
              << "   scale_level (slab) = " << pb.scale_level
              << "   outlevel (solvent) = " << pb.unilevel << '\n';

    if (pb.unilevel >= pb.scale_level || pb.scale_level >= pb.maxlevel)
      std::cerr << "  [ERROR] Level hierarchy violated: need outlevel < scale_level < maxlevel\n"
                << "          Check nlev_mem and nlev_sol.\n";
    else
      std::cout << "  Level hierarchy OK: " << pb.unilevel << " < " << pb.scale_level
                << " < " << pb.maxlevel << '\n';

    std::cout << "============================================\n\n";
  }

  pb.simple_conn_num_vertices = 8;
  pb.simple_conn_num_trees = 1;

  pb.simple_conn_p = std::make_unique<double[]> (pb.simple_conn_num_vertices * 3);
  pb.simple_conn_t = std::make_unique<p4est_topidx_t[]> (pb.simple_conn_num_vertices);

  auto tmp_p = {pb.ll[0], pb.ll[1], pb.ll[2],
                pb.rr[0], pb.ll[1], pb.ll[2],
                pb.ll[0], pb.rr[1], pb.ll[2],
                pb.rr[0], pb.rr[1], pb.ll[2],
                pb.ll[0], pb.ll[1], pb.rr[2],
                pb.rr[0], pb.ll[1], pb.rr[2],
                pb.ll[0], pb.rr[1], pb.rr[2],
                pb.rr[0], pb.rr[1], pb.rr[2]
               };
  auto tmp_t = {1, 2, 3, 4, 5, 6, 7, 8};

  std::copy (tmp_p.begin (), tmp_p.end (), pb.simple_conn_p.get ());
  std::copy (tmp_t.begin (), tmp_t.end (), pb.simple_conn_t.get ());

  for (int i = 0; i < 6; i++)
    pb.bcells.push_back (std::make_pair (0, i));
}

void
init_tmesh_mem_two_box (poisson_boltzmann &pb)
{
  for (auto i = 0; i < pb.outlevel; ++i) {
    pb.tmsh.set_refine_marker (poisson_boltzmann::uniform_refinement);
    pb.tmsh.refine (0, 1);
  }

  auto refinement = [&pb]
  (tmesh_3d::quadrant_iterator q) -> int {
    int currentlevel = static_cast<int> (q->the_quadrant->level);
    int retval = 0;
    double x1, y1, z1;

    if (currentlevel >= pb.maxlevel)
      retval = 0;
    else {
      for (int ii = 0; ii < 8; ++ii) {
        if (! q->is_hanging (ii)) {
          x1 = q -> p (0, ii);
          y1 = q -> p (1, ii);
          z1 = q -> p (2, ii);

          if ( (x1 >= pb.l_prot[0]) && (x1 <= pb.r_prot[0])
               && (y1 >= pb.l_prot[1]) && (y1 <= pb.r_prot[1])
               && (z1 >= pb.l_prot[2]) && (z1 <= pb.r_prot[2])) {
            retval = 1;
            break;
          } else if ((x1 >= pb.l_mem[0]) && (x1 <= pb.r_mem[0])
                     && (y1 >= pb.l_mem[1]) && (y1 <= pb.r_mem[1])
                     && (z1 >= pb.l_mem[2]) && (z1 <= pb.r_mem[2]) && (currentlevel <= pb.scale_level)) {
            retval = 1;
            break;
          }
        }
      }
    }

    return (retval);
  };

  for (auto i = 0; i < pb.maxlevel; ++i) {
    pb.tmsh.set_refine_marker (refinement);
    pb.tmsh.refine (0, 1);
  }
}
