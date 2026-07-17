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

/**
 * @file pb_membrane.h
 * @brief Membrane / lipid support for NextGenPB (new module — clean branch).
 *
 * Provides functions for:
 * - Reading and broadcasting lipid atom data from a separate PQR/PDB file.
 * - Building the NanoShaper atom list (protein + lipid) for surface generation.
 * - Trimming lipid atoms outside the xy domain and zeroing boundary residues.
 * - Merging lipid charges into the solute source.
 *
 * Activated when `membrane_enabled = true` in the solver options. PBC and
 * membrane-voltage (ΔV) are separate concerns, added in later phases
 * (see piano_cleanup.md).
 */

#ifndef PB_MEMBRANE_H
#define PB_MEMBRANE_H

#include "pb_class.h"
#include "raytracer.h"

// ─── Lipid I/O ───────────────────────────────────────────────────────────────

/**
 * @brief Read lipid atoms from the file specified in `pb.lipid_file`.
 *
 * Reads `pb.lipid_file` (PQR or PDB format, controlled by `pb.lipid_filetype`)
 * and populates `pb.lipid_atoms`, `pb.pos_lipid_atoms`, `pb.charge_lipid_atoms`,
 * and `pb.r_lipid_atoms`.  Should be called on rank 0 before broadcast_lipid_vectors().
 *
 * @param pb Solver instance with `lipid_file` and `lipid_filetype` already set.
 */
void read_lipids (poisson_boltzmann &pb);

/**
 * @brief Broadcast lipid atom data from rank 0 to all MPI ranks.
 *
 * Mirrors the behaviour of `poisson_boltzmann::broadcast_vectors()` for the
 * protein atoms, but operates on the lipid atom vectors.
 *
 * @param pb Solver instance (lipid vectors populated on rank 0).
 */
void broadcast_lipid_vectors (poisson_boltzmann &pb);

// ─── NanoShaper supercell ────────────────────────────────────────────────────

/**
 * @brief Build the NanoShaper atom list by combining protein and lipid atoms.
 *
 * Returns a single atom list (protein + lipid central cell) suitable for passing
 * to `ray_cache_t::init_analytical_surf_ns()`.  Periodicity in xy is handled by
 * applying periodic BCs on the potential, not by replicating the geometry.
 * The caller is responsible for extending the NanoShaper bounding box in xy by
 * 3 × probe_radius per side so that the membrane surface closes outside the
 * computational domain (see poisson_boltzmann.cpp).
 *
 * @param pb Solver instance with protein and lipid atoms populated.
 * @return Combined atom list for NanoShaper.
 */
std::vector<NS::Atom>
build_ns_supercell (const poisson_boltzmann &pb);

// ─── Post-NanoShaper lipid trimming ──────────────────────────────────────────

/**
 * @brief Remove lipid atoms that lie outside the computational domain in xy.
 *
 * Must be called after NanoShaper has finished (the atoms were included in the
 * NS supercell to allow surface closure at the periodic boundary).  Rebuilds
 * `pb.pos_lipid_atoms`, `pb.charge_lipid_atoms`, and `pb.r_lipid_atoms` to
 * match the trimmed `pb.lipid_atoms`.  Safe to call on all MPI ranks.
 *
 * @param pb Solver instance with `l_cr`/`r_cr` already set by create_mesh().
 */
void trim_lipid_atoms (poisson_boltzmann &pb);

/**
 * @brief Zero charges of membrane residues near the xy domain boundary.
 *
 * Any residue that has at least one atom within (max_vdw_radius + probe_radius)
 * of the xy domain boundary has all its atomic charges set to zero — both in
 * `pb.lipid_atoms` and in the flat `pb.charge_lipid_atoms` vector.
 * This prevents singularities near the periodic faces and ensures consistent
 * postprocessing.  Call after trim_lipid_atoms().  Safe on all MPI ranks.
 *
 * @param pb Solver instance with `l_cr`/`r_cr` and trimmed lipid_atoms.
 */
void zero_boundary_residue_charges (poisson_boltzmann &pb);

/**
 * @brief Append lipid charges to the solute charge set so they source the field.
 *
 * Appends the surviving charged lipid atoms (`pb.pos_lipid_atoms` /
 * `pb.charge_lipid_atoms`) to `pb.pos_atoms` / `pb.charge_atoms`, so that they
 * contribute to the density map (RHS of the potential) and to the energy /
 * flux postprocessing.  Atoms whose charge has been zeroed by
 * zero_boundary_residue_charges() (boundary residues) are skipped.
 * Recomputes `pb.net_charge` to include the merged lipid charges.
 *
 * Call after trim_lipid_atoms() and zero_boundary_residue_charges(), and before
 * create_density_map().  Safe on all MPI ranks.
 *
 * @param pb Solver instance with trimmed, boundary-zeroed lipid atoms.
 */
void merge_lipid_charges_into_solute (poisson_boltzmann &pb);

// ─── Membrane slab mesh (MESH_SHAPE_MEM) ─────────────────────────────────────

/**
 * @brief Build the membrane slab geometry, level structure and connectivity.
 *
 * Called from poisson_boltzmann::create_mesh() when `mesh_shape ==
 * MESH_SHAPE_MEM`.  Sets up a cubic outer box whose side is derived from the
 * lipid patch extent, a z-slab covering both the membrane and the protein, and
 * the three-tier level structure (far solvent / slab / molecular surface).
 *
 * Populates `l_prot`/`r_prot`, `l_mem`/`r_mem` (refinement boxes read back by
 * init_tmesh_mem_two_box()), `l_c`/`r_c`, `l_cr`/`r_cr`, `cc`, `ll`/`rr`, the
 * level fields (`maxlevel`, `unilevel`, `outlevel`, `minlevel`, `scale_level`,
 * `scale_level_min_box`), `scale`, and the `simple_conn_*` connectivity arrays.
 * The caller then reads the connectivity into `tmsh`.
 *
 * Recomputes the protein/lipid extents from `pb.pos_atoms` / `pb.pos_lipid_atoms`
 * rather than reusing create_mesh()'s preamble locals.
 *
 * @param pb Solver instance with protein and lipid atoms populated.
 */
void build_membrane_slab_mesh (poisson_boltzmann &pb);

/**
 * @brief Refine the membrane mesh using the protein and membrane boxes.
 *
 * Two-box refinement: cells inside the protein box are refined up to `maxlevel`;
 * cells inside the membrane box up to `scale_level`.  Everything else stays at
 * `outlevel` (far solvent).  Call after create_mesh() has read the connectivity —
 * this is the MESH_SHAPE_MEM case of the driver's "Building Grid" dispatch,
 * alongside poisson_boltzmann::init_tmesh() and friends.
 *
 * @param pb Solver instance whose mesh geometry was set by build_membrane_slab_mesh().
 */
void init_tmesh_mem_two_box (poisson_boltzmann &pb);

// NB: PBC mortar assembly (MortarFacePair / assemble_mortar_block) is a separate
// concern (PBC formulation) and is intentionally NOT part of this membrane module
// on the clean branch — it will be reintroduced with the PBC work. See piano_cleanup.md.

#endif // PB_MEMBRANE_H
