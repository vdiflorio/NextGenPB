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
 * @brief Membrane and periodic boundary condition support for NextGenPB.
 *
 * Provides functions for:
 * - Reading and broadcasting lipid atom data from a separate PQR/PDB file.
 * - Building a NanoShaper supercell by replicating lipid atoms in xy (PBC workaround).
 * - [PLACEHOLDER] Periodic BC / mixing mass contribution to the FEM assembly.
 *
 * This module is activated when `membrane_enabled = true` in the solver options.
 * The "mixing mass" PBC method for the assembly is TBD and will be added here.
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

// ─── Mortar assembly ─────────────────────────────────────────────────────────

/**
 * @brief Describes one periodic face pair for mortar assembly.
 *
 * face_left / face_right : p4est face indices of the two periodic faces
 * d0, d1                 : tangential directions (0=x, 1=y, 2=z)
 * d0_min, L0             : physical range in d0 is [d0_min, d0_min + L0]
 * d1_min, L1             : physical range in d1 is [d1_min, d1_min + L1]
 * mortar_offset          : first row/col of this pair's DOFs in the augmented A
 */
struct MortarFacePair {
  int    face_left;
  int    face_right;
  int    d0, d1;
  double d0_min, L0;
  double d1_min, L1;
  size_t mortar_offset;
};

/**
 * @brief Accumulate mortar C / C^T contributions for one quadrant on one face.
 *
 * Implements the per-quadrant mortar integration (bc.md, Punto 2, Steps 1-5):
 * tensor-product 1D hat-function mass matrices are computed and added to A with
 * the appropriate sign (+1 left face, -1 right face).  Nodes on z=0 or z=Lz
 * (Dirichlet) are skipped.
 *
 * @param A            Augmented system matrix (already resized to N_phys + 2*ndofm).
 * @param q            Iterator pointing to the current quadrant.
 * @param face_idx     Face index (0..3) being assembled.
 * @param pair         Pair descriptor with tangential directions and mortar offset.
 * @param sign         +1.0 for left face, -1.0 for right face.
 * @param nm           2^minlevel  (mortar cells per direction).
 * @param Lz           Domain extent in z (used for Dirichlet filter).
 */
void assemble_mortar_block (distributed_sparse_matrix        &A,
                            tmesh_3d::quadrant_iterator       q,
                            int                               face_idx,
                            const MortarFacePair             &pair,
                            double                            sign,
                            int                               nm);

#endif // PB_MEMBRANE_H
