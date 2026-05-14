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
 * mortar_row_offset      : first LOCAL row index of this pair's DOFs in A
 *                          (= N_own + cumulative_ndofm, used for identity/eps rows)
 * mortar_col_offset      : first GLOBAL column index for this pair's mortar DOFs
 *                          (= N_global + cumulative_ndofm, unambiguous across ranks)
 * mortar_dense_offset    : first row in the dense mortar_C array for this pair
 *                          (= cumulative_ndofm, i.e. mortar_col_offset - N_global)
 */
struct MortarFacePair {
  int    face_left;
  int    face_right;
  int    d0, d1;
  double d0_min, L0;
  double d1_min, L1;
  size_t mortar_row_offset;    ///< local A row: N_own + cumulative
  size_t mortar_col_offset;    ///< global col:  N_global + cumulative
  size_t mortar_dense_offset;  ///< row in mortar_C dense array: cumulative
};

/**
 * @brief Accumulate mortar C / C^T contributions for one quadrant on one face.
 *
 * C block  (physical rows, mortar columns): written into A using the globally
 * unambiguous column index pair.mortar_col_offset + k.
 *
 * C^T block (mortar rows, physical columns): written into the dense array
 * mortar_C[row * N_global + physical_col], where
 * row = pair.mortar_dense_offset + k.
 *
 * @param A            System matrix (physical block only; no mortar rows needed).
 * @param mortar_C     Dense C^T accumulator, size n_pairs*ndofm * N_global.
 * @param N_global     Total number of physical DOFs across all ranks.
 * @param q            Iterator pointing to the current quadrant.
 * @param face_idx     Face index (0..3) being assembled.
 * @param pair         Pair descriptor (mortar_col_offset and mortar_dense_offset used).
 * @param sign         +1.0 for left face, -1.0 for right face.
 * @param nm           2^minlevel  (mortar cells per direction).
 */
void assemble_mortar_block (distributed_sparse_matrix        &A,
                            std::vector<double>              &mortar_C,
                            size_t                            N_global,
                            tmesh_3d::quadrant_iterator       q,
                            int                               face_idx,
                            const MortarFacePair             &pair,
                            double                            sign,
                            int                               nm,
                            bool                              fill_mortar_rows = true);

#endif // PB_MEMBRANE_H
