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
void read_lipids (poisson_boltzmann& pb);

/**
 * @brief Broadcast lipid atom data from rank 0 to all MPI ranks.
 *
 * Mirrors the behaviour of `poisson_boltzmann::broadcast_vectors()` for the
 * protein atoms, but operates on the lipid atom vectors.
 *
 * @param pb Solver instance (lipid vectors populated on rank 0).
 */
void broadcast_lipid_vectors (poisson_boltzmann& pb);

// ─── NanoShaper supercell (PBC workaround) ───────────────────────────────────

/**
 * @brief Build a NanoShaper supercell by combining protein atoms with replicated lipids.
 *
 * Since NanoShaper does not natively support periodic boundary conditions,
 * this function replicates `pb.lipid_atoms` across the neighbouring periodic
 * images in xy and returns a combined atom list (protein + replicated lipids)
 * suitable for passing to `ray_cache_t::init_analytical_surf_ns()`.
 *
 * The replication strategy (number of images, domain extension) is
 * @b [PLACEHOLDER] — to be defined.
 *
 * @param pb Solver instance with protein and lipid atoms populated.
 * @return Combined atom list for the NanoShaper supercell.
 */
std::vector<NS::Atom>
build_ns_supercell (const poisson_boltzmann& pb);

// ─── Assembly placeholder ────────────────────────────────────────────────────

/**
 * @brief [PLACEHOLDER] Apply periodic BC / mixing mass contribution to the system matrix.
 *
 * This function will implement the "mixing mass" method for periodic boundary
 * conditions once the method is defined.  Currently a no-op stub.
 *
 * @param pb        Solver instance.
 * @param ray_cache Ray cache (may be needed for surface queries at the periodic faces).
 */
void apply_mixing_mass_bc (poisson_boltzmann& pb, ray_cache_t& ray_cache);

#endif // PB_MEMBRANE_H
