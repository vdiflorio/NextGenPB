# Membrane module — NextGenPB

This directory contains the input template for protein-in-membrane Poisson–Boltzmann calculations.

---

## What the membrane module does

NextGenPB supports atomistic lipid bilayer membranes by combining:

1. **Unified molecular surface** — protein and lipid atoms are passed together to NanoShaper so that a single, closed SES/SAS surface is computed across the protein–membrane interface.

2. **Slab-layered mesh (`MESH_SHAPE_MEM`)** — a three-zone octree mesh adapted to the membrane geometry:
   - `maxlevel` near the molecular surface (local refinement via `loc_refinement`)
   - `maxlevel − nlev_mem` in the membrane slab (covers the lipid bilayer and any protein extending above/below)
   - `maxlevel − nlev_sol` in the far solvent

3. **Periodic boundary conditions in xy** — the electrostatic potential satisfies periodic BCs on the lateral faces of the computational domain. The membrane patch extends to the domain boundary; NanoShaper's bounding box is extended by `3 × probe_radius` in xy so the surface can close outside the domain without artifacts.

4. **Boundary residue charge zeroing** — after NanoShaper runs, lipid atoms outside the square computational domain are removed, and the charges of all atoms belonging to membrane residues that straddle the xy domain boundary are set to zero. This prevents near-boundary singularities and ensures consistent postprocessing.

---

## New parameters

All parameters are set in the input `.prm` file.

### `[membrane]` section

| Parameter | Default | Description |
|-----------|---------|-------------|
| `enabled` | `0` | Set to `1` to activate membrane mode |
| `lipid_file` | `lipids.pqr` | Path to the lipid atom file (PQR or PDB) |
| `lipid_filetype` | `pqr` | Format of the lipid file: `pqr` or `pdb` |
| `periodic_x` | `0` | Apply periodic BC on the potential in x |
| `periodic_y` | `0` | Apply periodic BC on the potential in y |
| `membrane_dielectric` | `2.0` | Dielectric constant of the lipid bilayer. Default = protein dielectric (spatially varying membrane ε not yet implemented) |
| `stern_membrane` | `0` | Enable Stern (ion-exclusion) layer on membrane surface |
| `stern_membrane_d` | `0.0` | Stern layer thickness on the membrane [Å] |

### `[mesh]` section (membrane-specific)

When membrane mode is enabled, `mesh_shape` is **always forced to `MESH_SHAPE_MEM = 6`** regardless of any user value. The slab mesh is controlled by:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `scale_max` | `2.0` | Target finest resolution [cells/Å]. Actual scale = `2^maxlevel / box_side`, where `maxlevel = ceil(log2(box_side × scale_max))` |
| `nlev_mem` | `2` | Level difference between `maxlevel` and the slab level. Slab cells are `2^nlev_mem` times coarser than surface cells |
| `nlev_sol` | `4` | Level difference between `maxlevel` and the far-solvent level |

The **box side** is derived automatically from the lipid atom positions plus VdW radii (`min(lip_lx, lip_ly)`); `cell_length_x`/`cell_length_y` are not used.

The **slab z-bounds** cover both the lipid bilayer and the protein (the protein may extend above/below the membrane):
```
z_slab_bot = min(z_prot_min, z_lip_min) − max_radius − 2 × probe_radius
z_slab_top = max(z_prot_max, z_lip_max) + max_radius + 2 × probe_radius
```

---

## Level hierarchy

```
z_box_top  ─────────────────  outlevel  = maxlevel − nlev_sol   (far solvent)
z_slab_top ─────────────────  scale_level = maxlevel − nlev_mem  (membrane slab)
           [protein + membrane] → maxlevel via loc_refinement
z_slab_bot ─────────────────  scale_level
z_box_bot  ─────────────────  outlevel                           (far solvent)
     ←──── box_side = min(lip_lx, lip_ly) ────→  (membrane fills the xy boundary)
```

Constraint: `nlev_sol > nlev_mem > 0` is required for a valid level hierarchy (`outlevel < scale_level < maxlevel`).

---

## File formats

The lipid file follows the same conventions as the protein file:
- **PQR**: each `ATOM`/`HETATM` record contains position, charge, and radius in a single file.
- **PDB**: standard coordinate file; radii and charges are loaded from the same parameter files as the protein (`radius_file`, `charge_file` in `[input]`).

---

## Example usage

```bash
./ngpb --prmfile data/membrane/example.prm
```

See `example.prm` in this directory for a fully annotated template.
