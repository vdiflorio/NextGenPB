# scripts

## vtu2grid.py: potential map on a uniform grid (OpenDX or Gaussian cube)

NextGenPB solves on an adaptive octree and writes the potential on the mesh
nodes as VTU files. Many programs need the potential on a uniform grid instead:

- Brownian dynamics (SDA), docking and other tools that read APBS maps (`.dx`);
- viewers such as VMD and PyMOL (`.dx` or `.cube`).

`vtu2grid.py` writes that grid without resampling the solution:

- **Nodal value.** Where a grid point is a mesh node, the script writes the
  nodal value, with no interpolation.
- **Q1 value of the cell.** Elsewhere (far from the molecule, where the octree
  cells are larger) it writes the value of the NextGenPB solution in the cell
  that holds the point. The solution is piecewise trilinear (Q1), so this
  value comes from the 8 vertices of the cell. It is the FEM solution itself,
  not an extra approximation.

The script reports which points are which, and it stops if points of the
uniform region would need interpolation.

Requirements: Python 3 and numpy. The VTU files are read directly, so VTK is
not needed.

### 1. Run NextGenPB with a potential map

In the `.prm` file:

```
[mesh]
mesh_shape  = 0
rand_center = 0      # the mesh centre must stay at the bounding-box midpoint
scale       = 2      # spacing 1/scale = 0.5 A in the uniform region
perfil1     = 0.5    # size of the uniform region (see "Choosing the mesh")
[../]

[model]
map_type      = vtu
potential_map = 1    # writes potential_map_*.vtu and total_potential_map.pvtu
[../]
```

### 2. Convert

Run the script in the folder of the run:

```
python3 vtu2grid.py --pqr mol.pqr -o pot.dx                           # 129^3 points, h = 1 A
python3 vtu2grid.py --pqr mol.pqr --scale 2 --length 128 -o pot.dx    # 257^3 points, h = 0.5 A
python3 vtu2grid.py --pqr mol.pqr --scale 2 --length 80  -o pot.cube  # for VMD / PyMOL
```

The format comes from the extension of `-o` (`.dx` or `.cube`), or from
`--format`. Without VTU arguments the script reads `potential_map_*.vtu` in
the current folder. A `.pvtu` file, a list of `.vtu` pieces or a glob pattern
can also be given.

### Options

| option | meaning | default |
|---|---|---|
| `vtu ...` | `.pvtu`, `.vtu` pieces or glob pattern | `potential_map_*.vtu` |
| `--pqr FILE` | molecule of the run: the grid centre is the midpoint of the bounding box of the atom centres; the atoms also go in the cube file | |
| `--center X Y Z` | grid centre [A], instead of the PQR centre | |
| `--scale S` | grid points per A, spacing h = 1/S | 1 |
| `--length L` | box side [A]; points per side = L·S + 1 (L·S must be an integer) | 128 |
| `--dime NX NY NZ` | points per side, instead of `--length` | |
| `--format dx\|cube` | output format | from `-o` |
| `-o FILE` | output file | `potential.dx` |
| `--hmap FILE` | also write, on the same grid, the edge [A] of the cell used at each point (0 = nodal value) | |
| `--units kT\|kcal` | kT/e, or kcal/(mol e) at temperature `--T` | kT |
| `--T K` | temperature for `--units kcal` (use the T of the run) | 298.15 |
| `--allow-interp` | accept grid points of the uniform region that are not mesh nodes | off |
| `--tol` | geometric tolerance, fraction of h | 1e-6 |

`--pqr` or `--center` is required.

### When the grid points are mesh nodes

In the uniform region the mesh vertices lie at `centre + i / scale_mesh`, and
the grid points lie at `centre + i / scale_grid`. The grid points are
therefore mesh nodes when both conditions hold:

1. **Same centre.** The grid and the mesh have the same centre. This is the
   default with `--pqr`: NextGenPB with `mesh_shape = 0` and `rand_center = 0`
   centres the mesh on the midpoint of the bounding box of the atom centres.
   APBS uses the same point for `gcent mol 1`.
2. **Integer ratio.** `scale_mesh / scale_grid` is an integer. For example, a
   mesh at `scale = 2` works with a grid at `--scale 1` (1 A) or `--scale 2`
   (0.5 A), but not with `--scale 4`.

When both conditions hold, a grid at 1 A is exactly every second point of the
grid at 0.5 A, with identical values. If they do not hold, the script stops:

```
ERROR: N points of the uniform region (cells of 0.5 A) are not mesh nodes: the grid is not aligned with the mesh.
```

`--allow-interp` writes the grid anyway, with Q1 values in those points.

### Output summary

Example for barnase (1BRS, chain A), mesh at `scale 2`, `perfil1 0.8`, grid
129^3 at 1 A:

```
mesh nodes : 233709 / 2146689 points (10.9 %), nodal value, no interpolation
  in cells of  2.000 A : 86812 points (Q1 value of the cell)
  in cells of  4.000 A : 173436 points (Q1 value of the cell)
  in cells of  8.000 A : 363734 points (Q1 value of the cell)
  in cells of 16.000 A : 1190790 points (Q1 value of the cell)
  in cells of 32.000 A : 98208 points (Q1 value of the cell)
max |diff| between cells sharing a point: 5.866e-11 (Q1 continuity)
```

- `mesh nodes` is the fraction of grid points with nodal values.
- The `in cells of` lines tell where the other points come from.
- The last line compares the value obtained from different cells at the same
  point, which also includes the hanging nodes at the interface between cells
  of different size. It should stay at round-off level (~1e-11). A large value
  means a corrupted map.

`--hmap` writes the same information point by point, to view next to the
potential.

### Choosing the mesh

The size of the grid and the size of the uniform region of the mesh are
separate choices:

- The grid box is set by `--length` and is usually much larger than the
  molecule; SDA uses 128 A.
- The uniform region of the mesh is set by `perfil1` in `[mesh]`. Outside it
  the octree cells double in size until they reach the domain boundary.

Grid points beyond the uniform region therefore take their value from coarse
cells. With 50 mM salt the Debye length is about 13.6 A, so cells of 8 to 16 A
are too coarse for an accurate potential there.

Barnase, mesh at `scale 2`, grid 129^3 at 1 A (relative RMS difference from
the `perfil1 0.5` run, outside the 1.5 A ion-exclusion layer):

| `perfil1` | uniform region beyond the atoms | nodes | 10-20 A | 20-40 A | 40-80 A |
|---|---|---|---|---|---|
| 0.8 | about 9 A | 1.6 M | 2.1 % | 6.2 % | 10.7 % |
| 0.5 | about 25 A | 6.3 M | reference | reference | reference |

Lower `perfil1` values extend the uniform region. A value small enough to
cover the whole grid box (about 0.3 for barnase) gives nodal values at every
point. At `scale 2` that is about 26 M cells.

### Units and file formats

- **Units.** `phi` is written in kT/e, the same units as APBS. `--units kcal`
  multiplies by R·T.
- **`.dx`.** APBS layout: header with `origin` and `delta` in A, values with z
  fastest, 3 per line.
- **`.cube`.** Gaussian layout:
  - lengths in Bohr;
  - one line per atom of `--pqr`, with the atomic number guessed from the atom
    name; ions are recognised when the residue name equals the atom name, so
    `CA` in residue `CA` is calcium;
  - values with z fastest, 6 per line, each z-row on new lines.

  Without `--pqr` a single dummy atom is written.

### Errors

| message | cause | fix |
|---|---|---|
| `grid points are outside the NextGenPB domain` | the grid box is larger than the computational domain | lower `perfil2` in `[mesh]`, or use a smaller `--length` |
| `points of the uniform region ... are not mesh nodes` | grid not aligned with the mesh | use the PQR centre and an integer `scale_mesh / scale_grid`; check `rand_center = 0` |
| `cells are not 8-node voxels` / `no raw appended data` | the file is not a NextGenPB VTU | use the `potential_map_*.vtu` files written by NextGenPB |
