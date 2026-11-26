#!/usr/bin/env python3
"""NextGenPB potential map (VTU) -> uniform grid, OpenDX (.dx) or Gaussian cube (.cube).

NextGenPB solves on an adaptive octree and writes phi on the mesh nodes
(model/potential_map = 1, model/map_type = vtu). Programs that read APBS or
DelPhi maps (SDA, Brownian dynamics, docking, VMD, PyMOL) need a uniform grid.

The value written at each grid point is the NextGenPB solution itself, the
piecewise-trilinear (Q1) FEM function on the octree:
  - a grid point that is a mesh node gets the nodal value, with no interpolation.
    This holds in the whole uniform region of the mesh when the grid is centred
    on the mesh centre (default: midpoint of the bounding box of the atom
    centres, as NextGenPB with mesh_shape = 0 and rand_center = 0, and APBS
    gcent mol 1) and scale_mesh / scale_grid is an integer;
  - a grid point inside a coarser cell (far from the molecule) gets the Q1
    value of that cell, from its 8 vertices. Hanging vertices carry the
    average of their parents, i.e. the FEM constraint, so the field is
    continuous across cells of different size.
The script reports how many points are nodes and the cell size used for the
others; --hmap writes that cell size on the same grid (0 = nodal value). If
points of the uniform region are not nodes (grid not aligned with the mesh)
the script stops, unless --allow-interp.

The VTU files are read directly (raw appended binary, as written by
NextGenPB): only numpy is needed.

Examples:
  vtu2grid.py --pqr mol.pqr -o pot.dx                         # 129^3, h = 1 A (SDA default)
  vtu2grid.py --pqr mol.pqr --scale 2 --length 128 -o pot.dx  # 257^3, h = 0.5 A
  vtu2grid.py total_potential_map.pvtu --pqr mol.pqr --scale 2 --length 80 -o pot.cube

Units: phi in kT/e as written by NextGenPB (as APBS); --units kcal writes
kcal/(mol e) at the temperature --T (1 kT = R T). The cube file has lengths in
Bohr and the atoms of --pqr (atomic number guessed from the atom name).
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

R_KCAL = 1.987204259e-3      # gas constant, kcal/(mol K)
BOHR = 0.529177210903        # A
CHUNK = 500000


# --------------------------------------------------------------------------- input

def read_vtu (path):
    """Points, 8-node connectivity and the point array of a NextGenPB VTU (raw appended)."""
    with open (path, "rb") as f:
        data = f.read ()
    mark = data.find (b"<AppendedData")
    if mark < 0:
        sys.exit ("%s: no raw appended data (not a NextGenPB VTU?)" % path)
    head = data[:mark].decode ("ascii", "replace")
    start = data.index (b"_", mark) + 1
    htype = re.search (r'header_type="(\w+)"', head)
    hdt = np.dtype ("<u8") if htype and htype.group (1) == "UInt64" else np.dtype ("<u4")

    def array (name):
        m = re.search (r'<DataArray[^>]*Name="%s"[^>]*>' % name, head)
        if m is None:
            sys.exit ("%s: no DataArray named %s" % (path, name))
        tag = m.group (0)
        dt = {"Float64": "<f8", "Float32": "<f4", "Int32": "<i4", "Int64": "<i8"}[
            re.search (r'type="(\w+)"', tag).group (1)]
        off = start + int (re.search (r'offset="(\d+)"', tag).group (1))
        nbytes = int (np.frombuffer (data, hdt, 1, off)[0])
        return np.frombuffer (data, dt, nbytes // np.dtype (dt).itemsize, off + hdt.itemsize)

    pts = array ("Array").reshape (-1, 3)
    conn = array ("connectivity")
    offs = array ("offsets")
    if not np.all (np.diff (np.concatenate (([0], offs))) == 8):
        sys.exit ("%s: cells are not 8-node voxels" % path)
    m = re.search (r'<PointData>.*?Name="(\w+)"', head, re.S)
    phi = array (m.group (1))
    return pts, conn.reshape (-1, 8), phi


def pieces (args):
    """VTU pieces from .pvtu files, .vtu files or glob patterns."""
    out = []
    for p in args:
        found = sorted (glob.glob (p)) or [p]
        for q in found:
            if q.endswith (".pvtu"):
                d = os.path.dirname (q)
                out += [os.path.join (d, s) for s in re.findall (r'Source="([^"]+)"', open (q).read ())]
            else:
                out.append (q)
    missing = [p for p in out if not os.path.isfile (p)]
    if not out or missing:
        sys.exit ("VTU files not found: %s" % (" ".join (missing) or " ".join (args)))
    return out


ELEMENTS = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "NA": 11, "MG": 12, "P": 15, "S": 16,
            "CL": 17, "K": 19, "CA": 20, "MN": 25, "FE": 26, "CO": 27, "NI": 28, "CU": 29,
            "ZN": 30, "BR": 35, "I": 53}


def atomic_number (name, resname):
    """Element from a PQR atom name; ions are recognised by residue == atom name."""
    n = re.sub (r"[^A-Za-z]", "", name).upper ()
    r = resname.upper ()
    if n == r and n in ELEMENTS:          # NA, CL, ZN, CA (calcium) ... ions
        return ELEMENTS[n]
    if n[:2] in ("CL", "BR") and n[:2] in ELEMENTS:
        return ELEMENTS[n[:2]]
    return ELEMENTS.get (n[:1], 0)


def read_pqr (path):
    """Atom centres, charges and atomic numbers of a PQR file (chain column optional)."""
    xyz, q, z = [], [], []
    for line in open (path):
        if line.startswith (("ATOM", "HETATM")):
            t = line.split ()
            xyz.append ([float (t[-5]), float (t[-4]), float (t[-3])])
            q.append (float (t[-2]))
            z.append (atomic_number (t[2], t[3]))
    if not xyz:
        sys.exit ("%s: no ATOM/HETATM records" % path)
    return np.array (xyz), np.array (q), np.array (z)


# --------------------------------------------------------------------------- sampling

class Sampler:
    """For every grid point: Q1 value, whether it is a mesh node, size of the cell used."""

    def __init__ (self, n, origin, h, tol):
        self.n, self.o, self.h, self.tol = n, origin, h, tol
        size = int (np.prod (n))
        self.val = np.full (size, np.nan)
        # a point is a true node only if it is a vertex of every cell that holds it:
        # a hanging vertex of a fine cell lies inside a face of the coarse neighbour
        self.notnode = np.zeros (size, bool)
        self.hcell = np.zeros (size)            # edge of the largest cell holding the point
        self.maxdup = 0.0

    def add (self, pts, conn, phi):
        n, o, h, tol = self.n, self.o, self.h, self.tol
        X = pts[conn]                                   # (ncell, 8, 3)
        lo, hi = X.min (1), X.max (1)
        ilo = np.maximum (np.ceil ((lo - o) / h - tol), 0).astype (np.int64)
        ihi = np.minimum (np.floor ((hi - o) / h + tol), n - 1).astype (np.int64)
        ext = np.maximum (ihi - ilo + 1, 0)
        cnt = ext.prod (1)
        cells = np.nonzero (cnt)[0]
        for c0 in range (0, len (cells), CHUNK // 8):
            cs = cells[c0:c0 + CHUNK // 8]
            k = cnt[cs]
            cid = np.repeat (cs, k)
            r = np.arange (k.sum ()) - np.repeat (np.cumsum (k) - k, k)
            e = ext[cid]
            ijk = ilo[cid] + np.stack ((r // (e[:, 1] * e[:, 2]), (r // e[:, 2]) % e[:, 1], r % e[:, 2]), 1)
            x = o + ijk * h
            span = hi[cid] - lo[cid]
            t = np.clip ((x - lo[cid]) / span, 0.0, 1.0)            # (m, 3)
            bit = X[cid] > (lo[cid] + 0.5 * span)[:, None, :]       # corner of each vertex
            w = np.where (bit, t[:, None, :], 1.0 - t[:, None, :]).prod (2)   # (m, 8)
            v = (w * phi[conn[cid]]).sum (1)
            isnode = w.max (1) > 1.0 - 1e-9
            idx = np.ravel_multi_index (ijk.T, tuple (n))
            seen = ~np.isnan (self.val[idx])
            if seen.any ():
                self.maxdup = max (self.maxdup, float (np.max (np.abs (self.val[idx[seen]] - v[seen]))))
            # keep the value from the coarsest cell: at a hanging vertex it is the
            # Q1 value of that cell (identical up to round-off, Q1 is continuous)
            hc = span.max (1)
            better = hc > self.hcell[idx] + 1e-12
            np.maximum.at (self.hcell, idx, hc)
            self.val[idx[better | ~seen]] = v[better | ~seen]
            self.notnode[idx[~isnode]] = True


# --------------------------------------------------------------------------- output

def write_dx (path, grid, origin, h, comment):
    nx, ny, nz = grid.shape
    with open (path, "w") as f:
        f.write ("# %s\n" % comment)
        f.write ("object 1 class gridpositions counts %d %d %d\n" % (nx, ny, nz))
        f.write ("origin %.6e %.6e %.6e\n" % tuple (origin))
        f.write ("delta %.6e 0.000000e+00 0.000000e+00\n" % h[0])
        f.write ("delta 0.000000e+00 %.6e 0.000000e+00\n" % h[1])
        f.write ("delta 0.000000e+00 0.000000e+00 %.6e\n" % h[2])
        f.write ("object 2 class gridconnections counts %d %d %d\n" % (nx, ny, nz))
        f.write ("object 3 class array type double rank 0 items %d data follows\n" % grid.size)
        v = grid.ravel (order="C")          # z fastest, as in APBS
        full = (v.size // 3) * 3
        np.savetxt (f, v[:full].reshape (-1, 3), fmt="%.6e")
        if full < v.size:
            f.write (" ".join ("%.6e" % x for x in v[full:]) + "\n")
        f.write ('attribute "dep" string "positions"\n')
        f.write ('object "regular positions regular connections" class field\n')
        f.write ('component "positions" value 1\ncomponent "connections" value 2\n')
        f.write ('component "data" value 3\n')


def write_cube (path, grid, origin, h, comment, atoms):
    """Gaussian cube: lengths in Bohr, x slowest, z fastest, 6 values per line per z-row."""
    nx, ny, nz = grid.shape
    if atoms is None:                                   # cube needs at least one atom
        xyz, z = np.zeros ((1, 3)), np.zeros (1, int)
    else:
        xyz, z = atoms
    o = origin / BOHR
    with open (path, "w") as f:
        f.write ("NextGenPB potential (vtu2grid.py)\n")
        f.write ("%s\n" % comment)
        f.write ("%5d %12.6f %12.6f %12.6f\n" % (len (z), *o))
        f.write ("%5d %12.6f %12.6f %12.6f\n" % (nx, h[0] / BOHR, 0.0, 0.0))
        f.write ("%5d %12.6f %12.6f %12.6f\n" % (ny, 0.0, h[1] / BOHR, 0.0))
        f.write ("%5d %12.6f %12.6f %12.6f\n" % (nz, 0.0, 0.0, h[2] / BOHR))
        p = xyz / BOHR
        for zi, (x, y, w) in zip (z, p):
            f.write ("%5d %12.6f %12.6f %12.6f %12.6f\n" % (zi, float (zi), x, y, w))
        rows = grid.reshape (nx * ny, nz)
        full = (nz // 6) * 6
        for row in rows:
            if full:
                np.savetxt (f, row[:full].reshape (-1, 6), fmt="%13.5E", delimiter="")
            if full < nz:
                f.write ("".join ("%13.5E" % x for x in row[full:]) + "\n")


# --------------------------------------------------------------------------- main

def main ():
    ap = argparse.ArgumentParser (description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument ("vtu", nargs="*", default=["potential_map_*.vtu"],
                     help=".pvtu file, .vtu pieces or glob pattern (default: potential_map_*.vtu)")
    ap.add_argument ("--pqr", help="molecule: grid centre = midpoint of the atom bounding box; atoms of the cube")
    ap.add_argument ("--center", nargs=3, type=float, help="grid centre [A] (instead of the PQR centre)")
    ap.add_argument ("--scale", type=float, default=1.0, help="grid points per A: h = 1/scale (default 1)")
    ap.add_argument ("--length", type=float, default=128.0,
                     help="box side [A]: dime = length*scale + 1 (default 128)")
    ap.add_argument ("--dime", nargs=3, type=int, help="points per side (instead of --length)")
    ap.add_argument ("--allow-interp", action="store_true",
                     help="accept grid points that are not nodes inside the uniform region")
    ap.add_argument ("--units", choices=["kT", "kcal"], default="kT", help="kT/e (default) or kcal/(mol e)")
    ap.add_argument ("--T", type=float, default=298.15, help="temperature [K] for --units kcal")
    ap.add_argument ("--tol", type=float, default=1e-6, help="geometric tolerance, fraction of h")
    ap.add_argument ("--hmap", help="also write the cell size used at each point (0 = node), same format")
    ap.add_argument ("--format", choices=["dx", "cube"], help="output format (default: from the -o extension)")
    ap.add_argument ("-o", "--out", default="potential.dx", help="output file (default potential.dx)")
    a = ap.parse_args ()

    fmt = a.format or ("cube" if a.out.lower ().endswith ((".cube", ".cub")) else "dx")
    atoms = None
    if a.pqr:
        xyz, _, zat = read_pqr (a.pqr)
        atoms = (xyz, zat)
    if a.center is not None:
        c = np.array (a.center)
    elif a.pqr:
        c = 0.5 * (xyz.min (0) + xyz.max (0))
    else:
        ap.error ("give --pqr or --center")

    h = np.full (3, 1.0 / a.scale)
    if a.dime:
        n = np.array (a.dime)
    else:
        m = a.length * a.scale
        if abs (m - round (m)) > 1e-9:
            ap.error ("length * scale must be an integer")
        n = np.full (3, int (round (m)) + 1)
    origin = c - 0.5 * (n - 1) * h

    S = Sampler (n, origin, h, a.tol)
    files = pieces (a.vtu)
    for p in files:
        S.add (*read_vtu (p))
    print ("read %d VTU pieces" % len (files))

    miss = np.isnan (S.val)
    if miss.any ():
        sys.exit ("ERROR: %d grid points are outside the NextGenPB domain: lower perfil2 in [mesh] "
                  "or use a smaller box" % miss.sum ())

    node = ~S.notnode
    hc = S.hcell
    print ("centre  [A]: %.4f %.4f %.4f" % tuple (c))
    print ("origin  [A]: %.4f %.4f %.4f" % tuple (origin))
    print ("grid       : %d x %d x %d, h = %g %g %g A" % (*n, *h))
    print ("mesh nodes : %d / %d points (%.1f %%), nodal value, no interpolation"
           % (node.sum (), node.size, 100.0 * node.mean ()))
    hr = np.round (hc[~node], 6)
    for s in np.unique (hr):
        print ("  in cells of %6.3f A : %d points (Q1 value of the cell)" % (s, np.sum (hr == s)))
    finest = hc.min ()
    unaligned = np.sum (hr == np.round (finest, 6))
    if unaligned:
        msg = ("%d points of the uniform region (cells of %g A) are not mesh nodes: the grid is "
               "not aligned with the mesh. Use scale_mesh / scale_grid integer and the mesh centre."
               % (unaligned, finest))
        if not a.allow_interp:
            sys.exit ("ERROR: " + msg)
        print ("WARNING: " + msg)
    print ("max |diff| between cells sharing a point: %.3e (Q1 continuity)" % S.maxdup)

    grid = S.val.reshape (tuple (n))
    unit = "kT/e"
    if a.units == "kcal":
        grid = grid * R_KCAL * a.T
        unit = "kcal/(mol e) at T = %g K" % a.T
    comment = ("phi [%s], FEM Q1 solution, nodal values where the grid hits mesh nodes "
               "(%.1f %% of the points)" % (unit, 100.0 * node.mean ()))
    write = write_dx if fmt == "dx" else (lambda p, g, o_, h_, cm: write_cube (p, g, o_, h_, cm, atoms))
    write (a.out, grid, origin, h, comment)
    print ("phi range  : %.4g .. %.4g %s" % (grid.min (), grid.max (), unit))
    print ("written    : %s (%s)" % (a.out, fmt))
    if a.hmap:
        write (a.hmap, np.where (node, 0.0, hc).reshape (tuple (n)), origin, h,
               "edge [A] of the mesh cell used at each point; 0 = mesh node (nodal value)")
        print ("written    : %s" % a.hmap)


if __name__ == "__main__":
    main ()
