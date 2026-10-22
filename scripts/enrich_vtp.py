#!/usr/bin/env python3
"""
enrich_vtp.py - build a self-consistent, q_pol-enriched surface VTP for internal
potential reconstruction.

WHAT IT DOES
    For each structure it builds a binary (zlib-compressed) VTP whose nodes are
    the p4est edge crossings from `vertexdata.csv` (the EXACT source) and adds the
    per-node polarization charge `q_pol` (flux of D through the surface). Together
    with `phi` (surface potential) and `Normals`, this lets one reconstruct the
    electrostatic potential at any interior point (see the companion tool
    ngpb_potential.py), reproducing pot_field_fast.

    Per structure <CODE>:
      1. read   <ngpb>/<CODE>/vertexdata.csv          -> nodes, phi0, N, q_pol
      2. read   <ngpb>/<CODE>/triangulatedSurf_phi.vtp -> triangle connectivity
      3. re-index the connectivity onto the vertexdata nodes (nearest, used only
         for the per-node effective area in the small ionic geometric term)
      4. write  <out>/<code>.vtp  (Float64 points + q_pol, binary + zlib)
      5. round-trip validation of phi and q_pol

WHY THIS (and not the native VTP nodes)
    The published VTP is built on the NanoShaper mesh, whose vertices are NOT in
    exact one-to-one correspondence with the p4est crossings (they differ by up
    to ~1e-3 A and a few hundred fold-points are near-degenerate). phi_p is very
    sensitive to the (node, q_pol) pairing, so building the enriched VTP directly
    on the vertexdata nodes makes the reconstruction match the solver to machine
    precision.

PHYSICS
    q_pol(V) = -(phi2 - phi1) * wha(eps1,eps2,alpha) * flux_dir(eps1,eps2) * area_h
      wha(e1,e2,a)    = 1 / (a/e1 + (1-a)/e2)
      flux_dir(e1,e2) = +1 if e1<e2 else -1
      area_h          = 0.5  (= h_perp1*h_perp2/h_axis for a cubic h = 0.5 A grid;
                              the factor 4 from the edges shared by 4 cubes is
                              already folded in)
    (definitions: NextGenPB src/pb_class.cpp -> energy_fast / pot_field_fast)

USAGE
    # single structure
    python scripts/enrich_vtp.py --single 4LZI --jobs 1

    # whole dataset (parallel); paths can be overridden with --ngpb / --out
    python scripts/enrich_vtp.py --jobs 8

    # only structures whose output does not exist yet
    python scripts/enrich_vtp.py --jobs 8 --skip-done

REQUIREMENTS
    numpy, vtk, scipy
"""

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

# Default locations (override on the command line with --ngpb / --out).
DEFAULT_NGPB = Path("/mnt/vdiflorio/ngpb_runs")
DEFAULT_OUT  = Path("/mnt/vdiflorio/zenodo_staging/vtp")
SRC_NAME     = "triangulatedSurf_phi.vtp"
VDATA_NAME   = "vertexdata.csv"

# Per-node area (see header). h = 0.5 A fixed in this dataset.
AREA_H = 0.5

# Column indices in vertexdata.csv (known header).
C_PHI1, C_EPS1, C_ALPHA = 0, 1, 2
C_N0, C_N1, C_N2 = 3, 4, 5          # N_nu, N_nu1, N_nu2 (permuted by 'axis')
C_PHI2, C_EPS2 = 6, 7
C_X0, C_Y0, C_Z0, C_PHI0 = 24, 25, 26, 27
C_AXIS = 31


def parse_vertexdata(vd: np.ndarray):
    """From vertexdata -> (V, phi0, normals, q_pol) per node (vertexdata order)."""
    phi1, eps1, alpha = vd[:, C_PHI1], vd[:, C_EPS1], vd[:, C_ALPHA]
    phi2, eps2 = vd[:, C_PHI2], vd[:, C_EPS2]
    axis = vd[:, C_AXIS].astype(int)

    wha = 1.0 / (alpha / eps1 + (1.0 - alpha) / eps2)
    flux_dir = np.where(eps1 < eps2, 1.0, -1.0)
    q_pol = -(phi2 - phi1) * wha * flux_dir * AREA_H

    V = vd[:, [C_X0, C_Y0, C_Z0]].copy()
    phi0 = vd[:, C_PHI0].copy()

    # de-permute the normal: stored[jj] = trueN[(axis+jj)%3] -> trueN[k]=stored[(k-axis)%3]
    stored = vd[:, [C_N0, C_N1, C_N2]]
    n = np.zeros_like(stored)
    rows = np.arange(len(axis))
    for k in range(3):
        n[:, k] = stored[rows, (k - axis) % 3]
    return V, phi0, n, q_pol


def enrich_one(args: tuple) -> dict:
    code, src, vdata, dst = args
    import vtk
    try:
        from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk
    except ImportError:
        from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
    from scipy.spatial import cKDTree

    res = {"pdb_code": code, "status": "error", "valid": 0,
           "n_points": None, "n_tris": None, "error_msg": None}
    try:
        Path(dst).parent.mkdir(parents=True, exist_ok=True)

        vd = np.loadtxt(vdata, delimiter=",", skiprows=1)
        V, phi0, normals, q_pol = parse_vertexdata(vd)
        n_pts = V.shape[0]

        # triangle connectivity from the original VTP, re-indexed onto vertexdata nodes
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(str(src))
        reader.Update()
        pd_in = reader.GetOutput()
        if pd_in.GetNumberOfPoints() != n_pts:
            raise RuntimeError(
                f"node count mismatch: vertexdata={n_pts} vtp={pd_in.GetNumberOfPoints()}")
        vtp_pts = vtk_to_numpy(pd_in.GetPoints().GetData()).astype(np.float64)
        polys = vtk_to_numpy(pd_in.GetPolys().GetData()).reshape(-1, 4)[:, 1:]

        _, vtp2vd = cKDTree(V).query(vtp_pts)
        tris = vtp2vd[polys].astype(np.int64)            # connectivity in vertexdata order

        # --- build the new polydata ---
        pts = vtk.vtkPoints()
        pts.SetData(numpy_to_vtk(np.ascontiguousarray(V, np.float64), deep=1))
        out = vtk.vtkPolyData()
        out.SetPoints(pts)

        cells = vtk.vtkCellArray()
        conn = np.empty((tris.shape[0], 4), np.int64)
        conn[:, 0] = 3
        conn[:, 1:] = tris
        idarr = numpy_to_vtk(conn.ravel(), deep=1, array_type=vtk.VTK_ID_TYPE)
        cells.SetCells(tris.shape[0], idarr)
        out.SetPolys(cells)

        def add(name, a, dtype):
            v = numpy_to_vtk(np.ascontiguousarray(a, dtype), deep=1)
            v.SetName(name)
            out.GetPointData().AddArray(v)
            return v

        add("phi", phi0, np.float32)
        add("q_pol", q_pol, np.float64)
        nrm = add("Normals", normals, np.float32)
        out.GetPointData().SetNormals(nrm)

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(str(dst))
        writer.SetInputData(out)
        writer.SetDataModeToBinary()
        writer.SetCompressorTypeToZLib()
        if writer.Write() != 1:
            raise RuntimeError("vtkXMLPolyDataWriter.Write() returned 0")

        # round-trip check on q_pol and phi
        r2 = vtk.vtkXMLPolyDataReader()
        r2.SetFileName(str(dst))
        r2.Update()
        pd2 = r2.GetOutput()
        ok = (pd2.GetNumberOfPoints() == n_pts and
              pd2.GetNumberOfPolys() == tris.shape[0])
        ok = ok and np.array_equal(
            vtk_to_numpy(pd2.GetPointData().GetArray("q_pol")), q_pol)
        ok = ok and np.allclose(
            vtk_to_numpy(pd2.GetPointData().GetArray("phi")),
            phi0.astype(np.float32))

        res.update(status="enriched" if ok else "error",
                   valid=1 if ok else 0, n_points=n_pts, n_tris=tris.shape[0],
                   error_msg=None if ok else "round-trip failed")
        if not ok:
            Path(dst).unlink(missing_ok=True)
    except Exception as e:
        res["error_msg"] = str(e)[:500]
        Path(dst).unlink(missing_ok=True)
    return res


def get_pending(ngpb: Path, out: Path, only) -> list:
    items = []
    dirs = ([ngpb / only] if only else
            sorted(d for d in ngpb.iterdir() if d.is_dir()))
    for d in dirs:
        src, vdata = d / SRC_NAME, d / VDATA_NAME
        if src.exists() and vdata.exists():
            items.append((d.name, src, vdata, out / f"{d.name.lower()}.vtp"))
    return items


def main():
    ap = argparse.ArgumentParser(description="Build q_pol-enriched surface VTPs.")
    ap.add_argument("--ngpb", type=Path, default=DEFAULT_NGPB,
                    help="directory with <CODE>/ run folders")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT,
                    help="output directory for the enriched VTPs")
    ap.add_argument("--single", help="only this structure (e.g. 4LZI)")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--skip-done", action="store_true",
                    help="skip structures whose output already exists")
    ap.add_argument("--limit", type=int)
    args = ap.parse_args()

    items = get_pending(args.ngpb, args.out, args.single)
    if args.skip_done:
        items = [it for it in items if not it[3].exists()]
    if args.limit:
        items = items[:args.limit]
    if not items:
        print("Nothing to do.")
        return

    print(f"Structures to enrich: {len(items)}  (jobs={args.jobs})")
    ok = err = 0
    if args.jobs == 1:
        for it in items:
            r = enrich_one(it); _report(r)
            ok += r["valid"]; err += (r["valid"] == 0)
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futs = {ex.submit(enrich_one, it): it for it in items}
            for f in as_completed(futs):
                r = f.result(); _report(r)
                ok += r["valid"]; err += (r["valid"] == 0)
    print(f"\nDone: {ok} OK, {err} errors.")
    sys.exit(1 if err else 0)


def _report(r: dict):
    if r["valid"]:
        print(f"  [OK]  {r['pdb_code']:6s}  nodes={r['n_points']}  tris={r['n_tris']}")
    else:
        print(f"  [ERR] {r['pdb_code']:6s}  {r['error_msg']}")


if __name__ == "__main__":
    main()
