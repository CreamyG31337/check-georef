"""Best-fit one point cloud onto another (point-to-plane ICP), pure stdlib Python.

    python bestfit_clouds.py REFERENCE MOVING [--box minE minN minZ maxE maxN maxZ] [--shift dE dN dZ]

Inputs are .txt (x y z per line) or .las/.laz (read with the las2las64.exe that
ships with Trimble Business Center; both LAS files are cropped to their header
overlap at full density). Only surfaces present in both clouds are used.
Reports the offset of MOVING relative to REFERENCE: translation-only fit from
several starting shifts, how strongly each axis is constrained, per-surface-type
estimates, and a 6-DOF rotation check.

The same fit (shared_zone, sym3_eig, fit_planes, correspond, solve, median, icp) is copied
into the TBC macro tbc-macros/PointCloudTools/BestFitScans.py, which runs inside TBC's
IronPython and can't import this file. A fix to the fit here belongs there too, and back.
"""
import collections
import math
import os
import random
import sys

from config import CONFIG
from las_io import is_las, las_bbox, load_xyz as load

CELL = CONFIG["bestfit"]["cell"]                  # plane-fit cell size on the reference (m); raised to 0.10 for sparse clouds
MAX_SAMPLES = CONFIG["bestfit"]["max_samples"]    # moving-cloud points used for the fit
SEARCH = CONFIG["bestfit"]["search"]              # max distance from a moving point to a reference plane centroid (m)
SHARED_VOXEL = CONFIG["bestfit"]["shared_voxel"]  # voxel size used to find surfaces present in both clouds (m)


def bbox(pts):
    xs, ys, zs = zip(*pts)
    return min(xs), min(ys), min(zs), max(xs), max(ys), max(zs)


def overlap(a, b):
    o = tuple(max(a[i], b[i]) for i in range(3)) + tuple(min(a[i], b[i]) for i in range(3, 6))
    return None if any(o[i] >= o[i + 3] for i in range(3)) else o


def crop(pts, b, pad):
    return [p for p in pts if b[0] - pad <= p[0] <= b[3] + pad
            and b[1] - pad <= p[1] <= b[4] + pad and b[2] - pad <= p[2] <= b[5] + pad]


def shared_zone(ref, mov):
    """Keep only points near surfaces that both clouds contain (tolerant of offsets < ~SHARED_VOXEL)."""
    V = SHARED_VOXEL
    vk = lambda p: (math.floor(p[0] / V), math.floor(p[1] / V), math.floor(p[2] / V))
    nb = [(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]
    rv = collections.Counter(map(vk, ref))
    mc = collections.Counter(map(vk, mov))
    sm = {k for k, c in mc.items() if c >= 5 and any((k[0] + a, k[1] + b, k[2] + c) in rv for a, b, c in nb)}
    sr = {(k[0] + a, k[1] + b, k[2] + c) for k in sm for a, b, c in nb}
    return [p for p in ref if vk(p) in sr], [p for p in mov if vk(p) in sm]


def sym3_eig(c):
    """Eigenvalues (descending) and smallest eigenvector of a symmetric 3x3."""
    a, b, d, e, f, g = c  # xx, yy, zz, xy, xz, yz
    p1 = e * e + f * f + g * g
    q = (a + b + d) / 3
    p2 = (a - q) ** 2 + (b - q) ** 2 + (d - q) ** 2 + 2 * p1
    p = math.sqrt(p2 / 6) if p2 > 0 else 0.0
    if p < 1e-15:
        return (q, q, q), (0.0, 0.0, 1.0)
    B = [[(a - q) / p, e / p, f / p], [e / p, (b - q) / p, g / p], [f / p, g / p, (d - q) / p]]
    detB = (B[0][0] * (B[1][1] * B[2][2] - B[1][2] * B[2][1])
            - B[0][1] * (B[1][0] * B[2][2] - B[1][2] * B[2][0])
            + B[0][2] * (B[1][0] * B[2][1] - B[1][1] * B[2][0]))
    r = max(-1.0, min(1.0, detB / 2))
    phi = math.acos(r) / 3
    l1 = q + 2 * p * math.cos(phi)
    l3 = q + 2 * p * math.cos(phi + 2 * math.pi / 3)
    l2 = 3 * q - l1 - l3
    M = [[a - l3, e, f], [e, b - l3, g], [f, g, d - l3]]
    best, bn = None, -1
    for i, j in ((0, 1), (0, 2), (1, 2)):
        u, v = M[i], M[j]
        cr = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
        n = cr[0] ** 2 + cr[1] ** 2 + cr[2] ** 2
        if n > bn:
            best, bn = cr, n
    s = math.sqrt(bn) or 1.0
    return (l1, l2, l3), (best[0] / s, best[1] / s, best[2] / s)


def fit_planes(ref):
    cells = {}
    for p in ref:
        cells.setdefault((math.floor(p[0] / CELL), math.floor(p[1] / CELL), math.floor(p[2] / CELL)), []).append(p)
    planes = {}
    for k, cp in cells.items():
        if len(cp) < 10:
            continue
        n = len(cp)
        mx = sum(p[0] for p in cp) / n; my = sum(p[1] for p in cp) / n; mz = sum(p[2] for p in cp) / n
        xx = yy = zz = xy = xz = yz = 0.0
        for x, y, z in cp:
            dx, dy, dz = x - mx, y - my, z - mz
            xx += dx * dx; yy += dy * dy; zz += dz * dz; xy += dx * dy; xz += dx * dz; yz += dy * dz
        (l1, l2, l3), nv = sym3_eig((xx / n, yy / n, zz / n, xy / n, xz / n, yz / n))
        if l1 <= 0 or l3 / (l1 + l2 + l3) > 0.01 or l2 / l1 < 0.1 or math.sqrt(max(l3, 0)) > 0.006:
            continue  # not a clean plane (edge, vegetation, corner, noise)
        planes[k] = ((mx, my, mz), nv)
    return planes


def correspond(src, planes, t, maxres):
    out = []
    for p in src:
        q = (p[0] + t[0], p[1] + t[1], p[2] + t[2])
        ci, cj, ck = math.floor(q[0] / CELL), math.floor(q[1] / CELL), math.floor(q[2] / CELL)
        reach = max(1, math.ceil(SEARCH / CELL))
        best, bd = None, SEARCH * SEARCH
        for di in range(-reach, reach + 1):
            for dj in range(-reach, reach + 1):
                for dk in range(-reach, reach + 1):
                    pl = planes.get((ci + di, cj + dj, ck + dk))
                    if pl:
                        c = pl[0]
                        d = (q[0] - c[0]) ** 2 + (q[1] - c[1]) ** 2 + (q[2] - c[2]) ** 2
                        if d < bd:
                            best, bd = pl, d
        if best:
            c, n = best
            r = n[0] * (q[0] - c[0]) + n[1] * (q[1] - c[1]) + n[2] * (q[2] - c[2])
            # stay within the plane's footprint so a point isn't matched to an unrelated surface
            if abs(r) < maxres and bd - r * r < (CELL * 0.75) ** 2:
                out.append((q, n, r))
    return out


def solve(A, b):
    n = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for i in range(n):
        piv = max(range(i, n), key=lambda r: abs(M[r][i]))
        M[i], M[piv] = M[piv], M[i]
        if abs(M[i][i]) < 1e-12:
            raise ValueError("singular normal matrix (geometry does not constrain all axes)")
        for r in range(n):
            if r != i:
                f = M[r][i] / M[i][i]
                for c in range(i, n + 1):
                    M[r][c] -= f * M[i][c]
    return [M[i][n] / M[i][i] for i in range(n)]


def median(v):
    s = sorted(v)
    m = len(s) // 2
    return s[m] if len(s) % 2 else (s[m - 1] + s[m]) / 2


def icp(src, planes, dof6=False, subset=None, iters=30, start=(0.0, 0.0, 0.0)):
    t = list(start)
    rot = [0.0, 0.0, 0.0]
    maxres = SEARCH
    ctr = [sum(p[i] for p in src) / len(src) for i in range(3)]
    for it in range(iters):
        cor = correspond(src, planes, t, maxres)
        if subset:
            cor = [c for c in cor if subset(c[1])]
        if len(cor) < 50:
            raise ValueError(f"only {len(cor)} correspondences")
        res = [c[2] for c in cor]
        med = median(res)
        mad = median([abs(r - med) for r in res]) * 1.4826
        if it >= 2:
            maxres = max(0.01, 3 * mad)
        k = 6 if dof6 else 3
        AtA = [[0.0] * k for _ in range(k)]
        Atb = [0.0] * k
        for q, n, r in cor:
            if abs(r) > maxres:
                continue
            if dof6:
                px, py, pz = q[0] - ctr[0], q[1] - ctr[1], q[2] - ctr[2]
                row = (py * n[2] - pz * n[1], pz * n[0] - px * n[2], px * n[1] - py * n[0], n[0], n[1], n[2])
            else:
                row = n
            for i in range(k):
                Atb[i] -= row[i] * r
                for j in range(k):
                    AtA[i][j] += row[i] * row[j]
        x = solve(AtA, Atb)
        dt = x[3:] if dof6 else x
        if dof6:
            rot = [rot[i] + x[i] for i in range(3)]
            # apply the small rotation about the centroid to the source points
            wx, wy, wz = x[0], x[1], x[2]
            src = [(p[0] + (wy * (p[2] - ctr[2]) - wz * (p[1] - ctr[1])),
                    p[1] + (wz * (p[0] - ctr[0]) - wx * (p[2] - ctr[2])),
                    p[2] + (wx * (p[1] - ctr[1]) - wy * (p[0] - ctr[0]))) for p in src]
        t = [t[i] + dt[i] for i in range(3)]
        if math.sqrt(sum(d * d for d in dt)) < 5e-5 and it >= 3:
            break
    cor = correspond(src, planes, t, maxres)
    if subset:
        cor = [c for c in cor if subset(c[1])]
    inl = [c for c in cor if abs(c[2]) <= maxres]
    rms = math.sqrt(sum(c[2] ** 2 for c in inl) / len(inl))
    N = [[sum(c[1][i] * c[1][j] for c in inl) for j in range(3)] for i in range(3)]
    return t, rot, len(inl), rms, N, it + 1


def fmt_offset(t):
    # t moves MOVING onto REFERENCE, so MOVING sits at -t relative to REFERENCE
    e, n, z = -t[0], -t[1], -t[2]
    brg = (math.degrees(math.atan2(e, n)) + 360) % 360
    return (f"dE={e*1000:+7.1f}  dN={n*1000:+7.1f}  dZ={z*1000:+7.1f} mm   "
            f"horiz {math.hypot(e, n)*1000:5.1f} mm @ {brg:5.1f} deg")


def main():
    global CELL
    args = sys.argv[1:]
    shift = None
    if "--shift" in args:
        i = args.index("--shift")
        shift = tuple(float(v) for v in args[i + 1:i + 4])
        del args[i:i + 4]
    user_box = None
    if "--box" in args:
        i = args.index("--box")
        v = [float(x) for x in args[i + 1:i + 7]]
        user_box = (min(v[0], v[3]), min(v[1], v[4]), min(v[2], v[5]), max(v[0], v[3]), max(v[1], v[4]), max(v[2], v[5]))
        del args[i:i + 7]
    if len(args) != 2:
        sys.exit(__doc__)
    box = None
    if user_box:
        box = user_box
    elif all(is_las(a) for a in args):
        box = overlap(las_bbox(args[0]), las_bbox(args[1]))
        if not box:
            sys.exit("LAS headers do not overlap")
        box = tuple(v - 0.5 for v in box[:3]) + tuple(v + 0.5 for v in box[3:])
    ref, mov = load(args[0], box), load(args[1], box)
    if shift:  # self-test: apply a known offset to MOVING
        mov = [(p[0] + shift[0], p[1] + shift[1], p[2] + shift[2]) for p in mov]
    if user_box:  # .txt inputs aren't cropped by las2las
        ref, mov = crop(ref, user_box, 0), crop(mov, user_box, 0)
        if not ref or not mov:
            sys.exit(f"box is empty: reference {len(ref)} pts, moving {len(mov)} pts")
    print(f"reference {len(ref)} pts, moving {len(mov)} pts")
    ob = overlap(bbox(ref), bbox(mov))
    if not ob:
        sys.exit("clouds do not overlap")
    ref, mov = shared_zone(crop(ref, ob, 0.5), crop(mov, ob, 0.2))
    print(f"shared surfaces: reference {len(ref)} pts, moving {len(mov)} pts")
    if len(mov) < 500:
        sys.exit("too little shared surface to fit")
    planes = fit_planes(ref)
    if len(planes) < 200:
        CELL = 0.10
        planes = fit_planes(ref)
    hz = sum(1 for p in planes.values() if abs(p[1][2]) > 0.9)
    vt = sum(1 for p in planes.values() if abs(p[1][2]) < 0.3)
    print(f"reference planes ({CELL*100:.0f} cm cells): {len(planes)} "
          f"({hz} horizontal, {vt} vertical, {len(planes)-hz-vt} sloped)")
    random.seed(1)
    src = random.sample(mov, min(MAX_SAMPLES, len(mov)))

    print("\ntranslation-only, from several starting shifts:")
    fits = []
    for st in ((0, 0, 0), (0.06, 0.06, 0.12), (-0.06, -0.06, -0.12)):
        try:
            fits.append(icp(src, planes, start=st))
            print(f"  start {fmt_offset(st)[:38]} -> {fmt_offset(fits[-1][0])}")
        except ValueError as ex:
            print(f"  start {st}: {ex}")
    if not fits:
        sys.exit("no fit converged")
    t, _, n, rms, N, its = min(fits, key=lambda f: f[3])
    spread = max(math.dist(f[0], t) for f in fits)
    print(f"\nBEST: {fmt_offset(t)}")
    print(f"      {n} correspondences, rms {rms*1000:.1f} mm, start-to-start spread {spread*1000:.1f} mm")
    tot = N[0][0] + N[1][1] + N[2][2]
    for i, ax in enumerate("ENZ"):
        share = N[i][i] / tot
        sigma = rms / math.sqrt(N[i][i]) if N[i][i] > 0 else float("inf")
        note = "  <-- weakly constrained" if share < 0.02 else ""
        print(f"      {ax}: {share*100:4.1f}% of constraint, formal sd {sigma*1000:.2f} mm{note}")
    if spread > 0.002 or rms > 0.008:
        print("      WARNING: fits disagree or residuals are large; treat the result with caution")

    print()
    for label, f in (("horizontal surfaces", lambda nv: abs(nv[2]) > 0.9),
                     ("sloped/vertical surfaces", lambda nv: abs(nv[2]) <= 0.9)):
        try:
            ts, _, ns, rs, _, _ = icp(src, planes, subset=f, start=t)
            print(f"{label:26} {fmt_offset(ts)}   ({ns} pts, rms {rs*1000:.1f} mm)")
        except ValueError as ex:
            print(f"{label:26} {ex}")
    try:
        t6, r6, _, r6rms, _, _ = icp(src, planes, dof6=True, start=t)
        print(f"{'6-DOF check':26} {fmt_offset(t6)}")
        print(f"{'':26} rotation x {math.degrees(r6[0])*3600:+.1f}\"  y {math.degrees(r6[1])*3600:+.1f}\"  "
              f"z {math.degrees(r6[2])*3600:+.1f}\"   rms {r6rms*1000:.1f} mm")
    except ValueError as ex:
        print(f"6-DOF check: {ex}")


if __name__ == "__main__":
    main()
