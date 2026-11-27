"""Compare gridded surfaces exported as XYZ on a shared grid (e.g. repeat bathymetry matrices).

    python grid_compare.py OLDER.xyz NEWER.xyz [MORE.xyz ...] [options]

      --labels A,B,...              names for the files (default: file names)
      --box minE minN maxE maxN     only use cells inside this box
      --profile-y Y | --profile-x X section through every surface (default for --json: median Y)
      --ascii-map SIZE              print a coarse change map with SIZE m cells (e.g. 10)
      --json OUT.json               write 1 m change map, surface grid, histogram, profile,
                                    epoch table and spot checks for a web page

Cells are matched on identical XY (rounded to 1 mm), so the grids must share origin and
spacing; the tool warns when they don't. Volume = dz x cell area over common cells only:
no TIN, no interpolation, no boundary. With more than two files, consecutive pairs are
tabulated and the detailed report covers the last pair. Bands and thresholds are set in
config.toml [grid_compare].
"""
import collections
import json
import math
import os
import statistics
import sys

from config import CONFIG

BANDS = [tuple(b) for b in CONFIG["grid_compare"]["bands"]]
BIG = CONFIG["grid_compare"]["big_change"]              # m; cells changing more than this are summarised separately
SPOT_SEPARATION = CONFIG["grid_compare"]["spot_separation"]  # m between spot-check cells


def load_grid(path, box=None):
    d, dups = {}, 0
    with open(path) as f:
        for line in f:
            s = line.replace(",", " ").split()
            if len(s) < 3:
                continue
            try:
                x, y, z = float(s[0]), float(s[1]), float(s[2])
            except ValueError:
                continue
            if box and not (box[0] <= x <= box[2] and box[1] <= y <= box[3]):
                continue
            k = (round(x, 3), round(y, 3))
            dups += k in d
            d[k] = z
    return d, dups


def step(values):
    v = sorted(set(values))
    c = collections.Counter(round(v[i + 1] - v[i], 3) for i in range(len(v) - 1))
    return c.most_common(1)[0][0] if c else None


def pair(a, b, area):
    keys = [k for k in a if k in b]
    dz = [b[k] - a[k] for k in keys]
    st = {"cells": len(keys),
          "net": round(sum(dz) * area, 1),
          "fill": round(sum(v for v in dz if v > 0) * area, 1),
          "cut": round(sum(v for v in dz if v < 0) * area, 1),
          "mean": round(statistics.fmean(dz), 3) if dz else 0.0,
          "max": round(max(dz), 2) if dz else 0.0}
    return keys, dz, st


def spots(a, b, keys, dz):
    order = sorted(zip(dz, keys), key=lambda t: -abs(t[0]))
    picks = []
    for v, k in order:
        if all(abs(k[0] - p[0]) > SPOT_SEPARATION or abs(k[1] - p[1]) > SPOT_SEPARATION for p in picks):
            picks.append(k)
        if len(picks) == 6:
            break
    calm = [k for v, k in zip(dz, keys) if abs(v) < 0.02]
    picks += calm[::max(1, len(calm) // 3)][:3]
    return [{"x": k[0], "y": k[1], "a": a[k], "b": b[k], "dz": round(b[k] - a[k], 2)} for k in picks]


def profile(surfaces, labels, y=None, x=None):
    last = surfaces[-1]
    if y is not None:
        line = min({k[1] for k in last}, key=lambda v: abs(v - y))
        along = sorted({k[0] for s in surfaces for k in s if k[1] == line})
        return f"Y = {line}", [{"x": u, **{lab: s.get((u, line)) for lab, s in zip(labels, surfaces)}} for u in along]
    line = min({k[0] for k in last}, key=lambda v: abs(v - x))
    along = sorted({k[1] for s in surfaces for k in s if k[0] == line})
    return f"X = {line}", [{"y": u, **{lab: s.get((line, u)) for lab, s in zip(labels, surfaces)}} for u in along]


def ascii_map(keys, dz, size):
    acc = collections.defaultdict(list)
    for k, v in zip(keys, dz):
        acc[(int(k[0] // size), int(k[1] // size))].append(v)
    gx = sorted({g[0] for g in acc})
    gy = sorted({g[1] for g in acc}, reverse=True)

    def sym(v):
        return "#" if v > 2 else "+" if v > 0.5 else "o" if v > 0.1 else "." if v > -0.1 else "-" if v > -0.5 else "="
    print(f"\nmean dz per {size:g} m cell  (# >2  + >0.5  o >0.1  . flat  - <-0.1  = <-0.5 m)"
          f"   E {gx[0] * size:.0f} -> {(gx[-1] + 1) * size:.0f}, N top {(gy[0] + 1) * size:.0f}")
    for j in gy:
        row = "".join(sym(statistics.fmean(acc[(i, j)])) if (i, j) in acc else " " for i in gx)
        print(f"{j * size:.0f} |{row}|")


def map_data(b, keys, dz, box, cell=1.0):
    xs, ys = [k[0] for k in keys], [k[1] for k in keys]
    x0, y0 = math.floor((box[0] if box else min(xs)) / cell) * cell, math.floor((box[1] if box else min(ys)) / cell) * cell
    x1, y1 = math.ceil((box[2] if box else max(xs)) / cell) * cell, math.ceil((box[3] if box else max(ys)) / cell) * cell
    nx, ny = int(round((x1 - x0) / cell)), int(round((y1 - y0) / cell))
    acc, accz = collections.defaultdict(list), collections.defaultdict(list)
    for k, v in zip(keys, dz):
        acc[(int((k[0] - x0) // cell), int((k[1] - y0) // cell))].append(v)
    for k, z in b.items():
        accz[(int((k[0] - x0) // cell), int((k[1] - y0) // cell))].append(z)
    rows = range(ny - 1, -1, -1)  # north up
    grid = [[round(100 * statistics.fmean(acc[(i, j)])) if acc.get((i, j)) else None for i in range(nx)] for j in rows]
    zgrid = [[round(statistics.fmean(accz[(i, j)]), 2) if accz.get((i, j)) else None for i in range(nx)] for j in rows]
    return {"grid": grid, "zgrid": zgrid, "NX": nx, "NY": ny, "X0": x0, "Y0": y0, "X1": x1, "Y1": y1}


def main():
    args = sys.argv[1:]
    opt = {}
    for flag, n in (("--labels", 1), ("--box", 4), ("--profile-y", 1), ("--profile-x", 1), ("--ascii-map", 1), ("--json", 1)):
        if flag in args:
            i = args.index(flag)
            opt[flag] = args[i + 1:i + 1 + n]
            del args[i:i + 1 + n]
    if len(args) < 2:
        sys.exit(__doc__)
    box = tuple(float(v) for v in opt["--box"]) if "--box" in opt else None
    labels = opt["--labels"][0].split(",") if "--labels" in opt else [os.path.splitext(os.path.basename(p))[0] for p in args]
    if len(labels) != len(args):
        sys.exit("--labels needs one name per file")

    surfaces = []
    for path, lab in zip(args, labels):
        g, dups = load_grid(path, box)
        dx, dy = step(k[0] for k in g), step(k[1] for k in g)
        print(f"{lab}: {len(g):,} cells, spacing {dx} x {dy} m" + (f", {dups} duplicate XY lines" if dups else ""))
        surfaces.append(g)
    area = step(k[0] for k in surfaces[-1]) * step(k[1] for k in surfaces[-1])

    epochs = []
    if len(surfaces) > 2:
        print(f"\n{'interval':<30} {'cells':>9} {'net m3':>11} {'fill m3':>10} {'cut m3':>10} {'mean dz':>8} {'max dz':>7}")
    for i in range(len(surfaces) - 1):
        _, _, st = pair(surfaces[i], surfaces[i + 1], area)
        epochs.append({"from": labels[i], "to": labels[i + 1], **st})
        if len(surfaces) > 2:
            print(f"{labels[i] + ' -> ' + labels[i + 1]:<30} {st['cells']:>9,} {st['net']:>+11,.1f} {st['fill']:>+10,.1f} "
                  f"{st['cut']:>10,.1f} {st['mean']:>+8.3f} {st['max']:>+7.2f}")

    a, b = surfaces[-2], surfaces[-1]
    la, lb = labels[-2], labels[-1]
    keys, dz, st = pair(a, b, area)
    n = len(dz)
    print(f"\n=== {la} -> {lb}   (dz = {lb} - {la}, cell area {area:g} m2)")
    print(f"common cells {n:,} ({n * area:,.0f} m2); only in {la}: {len(a) - n:,}; only in {lb}: {len(b) - n:,}")
    if n < 0.5 * min(len(a), len(b)):
        print("WARNING: under half the cells co-locate; the grids may have different origins or spacing (resample first)")
    if not n:
        sys.exit(1)
    same = sum(1 for v in dz if v == 0)
    print(f"identical z: {same:,} ({100 * same / n:.1f}%)" + ("  <-- these may be the same data" if same > 0.9 * n else ""))
    print(f"dz mean {st['mean']:+.3f}  median {statistics.median(dz):+.3f}  sd {statistics.pstdev(dz):.3f}  "
          f"min {min(dz):+.2f}  max {max(dz):+.2f} m")
    calm = [v for v in dz if abs(v) < 0.10]
    if calm:
        print(f"stable cells (|dz| < 0.10 m): {100 * len(calm) / n:.0f}%, median dz {statistics.median(calm) * 1000:+.1f} mm"
              f"  <- vertical tie between the surveys")
    print(f"volume: net {st['net']:+,.1f} m3   fill {st['fill']:+,.1f}   cut {st['cut']:,.1f}")
    pct = {str(lim): round(100 * sum(1 for v in dz if abs(v) <= lim) / n, 1) for lim in (0.05, 0.1, 0.2, 0.5, 1.0)}
    print("share of common cells with |dz| <= " + "   ".join(f"{k} m: {v}%" for k, v in pct.items()))

    hist = []
    print(f"\n{'dz band (m)':<16} {'cells':>9} {'% cells':>8} {'net m3':>10}")
    for lo, hi in BANDS:
        sel = [v for v in dz if lo <= v < hi]
        hist.append({"lo": lo, "hi": hi, "n": len(sel), "vol": round(sum(sel) * area, 1)})
        if sel:
            label = f"< {hi:+g}" if lo <= -99 else f"> {lo:+g}" if hi >= 99 else f"{lo:+g} .. {hi:+g}"
            print(f"{label:<16} {len(sel):>9,} {100 * len(sel) / n:>7.1f}% {sum(sel) * area:>+10,.1f}")

    big = [(k, v) for k, v in zip(keys, dz) if abs(v) > BIG]
    if big:
        bx, by = [k[0] for k, _ in big], [k[1] for k, _ in big]
        print(f"\ncells with |dz| > {BIG:g} m: {len(big):,} ({len(big) * area:,.0f} m2, {100 * len(big) / n:.1f}% of area), "
              f"net {sum(v for _, v in big) * area:+,.1f} m3, within E {min(bx):.1f}..{max(bx):.1f} N {min(by):.1f}..{max(by):.1f}")

    sp = spots(a, b, keys, dz)
    print(f"\nspot checks (largest changes, >{SPOT_SEPARATION:g} m apart, then stable cells):")
    for s in sp:
        print(f"  E {s['x']:.2f}  N {s['y']:.2f}   {la} {s['a']:8.2f}   {lb} {s['b']:8.2f}   dz {s['dz']:+.2f}")

    if "--ascii-map" in opt:
        ascii_map(keys, dz, float(opt["--ascii-map"][0]))

    prof, prof_label = [], None
    if "--profile-y" in opt or "--profile-x" in opt or "--json" in opt:
        if "--profile-x" in opt:
            prof_label, prof = profile(surfaces, labels, x=float(opt["--profile-x"][0]))
        else:
            y = float(opt["--profile-y"][0]) if "--profile-y" in opt else statistics.median(k[1] for k in keys)
            prof_label, prof = profile(surfaces, labels, y=y)
        if "--profile-y" in opt or "--profile-x" in opt:
            print(f"\nprofile along {prof_label}")
            print(f"{'':>11} " + " ".join(f"{lab[:10]:>10}" for lab in labels))
            for row in prof[::max(1, len(prof) // 60)]:
                pos = row.get("x", row.get("y"))
                print(f"{pos:11.2f} " + " ".join(f"{row[lab]:10.2f}" if row[lab] is not None else f"{'--':>10}" for lab in labels))

    if "--json" in opt:
        out = {**map_data(b, keys, dz, box), "profile": prof, "profileLine": prof_label, "epochs": epochs, "hist": hist,
               "counts": dict(zip(labels, (len(s) for s in surfaces))), "common": n, "onlyA": len(a) - n, "onlyB": len(b) - n,
               "spot": sp, "pct": pct, "moundArea": len(big) * area, "moundVol": round(sum(v for _, v in big) * area)}
        with open(opt["--json"][0], "w") as f:
            json.dump(out, f)
        print(f"\nwrote {opt['--json'][0]}")


if __name__ == "__main__":
    main()
