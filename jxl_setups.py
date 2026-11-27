"""Summarise total-station setups in Trimble JobXML (.jxl) files.

    python jxl_setups.py JOB.jxl [JOB2.jxl ...] [--control FILE.csv ...] [--no-control]

For each job: coordinate system, instrument, and every station with its resected
coordinates, the scans shot from it, and each backsight's target (prism constant,
height) with observed-minus-keyed coordinates. Keyed backsight coordinates are also
compared with control CSVs (name,north,east,elevation,code) given with --control, or
listed under [control] files in config.toml. When several control sets exist for a
site, this shows at a glance which one a job was set up on.
"""
import csv
import os
import sys
import xml.etree.ElementTree as ET

from config import CONFIG


def text(el, path):
    v = el.find(path) if el is not None else None
    return (v.text or "") if v is not None else ""


def num(el, path):
    try:
        return float(text(el, path))
    except ValueError:
        return None


def coords(el, tag):
    g = el.find(tag)
    if g is None:
        return None
    e, n, z = num(g, "East"), num(g, "North"), num(g, "Elevation")
    return None if e is None or n is None else (e, n, z)


def load_control(path):
    pts = {}
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as f:
        for row in csv.reader(f):
            if len(row) < 4:
                continue
            try:
                n, e, z = float(row[1]), float(row[2]), float(row[3])
            except ValueError:
                continue
            name = row[0].strip()
            pts.setdefault(name, (e, n, z))
            for suffix in CONFIG["control"]["name_suffixes"]:
                if suffix and name.endswith(suffix):
                    pts.setdefault(name[:-len(suffix)], (e, n, z))
    return pts


def mm(a, b):
    if a is None or b is None:
        return "      n/a          "
    dz = "   n/a" if a[2] is None or b[2] is None else f"{(a[2] - b[2]) * 1000:+6.1f}"
    return f"{(a[0] - b[0]) * 1000:+6.1f} {(a[1] - b[1]) * 1000:+6.1f} {dz}"


def fmt_xyz(p):
    if p is None:
        return "n/a"
    return f"E {p[0]:.3f} N {p[1]:.3f}" + (f" Z {p[2]:.3f}" if p[2] is not None else "")


def summarise(path, controls):
    root = ET.parse(path).getroot()
    cs = root.find("Environment/CoordinateSystem")
    print(f"\n=== {os.path.basename(path)}   {root.get('product')} {root.get('productVersion')}, created {root.get('TimeStamp')}")
    if text(cs, "Projection/Type") == "ScaleOnly":
        print(f"coordinate system: Scale Only {text(cs, 'Projection/Scale')}, no datum/geoid")
    else:
        print(f"coordinate system: {text(cs, 'SystemName')} {text(cs, 'ZoneName')} / {text(cs, 'DatumName')} / "
              f"geoid {text(cs, 'VerticalAdjustment/GeoidName') or 'none'} / "
              f"project height {text(cs, 'LocalSite/ProjectLocationHeight') or '-'}")

    targets, keyed, instruments, stations, by_id = {}, {}, [], [], {}
    for r in root.find("FieldBook"):
        if r.tag == "InstrumentRecord":
            inst = f"{text(r, 'Model')} s/n {text(r, 'Serial')}"
            if inst not in instruments:
                instruments.append(inst)
        elif r.tag == "TargetRecord":
            targets[r.get("ID")] = (text(r, "Name"), num(r, "PrismConstant"), num(r, "TargetHeight"))
        elif r.tag == "StationRecord":
            st = {"name": text(r, "StationName"), "time": r.get("TimeStamp"), "type": text(r, "StationType"),
                  "coords": None, "obs": {}, "scans": []}
            stations.append(st)
            by_id[r.get("ID")] = st
        elif r.tag == "ScanningFrameRecord":
            st = by_id.get(text(r, "StationID"))
            if st:
                st["scans"].append(text(r, "ScanFrameName"))
        elif r.tag == "PointRecord" and text(r, "Deleted") != "true":
            name, cls = text(r, "Name"), text(r, "Classification")
            g = coords(r, "Grid")
            st = by_id.get(text(r, "StationID"))
            if cls in ("MTA", "BackSight") and st and coords(r, "ComputedGrid"):
                o = st["obs"].setdefault(name, {"MTA": [], "BackSight": [], "target": None, "keyed": keyed.get(name)})
                o[cls].append(coords(r, "ComputedGrid"))
                o["target"] = targets.get(text(r, "TargetID"), o["target"])
            elif g:
                if stations and name == stations[-1]["name"] and stations[-1]["coords"] is None:
                    stations[-1]["coords"] = g  # the resection / setup result
                else:
                    keyed[name] = g
    if instruments:
        print("instrument: " + "; ".join(instruments))

    for st in stations:
        print(f"\nstation {st['name']}  {st['time']}  {st['type']}  {fmt_xyz(st['coords'] or keyed.get(st['name']))}")
        if st["scans"]:
            print("  scans: " + " | ".join(st["scans"]))
        if not st["obs"]:
            print("  no backsight observations recorded (e.g. 'used last station setup')")
            continue
        print(f"  {'backsight':<12} {'target':<22} {'PC mm':>6} {'HT m':>6}  observed - keyed (mm)  dE     dN     dZ")
        for name, o in st["obs"].items():
            obs = o["MTA"] or o["BackSight"]
            mean = tuple(sum(p[i] for p in obs) / len(obs) for i in range(3))
            t = o["target"] or ("?", None, None)
            pc = f"{t[1] * 1000:+.0f}" if t[1] is not None else "?"
            ht = f"{t[2]:.3f}" if t[2] is not None else "?"
            print(f"  {name:<12} {t[0][:22]:<22} {pc:>6} {ht:>6}  {'':>21} {mm(mean, o['keyed'])}")

    if controls:
        names = sorted({n for st in stations for n in st["obs"]})
        print("\n  keyed backsight coordinates minus control file (mm)   dE     dN     dZ")
        for cpath, cpts in controls:
            print(f"  [{os.path.basename(cpath)}]")
            for n in names:
                print(f"    {n:<12} {mm(keyed.get(n), cpts.get(n)) if n in cpts else '   not in file'}")


def main():
    args = sys.argv[1:]
    paths = []
    while "--control" in args:
        i = args.index("--control")
        paths.append(args[i + 1])
        del args[i:i + 2]
    use_config = "--no-control" not in args
    args = [a for a in args if a != "--no-control"]
    if not args:
        sys.exit(__doc__)
    if not paths and use_config:
        paths = CONFIG["control"]["files"]
    controls = []
    for p in paths:
        if os.path.exists(p):
            controls.append((p, load_control(p)))
        else:
            print(f"control file not found, skipped: {p}")
    for path in args:
        summarise(path, controls)


if __name__ == "__main__":
    main()
