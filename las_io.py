"""Shared helpers for LAS/LAZ and XYZ point data, standard library only.

LAS/LAZ points are decoded with LAStools' las2las (Trimble Business Center installs a
copy), so no Python packages have to be installed. Set its path in config.toml.
"""
import hashlib
import os
import struct
import subprocess
import tempfile

from config import CONFIG

LAS2LAS = CONFIG["tools"]["las2las"]


def is_las(path):
    return os.path.splitext(path)[1].lower() in (".las", ".laz")


def las_header(path):
    """Public header fields of a LAS/LAZ file (LAZ compresses only the points, not the header)."""
    with open(path, "rb") as f:
        h = f.read(375)
    if h[:4] != b"LASF":
        raise ValueError(f"{path} is not a LAS/LAZ file")
    minor = h[25]
    fmt = h[104]
    n = struct.unpack("<I", h[107:111])[0]
    if minor >= 4 and len(h) >= 255:
        n = struct.unpack("<Q", h[247:255])[0] or n
    sx, sy, sz, ox, oy, oz, maxx, minx, maxy, miny, maxz, minz = struct.unpack("<12d", h[131:227])
    return {
        "version": f"1.{minor}",
        "point_format": fmt & 0x3F,
        "compressed": bool(fmt & 0xC0),
        "points": n,
        "vlrs": struct.unpack("<I", h[100:104])[0],
        "scale": (sx, sy, sz),
        "offset": (ox, oy, oz),
        "min": (minx, miny, minz),
        "max": (maxx, maxy, maxz),
    }


def las_bbox(path):
    h = las_header(path)
    return h["min"] + h["max"]


def run_las2las(args):
    if not os.path.exists(LAS2LAS):
        raise FileNotFoundError(f"las2las not found at {LAS2LAS}; set [tools] las2las in config.toml")
    r = subprocess.run([LAS2LAS] + [str(a) for a in args], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"las2las failed ({r.returncode}): {(r.stderr or r.stdout).strip()}")
    return r


def las_to_xyz(path, box=None, thin=None, every_nth=None):
    """Decode LAS/LAZ to a cached 'x y z' text file in the temp folder and return its path.

    box = (minx, miny, minz, maxx, maxy, maxz) crops, thin = grid thinning size (m),
    every_nth keeps one point in n.
    """
    opts = []
    if box:
        opts += ["-keep_xyz"] + [f"{v:.4f}" for v in box]
    if thin:
        opts += ["-thin_with_grid", str(thin)]
    if every_nth and every_nth > 1:
        opts += ["-keep_every_nth", str(int(every_nth))]
    key = (os.path.abspath(path), os.path.getmtime(path), os.path.getsize(path), opts)
    tag = hashlib.md5(repr(key).encode()).hexdigest()[:12]
    out = os.path.join(tempfile.gettempdir(), f"las_io_{tag}.txt")
    if not os.path.exists(out):
        part = out[:-4] + ".part.txt"  # renamed when complete so an interrupted run is never reused
        run_las2las(["-i", path] + opts + ["-otxt", "-oparse", "xyz", "-o", part])
        os.replace(part, out)
    return out


def read_xyz(path):
    pts = []
    with open(path) as f:
        for line in f:
            parts = line.replace(",", " ").split()
            if len(parts) >= 3:
                try:
                    pts.append((float(parts[0]), float(parts[1]), float(parts[2])))
                except ValueError:
                    pass  # header or comment line
    return pts


def load_xyz(path, box=None, thin=None):
    """Points from a text file (x y z first on each line) or LAS/LAZ, as a list of (x, y, z)."""
    if is_las(path):
        return read_xyz(las_to_xyz(path, box=box, thin=thin))
    return read_xyz(path)
