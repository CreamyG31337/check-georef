"""Synthetic tests of the best-fit (point-to-plane ICP) in bestfit_clouds.py.

The same fit is copied into tbc-macros/PointCloudTools/BestFitScans.py, so these also guard that copy's method.
Run from the repo folder:  python -m pytest
"""
import math
import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import bestfit_clouds as bf  # noqa: E402

X0, Y0, Z0 = 314000.0, 5541000.0, 200.0    # site-like coordinates


def corner(step, rng, jitter=0.0, noise=0.001):
    """Floor + two walls of a 1 m room corner: constrains E, N and Z."""
    pts = []
    n = int(round(1 / step))
    for i in range(n):
        for j in range(n):
            a = (i + 0.5 + rng.uniform(-jitter, jitter)) * step
            b = (j + 0.5 + rng.uniform(-jitter, jitter)) * step
            e = lambda: rng.gauss(0, noise)
            pts.append((X0 + a, Y0 + b, Z0 + e()))           # floor
            pts.append((X0 + e(), Y0 + a, Z0 + b))           # wall on x = X0
            pts.append((X0 + a, Y0 + e(), Z0 + b))           # wall on y = Y0
    return pts


@pytest.fixture(scope="module")
def planes():
    return bf.fit_planes(corner(0.0125, random.Random(1)))


def test_sym3_eig_diagonal():
    vals, vec = bf.sym3_eig((3.0, 2.0, 1.0, 0.0, 0.0, 0.0))
    assert vals == pytest.approx((3.0, 2.0, 1.0))
    assert [abs(v) for v in vec] == pytest.approx([0.0, 0.0, 1.0])


def test_fit_planes_finds_the_three_surfaces(planes):
    normals = {tuple(round(abs(c)) for c in nv) for _, nv in planes.values()}
    assert normals == {(0, 0, 1), (1, 0, 0), (0, 1, 0)}


def test_icp_recovers_a_known_shift(planes):
    rng = random.Random(2)
    shift = (0.012, -0.008, 0.005)
    mov = [(x + shift[0], y + shift[1], z + shift[2])
           for x, y, z in rng.sample(corner(0.0125, rng, jitter=0.4), 3000)]
    t, rot, n, rms, N, iters = bf.icp(mov, planes)
    assert t == pytest.approx([-s for s in shift], abs=0.0015)   # t moves MOVING back onto REFERENCE
    assert n > 1000 and rms < 0.003


def test_icp_refuses_without_correspondences(planes):
    far = [(X0 + 50 + i * 0.01, Y0 + 50, Z0) for i in range(200)]
    with pytest.raises(ValueError):
        bf.icp(far, planes)


def test_shared_zone_drops_surfaces_only_one_cloud_has():
    rng = random.Random(3)
    ref = corner(0.05, rng)
    extra = [(X0 + 5 + rng.random(), Y0 + 5 + rng.random(), Z0) for _ in range(500)]   # only in MOVING
    r, m = bf.shared_zone(ref, ref + extra)
    assert not set(extra) & set(m)                           # the surface only MOVING has is gone
    assert len(m) > 0.95 * len(ref) and len(r) > 0.95 * len(ref)   # (voxels need >= 5 moving points)


def test_fmt_offset_reports_moving_relative_to_reference():
    s = bf.fmt_offset((-0.003, 0.004, 0.0))      # MOVING sits 3 mm east, 4 mm south
    assert "dE=   +3.0" in s and "dN=   -4.0" in s and "horiz   5.0 mm" in s
    assert math.isclose(float(s.split("@")[1].split()[0]), 143.1, abs_tol=0.1)
