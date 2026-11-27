"""Translate LAS/LAZ files by a fixed offset, write a provenance note, and verify the result.

    python shift_laz.py dE dN dZ FILE [FILE ...] [--suffix TEXT] [--reason "why"]

Writes "<name><suffix>.laz" next to each input (suffix default from config.toml) and never
overwrites anything. The note "<output>.txt" records the shift; if "<input>.txt" export
notes exist (Trimble Business Center writes them) they are copied below it. Each output
is checked: point count, VLRs, format and scale unchanged; bounds moved by exactly the
shift; a sample of decoded points moved by exactly the shift.
"""
import datetime
import os
import sys

from config import CONFIG
from las_io import las_header, las_to_xyz, read_xyz, run_las2las


def verify(src, dst, shift):
    h0, h1 = las_header(src), las_header(dst)
    tol = max(h1["scale"]) * 1.01  # coordinates are stored as integers of the scale
    problems = []
    for key in ("points", "vlrs", "point_format", "scale"):
        if h0[key] != h1[key]:
            problems.append(f"{key} changed: {h0[key]} -> {h1[key]}")
    for i, ax in enumerate("ENZ"):
        for end in ("min", "max"):
            moved = h1[end][i] - h0[end][i]
            if abs(moved - shift[i]) > tol:
                problems.append(f"{end} {ax} moved {moved:+.4f} m, expected {shift[i]:+.4f}")
    nth = max(1, h0["points"] // 500)
    a = read_xyz(las_to_xyz(src, every_nth=nth))
    b = read_xyz(las_to_xyz(dst, every_nth=nth))
    if len(a) != len(b):
        problems.append(f"sampled {len(a)} input points but {len(b)} output points")
    worst = max((abs((q[i] - p[i]) - shift[i]) for p, q in zip(a, b) for i in range(3)), default=0.0)
    if worst > tol:
        problems.append(f"sampled points deviate from the shift by up to {worst * 1000:.2f} mm")
    return problems, len(a)


def main():
    args = sys.argv[1:]
    opts = {"--suffix": CONFIG["shift"]["suffix"], "--reason": ""}
    for k in opts:
        if k in args:
            i = args.index(k)
            opts[k] = args[i + 1]
            del args[i:i + 2]
    try:
        shift = tuple(float(v) for v in args[:3])
    except ValueError:
        sys.exit(__doc__)
    if len(shift) != 3 or len(args) < 4:
        sys.exit(__doc__)

    failed = False
    for src in args[3:]:
        base, ext = os.path.splitext(src)
        dst = base + opts["--suffix"] + ext
        if os.path.exists(dst):
            print(f"SKIP {os.path.basename(src)}: {os.path.basename(dst)} already exists")
            failed = True
            continue
        run_las2las(["-i", src, "-translate_xyz"] + [f"{v:.6f}" for v in shift]
                    + ["-olaz" if ext.lower() == ".laz" else "-olas", "-o", dst])

        note = [f"SHIFTED {datetime.date.today()}: translate dE {shift[0]:+.4f} dN {shift[1]:+.4f} "
                f"dZ {shift[2]:+.4f} m (las2las -translate_xyz)",
                f"Source: {os.path.basename(src)}"]
        if opts["--reason"]:
            note.append(f"Reason: {opts['--reason']}")
        note.append("Do NOT apply this shift again.")
        with open(dst + ".txt", "wb") as f:
            f.write(("\n".join(note) + "\n\n").encode("utf-8"))
            if os.path.exists(src + ".txt"):
                with open(src + ".txt", "rb") as notes:
                    f.write(notes.read())  # raw copy: export notes aren't necessarily UTF-8

        problems, sampled = verify(src, dst, shift)
        status = "FAIL" if problems else "OK  "
        print(f"{status} {os.path.basename(dst)}  ({las_header(dst)['points']:,} pts, {sampled} sampled)")
        for p in problems:
            print("     " + p)
        failed = failed or bool(problems)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
