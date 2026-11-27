# check georef

Small command-line tools for checking whether survey point clouds, scans and surfaces line up,
and finding out why when they don't.

Plain Python with no third-party packages, so they run on locked-down machines where installing
compiled packages isn't an option. LAS/LAZ files are read with LAStools' `las2las`; Trimble
Business Center installs a copy.

Requires Python 3.11+.

## Setup

Copy `config.example.toml` to `config.toml` and adjust it. Every setting is optional:
- the `las2las` path
- best-fit sizes
- default output suffix
- control CSV files
- change bands

`config.toml` is gitignored, as are `notes/`, `archive/` and `private/`, so project paths and
site details stay local. Set `CHECK_GEOREF_CONFIG` to use a config file somewhere else.

## Tools

| Script | What it does |
|---|---|
| `bestfit_clouds.py` | Point-to-plane best fit of one cloud onto another; reports the offset, how well each axis is constrained, and a rotation check |
| `shift_laz.py` | Applies a fixed dE dN dZ to LAS/LAZ files, writes a note file, and verifies the result |
| `jxl_setups.py` | Reads Trimble JobXML: coordinate system, instrument, each station's backsights (target, prism constant, height, observed − keyed), scans per station; compares keyed backsights with control CSVs |
| `grid_compare.py` | Differences XYZ surfaces on a shared grid: volumes, change bands, stable-area vertical tie, spot checks, profiles, ASCII map, JSON for a web page |
| `las_io.py` | Shared helpers: LAS header reader, las2las wrapper (crop / thin / sample to cached text), XYZ loading |
| `config.py` | Loads `config.toml` over the built-in defaults |

The Best-Fit Scans macro in [tbc-macros](https://github.com/CreamyG31337/tbc-macros) brings the
`bestfit_clouds.py` best-fit workflow into Trimble Business Center: it fits selected point cloud
regions and shifts the moving scans in TBC. See that repository for installation, usage and
current limitations.

## When two scans don't line up

1. **Check which control each job used.**
   ```
   python jxl_setups.py old-job.jxl new-job.jxl --control control-local.csv --control control-grid.csv
   ```
   Control CSVs are `name,north,east,elevation,code`. List them under `[control] files` in
   `config.toml` to skip the options.
   - The usual cause is backsights keyed from different control sets, or different monuments.
   - A job's coordinate system setting (scale-only vs a projection) only changes distances by
     a few ppm; it doesn't move heights.
2. **Measure the offset.** Export both clouds as LAS/LAZ, then:
   ```
   python bestfit_clouds.py reference.laz moving.laz [--box minE minN minZ maxE maxN maxZ]
   ```
   - Pick a box that both scans saw, with flat surfaces *and* walls facing at least two directions.
   - The output flags weakly constrained axes.
   - Expect a few mm difference between areas scanned from different setups.
3. **Apply a shift.** This writes new files and never overwrites:
   ```
   python shift_laz.py -0.010 -0.009 -0.078 a.laz b.laz --reason "moving scans were set up on the other control set"
   ```
4. **Better still, fix the control in the processing software and reprocess**, so setups, scans
   and staked points all move together.

## Surface differencing

```
python grid_compare.py week1.xyz week2.xyz week3.xyz --ascii-map 10 --profile-y <northing> --json change.json
```
- Cells match on identical XY, so the grids need the same origin and spacing.
- Volumes cover common cells only.
- A stable-area median dz well away from zero means the surveys aren't tied vertically.
- Nearly all cells identical means both files hold the same data.

## Limits and ideas

- `bestfit_clouds.py` is pure Python. It handles crops of a few hundred thousand points in about
  1–2 minutes but is too slow for whole projects. It solves a translation only; rotation is
  reported as a check.
- Per-setup fits need per-setup exports, since exported LAS may not carry a source ID per point.
- `grid_compare.py` doesn't resample grids with different origins.
- The TBC macro is still experimental and runs the fit on TBC's UI thread, so use it on a copy of
  a project and keep regions to the overlap.

## Tests

`tests/test_bestfit.py` checks the fit on a synthetic room corner (floor and two walls). It must
recover a known 12 / -8 / 5 mm shift to within 1.5 mm, find the three planes, and drop surfaces that
only one cloud has. The TBC macro's BestFitScans.py carries a copy of this fit, so run these tests
after changing either copy: `python -m pytest` (needs pytest).
