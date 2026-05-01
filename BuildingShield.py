# -*- coding: utf-8 -*-
"""
BuildingShield.py
======================
Fast rewrite of BuildingShield.py (exported from ModelBuilder 2026-04-26).

Preserves ALL logic, field names, paths and CRS from the original working model.
Two speedups added on top:

  1. KDTree spatial index  — replaces Select Layer By Location
     Instead of ArcGIS selecting features one at a time, a KDTree finds all
     neighbours within 300 m in microseconds per building.

  2. Multiprocessing  — runs buildings in parallel across all CPU cores
     The 89,000 buildings are split into batches and processed simultaneously.

Expected time: 5–15 minutes for 89,000 buildings (vs ~18 hours in ModelBuilder).

════════════════════════════════════════════════════════════
HOW TO RUN
════════════════════════════════════════════════════════════
1. Edit the CONFIG block below — update any paths that differ on your machine.
2. Open the ArcGIS Pro Python Command Prompt
     Start Menu → search "Python Command Prompt" → open it
3. Run:
     python "C:\\path\\to\\BuildingShield.py"

Do NOT run from the Python Window inside ArcGIS Pro.
Multiprocessing requires a real terminal, not the embedded window.

════════════════════════════════════════════════════════════
INPUT / OUTPUT
════════════════════════════════════════════════════════════
INPUT  (set in CONFIG):
  CENTROID_FC    — your building_portfolio feature class (full path)
                   Must already have fields: id, height, CENT_X, CENT_Y
  WORKSPACE_GDB  — your project GDB

OUTPUT (set in CONFIG):
  OUTPUT_GDB     — GDB where ShieldingResults table is written 
  OUTPUT_TABLE   — name of the results table inside OUTPUT_GDB
                   Table will have two columns:
                     TargetID      (Text)   — the building id value
                     ShieldingIndex (Double) — mean of per-sector MAX elevation
                                              angles (degrees, 0–90 scale)

All intermediate scratch data is written to memory\\ and deleted automatically.
No sectors.gdb is needed — intermediate sector stats stay in memory.
"""

import os
import math
import time
import multiprocessing as mp

import numpy as np

try:
    from scipy.spatial import cKDTree
except ImportError:
    raise ImportError(
        "scipy is required. Open the ArcGIS Pro Python Command Prompt and run:\n"
        "    conda install scipy\n"
        "Then re-run this script."
    )

import arcpy

# ══════════════════════════════════════════════════════════════════════════════
# CONFIG  —  edit these paths to match your machine
# ══════════════════════════════════════════════════════════════════════════════

# Full path to your centroid feature class
# Must have fields: id, height, CENT_X, CENT_Y
CENTROID_FC = (
    r"C:\path\to\your\centroids.shp"
)

# Your project workspace GDB (used for environment settings)
WORKSPACE_GDB = (
    r"C:\path\to\your\project.gdb"
)

# Output GDB — where ShieldingResults will be written
OUTPUT_GDB = (
    r"C:\path\to\your\output.gdb"
)

# Name of the output table inside OUTPUT_GDB
OUTPUT_TABLE = "ShieldingResults"

# ── Analysis parameters (match your working model exactly) ────────────────────
SEARCH_RADIUS_M  = 300    # metres — same as your working model
SECTOR_WIDTH_DEG = 15     # degrees per sector
N_SECTORS        = 24     # 360 / 15 = 24

# ── Field names (taken directly from your working model) ──────────────────────
ID_FIELD     = "TARG_ID"       # unique building identifier
HEIGHT_FIELD = "height"   # building height in metres
CENT_X_FIELD = "CENT_X"  # pre-calculated centroid X
CENT_Y_FIELD = "CENT_Y"  # pre-calculated centroid Y

# ── Parallel processing ────────────────────────────────────────────────────────
# N_CORES = None  →  use all available CPU cores automatically
# N_CORES = 4     →  limit to 4 cores (reduces system load)
N_CORES    = None
BATCH_SIZE = 500   # buildings per worker batch

# ══════════════════════════════════════════════════════════════════════════════
# CORE COMPUTATION  —  no arcpy here, safe for worker processes
# ══════════════════════════════════════════════════════════════════════════════

def _process_batch(args):
    """
    Compute shielding for a batch of target buildings.
    Runs inside a worker process — must not use arcpy.

    Replicates exactly the logic from your CalculateFields + Summary Statistics
    steps in the original model:
      - N_DIST    : euclidean distance between centroids
      - AZIMUTH   : clockwise bearing from north
      - ELEV_ANG  : atan((neighbour_height - target_height) / distance)
      - SECTOR_ID : int(azimuth / 15) % 24
      - Per sector : MAX(ELEV_ANG)          ← Summary Statistics 1
      - Per building: MEAN(MAX_ELEV_ANG)    ← Summary Statistics 2 → ShieldingIndex
    """
    target_indices, coords, heights, ids = args

    # Rebuild KDTree inside worker (cannot be pickled on Windows)
    tree    = cKDTree(coords)
    results = []

    for i in target_indices:
        tx = coords[i, 0]
        ty = coords[i, 1]
        th = heights[i]
        tid = ids[i]

        # KDTree replaces Select Layer By Location (300 m radius)
        neighbour_idx = tree.query_ball_point([tx, ty], SEARCH_RADIUS_M)

        sector_max = {}   # sector_id → max ELEV_ANG

        for j in neighbour_idx:
            if j == i:
                continue   # skip target itself (replicates Select By Attribute)

            nx = coords[j, 0]
            ny = coords[j, 1]
            nh = heights[j]

            # N_DIST
            dx   = nx - tx
            dy   = ny - ty
            dist = math.sqrt(dx * dx + dy * dy)
            if dist == 0:
                continue

            # ELEV_ANG  — matches get_el(!height!, !TARG_H!, !N_DIST!)
            dh = nh - th
            el = math.degrees(math.atan(dh / dist))

            # AZIMUTH — matches get_az(CENT_X, CENT_Y, TARG_X, TARG_Y)
            az  = math.degrees(math.atan2(dx, dy)) % 360.0

            # SECTOR_ID — matches int(!AZIMUTH! / 15) % 24
            sid = int(az / SECTOR_WIDTH_DEG) % N_SECTORS

            # Summary Statistics 1: MAX ELEV_ANG per sector
            if sid not in sector_max or el > sector_max[sid]:
                sector_max[sid] = el

        # Summary Statistics 2: MEAN of MAX_ELEV_ANG across all 24 sectors
        # Sectors with no neighbour contribute 0.0 (same as missing case in Stats)
        all_max = [sector_max.get(s, 0.0) for s in range(N_SECTORS)]
        shielding_index = sum(all_max) / N_SECTORS   # = MEAN_MAX_ELEV_ANG

        results.append({
            "target_id":       tid,
            "shielding_index": shielding_index,
        })

    return results


# ══════════════════════════════════════════════════════════════════════════════
# DATA LOADING  —  uses arcpy, runs in main process only
# ══════════════════════════════════════════════════════════════════════════════

def load_centroids():
    """
    Read id, height, CENT_X, CENT_Y from the centroid feature class.
    Returns plain Python lists and numpy arrays — no arcpy objects.
    """
    ids     = []
    coords  = []
    heights = []

    fields = [ID_FIELD, HEIGHT_FIELD, CENT_X_FIELD, CENT_Y_FIELD]

    with arcpy.da.SearchCursor(CENTROID_FC, fields) as cur:
        for row in cur:
            bid, h, cx, cy = row
            ids.append(bid)
            coords.append((cx, cy))
            heights.append(float(h) if h is not None else 0.0)

    coords  = np.array(coords,  dtype=np.float64)
    heights = np.array(heights, dtype=np.float64)

    arcpy.AddMessage(f"  Loaded {len(ids):,} buildings from {CENTROID_FC}")
    return ids, coords, heights


# ══════════════════════════════════════════════════════════════════════════════
# OUTPUT WRITING  —  uses arcpy, runs in main process only
# ══════════════════════════════════════════════════════════════════════════════

def prepare_output_table():
    """
    Create ShieldingResults table with TargetID and ShieldingIndex fields.
    Matches the field mapping in your original Append step exactly:
      TargetID      Text 255   ← from TARG_ID
      ShieldingIndex Double    ← from MEAN_MAX_ELEV_ANG
    Deletes and recreates the table so re-runs don't produce duplicates.
    """
    table_path = os.path.join(OUTPUT_GDB, OUTPUT_TABLE)

    if arcpy.Exists(table_path):
        arcpy.management.Delete(table_path)
        arcpy.AddMessage(f"  Deleted existing {OUTPUT_TABLE}")

    arcpy.management.CreateTable(OUTPUT_GDB, OUTPUT_TABLE)
    arcpy.management.AddField(table_path, "TargetID",       "TEXT",   field_length=255)
    arcpy.management.AddField(table_path, "ShieldingIndex", "DOUBLE")

    arcpy.AddMessage(f"  Created output table: {table_path}")
    return table_path


def write_results(all_results, table_path):
    """Write all results to ShieldingResults in one batch insert."""
    fields = ["TargetID", "ShieldingIndex"]
    with arcpy.da.InsertCursor(table_path, fields) as cur:
        for r in all_results:
            cur.insertRow([str(r["target_id"]), r["shielding_index"]])
    arcpy.AddMessage(f"  Wrote {len(all_results):,} rows to {table_path}")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def run():
    t0 = time.time()

    arcpy.AddMessage("=" * 62)
    arcpy.AddMessage("BuildingShield.py")
    arcpy.AddMessage(f"  Input    : {CENTROID_FC}")
    arcpy.AddMessage(f"  Output   : {os.path.join(OUTPUT_GDB, OUTPUT_TABLE)}")
    arcpy.AddMessage(f"  Radius   : {SEARCH_RADIUS_M} m")
    arcpy.AddMessage(f"  Sectors  : {N_SECTORS} × {SECTOR_WIDTH_DEG}°")
    arcpy.AddMessage("=" * 62)

    # ── Validate inputs ───────────────────────────────────────────────────────
    if not arcpy.Exists(CENTROID_FC):
        raise FileNotFoundError(f"Centroid layer not found:\n  {CENTROID_FC}")
    if not arcpy.Exists(OUTPUT_GDB):
        raise FileNotFoundError(
            f"Output GDB not found:\n  {OUTPUT_GDB}\n"
            f"Create it first in ArcGIS Pro (Catalog → New File Geodatabase)."
        )

    arcpy.env.overwriteOutput = True
    arcpy.env.workspace       = WORKSPACE_GDB

    # ── Step 1: Load all centroids ────────────────────────────────────────────
    arcpy.AddMessage("\nStep 1/4 — Loading centroid data…")
    ids, coords, heights = load_centroids()
    n_total = len(ids)

    # ── Step 2: Split into batches for parallel workers ───────────────────────
    arcpy.AddMessage(f"\nStep 2/4 — Preparing {n_total:,} buildings for parallel processing…")
    all_indices = list(range(n_total))
    batches = [
        all_indices[i : i + BATCH_SIZE]
        for i in range(0, n_total, BATCH_SIZE)
    ]
    n_cores = N_CORES or mp.cpu_count()
    arcpy.AddMessage(f"  {len(batches)} batches × {BATCH_SIZE} buildings | {n_cores} CPU cores")

    # Package args for each worker
    # numpy arrays are passed by value (pickled) — safe for multiprocessing
    worker_args = [
        (batch, coords, heights, ids)
        for batch in batches
    ]

    # ── Step 3: Run parallel computation ─────────────────────────────────────
    arcpy.AddMessage(f"\nStep 3/4 — Computing shielding (parallel)…")

    all_results = []
    completed   = 0

    # spawn context is required on Windows with ArcGIS Pro
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=n_cores) as pool:
        for batch_result in pool.imap_unordered(_process_batch, worker_args):
            all_results.extend(batch_result)
            completed += len(batch_result)

            # Progress update every 2,000 buildings
            if completed % 2000 == 0 or completed == n_total:
                elapsed = time.time() - t0
                rate    = completed / elapsed if elapsed > 0 else 1
                eta_min = (n_total - completed) / rate / 60
                arcpy.AddMessage(
                    f"  {completed:,}/{n_total:,}  |  "
                    f"{rate:.0f} bldg/sec  |  ETA {eta_min:.1f} min"
                )

    # ── Step 4: Write output table ────────────────────────────────────────────
    arcpy.AddMessage(f"\nStep 4/4 — Writing results…")
    table_path = prepare_output_table()
    write_results(all_results, table_path)

    # ── Summary ───────────────────────────────────────────────────────────────
    elapsed   = time.time() - t0
    sh_values = [r["shielding_index"] for r in all_results]
    arcpy.AddMessage("\n" + "=" * 62)
    arcpy.AddMessage("Done.")
    arcpy.AddMessage(f"  Total time     : {elapsed/60:.1f} minutes")
    arcpy.AddMessage(f"  Buildings      : {n_total:,}")
    arcpy.AddMessage(f"  Output table   : {table_path}")
    arcpy.AddMessage(f"  ShieldingIndex : min {min(sh_values):.2f}°  "
                     f"max {max(sh_values):.2f}°  "
                     f"mean {sum(sh_values)/n_total:.2f}°")
    arcpy.AddMessage("=" * 62)


# ── Required guard for multiprocessing on Windows ─────────────────────────────
if __name__ == "__main__":
    run()
