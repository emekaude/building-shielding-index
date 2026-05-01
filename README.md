# Building Shielding Index for Wind Load Characterization

[![ArcGIS Pro](https://img.shields.io/badge/ArcGIS_Pro-3.x-blue)](https://www.esri.com/en-us/arcgis/products/arcgis-pro/overview)
[![Python](https://img.shields.io/badge/Python-3.x-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

---

## Overview

This project computes a **directional building shielding index** for every building in an urban dataset. The shielding index quantifies how much a given building is sheltered from wind by its surrounding neighbours, expressed as the mean of per-sector maximum elevation angles across 24 azimuthal sectors.

The index is designed as one factor in a **wind load characterization framework**, providing a spatially explicit shielding modifier that can be applied to reference wind speeds. It is intended as a first step toward supporting research and validating the shielding correction factors prescribed in structural engineering codes (e.g. Eurocode 1, ASCE 7).

> **Note:** Building shielding is one of several factors relevant to wind load assessment. This tool isolates and quantifies the geometric shielding component from surrounding buildings. It does not replace a full wind engineering analysis.

---

## Background and Motivation

Wind load design in structural engineering relies on code-based exposure and shielding modifiers that account for surrounding obstructions. These modifiers are typically applied in a generalised way, using broad terrain categories rather than building-specific geometric context.

In dense urban environments, the actual shielding experienced by a building depends strongly on:

- the height of surrounding buildings relative to the target building,
- the proximity and directional distribution of those neighbours,
- the presence of gaps or sheltered corridors in the urban fabric.

This tool operationalises those factors into a single per-building index derived from real building footprint and height data, enabling spatial analysis of wind exposure across entire urban portfolios.

---

## Methodology

### Input data

- Building footprint polygons with a unique ID field and a height field (metres)
- Pre-computed centroid coordinates (`CENT_X`, `CENT_Y`) in a projected CRS (UTM or equivalent)

### Workflow

For each target building the model:

1. Identifies all neighbouring building centroids within a configurable search radius (default 300 m)
2. Removes the target building from its own neighbour set
3. Computes for each neighbour:
   - **Distance** — Euclidean distance between centroids
   - **Azimuth** — clockwise bearing from north from target to neighbour
   - **Elevation angle** — `arctan((H_neighbour − H_target) / distance)`
4. Bins neighbours into 24 azimuthal sectors of 15° each
5. Computes the **maximum elevation angle** per sector (the controlling obstruction in each direction)
6. Computes the **mean of sector maxima** as the final shielding index

### Shielding index interpretation

| Value | Meaning |
|---|---|
| Negative (°) | Target building is taller than most neighbours — exposed |
| Near 0° | Target and neighbours are of similar height — neutral |
| Positive (°) | Target building is shorter than surrounding neighbours — shielded |
| Higher positive (°) | Strongly shielded — tall dense neighbours in many directions |

The index is in degrees (0–90 scale) and can be normalised to 0–1 by dividing by 90.

### Coordinate reference system

The tool requires a **projected CRS in metres**. The reference implementation uses WGS 1984 UTM Zone 32N (EPSG:32632) for the Padova municipality dataset.

---

## Repository Structure

```
building-shielding-index/
│
├── ModelBuilder/
│   └── BuildingShield.py          # ArcGIS ModelBuilder export (original workflow)
│
├── BuildingShield_Fast.py         # Optimised Python implementation
│
├── README.md                      # This file
└── LICENSE
```

---

## Implementations

This repository contains two implementations of the same analytical workflow.

### 1. ModelBuilder (`ModelBuilder/BuildingShield.py`)

The original workflow built in ArcGIS Pro ModelBuilder and exported as a Python script. This version uses native ArcGIS geoprocessing tools chained sequentially:

- Iterate Feature Selection → Make Feature Layer → Select Layer By Location → Select Layer By Attribute → Copy Features → Calculate Fields → Summary Statistics × 2 → Append

Suitable for smaller datasets or for users who want to inspect and modify the workflow visually inside ModelBuilder.

**Limitation:** Processes one building at a time. At approximately 1 building per second, processing time scales linearly with dataset size (~18 hours for 89,000 buildings).

---

### 2. Fast Python Script (`BuildingShield_Fast.py`)

A ground-up rewrite of the same analytical logic with two major performance optimisations:

**KDTree spatial index (`scipy.spatial.cKDTree`)**
Replaces ArcGIS Select Layer By Location. All building centroids are indexed once at startup. Finding all neighbours within the search radius becomes a single tree query per building rather than a full spatial scan, reducing the core computation by approximately 99%.

**Multiprocessing (`multiprocessing` — spawn context)**
Buildings are split into batches and processed simultaneously across all available CPU cores. On a 16-core machine this adds a further ~15× speedup on top of the KDTree gains.

**Benchmark — Padova municipality dataset (89,466 buildings)**

| Implementation | Time |
|---|---|
| ModelBuilder | ~18 hours (estimated) |
| Fast Python (16 cores) | **0.7 minutes** |

The fast script produces identical output to the ModelBuilder workflow.

---

## Requirements

- ArcGIS Pro 3.x (provides `arcpy`, `numpy`, `scipy`)
- Python 3.x (bundled with ArcGIS Pro)
- `scipy` — ships with ArcGIS Pro's conda environment

To verify scipy is available, open the ArcGIS Pro Python Command Prompt and run:
```
python -c "from scipy.spatial import cKDTree; print('OK')"
```

If it is missing:
```
conda install scipy
```

---

## Usage

### Prepare your data

Before running the script, your building centroid layer must have the following fields pre-calculated:

| Field | Type | Description |
|---|---|---|
| `id` | Text / Integer | Unique building identifier |
| `height` | Double | Building height in metres |
| `CENT_X` | Double | Centroid X coordinate (projected CRS, metres) |
| `CENT_Y` | Double | Centroid Y coordinate (projected CRS, metres) |

`CENT_X` and `CENT_Y` can be calculated in ArcGIS Pro using **Calculate Geometry Attributes** on a centroid feature class derived from your building footprints using **Feature To Point**.

---

### Running the fast script

**Step 1 — Edit the CONFIG block** at the top of `BuildingShield_Fast.py`:

```python
CENTROID_FC   = r"C:\path\to\your\centroids.shp"   # or .gdb feature class
WORKSPACE_GDB = r"C:\path\to\your\project.gdb"
OUTPUT_GDB    = r"C:\path\to\your\output.gdb"
OUTPUT_TABLE  = "ShieldingResults"

SEARCH_RADIUS_M  = 300    # neighbour search radius in metres
SECTOR_WIDTH_DEG = 15     # degrees per azimuthal sector
N_SECTORS        = 24     # must equal 360 / SECTOR_WIDTH_DEG

ID_FIELD     = "id"
HEIGHT_FIELD = "height"
CENT_X_FIELD = "CENT_X"
CENT_Y_FIELD = "CENT_Y"

N_CORES    = None   # None = use all available CPU cores
BATCH_SIZE = 500
```

**Step 2 — Open the ArcGIS Pro Python Command Prompt**

Start Menu → search **Python Command Prompt** → open it.

> Do not use the Python Window embedded inside ArcGIS Pro. Multiprocessing requires a real terminal.

**Step 3 — Run**

```bash
python "C:\path\to\BuildingShield_Fast.py"
```

---

### Output

The script writes a table called `ShieldingResults` inside your output GDB with one row per building:

| Field | Type | Description |
|---|---|---|
| `TargetID` | Text | Building unique ID |
| `ShieldingIndex` | Double | Mean of sector-maximum elevation angles (degrees) |

This table can be joined back to your building footprint layer in ArcGIS Pro using `TargetID` = your building ID field for spatial visualisation and further analysis.

---

## Configuration reference

| Parameter | Default | Description |
|---|---|---|
| `SEARCH_RADIUS_M` | 300 | Radius in metres within which neighbours are considered |
| `SECTOR_WIDTH_DEG` | 15 | Angular width of each azimuthal sector |
| `N_SECTORS` | 24 | Number of sectors (must equal 360 / SECTOR_WIDTH_DEG) |
| `N_CORES` | None (all) | Number of CPU cores to use |
| `BATCH_SIZE` | 500 | Buildings per worker batch |

Increasing `SEARCH_RADIUS_M` increases computation time because more neighbours are processed per building. The relationship is roughly quadratic — doubling the radius approximately quadruples the number of neighbours.

---

## Data source

The building footprint and height data used in this project are sourced from the **GlobalBuildingAtlas** dataset developed by the [EO Data Science group](https://www.asg.ed.tum.de/en/sipeo/home/) at the Technical University of Munich (TUM).

### Citation

If you use this project with the GlobalBuildingAtlas dataset, please cite the original data publication:

```bibtex
@Article{essd-17-6647-2025,
  AUTHOR  = {Zhu, X. X. and Chen, S. and Zhang, F. and Shi, Y. and Wang, Y.},
  TITLE   = {GlobalBuildingAtlas: an open global and complete dataset of
             building polygons, heights and LoD1 3D models},
  JOURNAL = {Earth System Science Data},
  YEAR    = {2025},
  DOI     = {10.5194/essd-17-6647-2025}
}
```

**Data repository:**
Zhu, X.X., Chen, S., Zhang, F., Shi, Y., Wang, Y. (2025). GlobalBuildingAtlas. Technical University of Munich (TUM). [https://doi.org/10.14459/2025mp1782307](https://doi.org/10.14459/2025mp1782307)

**GitHub:** [https://github.com/zhu-xlab/GlobalBuildingAtlas](https://github.com/zhu-xlab/GlobalBuildingAtlas)

### Data license notice

The GlobalBuildingAtlas dataset is released under a **Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)** license.

> ⚠️ **This means any work derived from or validated against the GlobalBuildingAtlas data — including results produced by this tool — is restricted to non-commercial use.** If you intend to use outputs for commercial purposes, you must obtain data from an alternatively licensed source.

The ODbL-licensed polygon component (GBA.ODbLPolygon) is separately available under the Open Data Commons Open Database License (ODbL) v1.0, which permits broader reuse subject to share-alike conditions.

For full licensing details refer to the [GlobalBuildingAtlas GitHub repository](https://github.com/zhu-xlab/GlobalBuildingAtlas).

---

## Tested dataset

- **Area:** Padova municipality, Italy
- **Buildings:** 89,466
- **CRS:** WGS 1984 UTM Zone 32N (EPSG:32632)
- **Search radius:** 300 m / 800 m
- **Hardware:** 16-core CPU
- **Result:** ShieldingIndex range −10.15° to 14.05°, mean 1.72° (300 m radius)

---

## Limitations and future work

- The model uses **centroid-to-centroid distance** rather than building envelope geometry. This is a simplification that introduces some positional error for large footprints.
- **Terrain is assumed flat.** Buildings on sloping ground may produce incorrect elevation angles if heights are above-ground rather than absolute elevation.
- The shielding index is a **geometric proxy** for wind sheltering. It does not account for aerodynamic effects such as channelling, wake turbulence, or pressure coefficients.
- Future work could incorporate directional wind rose weighting to produce a wind-climate-adjusted shielding index, combining the geometric modifier with local prevailing wind frequency distributions.

---

## Citation

If you use this tool in your research, please cite this repository and reference the methodology described in the README.

---

## License

MIT License — see [LICENSE](LICENSE) for details.
