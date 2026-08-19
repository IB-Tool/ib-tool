# IB-Tool

**Automatic delineation of the *Innenbereich* — the coherently built-up part of a municipality (§ 34 BauGB) — from building footprints, road network and land use data.** The international publication of the method describes the same delineation as an *Urban Growth Boundary (UGB)*.

The method combines semantic and spatial building filtering, block-based building coverage ratio (BCR) analysis, MST-based aggregation of building clusters, and boundary refinement against the road network. It produces a fine-grained boundary that follows individual buildings rather than administrative units.

> Harig, O.; Hecht, R.; Burghardt, D.; Meinel, G. **Automatic Delineation of Urban Growth Boundaries Based on Topographic Data Using Germany as a Case Study.** *ISPRS Int. J. Geo-Inf.* **2021**, *10*(5), 353. <https://doi.org/10.3390/ijgi10050353>

---

## Repositories

**Start here → [IB-Tool 3](https://github.com/IB-Tool/IB-Tool-3)** — the current implementation.

| Repository | What it is | Platform | Status |
|---|---|---|---|
| **[IB-Tool 3](https://github.com/IB-Tool/IB-Tool-3)** | The delineation itself: filtering, density analysis, MST aggregation, boundary refinement. Guided 4-step dialog, GeoPackage output. | QGIS 3.40+ | active · v0.2.2 |
| **[Data Wizard](https://github.com/IB-Tool/data_wizard)** | Turns raw ATKIS Basis-DLM shapefiles into the `HU` / `RN` / `Aux` GeoPackages IB-Tool 3 expects — CRS handling, study-area clipping, layer mapping. | QGIS 3.40+ | experimental · v0.1 |
| **[Partitioning](https://github.com/IB-Tool/Partitioning)** | Divides the study area into settlement units (`PART_<id>`) via kernel density estimation and Voronoi tessellation, so large datasets can be processed partition by partition. | QGIS 3.40+ | experimental · v0.1 |
| **[ib-tool](https://github.com/IB-Tool/ib-tool)** *(this repository)* | The original ArcPy toolset from the 2021 publication. Superseded by IB-Tool 3, kept for reference and for ArcGIS users. | ArcGIS Advanced + Spatial Analyst | frozen · 2021 |

---

## Typical workflow

```
raw ATKIS Basis-DLM  ──►  Data Wizard  ──►  HU / RN / Aux
                                              │
building footprints (HU)  ──►  Partitioning  ──►  Part
                                              │
                                              ▼
                                        IB-Tool 3  ──►  Innenbereich / UGB (GeoPackage)
```

Neither companion plugin is required: any data meeting the specification in [IB-Tool 3 → `docs/input-data.md`](https://github.com/IB-Tool/IB-Tool-3/blob/master/docs/input-data.md) works, regardless of how it was produced. All three QGIS plugins install side by side and appear in the same **Plugins → IB-Tool** menu.

## Which version do I need?

- **You work in QGIS** → [IB-Tool 3](https://github.com/IB-Tool/IB-Tool-3). Start with [`docs/quickstart.md`](https://github.com/IB-Tool/IB-Tool-3/blob/master/docs/quickstart.md); a sample dataset is included.
- **You work in ArcGIS** → this repository. Requires ArcGIS Advanced with a Spatial Analyst license; see [ARCGIS-TOOLSET.md](https://github.com/IB-Tool/ib-tool/blob/master/ARCGIS-TOOLSET.md).

---

## Documentation

| Topic | Where |
|---|---|
| Installation, first run, sample data | [IB-Tool 3 → `docs/quickstart.md`](https://github.com/IB-Tool/IB-Tool-3/blob/master/docs/quickstart.md) |
| Input layer specifications and validation checks | [IB-Tool 3 → `docs/input-data.md`](https://github.com/IB-Tool/IB-Tool-3/blob/master/docs/input-data.md) |
| Algorithm, pseudocode, accuracy results | [IB-Tool 3 → `docs/how-it-works.md`](https://github.com/IB-Tool/IB-Tool-3/blob/master/docs/how-it-works.md) |
| Parameters and how to tune them | [IB-Tool 3 → `docs/parameterization.md`](https://github.com/IB-Tool/IB-Tool-3/blob/master/docs/parameterization.md) |
| Terminology (*Innenbereich* vs. UGB) and which source covers what | [IB-Tool 3 → `docs/terminology.md`](https://github.com/IB-Tool/IB-Tool-3/blob/master/docs/terminology.md) |
| Development, CI and release conventions (shared across all three plugins) | [IB-Tool 3 → `docs/contributing.md`](https://github.com/IB-Tool/IB-Tool-3/blob/master/docs/contributing.md) |
| Method documentation in German (Harig 2024) | [`DOKU.md`](https://github.com/IB-Tool/ib-tool/blob/master/DOKU.md) |
| ArcGIS toolset (this repository) | [`ARCGIS-TOOLSET.md`](https://github.com/IB-Tool/ib-tool/blob/master/ARCGIS-TOOLSET.md) |

---

## This repository: the ArcGIS toolset (2021)

Three independent scripts for ArcGIS Advanced (ArcMap or ArcGIS Pro) with a Spatial Analyst license:

- `Partitioning.pyt` — creates the partitioning from building footprints
- `IB-Tool2.py` + `IB-Tool2_Config.txt` + `IB-Tool2_Filter.txt` — creates the delineation
- `Error_Classification.py` — assesses the result against an expert delineation

Full prerequisites, data preparation, filter and parameter reference: **[ARCGIS-TOOLSET.md](https://github.com/IB-Tool/ib-tool/blob/master/ARCGIS-TOOLSET.md)**. The method behind it is documented in German in [`DOKU.md`](https://github.com/IB-Tool/ib-tool/blob/master/DOKU.md).

This toolset is no longer developed. New work happens in [IB-Tool 3](https://github.com/IB-Tool/IB-Tool-3).

---

## Citation

**Method:**

> Harig, O.; Hecht, R.; Burghardt, D.; Meinel, G. Automatic Delineation of Urban Growth Boundaries Based on Topographic Data Using Germany as a Case Study. *ISPRS Int. J. Geo-Inf.* **2021**, *10*(5), 353. <https://doi.org/10.3390/ijgi10050353>

**Original toolset (this repository):**

> Harig, O. (2021). *Toolset for the delineation of settlements on the basis building footprints, road network and land use data* (v1.0). <https://doi.org/10.26084/IOERFDZ-SOFT-001>

The parameterisation and the grounding in § 34 BauGB are covered by Harig (2024); further sources and what each one covers are listed in [IB-Tool 3 → `docs/terminology.md`](https://github.com/IB-Tool/IB-Tool-3/blob/master/docs/terminology.md#references--which-source-covers-what).

---

## License

- **ib-tool** (ArcGIS toolset): MIT — see [`LICENSE.md`](https://github.com/IB-Tool/ib-tool/blob/master/LICENSE.md)
- **IB-Tool 3, Data Wizard, Partitioning** (QGIS plugins): GPL-2.0-or-later

## Author

**Oliver Harig** — conception and implementation
Conception (2021 publication): Robert Hecht, Dirk Burghardt, Gotthard Meinel
