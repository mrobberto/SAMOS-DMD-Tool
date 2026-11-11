# SAMOS Target Tool — Chapter 1 preparation (v1.4.7)

**Working documentation, 11 November 2026**
**Scope:** Standalone, offline preparation of SAMOS target catalogs, slit configurations, finding charts, and region files. This tool does **not** command the telescope, DMD, motors, or detectors.

This README consolidates `README_baseline.md`, `README_FIRST.txt`, and `README_1_4_7.txt` and replaces their superseded instructions. Historical READMEs may be retained in `Archive/`, but this file is the current entry point. The SAMOS Observer Manual, Chapter 1, provides the detailed observing rationale.

## 1. Package layout

From the `SAMOS_Target_Tool_v1_4_7/` directory:

```text
samos_prepare.py          Main preparation command; partial-run recovery
samos_target_tool/         Optimizer, catalog access, mask finalizer, survey fetch
examples/                  Example virgin input catalogs, including Dolidze25
output/                    Generated runs (each field directory is immutable)
test_samos_prepare.py      Six regression tests with built-in unittest
pyproject.toml             Python package/dependency metadata
README.md                  This consolidated guide
Archive/                   Historical scripts and documentation, if retained
```

Generated science catalogs and region files are *not* input to the instrument until validated against a real WCS image and the live SAMOS GUI.

## 2. Python environment

Python 3.10+ is required; the working installation uses Python 3.11 in the existing `samos` conda environment. With that environment active, from this directory:

```bash
python -m pip install -e .
python test_samos_prepare.py
```

`test_samos_prepare.py` uses built-in `unittest`; **pytest is not required**. SkyMapper catalog retrieval uses `pandas.read_html()` and requires `lxml`. If missing:

```bash
conda install -c conda-forge lxml
python -c "import lxml.etree; print('lxml OK')"
```

A `pandas` warning about parsing XML as HTML can occur during SkyMapper retrieval; this warning alone is not a fatal error. Network access is needed for online Gaia / photometry / image downloads.

## 3. Required observer inputs

Supply a **virgin**, never-overwritten CSV with one science target per row:

```text
target_id,RA,DEC,priority,Gmag
```

`target_id` is persistent; `RA` and `DEC` are decimal degrees; priority 1 is highest. `Gmag` is recommended for the ordering of repeat targets. The field name (`--name`) and pointing-reference star (`--pointing-star RA_DEG DEC_DEG`) are mandatory. The star defines the center of the nominal 3' × 3' planning field and its coordinates appear in exported region filenames.

Do not enter preassigned T-mask numbers or slit IDs in the virgin catalog. Keep the original CSV unchanged.

## 4. Standard Dolidze 25 run

```bash
python samos_prepare.py \
  --input examples/Dolidze25/input.csv \
  --name Dolidze25 \
  --pointing-star 101.26963621981825 0.22980340015201453 \
  --slit-width 3 \
  --slit-length 11 \
  --trace-gap 1 \
  --n-masks 5 \
  --n-sky 5 \
  --survey auto \
  --band i \
  --outdir output
```

The geometry flags are optional; their defaults are **3**, **11**, and **1** planning mirrors respectively. Width and length can also be entered using the aliases `--slit-width-pix` and `--slit-length-pix`; `--trace-gap` also accepts `--min-trace-gap-mirrors`. The generated field name becomes `Dolidze25_3x11`; an old geometry suffix supplied through `--name` is replaced rather than duplicated.

At the current planning scale of **0.18 arcsec per mirror**, the 3 × 11 slit has projected dimensions **0.540 arcsec × 1.980 arcsec**. A one-mirror gap corresponds to **0.180 arcsec** and sky/filler trace-center separation **(11 + 1) × 0.18 = 2.160 arcsec**. The optimizer receives integer geometry in planning pixels; the finalizer receives the corresponding angular geometry.

**Important limitation:** The verified v1.4.7 wrapper passes the configurable trace gap to the **sky/filler finalizer**, but the science-target optimizer has no independently verified gap argument. Its spectral-collision model still needs examination. A 1-mirror gap must not be interpreted as measured clearance for all science spectra on the detector.

### Useful options

| Argument | Purpose / default |
| --- | --- |
| `--n-masks` | Number of T patterns; default 5 |
| `--n-sky` | Candidate blank-sky apertures per ALL mask; default 5 |
| `--max-fillers` | Maximum ancillary Gaia filler apertures; default 10 |
| `--max-repeats` | Maximum repeat apertures on a later mask; default 100 |
| `--no-repeat-targets` | Disable repeat-science apertures |
| `--n-trials`, `--seed` | Optimizer trials (20) and random seed (42) |
| `--priority-column` | Input priority column; default `priority` |
| `--filler-min-ra-sep` | Separate filler selection parameter; default 1.5 |
| `--spectral-length-pix` | Assumed spectral footprint length; default 2883 |
| `--survey` | `auto` (SkyMapper first), `skymapper`, or `panstarrs` |
| `--band` | Survey FITS passband; default `i` |
| `--epoch` | Optional epoch passed to Gaia query |
| `--reference-csv` | Existing Gaia reference cache for offline/cached runs |
| `--fits-image` | Existing real survey FITS with celestial WCS |
| `--skip-photometry` | Skip live survey photometry query |
| `--skip-survey-image` | Skip survey FITS image download |

With `--survey auto`, photometric catalog retrieval initially uses SkyMapper; the **image** download tries SkyMapper and then Pan-STARRS if SkyMapper fails. The two paths are not identical fallbacks.

### Cached/offline preparation

To reuse a Gaia reference and WCS FITS without downloading a new survey image:

```bash
python samos_prepare.py \
  --input examples/Dolidze25/input.csv \
  --name Dolidze25 \
  --pointing-star 101.26963621981825 0.22980340015201453 \
  --reference-csv /path/to/gaia_reference.csv \
  --fits-image /path/to/real_wcs_image.fits \
  --skip-photometry \
  --outdir output_offline
```

A fully offline run without a FITS image can use `--skip-survey-image`, but no on-image astrometric validation is then possible. A cached catalog without epoch information cannot automatically correct proper motions.

## 5. Outputs and naming

All generated products live under `output/<field>_<width>x<length>/`, for example:

```text
output/Dolidze25_3x11/
  catalogs/
    gaia_reference.csv
    gaia_dr3_full.csv                     # live Gaia queries only
    Dolidze25_3x11-GAIA_in_field_101.2696362_+0.2298034.csv
    Dolidze25_3x11-SkyMapper_in_field_101.2696362_+0.2298034.csv
    Dolidze25_3x11-SkyMapper_in_field_101.2696362_+0.2298034_i.fits
  regions/
    Dolidze25_3x11-T00-all_RADEC=101.2696362+0.2298034.reg
    Dolidze25_3x11-T00-EVEN_RADEC=101.2696362+0.2298034.reg
    Dolidze25_3x11-T00-ODD_RADEC=101.2696362+0.2298034.reg
    ... (T01, T02, ...)
  figures/                                # when figures are generated
  manifests/
    target_assignments.csv
    slit_manifest.csv
    source_observation_index.csv
    mask_summary.csv
    multiplex_summary_3x11.csv
    region_wcs_validation_3x11.csv       # if real WCS FITS supplied
```

The listing gives an illustrative **online SkyMapper** run: cached/offline runs and a Pan-STARRS fallback have different subsets and survey names. The `in_field` convention applies to Gaia, SkyMapper, and Pan-STARRS GUI/survey products.

The GUI-format Gaia CSV columns are exactly **`MainID,ra_now,dec_now,Gmag`**. Gaia source IDs are preserved as complete strings and records sorted from brightest to faintest (increasing `Gmag`). The source catalog and the manifest use stable target identities, which must not be replaced with CSV row numbers.

The original `target_assignments.csv` records primary target allocation. `slit_manifest.csv` is the authoritative aperture inventory (SCIENCE / REPEAT / FILLER / SKY, T pattern, SLITID, source identity and cross-mask history); `source_observation_index.csv` records planned patterns for each non-sky source. `mask_summary.csv` lists per-mask ALL/EVEN/ODD counts and filenames; `multiplex_summary_3x11.csv` summarizes finalized ALL masks. Planned exposures do not mean executed exposures.

### DS9 region conventions

All sky-region box positions and sizes are **decimal degrees**, including explicit `d` units; angles are in degrees. A 3 × 11 box, for example, is:

```text
fk5
box(101.245996945000d,0.239454486000d,0.000550000000d,0.000150000000d,0.00000000) # color=green text={SLIT000}
```

The horizontal length is 0.000550 deg = 1.980 arcsec (11 planning mirrors); vertical width is 0.000150 deg = 0.540 arcsec (3 planning mirrors). Text labels, DS9 headers and comments are preserved; legal quoted strings such as `font="helvetica 10 normal roman"` are not angular units. Pixel-coordinate region files retain pixel units.

T patterns (T00, T01, ...) are independent observations. Each has one authoritative ALL list. EVEN and ODD are **exact partitions** of ALL, with global slit IDs within that T pattern assigned by increasing RA: `SLIT000, SLIT001, ...`; EVEN IDs 0,2,...; ODD 1,3,... . Slits are not renumbered when split. A unique aperture key is `(TMASK, SLITID)`, not SLITID alone. The source ID (e.g., Gaia ID) stays unchanged when a star is observed in a later mask.

For successive masks, the optimizer allocates new science targets first; the finalizer adds requested sky candidates and, for T01 onward, eligible repeat science targets, followed by Gaia fillers. Repeat apertures retain their original source identity. Color coding and full source association should be checked against the manifest, not guessed from a region label.

## 6. Validation before transfer

1. Run `python test_samos_prepare.py` (six regression tests; no `pytest`).
2. Check `slit_manifest.csv` and `multiplex_summary_<geometry>.csv` for source IDs, counts and required sky candidates.
3. Verify every celestial DS9 box uses degree units, the selected slit dimensions and correct orientation.
4. Confirm ALL equals the union of EVEN and ODD, preserving slit numbers and geometry.
5. Inspect the real survey image: target alignment, guide star, contaminants, and sky candidates. The WCS roundtrip report checks coordinate geometry, **not** whether boxes coincide with the intended stars.
6. Review schematic spectral-stripe occupancy; it is not a calibrated SAMI detector trace/collision model.
7. Only after inspection, transfer final `regions/` plus needed Gaia/survey catalogs to the SOAR workstation. Keep manifests, plots and original input with run documentation.

Do **not** send masks to the DMD solely on the basis of an offline WCS roundtrip. The GUI must apply and validate the current image WCS and preserve a new version when observers adjust slit positions.

### Dolidze 25 observed preparation result (11 Nov 2026)

With 3 × 11 geometry, 1-mirror sky/filler gap, five masks and five candidate sky slits per ALL mask, the **actual** run reported 100 input targets, 96 in field and **93 newly allocated**, distributed **T00 39, T01 25, T02 16, T03 10, T04 3**. The finalizer output reported **81 repeat apertures**, **19 Gaia fillers** and **25 candidate sky apertures** across all five ALL masks. These are planned apertures, not measured detector spectra or executed exposures. The T00 ALL mask was checked for degree-only geometry and EVEN/ODD partition consistency, but live astrometric registration still requires inspection.

## 7. Partial-run recovery and no-overwrite policy

Fresh runs **refuse to overwrite** a generated field directory, even if the earlier run failed. Choose a new `--outdir` for a clean rerun, or use the dedicated recovery switch only when all finalizer outputs already exist and the post-processing step failed:

```bash
python samos_prepare.py \
  --resume-existing \
  --input examples/Dolidze25/input.csv \
  --name Dolidze25 \
  --pointing-star 101.26963621981825 0.22980340015201453 \
  --slit-width 3 --slit-length 11 --trace-gap 1 \
  --n-masks 5 --n-sky 5 \
  --outdir output_retry
```

Recovery **does not rerun** the optimizer, Gaia, survey photometry/image download, or finalizer; it only converts/checks existing regions and generates summaries. It refuses to resume an output that already has a finalized `multiplex_summary_<geometry>.csv`. The distributed v1.4.7 recovery branch sets `backup_original=True`, creating `regions_original_before_degree_conversion/`; a local change to `backup_original=False` disables this extra copy, if preferred. That directory is not created during ordinary fresh runs.

## 8. Version notes and known limitations

- **v1.4.5:** Primary allocation, repeats, blank-sky candidates, fillers, and ALL/EVEN/ODD products; historical baseline.
- **v1.4.6:** Parameterized geometry and sky/filler trace gap; `_3x11` and `_in_field_` naming; Gaia ID/brightness export; degree-only celestial region geometry; WCS/source-ID audits. First real run exposed an overly broad double-quote check in DS9 processing.
- **v1.4.7:** Limits quote validation to numerical region geometry, keeps header/comment quotes, adds preflight of all region files and `--resume-existing`. Six `unittest` regression checks; real Dolidze 25 recovery and T00 region consistency subsequently completed on the working system.

**Pending acceptance checks:** (a) exact science-target spectral trace clearance under the specified gap; (b) intended-star centering and blank-sky integrity using a genuine, verified image; (c) live GUI WCS and mask-loading behavior. These are not demonstrated by the unit tests.
