from __future__ import annotations

from pathlib import Path
import pandas as pd

from .config import SAMOSGeometry


def write_outputs(
    assigned_table: pd.DataFrame,
    outdir: str | Path,
    target_name: str,
    ra_center_deg: float,
    dec_center_deg: float,
    geometry: SAMOSGeometry,
) -> None:
    """Write optimizer assignments; slit regions are generated only at finalization."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    assigned_table.to_csv(outdir / f"{target_name}_target_tool_assignments.csv", index=False)

    # Intermediate region files are intentionally not emitted.
    # The finalized ALL/EVEN/ODD masks are the only instrument regions.


def write_ds9_slit_regions(df: pd.DataFrame, path: str | Path, geometry: SAMOSGeometry) -> None:
    """Write a DS9 ICRS region file using SAMOS slit dimensions.

    Width is along RA on the sky, matching the orientation change in the
    prototype's Save_table_2_Astropy_slitregions method.
    """
    path = Path(path)
    slit_ra_arcsec = geometry.slit_length_pix * geometry.image_scale_arcsec_per_pix
    slit_dec_arcsec = geometry.slit_width_pix * geometry.image_scale_arcsec_per_pix

    with path.open("w") as f:
        f.write("# Region file format: DS9 version 4.1\n")
        f.write("global color=green dashlist=8 3 width=1 font='helvetica 10 normal roman' select=1 highlite=1\n")
        f.write("icrs\n")
        for _, row in df.sort_values(["RA", "DEC"], ascending=[True, True], kind="stable").iterrows():
            label = row.get("target_id", "")
            f.write(
                f"box({float(row['RA']):.8f},{float(row['DEC']):.8f},"
                f"{slit_ra_arcsec:.3f}\",{slit_dec_arcsec:.3f}\",0) # text={{{label}}}\n"
            )
