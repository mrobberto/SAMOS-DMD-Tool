from __future__ import annotations

import numpy as np
import pandas as pd

from .config import SAMOSGeometry


def standardize_radec_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with RA and DEC columns present.

    Common variants from Gaia, SkyMapper, Pan-STARRS, and user tables are
    accepted. RA/DEC are always decimal degrees.
    """
    out = df.copy()
    aliases = {
        "ra": "RA",
        "RA_ICRS": "RA",
        "raj2000": "RA",
        "dec": "DEC",
        "DE_ICRS": "DEC",
        "dej2000": "DEC",
    }
    for old, new in aliases.items():
        if old in out.columns and new not in out.columns:
            out = out.rename(columns={old: new})
    if "RA" not in out.columns or "DEC" not in out.columns:
        raise ValueError("Input catalog must contain RA and DEC columns in decimal degrees.")
    out["RA"] = pd.to_numeric(out["RA"], errors="coerce")
    out["DEC"] = pd.to_numeric(out["DEC"], errors="coerce")
    out = out.dropna(subset=["RA", "DEC"]).reset_index(drop=True)
    return out


def add_field_coordinates(
    df: pd.DataFrame,
    ra_center_deg: float,
    dec_center_deg: float,
    geometry: SAMOSGeometry,
) -> pd.DataFrame:
    """Add tangent-plane offsets and SAMOS pixel coordinates.

    This replaces the prototype's dependency on a downloaded 2MASS WCS for
    deciding whether a source lies inside the SAMOS 3 arcmin by 3 arcmin field.
    For this small field, the small-angle tangent-plane approximation is
    sufficient for target-mask preselection.
    """
    out = standardize_radec_columns(df)
    cosdec = np.cos(np.deg2rad(dec_center_deg))

    # Positive x_arcsec means increasing RA on the sky. The optimizer uses the
    # same sign convention as the prototype; a later hardware export layer can
    # apply the DMD orientation transform if needed.
    out["x_arcsec"] = (out["RA"] - ra_center_deg) * 3600.0 * cosdec
    out["y_arcsec"] = (out["DEC"] - dec_center_deg) * 3600.0
    out["x_pix"] = out["x_arcsec"] / geometry.image_scale_arcsec_per_pix
    out["y_pix"] = out["y_arcsec"] / geometry.image_scale_arcsec_per_pix

    out["in_field"] = (
        out["x_pix"].abs().le(geometry.half_field_x_pix)
        & out["y_pix"].abs().le(geometry.half_field_y_pix)
    )
    return out
