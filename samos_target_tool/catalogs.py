from __future__ import annotations

from pathlib import Path
import pandas as pd


def read_user_catalog(path: str | Path) -> pd.DataFrame:
    """Read a user-supplied CSV target catalog."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def resolve_simbad(name: str):
    """Resolve a target name through SIMBAD and return an astropy SkyCoord.

    Astroquery/SIMBAD changed its default coordinate column names across
    versions.  Older versions commonly returned ``RA``/``DEC`` sexagesimal
    strings, while newer versions may return lower-case ``ra``/``dec`` in
    decimal degrees.  This resolver accepts both forms.
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.simbad import Simbad

    table = Simbad.query_object(name)
    if table is None or len(table) == 0:
        raise ValueError(f"SIMBAD could not resolve {name!r}")

    colmap = {c.lower(): c for c in table.colnames}

    if "ra" not in colmap or "dec" not in colmap:
        raise KeyError(
            "SIMBAD result does not contain recognizable RA/DEC columns. "
            f"Columns returned were: {table.colnames}"
        )

    ra_col = colmap["ra"]
    dec_col = colmap["dec"]
    ra_val = table[ra_col][0]
    dec_val = table[dec_col][0]

    ra_unit = getattr(table[ra_col], "unit", None)
    dec_unit = getattr(table[dec_col], "unit", None)

    # Newer astroquery often returns RA/DEC as floating-point degrees.
    if ra_unit is not None and dec_unit is not None:
        try:
            return SkyCoord(ra_val * ra_unit, dec_val * dec_unit, frame="icrs")
        except Exception:
            pass

    # Numeric coordinates without explicit units are safest as decimal degrees.
    try:
        return SkyCoord(float(ra_val) * u.deg, float(dec_val) * u.deg, frame="icrs")
    except Exception:
        pass

    # Older astroquery commonly returns sexagesimal RA and DEC strings.
    return SkyCoord(str(ra_val), str(dec_val), unit=(u.hourangle, u.deg), frame="icrs")


def query_gaia_dr3(ra_deg: float, dec_deg: float, radius_arcmin: float = 2.0) -> pd.DataFrame:
    """Query Gaia DR3 around a position and return a pandas table."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.gaia import Gaia

    Gaia.ROW_LIMIT = -1
    center = SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.deg, u.deg), frame="icrs")
    job = Gaia.cone_search_async(center, radius=radius_arcmin * u.arcmin)
    tab = job.get_results()
    df = tab.to_pandas()
    df = df.rename(columns={"ra": "RA", "dec": "DEC", "phot_g_mean_mag": "Gmag", "source_id": "target_id"})
    return df
