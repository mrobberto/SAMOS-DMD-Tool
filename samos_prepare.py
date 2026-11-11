#!/usr/bin/env python3
"""SAMOS Chapter 1 preparation v1.4.7 — configurable slit geometry and degree-only masks.

A revision of the Chapter 1 ``samos_prepare_3x11_fixed`` wrapper.  No
installed samos_target_tool module is modified.  The science optimizer and
mask finalizer receive matching slit dimensions, while the finalizer's sky and
filler collision separations also include the requested trace gap.

Example (default geometry 3 x 11, gap 1 planning mirror)::

  python samos_prepare.py \
  --resume-existing \
  --input examples/Dolidze25/input.csv \
  --name Dolidze25 \
  --pointing-star 101.26963621981825 0.22980340015201453 \
  --slit-width 3 \
  --slit-length 11 \
  --trace-gap 1 \
  --n-masks 5 \
  --n-sky 5 \
  --outdir output

Alternative geometry::

  python samos_prepare.py ... --slit-width 4 --slit-length 13 --trace-gap 2

Outputs are immutable under <outdir>/<target>_<width>x<length>/, with
<target>_<width>x<length> in mask names.  Gaia and photometry catalogs use
<target>_<width>x<length>-<survey>_in_field_<RA>_<Dec>.csv.

Limitations: the existing optimizer CLI only exposes width and length, not a
separate trace-gap parameter.  --trace-gap controls finalizer sky/filler
separations; verify science-trace spacing against your optimizer implementation.
"""

import argparse
from collections import Counter
from pathlib import Path
import math
import re
import shutil
import sys

import pandas as pd

# The imports below are intentionally deferred to main(): pure formatting
# tests and --help work even without the installed telescope target tool.

SLIT_WIDTH = 3              # default, integer planning mirrors (dispersion)
SLIT_LENGTH = 11            # default, integer planning mirrors (cross-dispersion)
TRACE_GAP = 1               # default, clear planning mirrors between traces
PLANNING_SCALE = 0.18       # arcsec / planning mirror; matches v1.4.5 geometry

SKY_FRAMES = {'fk5', 'icrs', 'j2000'}
PIXEL_FRAMES = {'image', 'physical'}
REGION_SHAPES = {'box', 'circle', 'ellipse'}
BOX_RE = re.compile(r'(?i)(?<![\w])(?P<shape>box|circle|ellipse)\s*\((?P<args>[^()]*)\)')
OLD_GEOMETRY_RE = re.compile(r'_(\d+)x(\d+)$', re.IGNORECASE)


def geometry_name(field: str, slit_width: int, slit_length: int) -> tuple[str, str]:
    """Accept either 'Dolidze25' or legacy 'Dolidze25_3x11' as --name."""
    base = OLD_GEOMETRY_RE.sub('', field.strip())
    if not base or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.+-]*', base):
        raise ValueError('--name must contain only path-safe characters and no spaces')
    return base, f'{base}_{slit_width}x{slit_length}'


def field_tag(ra_deg: float, dec_deg: float) -> str:
    """Filename-safe field coordinates (no spaces or colon delimiters)."""
    return f'{ra_deg:.7f}_{dec_deg:+.7f}'


def catalog_filename(field_with_geometry: str, survey: str,
                     ra_deg: float, dec_deg: float) -> str:
    canonical = {'gaia': 'GAIA', 'skymapper': 'SkyMapper',
                 'panstarrs': 'PanStarrs'}[survey.lower()]
    return f'{field_with_geometry}-{canonical}_in_field_{field_tag(ra_deg, dec_deg)}.csv'


def _canonical_source_ids(series: pd.Series) -> pd.Series:
    """Protect 19-digit Gaia IDs against conversion to floating point."""
    # Once a 19-digit ID has gone through float64 the lost bits are
    # irrecoverable.  Refuse to export a seemingly valid but corrupted ID.
    if pd.api.types.is_float_dtype(series.dtype):
        raise ValueError('Gaia SOURCE_ID is floating-point; full 64-bit identifiers '
                         'may be corrupted. Re-read the reference CSV with ID as string.')
    ids = series.astype('string').str.strip()
    if ids.isna().any() or not ids.str.fullmatch(r'\d+').all():
        raise ValueError('Gaia SOURCE_ID must be a nonempty integer identifier')
    return ids


def export_gaia_in_field(g: pd.DataFrame, destination: Path) -> pd.DataFrame:
    """Original GUI-compatible Gaia schema, complete IDs, bright stars first."""
    mandatory = ('SOURCE_ID', 'ra_now', 'dec_now', 'phot_g_mean_mag')
    missing = [key for key in mandatory if key not in g.columns]
    if missing:
        raise ValueError(f'gaia_local returned a catalog missing {missing}')
    out = pd.DataFrame({
        'MainID': _canonical_source_ids(g['SOURCE_ID']),
        'ra_now': pd.to_numeric(g['ra_now'], errors='raise'),
        'dec_now': pd.to_numeric(g['dec_now'], errors='raise'),
        'Gmag': pd.to_numeric(g['phot_g_mean_mag'], errors='coerce'),
    })
    if (not out['ra_now'].map(math.isfinite).all()
            or not out['dec_now'].map(math.isfinite).all()
            or not out['ra_now'].between(0, 360, inclusive='left').all()
            or not out['dec_now'].between(-90, 90).all()):
        raise ValueError('Gaia catalog contains invalid coordinates')
    if out['MainID'].duplicated().any():
        raise ValueError('Gaia catalog contains duplicated source identifiers')
    out = out.sort_values('Gmag', kind='stable', na_position='last').reset_index(drop=True)
    destination.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(destination, index=False, float_format='%.12f')
    # Explicit roundtrip check against pandas' float-based ID inference.
    chk = pd.read_csv(destination, dtype={'MainID': 'string'})
    if chk['MainID'].tolist() != out['MainID'].tolist():
        raise RuntimeError('Gaia identifier changed during CSV serialization')
    return out


def _sexagesimal(token: str, *, ra: bool) -> float:
    """Read an existing DS9 sexagesimal coordinate and emit degrees."""
    sign = -1 if token.startswith('-') else 1
    parts = token.lstrip('+-').split(':')
    if len(parts) != 3:
        raise ValueError(f'Unsupported sexagesimal coordinate: {token!r}')
    components = [float(p) for p in parts]
    if components[1] >= 60 or components[2] >= 60:
        raise ValueError(f'Invalid sexagesimal coordinate: {token!r}')
    return sign * (components[0] + components[1] / 60 + components[2] / 3600) * (15 if ra else 1)


def degree_value(raw: str, *, coordinate: str = '', size: bool = False) -> float:
    """Read a DS9 celestial token, handling degrees, arcsec, arcmin and HMS."""
    token = raw.strip()
    if ':' in token:
        if size or coordinate not in ('ra', 'dec'):
            raise ValueError(f'Unexpected sexagesimal value {token!r}')
        return _sexagesimal(token, ra=coordinate == 'ra')
    if token.endswith('"'):
        multiplier, token = 1 / 3600, token[:-1]
    elif token.endswith("'"):
        multiplier, token = 1 / 60, token[:-1]
    elif token.lower().endswith('deg'):
        multiplier, token = 1, token[:-3]
    elif token.lower().endswith('d'):
        multiplier, token = 1, token[:-1]
    else:
        multiplier = 1  # DS9 celestial unqualified numbers are degrees
    value = float(token) * multiplier
    if not math.isfinite(value):
        raise ValueError(f'Nonfinite region value: {raw!r}')
    return value


def _d(value: float) -> str:
    return f'{value:.12f}d'


def _shape_line_to_degrees(line: str, *, expected_width: float,
                           expected_length: float, source: str) -> tuple[str, list[tuple[float, float]]]:
    """Convert sky-region parameters to *explicit* decimal degrees.

    Preserve trailing DS9 comments/tags verbatim; only numeric args change.
    Boxes follow the finalizer's length-by-width convention.  Reject unexpected
    geometries instead of silently changing scientific apertures.
    """
    positions = []

    def sub(match: re.Match) -> str:
        shape = match.group('shape').lower()
        items = [x.strip() for x in match.group('args').split(',')]
        counts = {'box': 5, 'circle': 3, 'ellipse': 5}
        if len(items) != counts[shape]:
            raise ValueError(f'{source}: expected {counts[shape]} {shape} args: {line}')
        ra = degree_value(items[0], coordinate='ra')
        dec = degree_value(items[1], coordinate='dec')
        if not (0 <= ra < 360 and -90 <= dec <= 90):
            raise ValueError(f'{source}: out-of-range RA/Dec: {ra}, {dec}')
        positions.append((ra, dec))
        sizes = [degree_value(raw, size=True) for raw in items[2:-1 if shape != 'circle' else None]]
        if shape == 'circle':
            if sizes[0] <= 0:
                raise ValueError(f'{source}: invalid circle radius')
            return f'circle({_d(ra)},{_d(dec)},{_d(sizes[0])})'
        if any(value <= 0 for value in sizes):
            raise ValueError(f'{source}: nonpositive region dimension')
        angle = degree_value(items[-1])
        if shape == 'box':
            # Original finalizer geometry: DS9 box major axis = slit length,
            # minor axis = slit width.  Never assume the swapped geometry.
            tol = 0.00051 / 3600  # legacy validation tolerance, in degrees
            if abs(sizes[0] - expected_length) > tol or abs(sizes[1] - expected_width) > tol:
                raise ValueError(f'{source}: box dimensions '
                                 f'{sizes[0]*3600:.5f}" x {sizes[1]*3600:.5f}" '
                                 'do not match the requested slit geometry')
        result = ','.join([_d(ra), _d(dec)] + [_d(v) for v in sizes]
                          + [f'{angle:.8f}'])
        return f'{shape}({result})'

    new_line, count = BOX_RE.subn(sub, line)
    if count == 0 and re.search(r'(?i)(?<!\w)[a-z_]+\s*\(', line.split('#', 1)[0]):
        raise ValueError(f'{source}: unsupported region shape in {line!r}')
    return new_line, positions


def rewrite_ds9_as_degrees(reg_path: Path, *, width_arcsec: float,
                           length_arcsec: float, dry_run: bool = False) -> list[tuple[int, float, float]]:
    """Convert a celestial DS9 region file in place without changing identities."""
    lines = reg_path.read_text(encoding='utf-8').splitlines()
    active_frame = None
    saw_pixel = False
    positions = []
    output = []
    for lineno, line in enumerate(lines, start=1):
        core = line.split('#', 1)[0].strip()
        lowered = core.lower()
        if not core or core.startswith('#') or lowered.startswith('global '):
            output.append(line)
            continue
        # DS9 supports standalone frame declarations and inline 'fk5;box(...)'.
        section = core.split(';', 1)
        frame = section[0].strip().lower()
        if frame in SKY_FRAMES | PIXEL_FRAMES:
            active_frame = frame
            if len(section) == 1:
                output.append(line)
                continue
        if active_frame in PIXEL_FRAMES:
            saw_pixel = True
            output.append(line)  # pixel region dimensions must remain pixels
            continue
        if active_frame not in SKY_FRAMES:
            raise ValueError(f'{reg_path}:{lineno}: unsupported or missing DS9 sky frame: {core!r}')
        # DS9 annotations and global settings can legitimately contain quotes
        # (e.g. font="helvetica..." or # text={"bright star"}).  Work ONLY
        # on the region expression before #, retaining the metadata verbatim.
        region_core, separator, metadata = line.partition('#')
        converted_core, centers = _shape_line_to_degrees(
            region_core, expected_width=width_arcsec / 3600,
            expected_length=length_arcsec / 3600,
            source=f'{reg_path.name}:{lineno}')
        for match in BOX_RE.finditer(converted_core):
            if '"' in match.group('args') or "'" in match.group('args'):
                raise ValueError(f'{reg_path.name}:{lineno}: a celestial shape '
                                 'still contains arcsecond/arcminute units')
        output.append(converted_core + (separator + metadata if separator else ''))
        positions.extend((lineno, ra, dec) for ra, dec in centers)
    if saw_pixel and positions:
        raise ValueError(f'{reg_path}: mixed sky and pixel coordinate frames not supported')
    if not saw_pixel and not positions:
        raise ValueError(f'{reg_path}: no sky-region shapes found')
    if not saw_pixel:
        new_text = '\n'.join(output) + '\n'
        # Only write after the entire file has passed numerical validation.
        # DS9 header and comment quotes are legal and must be preserved.
        if not dry_run and new_text != reg_path.read_text(encoding='utf-8'):
            temp_path = reg_path.with_suffix('.reg.pending')
            temp_path.write_text(new_text, encoding='utf-8')
            temp_path.replace(reg_path)
    return positions


def validate_wcs_projection(fitsfile: str | None,
                            centers: list[tuple[str, int, float, float]],
                            output: Path) -> None:
    """Check that each sky slit center projects and roundtrips with the FITS WCS."""
    if not fitsfile:
        print('NOTICE: No WCS FITS supplied; sky-to-detector positioning NOT validated.')
        return
    from astropy.io import fits
    from astropy.wcs import WCS
    with fits.open(fitsfile) as hdul:
        hdr = hdul[0].header
        data = hdul[0].data
        if data is None or data.ndim != 2:
            raise ValueError(f'{fitsfile}: WCS validation requires a 2-D FITS image')
        wcs = WCS(hdr).celestial
        if not wcs.has_celestial:
            raise ValueError(f'{fitsfile}: no celestial WCS; cannot validate masks')
        height, width = data.shape
    records = []
    for fname, lineno, ra, dec in centers:
        x, y = wcs.all_world2pix(ra, dec, 0)
        if not math.isfinite(x) or not math.isfinite(y):
            raise ValueError(f'{fname}:{lineno}: WCS projection is nonfinite')
        ra2, dec2 = wcs.all_pix2world(x, y, 0)
        dra = ((ra2 - ra + 180) % 360 - 180) * math.cos(math.radians(dec))
        delta_arcsec = 3600 * math.hypot(dra, dec2 - dec)
        if not math.isfinite(delta_arcsec) or delta_arcsec > 0.01:
            raise ValueError(f'{fname}:{lineno}: WCS roundtrip failed: {delta_arcsec:.5g} arcsec')
        records.append(dict(region_file=fname, line=lineno, ra_deg=ra, dec_deg=dec,
                            fits_x=x, fits_y=y, in_image=0 <= x < width and 0 <= y < height,
                            wcs_roundtrip_arcsec=delta_arcsec))
    if records:
        table = pd.DataFrame(records)
        table.to_csv(output, index=False)
        print(f'WCS positions checked: {len(table)} in {fitsfile}')
        print(f'  Inside FITS image: {int(table.in_image.sum())}/{len(table)}')
        print(f'  Pixel-projection audit: {output}')
        if not table.in_image.any():
            raise ValueError('No region centers project inside the supplied FITS image. '
                             'Likely wrong image, target coordinates, or WCS.')


def _summarize(folder: Path, expected_sky: int, expected_masks: int,
               label: str, virgin_catalog: Path | None = None) -> None:
    path = folder / 'manifests' / 'slit_manifest.csv'
    if not path.is_file():
        raise RuntimeError(f'Cannot calculate finalized counts: missing {path}')
    df = pd.read_csv(path, dtype=str).fillna('')
    required = {'TMASK', 'SLITTYPE'}
    if not required.issubset(df.columns):
        raise RuntimeError(f'Manifest must contain {sorted(required)}; found {list(df.columns)}')
    df['TMASK'] = df['TMASK'].str.strip()
    df['SLITTYPE'] = df['SLITTYPE'].str.strip().str.upper()
    unexpected_types = sorted(set(df.SLITTYPE) - {'SCIENCE', 'REPEAT', 'FILLER', 'SKY'})
    if unexpected_types:
        print('WARNING: unexpected SLITTYPE classifications:', unexpected_types)
    rows = []
    for mask, group in sorted(df.groupby('TMASK'), key=lambda x: x[0]):
        counts = Counter(group['SLITTYPE'])
        rows.append({'TMASK': mask, 'new_science': counts['SCIENCE'],
                     'repeat_science': counts['REPEAT'], 'gaia_fillers': counts['FILLER'],
                     'blank_sky': counts['SKY'], 'total_apertures': len(group)})
    summary = pd.DataFrame(rows)
    if summary.empty:
        raise RuntimeError('Finalized slit manifest contains no apertures')
    names = {f'T{i:02d}' for i in range(expected_masks)}
    seen = set(summary.TMASK)
    if seen != names:
        print(f'NOTICE: Found {len(seen)} nonempty masks, not {expected_masks}. '
              f'Missing: {sorted(names - seen)}')
    dest = folder / 'manifests' / f'multiplex_summary_{label}.csv'
    summary.to_csv(dest, index=False)
    print(f'\nFINALIZED ALL-MASK APERTURE COUNTS ({label}; SKY separate)')
    print(summary.to_string(index=False))
    if 'target_id' in df.columns:
        science = df[df.SLITTYPE.isin({'SCIENCE', 'REPEAT'})]
        science_ids = set(science.target_id.str.strip()) - {''}
        unique = len(science_ids)
        print(f'Unique science target IDs across all T masks: {unique}')
        if (science.target_id.str.strip() == '').any():
            raise RuntimeError('Finalizer lost one or more SCIENCE/REPEAT target IDs')
        if virgin_catalog is not None:
            original = pd.read_csv(virgin_catalog, dtype=str)
            if 'target_id' in original.columns:
                original_ids = set(original.target_id.dropna().str.strip()) - {''}
                missing_ids = science_ids - original_ids
                if missing_ids:
                    raise RuntimeError('Science IDs in final manifest do not match '
                                       f'virgin input catalog: {sorted(missing_ids)[:10]}')
                print(f'Science target identity checked against virgin input: {unique} IDs')
            else:
                print('NOTICE: virgin catalog has no target_id column; '
                      'science target-ID audit skipped')
    else:
        print('NOTICE: manifest has no target_id column; source-ID audit skipped')
    for k in ('new_science', 'repeat_science', 'gaia_fillers', 'blank_sky'):
        print(f'Total {k}: {int(summary[k].sum())}')
    for row in rows:
        if row['blank_sky'] != expected_sky:
            raise RuntimeError(f'{row["TMASK"]} has {row["blank_sky"]} sky apertures, '
                               f'expected {expected_sky}')
    print('Saved:', dest)
    print('NOTE: EVEN/ODD are subsets of ALL; do not add them to these totals.')
    print('NOTE: SKY positions are candidates; verify against real survey images.')



def finalize_region_products(folder: Path, *, slit_width: int, slit_length: int,
                             n_sky: int, n_masks: int, fitsfile: str | None,
                             virgin_catalog: Path, backup_original: bool = False) -> None:
    """Normalize and audit existing finalizer outputs, without re-optimizing.

    A preflight conversion of all files runs before any mutation. Recovery
    from a failed prior run additionally saves every original region file.
    """
    mask = folder / 'regions'
    man = folder / 'manifests'
    label = f'{slit_width}x{slit_length}'
    region_files = sorted(mask.rglob('*.reg'))
    if not region_files:
        raise RuntimeError(f'No DS9 region files found under {mask}')
    width_as = slit_width * PLANNING_SCALE
    length_as = slit_length * PLANNING_SCALE
    for region_path in region_files:
        if label not in region_path.name:
            raise ValueError(f'{region_path.name}: missing {label} in region filename')
        rewrite_ds9_as_degrees(region_path, width_arcsec=width_as,
                               length_arcsec=length_as, dry_run=True)
    if backup_original:
        original_dir = folder / 'regions_original_before_degree_conversion'
        if original_dir.exists():
            raise FileExistsError(f'Original region backup already exists: {original_dir}. '
                                  'Aborting to avoid overwriting the backup.')
        for region_path in region_files:
            relative = region_path.relative_to(mask)
            copy_to = original_dir / relative
            copy_to.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(region_path, copy_to)
        print(f'Original region files preserved in: {original_dir}')
    all_centers = []
    for region_path in region_files:
        centers = rewrite_ds9_as_degrees(region_path, width_arcsec=width_as,
                                         length_arcsec=length_as)
        all_centers.extend((region_path.name, lineno, ra, dec)
                           for lineno, ra, dec in centers)
    print(f'DS9 region validation: {len(all_centers)} celestial regions in '
          f'{len(region_files)} files; numeric sky sizes in decimal degrees '
          '(DS9 header and comments preserved)')
    validate_wcs_projection(fitsfile, all_centers,
                            man / f'region_wcs_validation_{label}.csv')
    _summarize(folder, n_sky, n_masks, label, virgin_catalog)

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input', required=True, help='Virgin science target CSV, never modified')
    p.add_argument('--name', required=True, help='Field name (Dolidze25 or Dolidze25_3x11)')
    p.add_argument('--pointing-star', nargs=2, type=float, required=True,
                   metavar=('RA_DEG', 'DEC_DEG'))
    p.add_argument('--slit-width', '--slit-width-pix', type=int, default=SLIT_WIDTH,
                   dest='slit_width', metavar='MIRRORS', help='Slit width [default 3]')
    p.add_argument('--slit-length', '--slit-length-pix', type=int, default=SLIT_LENGTH,
                   dest='slit_length', metavar='MIRRORS', help='Slit length [default 11]')
    p.add_argument('--trace-gap', '--min-trace-gap-mirrors', type=int, default=TRACE_GAP,
                   dest='trace_gap', metavar='MIRRORS',
                   help='Minimum clear trace gap for finalizer sky/filler placements [default 1]')
    p.add_argument('--outdir', default='output',
                   help='Root output directory; generated field folder is never overwritten')
    p.add_argument('--resume-existing', action='store_true',
                   help='Finish a PARTIAL existing output after a conversion error; '
                        'does not rerun optimizer, Gaia, SkyMapper or finalizer')
    p.add_argument('--n-masks', type=int, default=5)
    p.add_argument('--n-trials', type=int, default=20)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--priority-column', default='priority')
    p.add_argument('--n-sky', type=int, default=5)
    p.add_argument('--max-fillers', type=int, default=10)
    p.add_argument('--max-repeats', type=int, default=100)
    p.add_argument('--no-repeat-targets', action='store_true')
    p.add_argument('--filler-min-ra-sep', type=float, default=1.5)
    p.add_argument('--spectral-length-pix', type=int, default=2883)
    p.add_argument('--reference-csv', help='Cached Gaia reference for offline run')
    p.add_argument('--fits-image', help='Real FITS image with celestial WCS')
    p.add_argument('--survey', choices=['auto', 'skymapper', 'panstarrs'], default='auto')
    p.add_argument('--band', default='i')
    p.add_argument('--skip-survey-image', action='store_true')
    p.add_argument('--skip-photometry', action='store_true')
    p.add_argument('--epoch', type=float)
    a = p.parse_args(argv)
    try:
        field, full_name = geometry_name(a.name, a.slit_width, a.slit_length)
    except ValueError as error:
        p.error(str(error))
    if a.slit_width < 1 or a.slit_length < 1 or a.trace_gap < 0:
        p.error('Slit width/length must be positive; trace gap must be nonnegative')
    if a.n_masks < 1 or a.n_trials < 1 or a.n_sky < 0:
        p.error('--n-masks and --n-trials must be positive; --n-sky nonnegative')
    if not Path(a.input).is_file():
        p.error(f'Input catalog does not exist: {a.input}')
    if a.reference_csv and not Path(a.reference_csv).is_file():
        p.error(f'Reference catalog does not exist: {a.reference_csv}')
    if a.fits_image and not Path(a.fits_image).is_file():
        p.error(f'FITS image does not exist: {a.fits_image}')
    ra, dec = a.pointing_star
    if not (math.isfinite(ra) and math.isfinite(dec) and 0 <= ra < 360 and -90 <= dec <= 90):
        p.error('Invalid RA/Dec in --pointing-star (degrees)')

    width_as = a.slit_width * PLANNING_SCALE
    length_as = a.slit_length * PLANNING_SCALE
    trace_sep_as = (a.slit_length + a.trace_gap) * PLANNING_SCALE

    folder = Path(a.outdir) / full_name
    cat, mask, fig, man = (folder / item for item in
                           ('catalogs', 'regions', 'figures', 'manifests'))
    label = f'{a.slit_width}x{a.slit_length}'
    if a.resume_existing:
        if not folder.is_dir() or not mask.is_dir() or not man.is_dir():
            raise FileNotFoundError(f'Cannot resume: missing existing partial output: {folder}')
        if (man / f'multiplex_summary_{label}.csv').exists():
            raise RuntimeError('This output already has a finalized summary; '
                               'refusing to alter a completed observation plan.')
        if not (man / 'slit_manifest.csv').is_file():
            raise RuntimeError('Cannot resume: finalizer slit_manifest.csv is missing.')
        fitsfile = a.fits_image
        if not fitsfile:
            downloaded = sorted(cat.glob(f'{full_name}-*_in_field_*_{a.band}.fits'))
            if len(downloaded) == 1:
                fitsfile = str(downloaded[0])
                print('Reusing existing survey image for WCS validation:', fitsfile)
            elif len(downloaded) > 1:
                raise RuntimeError('Multiple downloaded FITS found; specify --fits-image.')
            else:
                print('NOTICE: No existing FITS WCS image found for audit.')
        print('Resuming post-finalizer region processing ONLY:', folder)
        finalize_region_products(folder, slit_width=a.slit_width,
                                 slit_length=a.slit_length, n_sky=a.n_sky,
                                 n_masks=a.n_masks, fitsfile=fitsfile,
                                 virgin_catalog=Path(a.input), backup_original=False)
        print('OUTPUT:', folder)
        return 0

    from samos_target_tool.cli import main as optimize_main
    from samos_target_tool.mask_finalize import run as finalize_run, parser as finalize_parser
    from samos_target_tool import reference_catalogs, survey_fetch

    final_parser = finalize_parser()
    supported = {option for action in final_parser._actions for option in action.option_strings}
    required = {'--slit-width-arcsec', '--slit-length-arcsec',
                '--sky-trace-separation-arcsec', '--filler-trace-sep'}
    missing = sorted(required - supported)
    if missing:
        raise RuntimeError('STOP: Installed mask_finalizer lacks expected options: '
                           + ', '.join(missing) + '. No files written.')

    if folder.exists():
        raise FileExistsError(f'Output already exists (never overwrite): {folder}')
    print(f'Field: {full_name}')
    print(f'Slits: {a.slit_width} x {a.slit_length} planning mirrors '
          f'({width_as:.3f}" x {length_as:.3f}")')
    print(f'Min trace gap: {a.trace_gap} mirror(s); sky/filler separation '
          f'{trace_sep_as:.3f} arcsec')
    print('NOTE: Science-mask optimizer accepts dimensions, but has no verified '
          'separate trace-gap option in this wrapper.')
    for d in (cat, mask, fig, man):
        d.mkdir(parents=True, exist_ok=False)

    args = ['--target-name', full_name, '--catalog-csv', a.input,
            '--center-radec', str(ra), str(dec),
            '--n-masks', str(a.n_masks), '--n-trials', str(a.n_trials),
            '--random-seed', str(a.seed), '--outdir', str(man),
            '--slit-width-pix', str(a.slit_width),
            '--slit-length-pix', str(a.slit_length)]
    if a.priority_column:
        args += ['--priority-column', a.priority_column]
    code = optimize_main(args)
    if code:
        raise RuntimeError(f'Optimizer returned error status {code}')
    assignments = man / 'target_assignments.csv'
    legacy = man / f'{full_name}_target_tool_assignments.csv'
    if not legacy.is_file():
        raise FileNotFoundError(f'Missing optimizer assignments: {legacy}')
    legacy.replace(assignments)

    if a.reference_csv:
        reference = pd.read_csv(a.reference_csv,
                                dtype={'target_id': 'string', 'SOURCE_ID': 'string',
                                       'source_id': 'string', 'MainID': 'string'})
        reference_file = cat / 'gaia_reference.csv'
        reference.to_csv(reference_file, index=False)
    else:
        gaia = reference_catalogs.gaia_online(ra, dec, a.epoch)
        reference_file = cat / 'gaia_reference.csv'
        gaia.to_csv(cat / 'gaia_dr3_full.csv', index=False)
        reference = pd.DataFrame({
            'target_id': gaia.source_id.astype('string'),
            'RA': gaia.ra_now, 'DEC': gaia.dec_now,
            'Gmag': gaia.phot_g_mean_mag,
        })
        reference.to_csv(reference_file, index=False)
    g = reference_catalogs.gaia_local(reference_file, ra, dec)
    gaia_output = cat / catalog_filename(full_name, 'gaia', ra, dec)
    clean_gaia = export_gaia_in_field(g, gaia_output)
    print(f'GAIA: {gaia_output} — {len(clean_gaia)} sources (ordered by Gmag)')
    if not a.skip_photometry:
        choice = 'skymapper' if a.survey == 'auto' else a.survey
        phot = (reference_catalogs.skymapper if choice == 'skymapper'
                else reference_catalogs.panstarrs)(ra, dec)
        path = cat / catalog_filename(full_name, choice, ra, dec)
        phot.to_csv(path, index=False)
        print(f'{choice}: {path}')

    fitsfile = a.fits_image
    if not fitsfile and not a.skip_survey_image:
        failures = []
        choices = ['skymapper', 'panstarrs'] if a.survey == 'auto' else [a.survey]
        for choice in choices:
            candidate = cat / (catalog_filename(full_name, choice, ra, dec)
                               .replace('.csv', f'_{a.band}.fits'))
            try:
                (survey_fetch.skymapper if choice == 'skymapper'
                 else survey_fetch.panstarrs)(ra, dec, str(candidate), a.band, 3.5)
                fitsfile = str(candidate)
                break
            except Exception as error:
                failures.append(f'{choice}: {error}')
                if candidate.exists():
                    candidate.unlink()  # remove incomplete failed survey download
        else:
            raise RuntimeError('Real survey image unavailable: ' + '; '.join(failures))

    cmd = ['--assignments', str(assignments), '--reference', str(reference_file),
           '--pointing-star', str(ra), str(dec), '--name', full_name,
           '--outdir', str(mask), '--n-sky', str(a.n_sky),
           '--max-fillers', str(a.max_fillers),
           '--max-repeats', str(a.max_repeats),
           '--filler-min-ra-sep', str(a.filler_min_ra_sep),
           '--spectral-length-pix', str(a.spectral_length_pix),
           '--slit-width-arcsec', str(width_as),
           '--slit-length-arcsec', str(length_as),
           '--sky-trace-separation-arcsec', str(trace_sep_as),
           '--filler-trace-sep', str(trace_sep_as)]
    if a.no_repeat_targets:
        cmd += ['--no-repeat-targets']
    if fitsfile:
        cmd += ['--fits-image', str(fitsfile)]
    finalize_run(final_parser.parse_args(cmd))

    for obsolete in (man / 'slit_manifest_v000.csv', man / 'slit_summary.csv',
                     man / f'{full_name}_target_tool_assignments.csv'):
        if obsolete.exists():
            obsolete.unlink()
    for file in list(mask.glob('*.png')):
        shutil.move(str(file), str(fig / file.name))
    for filename in ('slit_manifest.csv', 'mask_summary.csv',
                     'source_observation_index.csv'):
        file = mask / filename
        if file.exists():
            shutil.move(str(file), str(man / filename))

    finalize_region_products(folder, slit_width=a.slit_width, slit_length=a.slit_length,
                             n_sky=a.n_sky, n_masks=a.n_masks, fitsfile=fitsfile,
                             virgin_catalog=Path(a.input))
    print('\nOUTPUT:', folder)
    return 0


if __name__ == '__main__':
    sys.exit(main())
