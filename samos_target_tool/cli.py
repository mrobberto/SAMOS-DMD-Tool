from __future__ import annotations

import argparse
from pathlib import Path

from .catalogs import query_gaia_dr3, read_user_catalog, resolve_simbad
from .config import SAMOSGeometry, TargetToolConfig
from .geometry import add_field_coordinates
from .optimizer import MaskOptimizer
from .ranking import assign_priority_bins
from .regions import write_outputs


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Prepare SAMOS DMD slit masks from a target catalog.")
    p.add_argument("--target-name", required=True)
    center = p.add_mutually_exclusive_group(required=True)
    center.add_argument("--center-simbad", help="Resolve field center with SIMBAD")
    center.add_argument("--center-radec", nargs=2, type=float, metavar=("RA_DEG", "DEC_DEG"))

    cat = p.add_mutually_exclusive_group(required=True)
    cat.add_argument("--catalog-csv", help="CSV catalog with RA and DEC columns in decimal degrees")
    cat.add_argument("--catalog", choices=["gaia"], help="Online catalog to query")

    p.add_argument("--query-radius-arcmin", type=float, default=2.0)
    p.add_argument("--outdir", default=None)
    p.add_argument("--n-masks", type=int, default=5)
    p.add_argument("--n-trials", type=int, default=20)
    p.add_argument("--random-seed", type=int, default=42)
    p.add_argument("--max-targets", type=int, default=1000)
    p.add_argument("--sort-column", default="Gmag", help="Column used to auto-create priorities")
    p.add_argument("--priority-column", default=None, help="Existing priority column; 1 is highest")
    p.add_argument("--reverse-sort", action="store_true", help="Use when larger sort values are higher priority")
    p.add_argument("--field-size-arcsec", type=float, default=180.0)
    p.add_argument("--slit-width-pix", type=int, default=3)
    p.add_argument("--slit-length-pix", type=int, default=9)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.center_simbad:
        center = resolve_simbad(args.center_simbad)
        ra_center = float(center.ra.deg)
        dec_center = float(center.dec.deg)
    else:
        ra_center, dec_center = map(float, args.center_radec)

    if args.catalog_csv:
        raw = read_user_catalog(args.catalog_csv)
    else:
        raw = query_gaia_dr3(ra_center, dec_center, radius_arcmin=args.query_radius_arcmin)

    geometry = SAMOSGeometry(
        field_width_arcsec=args.field_size_arcsec,
        field_height_arcsec=args.field_size_arcsec,
        slit_width_pix=args.slit_width_pix,
        slit_length_pix=args.slit_length_pix,
    )
    tool_config = TargetToolConfig(
        n_masks=args.n_masks,
        n_trials=args.n_trials,
        random_seed=args.random_seed,
        max_targets=args.max_targets,
        sort_column=args.sort_column,
        priority_column=args.priority_column,
        lower_sort_is_better=not args.reverse_sort,
    )

    table = add_field_coordinates(raw, ra_center, dec_center, geometry)
    if args.priority_column:
        table = table.rename(columns={args.priority_column: "priority"})
    else:
        table = assign_priority_bins(
            table,
            sort_column=args.sort_column,
            lower_is_better=not args.reverse_sort,
            output_column="priority",
        )

    optimizer = MaskOptimizer(geometry, tool_config)
    assigned, summary = optimizer.optimize(table)

    outdir = Path(args.outdir or args.target_name)
    write_outputs(assigned, outdir, args.target_name, ra_center, dec_center, geometry)

    print("SAMOS Target Tool summary")
    print("-------------------------")
    print(f"Target       : {args.target_name}")
    print(f"Center       : RA={ra_center:.8f} deg  DEC={dec_center:+.8f} deg")
    print(f"Input rows   : {summary.n_input}")
    print(f"In field     : {summary.n_in_field}")
    print(f"Assigned     : {summary.n_assigned}")
    print(f"By mask      : {summary.assigned_by_mask}")
    print(f"By priority  : {summary.assigned_by_priority}")
    print(f"Output dir   : {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
