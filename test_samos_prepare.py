#!/usr/bin/env python3
"""Offline regression tests for samos_prepare.py, using standard-library unittest.

No pytest, telescope connection, or installed samos_target_tool required.
Pandas is required by samos_prepare.py itself.

Run: python test_samos_prepare.py
"""
import argparse
import importlib.util
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

import pandas as pd

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('ch1_samos_prepare', ROOT / 'samos_prepare.py')
app = importlib.util.module_from_spec(spec)
spec.loader.exec_module(app)


class SAMOSPrepareTests(unittest.TestCase):
    def test_names(self):
        self.assertEqual(app.geometry_name('Dolidze25', 3, 11), ('Dolidze25', 'Dolidze25_3x11'))
        self.assertEqual(app.geometry_name('Dolidze25_3x11', 5, 13), ('Dolidze25', 'Dolidze25_5x13'))
        self.assertEqual(
            app.catalog_filename('Dolidze25_3x11', 'gaia', 101.2696362, 0.2298034),
            'Dolidze25_3x11-GAIA_in_field_101.2696362_+0.2298034.csv')
        self.assertIn('PanStarrs_in_field', app.catalog_filename('Dolidze25_3x11', 'panstarrs', 0, -36))

    def test_gaia(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'catalog.csv'
            table = pd.DataFrame({
                'SOURCE_ID': pd.Series(['3125505670478616192', '3125505670478616193'], dtype='string'),
                'ra_now': [101.27, 101.28], 'dec_now': [0.23, 0.24],
                'phot_g_mean_mag': [19, 11],
            })
            result = app.export_gaia_in_field(table, path)
            self.assertEqual(list(result.columns), ['MainID', 'ra_now', 'dec_now', 'Gmag'])
            self.assertEqual(result['MainID'].tolist(),
                             ['3125505670478616193', '3125505670478616192'])
            self.assertEqual(pd.read_csv(path, dtype={'MainID': str})['MainID'].tolist(),
                             result['MainID'].tolist())
            with self.assertRaisesRegex(ValueError, 'floating-point'):
                app.export_gaia_in_field(
                    table.assign(SOURCE_ID=[3.125505670478616e18, 3.125505670478617e18]), path)

    def test_degrees(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'Dolidze25_3x11-T00-all_RADEC=101.2+0.2.reg'
            path.write_text('# Region file format: DS9\n'
                            'global color=green font="helvetica 10 normal roman"\n'
                            'fk5\n'
                            'box(101.2696362198,0.2298034,1.980",0.540",35) # text={"guide star, box(0,0,0)"} tag={SCIENCE}\n'
                            'circle(101.270,0.230,1.00") # color=cyan\n')
            centers = app.rewrite_ds9_as_degrees(path, width_arcsec=0.54, length_arcsec=1.98)
            result = path.read_text()
            self.assertEqual(len(centers), 2)
            # A quotation mark in a DS9 font definition or annotation is NOT
            # an arcsecond dimension and must survive normalization.
            self.assertIn('font="helvetica 10 normal roman"', result)
            self.assertIn('text={"guide star, box(0,0,0)"}', result)
            for line in result.splitlines():
                numeric_core = line.split('#', 1)[0]
                if 'box(' in numeric_core or 'circle(' in numeric_core:
                    self.assertNotIn('"', numeric_core)
            self.assertIn('0.000550000000d', result)
            self.assertIn('0.000150000000d', result)
            self.assertIn('tag={SCIENCE}', result)
            self.assertIn('35.00000000', result)
            # Rewriting an already-canonical file should not change it.
            app.rewrite_ds9_as_degrees(path, width_arcsec=0.54, length_arcsec=1.98)
            self.assertEqual(path.read_text(), result)
            pixel = Path(directory) / 'Dolidze25_3x11-T00-all_PIXELS.reg'
            pixel.write_text('image\nbox(10,20,11,3,35) # tag={SCIENCE}\n')
            self.assertEqual(app.rewrite_ds9_as_degrees(pixel,
                             width_arcsec=0.54, length_arcsec=1.98), [])
            self.assertIn('11,3', pixel.read_text())

    def test_bad_geometry_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'Dolidze25_3x11-T00-all.reg'
            body = 'fk5\nbox(101d,0d,2.98",0.54",0)\n'
            path.write_text(body)
            with self.assertRaisesRegex(ValueError, 'do not match'):
                app.rewrite_ds9_as_degrees(path, width_arcsec=0.54, length_arcsec=1.98)
            self.assertEqual(path.read_text(), body)

    def test_resume_existing_partial_output_without_network(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            field = root / 'output_retry' / 'Dolidze25_3x11'
            region_dir = field / 'regions'
            man = field / 'manifests'
            region_dir.mkdir(parents=True)
            man.mkdir(parents=True)
            catalog = root / 'input.csv'
            catalog.write_text('target_id,RA,DEC,priority\nA,101.27,0.23,1\n')
            pd.DataFrame([
                {'TMASK': 'T00', 'SLITTYPE': 'SCIENCE', 'target_id': 'A'},
                {'TMASK': 'T00', 'SLITTYPE': 'SKY', 'target_id': ''},
            ]).to_csv(man / 'slit_manifest.csv', index=False)
            region = region_dir / 'Dolidze25_3x11-T00-EVEN_RADEC=101.269636+0.229803.reg'
            original = ('# Region file format: DS9\n'
                        'global color=green font="helvetica 10 normal roman"\n'
                        'fk5\n'
                        'box(101.269636,0.229803,1.980",0.540",0) '
                        '# text={"SCIENCE target"} tag={SCIENCE}\n')
            region.write_text(original)
            argv = ['--resume-existing', '--input', str(catalog),
                    '--name', 'Dolidze25', '--pointing-star', '101.269636', '0.229803',
                    '--n-masks', '1', '--n-sky', '1',
                    '--outdir', str(root / 'output_retry')]
            self.assertEqual(app.main(argv), 0)
            backup = field / 'regions_original_before_degree_conversion' / region.name
            self.assertFalse(backup.exists())
            self.assertIn('0.000550000000d', region.read_text())
            self.assertIn('font="helvetica 10 normal roman"', region.read_text())
            self.assertTrue((man / 'multiplex_summary_3x11.csv').exists())
            with self.assertRaisesRegex(RuntimeError, 'finalized summary'):
                app.main(argv)

    def test_full_wrapper_mock(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp_path = Path(directory)
            folder = tmp_path / 'out'
            catalog = tmp_path / 'input.csv'
            catalog.write_text('target_id,RA,DEC,priority\nA,101.27,0.23,1\n')
            calls = []

            def cli_main(args):
                calls.append(('optimizer', args))
                manifest_dir = Path(args[args.index('--outdir') + 1])
                manifest_dir.mkdir(parents=True, exist_ok=True)
                name = args[args.index('--target-name') + 1]
                (manifest_dir / f'{name}_target_tool_assignments.csv').write_text('target_id,TMASK\nA,T00\n')
                return 0

            def fin_parser():
                parser = argparse.ArgumentParser()
                for key in ('assignments', 'reference', 'pointing-star', 'name', 'outdir',
                            'n-sky', 'max-fillers', 'max-repeats', 'filler-min-ra-sep',
                            'spectral-length-pix', 'slit-width-arcsec', 'slit-length-arcsec',
                            'sky-trace-separation-arcsec', 'filler-trace-sep', 'fits-image'):
                    parser.add_argument('--' + key, nargs=2 if key == 'pointing-star' else None)
                parser.add_argument('--no-repeat-targets', action='store_true')
                return parser

            def fin_run(args):
                calls.append(('finalizer', args))
                region_dir = Path(args.outdir)
                region_dir.mkdir(parents=True, exist_ok=True)
                pd.DataFrame([
                    {'TMASK': 'T00', 'SLITTYPE': 'SCIENCE', 'target_id': 'A'},
                    {'TMASK': 'T00', 'SLITTYPE': 'SKY', 'target_id': ''},
                ]).to_csv(region_dir / 'slit_manifest.csv', index=False)
                (region_dir / f'{args.name}-T00-all_RADEC=101.269636+0.229803.reg').write_text(
                    'global color=green font="helvetica 10 normal roman"\nfk5\n'
                    'box(101.269636,0.229803,2.340",0.720",0) # text={"target A"} tag={SCIENCE}\n')

            def gaia_local(path, ra, dec):
                return pd.DataFrame({
                    'SOURCE_ID': pd.Series(['3125505670478616192'], dtype='string'),
                    'ra_now': [ra], 'dec_now': [dec], 'phot_g_mean_mag': [11.2],
                })

            root = types.ModuleType('samos_target_tool')
            cli = types.ModuleType('samos_target_tool.cli')
            cli.main = cli_main
            fin = types.ModuleType('samos_target_tool.mask_finalize')
            fin.run = fin_run
            fin.parser = fin_parser
            refs = types.ModuleType('samos_target_tool.reference_catalogs')
            refs.gaia_local = gaia_local
            survey = types.ModuleType('samos_target_tool.survey_fetch')
            root.cli = cli
            root.mask_finalize = fin
            root.reference_catalogs = refs
            root.survey_fetch = survey
            fake_modules = {
                'samos_target_tool': root,
                'samos_target_tool.cli': cli,
                'samos_target_tool.mask_finalize': fin,
                'samos_target_tool.reference_catalogs': refs,
                'samos_target_tool.survey_fetch': survey,
            }
            cache = tmp_path / 'reference.csv'
            cache.write_text('target_id,RA,DEC,Gmag\n3125505670478616192,101.27,0.23,11.2\n')
            argv = [
                '--input', str(catalog), '--reference-csv', str(cache),
                '--name', 'Dolidze25_3x11', '--pointing-star', '101.2696362', '0.2298034',
                '--slit-width', '4', '--slit-length', '13', '--trace-gap', '2',
                '--n-masks', '1', '--n-sky', '1', '--skip-photometry', '--skip-survey-image',
                '--outdir', str(folder),
            ]
            with patch.dict(sys.modules, fake_modules):
                self.assertEqual(app.main(argv), 0)
                dest = folder / 'Dolidze25_4x13'
                self.assertTrue((dest / 'catalogs' /
                    'Dolidze25_4x13-GAIA_in_field_101.2696362_+0.2298034.csv').exists())
                regions = list((dest / 'regions').glob('*.reg'))
                self.assertEqual(len(regions), 1)
                self.assertIn('4x13', regions[0].name)
                self.assertIn('font="helvetica 10 normal roman"', regions[0].read_text())
                self.assertIn('0.000650000000d', regions[0].read_text())
                self.assertIn('0.000200000000d', regions[0].read_text())
                self.assertTrue((dest / 'manifests' / 'multiplex_summary_4x13.csv').exists())
                opt_args = calls[0][1]
                self.assertEqual(opt_args[opt_args.index('--slit-width-pix') + 1], '4')
                self.assertEqual(opt_args[opt_args.index('--slit-length-pix') + 1], '13')
                fin_args = calls[1][1]
                self.assertAlmostEqual(float(fin_args.sky_trace_separation_arcsec), (13 + 2) * .18)
                with self.assertRaises(FileExistsError):
                    app.main(argv)


if __name__ == '__main__':
    unittest.main(verbosity=2)
