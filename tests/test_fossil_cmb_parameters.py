#!/usr/bin/env python3
"""Validate the explicit primary-CMB/Draine parameter path and legacy baseline."""
import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from fossil_parameter_tests import fixture_manifest, generator


def physical_manifest():
    manifest = fixture_manifest()
    manifest['primary_cmb_anisotropy'] = True
    manifest['metadata']['primary_cmb_anisotropy'] = True
    manifest['cmb'] = {'unit_convention': 'SI', 'lmax_master': 383}
    manifest['freefree_model'] = 'draine_gaunt_amplitude_plus_Te'
    manifest['model'] = 'fixture_primary_cmb_draine'
    manifest['signature'] = 'fixture_signature'
    manifest['analysis_signature'] = 'fixture_fullsky_analysis_signature'
    manifest['analysis_mask'] = 'fullsky'
    for nside in ('64', '128'):
        # The original fixture shares a resolution object; detach it first.
        manifest['resolutions'][nside] = copy.deepcopy(manifest['resolutions'][nside])
        resolution = manifest['resolutions'][nside]
        resolution['postprocess_mask'] = resolution['mask']
        resolution['mask'] = resolution['fullsky_mask']
        ff = resolution['components'][2]
        ff['type'] = 'freefree'
        ff['parameters'] = [{'map': f'initial/freefree_temperature_{nside}.fits',
                             'truth_map': f'truth/freefree_temperature_{nside}.fits'}]
        resolution['components'].append(dict(name='cmb', type='cmb', unit='MJy/sr',
            reference_frequency_GHz=50., reference_band=1,
            amplitude=f'initial/cmb_amplitude_{nside}.fits',
            truth_amplitude=f'truth/cmb_amplitude_{nside}.fits', parameters=[]))
        resolution['cmb_prior'] = f'nside{nside}/cmb/prior_Dl_TT.txt'
    return manifest


class CMBParameterTests(unittest.TestCase):
    def test_free_and_fixed_cmb_modes_follow_stage(self):
        manifest = physical_manifest()
        for nside in (64, 128):
            for count in (64, 131):
                for stage, settings in generator.STAGES.items():
                    with self.subTest(nside=nside, count=count, stage=stage):
                        p = generator.build_parameters(manifest, nside, stage, '/inputs', '/chains/test', band_count=count)
                        free = settings['local_amplitudes']
                        self.assertEqual(p['NUM_FG_SIGNAL_COMPONENTS'], 6)
                        self.assertEqual([p[f'FG_LABEL{i:02d}'] for i in range(1, 7)],
                                         ['dust', 'synch', 'freefree', 'ame1', 'ame2', 'cmb'])
                        self.assertEqual(p['COMP_TYPE06'], 'cmb')
                        self.assertEqual(p['SAMPLE_AMPLITUDES06'], free)
                        self.assertEqual(p['ENFORCE_ZERO_CL'], not free)
                        self.assertEqual(p['FOSSIL_CMB_MIN_L'], 2)
                        self.assertTrue(p['FOSSIL_SI_CMB'])
                        self.assertTrue(p['CONDITION_ON_CL'])
                        self.assertTrue(p['SAMPLE_TT_MODES'])
                        self.assertEqual(p['CMB_MONO_DIPOLE_MASK'], 'none')
                        self.assertEqual(p['NUM_JOINT_ALMS_CL_STEPS'], 0)
                        self.assertEqual(p['IN_POWSPECFILE'], f'/inputs/nside{nside}/cmb/prior_Dl_TT.txt')
                        expected = 'initial' if free else 'truth'
                        self.assertTrue(p['INITIAL_AMPLITUDE_MAP06'].startswith(expected + '/'))
                        self.assertEqual(p['FG_UNIT06'], 'MJy/sr')
                        self.assertEqual(p['REFERENCE_FREQUENCY06'], 50)
                        self.assertEqual(p['FREQ_C' + generator.band_suffix(p['REFERENCE_BAND06'])], 50)
                        self.assertFalse(any(key.startswith(('INIT_INDEX_MAP06_', 'FWHM_PAR06_', 'REGION_DEFINITION06_'))
                                             for key in p))
                        # The readable renderer must account for every new key.
                        rendered = generator.render(p, stage, nside)
                        self.assertIn('Primary CMB anisotropy: included', rendered)
                        self.assertIn('Draine free-free', rendered)

    def test_freefree_has_one_temperature_and_no_powerlaw_shape(self):
        manifest = physical_manifest()
        for stage, settings in generator.STAGES.items():
            p = generator.build_parameters(manifest, 64, stage, '/inputs', '/chains/test')
            self.assertEqual(p['COMP_TYPE03'], 'freefree')
            self.assertEqual(p['DEFAULT_T_E03'], 7000)
            self.assertEqual(p['T_E_PRIOR_UNIFORM_LOW03'], 3000)
            self.assertEqual(p['T_E_PRIOR_UNIFORM_HIGH03'], 15000)
            self.assertEqual(p['T_E_PRIOR_GAUSSIAN_STDDEV03'], 500 if settings['local_shapes'] else 0)
            self.assertNotIn('INDEX_PRIOR_GAUSSIAN_STDDEV03', p)
            self.assertNotIn('CURVATURE_PRIOR_GAUSSIAN_STDDEV03', p)
            self.assertNotIn('INIT_INDEX_MAP03_02', p)
            prefix = 'initial/' if settings['local_shapes'] else 'truth/'
            self.assertTrue(p['INIT_INDEX_MAP03_01'].startswith(prefix))

    def test_instrument_and_units_stay_fixed_with_cmb(self):
        p = generator.build_parameters(physical_manifest(), 128, 'forecast_sample', '/inputs', '/chains/test')
        self.assertFalse(p['POLARIZATION'])
        self.assertFalse(p['SAMPLE_INSIDE_MASK'])
        self.assertFalse(p['OUTPUT_MAPS'])  # Standard harmonic diagnostic would use uK_CMB.
        self.assertTrue(p['OUTPUT_FOREGROUNDS'])  # CMB reference output retains MJy/sr.
        for i in range(1, 132):
            suffix = generator.band_suffix(i)
            self.assertEqual(p['FREQ_UNIT' + suffix], 'MJy/sr')
            self.assertEqual(p['BANDPASS_TYPE' + suffix], 'delta')
            self.assertTrue(p['FIX_MONOPOLE' + suffix])
            self.assertTrue(p['FIX_DIPOLE' + suffix])
            for key in ('GAIN_PRIOR_RMS', 'BP_RMS', 'NOISEAMP_PRIOR_RMS', 'A2T', 'F2T'):
                self.assertEqual(p[key + suffix], 0)

    def test_pdf_only_baseline_does_not_enable_optin_cmb_flags(self):
        baseline = fixture_manifest()
        for stage in generator.STAGES:
            p = generator.build_parameters(baseline, 64, stage, '/inputs', '/chains/test')
            self.assertEqual(p['NUM_FG_SIGNAL_COMPONENTS'], 5)
            self.assertEqual(p['COMP_TYPE03'], 'power_law')
            self.assertNotIn('FOSSIL_SI_CMB', p)
            self.assertNotIn('FOSSIL_CMB_MIN_L', p)
            self.assertTrue(p['ENFORCE_ZERO_CL'])
            self.assertEqual(p['IN_POWSPECFILE'], 'none')

    def test_fullsky_likelihood_processing_and_component_masks(self):
        manifest = physical_manifest()
        for nside in (64, 128):
            for count in (64, 131):
                for stage in generator.STAGES:
                    p = generator.build_parameters(manifest, nside, stage, '/inputs', '/chains/test', band_count=count)
                    fullsky = manifest['resolutions'][str(nside)]['fullsky_mask']
                    for key in ('MASKFILE', 'MASKFILE_CORR', 'PROCESSING_MASK', 'OVERVIEW_MASK_HIGHFG',
                                'OVERVIEW_MASK_LOWFG', 'GAIN_APODIZATION_MASK'):
                        self.assertEqual(p[key], fullsky)
                    for component in range(1, 7):
                        for key in ('COMPONENT_MASK', 'INDEX_SAMPLING_MASK', 'PRIOR_OPTIMIZATION_MASK'):
                            self.assertEqual(p[key + f'{component:02d}'], 'fullsky')
                    for band in range(1, count+1):
                        self.assertEqual(p['MASKFILE_CALIB' + generator.band_suffix(band)], fullsky)
                    self.assertIn('Full-sky likelihood mask', generator.render(p, stage, nside))

    def test_cmb_order_is_validated(self):
        manifest = physical_manifest()
        manifest['resolutions']['64']['components'][0], manifest['resolutions']['64']['components'][5] = (
            manifest['resolutions']['64']['components'][5], manifest['resolutions']['64']['components'][0])
        with self.assertRaises(ValueError):
            generator.build_parameters(manifest, 64, 'forecast_sample', '/inputs', '/chains/test')

    def test_relative_cwd_prior_is_not_prefixed_twice(self):
        with tempfile.TemporaryDirectory() as directory:
            previous = Path.cwd()
            try:
                os.chdir(directory)
                Path('inputs/cmb').mkdir(parents=True)
                Path('inputs/cmb/prior.txt').write_text('2 100 0\n')
                generator.preflight_paths({'INPUT_DIRECTORY': 'inputs', 'CHAIN_DIRECTORY': 'chains',
                                          'IN_POWSPECFILE': 'inputs/cmb/prior.txt'}, verify_inputs=True)
            finally:
                os.chdir(previous)

    def test_static_suite_count_paths_hashes_and_restart_dependencies(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / 'manifest.json'
            manifest.write_text(json.dumps(physical_manifest()))
            script = Path(__file__).resolve().parents[1] / 'examples/fossil/generate_cmb_draine_suite.py'
            subprocess.run([sys.executable, str(script), '--manifest', str(manifest),
                            '--output-root', str(root / 'suite'), '--cluster-root', '/remote/fossil'],
                           check=True, capture_output=True, text=True)
            selected = root / 'suite/cmb_draine64'
            records = json.loads((selected / 'configurations.json').read_text())
            self.assertEqual(len(records), 26)
            self.assertEqual(len(list(selected.glob('*.txt'))), 26)
            self.assertFalse(list(selected.glob('*.par')))
            self.assertEqual({r['seed'] for r in records if r['stage'] == 'forecast_sample'},
                             {20261102, 20261103, 20261104, 20261105})
            readme = (selected / 'README.md').read_text()
            for record in records:
                self.assertEqual(record['input_directory'], '/remote/fossil/fossil_inputs_cmb_draine')
                self.assertTrue(record['chain_directory'].startswith('/remote/fossil/tests_64channels_cmb_draine_fullsky/'))
                self.assertEqual(record['analysis_signature'], 'fixture_fullsky_analysis_signature')
                self.assertEqual(record['analysis_mask'], 'fullsky')
                self.assertEqual(record['postprocess_mask'], 'mask70.fits')
                self.assertIn('# Analysis signature: fixture_fullsky_analysis_signature',
                              (selected / record['file']).read_text())
                self.assertEqual(record['mpi_ranks'], 64)
                self.assertEqual(record['sha256'], hashlib.sha256((selected / record['file']).read_bytes()).hexdigest())
                self.assertIn('](' + record['file'] + ')', readme)
                if record['case'] == '12':
                    self.assertEqual(record['first_iteration'], 1001)
                    self.assertEqual(record['final_iteration'], 2000)
                    self.assertEqual(record['expected_iterations'], 1000)
                    self.assertEqual(len(record['depends_on']), 1)
                    self.assertTrue((selected / record['depends_on'][0]).is_file())
                    self.assertIn('globals_sample_split1000_64bands_seed20261102/global_restart_c0001.txt',
                                  record['global_restart_file'])
            full = root / 'suite/cmb_draine131'
            records = json.loads((full / 'configurations.json').read_text())
            self.assertEqual(len(records), 2)
            self.assertEqual({r['nside'] for r in records}, {64, 128})
            self.assertTrue(all(r['mpi_ranks'] == 131 for r in records))
            self.assertTrue(all('/tests_131channels_cmb_draine_fullsky/' in r['chain_directory'] for r in records))


if __name__ == '__main__':
    unittest.main()
