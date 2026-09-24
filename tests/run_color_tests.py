#!/usr/bin/env python3
"""Compile colour-correction tests without the cluster's HEALPix/MPI stack.

Integration fixtures compile the actual production procedures, extracted verbatim.
Only unrelated SEDs, splines, parameter I/O and single-rank MPI are stubbed.
They do not replace a full Commander build or a multi-rank cluster test.
"""
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'src' / 'commander'
SUPPORT = ROOT / 'tests' / 'color_support'


def procedure(text, name):
    pattern = rf'^  *(?:(?:logical|pure) +)?(function|subroutine) +{name}\b.*?^  *end +\1 +{name}\b[^\n]*\n'
    match = re.search(pattern, text, re.MULTILINE | re.DOTALL | re.IGNORECASE)
    if match is None:
        raise RuntimeError(f'Cannot extract production procedure {name}')
    return match.group(0)


def run():
    fg = (SRC / 'comm_fg_component_mod.f90').read_text()
    bp = (SRC / 'comm_bp_mod.f90').read_text()
    sampler = (SRC / 'comm_fg_mod.f90').read_text()
    types = fg[fg.index('  type fg_par_comp'):fg.index('  integer(i4b)                     :: num_fg_comp')]
    bandtype = bp[bp.index('  type bandinfo'):bp.index('  real(dp)                                      :: ind_iras')]
    funcs = ['color_corrected_component', 'monochromatic_fg_spectrum', 'reconstruct_color_sky',
             'foreground_color_factor', 'update_sky_color_corrections', 'get_effective_fg_spectrum',
             'get_effective_deriv_fg_spectrum', 'reorder_fg_params', 'deallocate_fg_params',
             'compute_power_law_spectrum', 'compute_faraday_rotation', 'update_fg_pix_response_map']
    fixtures = (SUPPORT / 'context.f90.in').read_text().replace('@TYPES@', types + bandtype)
    fixtures = fixtures.replace('@PROCEDURES@', '\n'.join(procedure(fg, f) for f in funcs)
                                + procedure(bp, 'read_color_correction_parameters')
                                + procedure(sampler, 'sample_QU_block_region'))
    compiler = shlex.split(os.environ.get('FC', 'gfortran'))
    flags = shlex.split(os.environ.get('FFLAGS', ''))
    with tempfile.TemporaryDirectory(prefix='commander-color-tests-') as tmp:
        work = Path(tmp)
        (work / 'context.f90').write_text(fixtures)
        commander = (SRC / 'commander.f90').read_text()
        start = commander.index('       if (sample_fg_pix .and. any(bp%use_color_corr) .and. iter > first_iteration) then')
        end = commander.index('       !call set_sample_temp_coeffs', start)
        schedule = commander[start:end]
        (work / 'cadence.f90').write_text((SUPPORT/'cadence.f90.in').read_text().replace('@SCHEDULE@', schedule))
        common = compiler + flags + ['-std=f2008', '-ffree-line-length-none', '-fcheck=all', '-O0',
                                     '-J', tmp, '-I', tmp]
        for name, sources in [
            ('numeric', [SRC / 'comm_color_correction_mod.f90', ROOT / 'tests/test_color_corrections.f90']),
            ('cadence', [work / 'cadence.f90']),
            ('integration', [SRC / 'comm_color_correction_mod.f90', work / 'context.f90',
                             SUPPORT / 'integration.f90']),
        ]:
            exe = work / name
            subprocess.run(common + list(map(str, sources)) + ['-o', str(exe)], cwd=tmp, check=True)
            subprocess.run([str(exe)], cwd=tmp, check=True)
        # The same executable verifies invalid configurations fail, not silently disable CC.
        for mode in ['missing_limit', 'delta_without_limits', 'invalid_limits', 'nonfinite_coeff', 'moving_bandpass']:
            result = subprocess.run([str(work/'integration'), mode], cwd=tmp, capture_output=True, text=True)
            if result.returncode == 0:
                raise RuntimeError(f'Invalid configuration accepted: {mode}')
            if 'EXPECTED_CONFIG_ABORT' not in result.stdout:
                raise RuntimeError(f'Unexpected failure for {mode}: {result.stdout}\n{result.stderr}')
            print('PASS: rejected', mode)


if __name__ == '__main__':
    run()
