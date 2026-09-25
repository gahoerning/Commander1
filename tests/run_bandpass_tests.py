#!/usr/bin/env python3
"""Compile the full production bandpass module with lightweight dependency fixtures.

Uses real bandpass file I/O, initialization and production trapezoidal quadrature.
Parameter parsing, random numbers, sorting and MPI use test stand-ins. This does
not replace a complete cluster build or multi-rank Commander run.
"""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

from run_color_tests import procedure

ROOT = Path(__file__).resolve().parents[1]
SUPPORT = ROOT / 'tests/bandpass_support'


def run():
    compiler = shlex.split(os.environ.get('FC', 'gfortran'))
    flags = shlex.split(os.environ.get('FFLAGS', ''))
    math = (ROOT/'src/include/math_tools.f90').read_text()
    deps = (SUPPORT/'dependencies.f90.in').read_text().replace('@TSUM@', procedure(math, 'tsum'))
    with tempfile.TemporaryDirectory(prefix='commander-bandpass-tests-') as tmp:
        work = Path(tmp)
        (work/'dependencies.f90').write_text(deps)
        exe = work/'bandpass-tests'
        command = compiler + flags + ['-std=f2008', '-ffree-line-length-none', '-fcheck=all',
                                      '-O0', '-J', tmp, '-I', tmp]
        sources = [work/'dependencies.f90', ROOT/'src/commander/comm_bp_mod.f90', SUPPORT/'integration.f90']
        subprocess.run(command + list(map(str, sources)) + ['-o', str(exe)], cwd=tmp, check=True)
        base = dict(NUMBAND='1', T_CMB='2.7255', CHAIN_DIRECTORY=f"'{tmp}'",
                    GAIN_INIT="'unused'", BANDPASS_INIT="'unused'", MJYSR_CONVENTION="'IRAS'",
                    APPLY_BP_CORRECTIONS='F', APPLY_GAIN_CORRECTIONS='F',
                    APPLY_COLOR_CORRECTIONS='T', POLARIZATION='F', GAIN_INIT_RMS='0', BP_INIT_RMS='0',
                    BANDPASS_TYPE01="'CBASS'", FREQ_C01='4.76', FREQ_UNIT01="'uK_ant'",
                    FREQ_LABEL01="'C-BASS I'", GAIN_CALIB_COMPONENT01='0', GAIN_PRIOR_RMS01='0',
                    A2T01='-1', F2T01='-1', CO2T01='-1', A2SZ01='-1',
                    USE_COLOR_CORRECTION01='F', BANDPASS_CAL_INDEX01='-0.299',
                    BANDPASS01=f"'{ROOT/'examples/bandpasses/cbass_I_20130625.dat'}'")
        cases = [('valid', {}, None),
                 ('polarization', {'POLARIZATION': 'T'}, 'Only Stokes I'),
                 ('wrong_unit', {'FREQ_UNIT01': "'uK_cmb'"}, 'FREQ_UNIT must be uK_ant'),
                 ('polynomial', {'USE_COLOR_CORRECTION01': 'T'}, 'Set USE_COLOR_CORRECTION=F'),
                 ('missing_index', {'BANDPASS_CAL_INDEX01': None}, 'MISSING_PARAMETER: BANDPASS_CAL_INDEX01'),
                 ('nan_index', {'BANDPASS_CAL_INDEX01': 'NaN'}, 'BANDPASS_CAL_INDEX must be finite'),
                 ('bad_reference', {'FREQ_C01': '0'}, 'FREQ_C must be finite and positive'),
                 ('nan_reference', {'FREQ_C01': 'NaN'}, 'FREQ_C must be finite and positive'),
                 ('override', {'A2T01': '1'}, 'automatic bandpass conversions'),
                 ('nan_override', {'F2T01': 'NaN'}, 'automatic bandpass conversions'),
                 ('negative_gain', {}, 'Use nonnegative linear responses'),
                 ('zero_gain', {}, 'Use nonnegative linear responses'),
                 ('nan_gain', {}, 'responses must be finite'),
                 ('nan_frequency', {}, 'responses must be finite'),
                 ('descending_frequency', {}, 'strictly increasing'),
                 ('zero_frequency', {}, 'strictly increasing'),
                 ('one_sample', {}, 'at least two samples'),
                 ('overflow_calibration', {}, 'Non-finite calibration spectrum')]
        for mode, changes, error in cases:
            config = base | changes
            params = work/'params.txt'
            params.write_text('\n'.join(f'{k} = {v}' for k,v in config.items() if v is not None)+'\n')
            result = subprocess.run([str(exe), str(params), mode], cwd=tmp, capture_output=True, text=True)
            if error is None:
                if result.returncode:
                    raise RuntimeError(result.stdout + result.stderr)
                print(result.stdout.strip())
            else:
                if result.returncode == 0 or error not in result.stdout:
                    raise RuntimeError(f'{mode}: expected rejection containing {error!r}\n{result.stdout}\n{result.stderr}')
                print('PASS: rejected', mode)

        # SPASS_RJ uses the same production module but a distinct RJ normalization.
        sources[-1] = SUPPORT/'spass_integration.f90'
        subprocess.run(command + list(map(str, sources)) + ['-o', str(exe)], cwd=tmp, check=True)
        base.update(BANDPASS_TYPE01="'SPASS_RJ'", FREQ_C01='2.303', FREQ_LABEL01="'S-PASS I'",
                    BANDPASS01=f"'{ROOT/'examples/bandpasses/spass_dr1_two_windows.dat'}'")
        del base['BANDPASS_CAL_INDEX01']  # This convention has no calibrator spectral index.
        for mode, changes, error in cases:
            if mode in {'missing_index', 'nan_index', 'overflow_calibration'}:
                continue
            config = base | changes
            params.write_text('\n'.join(f'{k} = {v}' for k,v in config.items() if v is not None)+'\n')
            result = subprocess.run([str(exe), str(params), mode], cwd=tmp, capture_output=True, text=True)
            if error is None:
                if result.returncode:
                    raise RuntimeError(result.stdout + result.stderr)
                print(result.stdout.strip())
            else:
                if result.returncode == 0 or error not in result.stdout:
                    raise RuntimeError(f'SPASS_RJ {mode}: expected {error!r}\n{result.stdout}\n{result.stderr}')
                print('PASS: SPASS_RJ rejected', mode)


if __name__ == '__main__':
    run()
