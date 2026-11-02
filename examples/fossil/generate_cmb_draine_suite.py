#!/usr/bin/env python3
"""Generate the versioned Draine/free-free plus primary-CMB static suites.

Creates 26 selected-channel configurations and two full-band forecasts in new
directories, preserving the original PDF-only parameter suite. This script does
not launch Commander, write input maps, or transfer files to the cluster.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil

import generate_parameters as generator

DEFAULT_CLUSTER_ROOT = '/mn/stornext/d23/cmbco/AST9240/2024/gabrieah/fossil_component_separation'
SEEDS = (20261102, 20261103, 20261104, 20261105)
CASES = (
    ('01', 'globals_optimize'), ('02', 'globals_sample'),
    ('03', 'amplitudes_optimize'), ('04', 'global_shapes_sample'),
    ('05', 'amplitudes_sample'), ('06', 'forecast_optimize'),
    ('07', 'forecast_sample'), ('08', 'no_globals_regression'),
)
DESCRIPTIONS = {
    '01': ('Recover the six physical monopole amplitudes with a completely known local sky',
           'Noiseless; 12 optimize', 'Global amplitudes; Galactic maps/shapes and primary CMB fixed at truth'),
    '02': ('Validate the exact Gaussian conditional for physical monopole amplitudes',
           'Noisy; 2000 sample', 'Global amplitudes; all local quantities including primary CMB fixed at truth'),
    '03': ('Exercise joint local-amplitude/CMB CG coupled to the global block',
           'Noiseless; 12 optimize', 'Galactic amplitudes, harmonic CMB realization, global amplitudes; all spectral shapes fixed'),
    '04': ('Explore global spectral shapes with a completely known local sky',
           'Noisy; 5000 sample', 'Global amplitudes and tSZ/CIB shapes; Galactic quantities and primary CMB fixed at truth'),
    '05': ('Study amplitude degeneracies with a free primary-CMB realization',
           'Noisy; 5000 sample', 'Galactic amplitudes, harmonic CMB, global amplitudes; all spectral shapes fixed'),
    '06': ('Optimize the complete matched component-separation model',
           'Noiseless; 50 optimize', 'Galactic amplitudes/shapes, harmonic CMB realization, global amplitudes/shapes'),
    '07': ('Run the first full posterior forecast chain',
           'Noisy; 10000 sample', 'Galactic amplitudes/shapes, harmonic CMB realization, global amplitudes/shapes'),
    '08': ('Optional disabled-global software-path control',
           'Noiseless; 2 optimize', 'No globals fitted; all local quantities including primary CMB fixed'),
    '09': ('Run an additional independent full forecast chain',
           'Noisy; 10000 sample', 'Same complete model as case 07; independent seed'),
    '11': ('Run the first half of the global-state/RNG restart check',
           'Noisy; 1000 sample', 'Same fixed-local conditional as case 02'),
    '12': ('Resume the global-state/RNG check to iteration 2000',
           'Noisy; iterations 1001–2000', 'Restore the global state and RNG from case 11; fixed local sky'),
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_inputs(parameters, local_export, remote_input):
    """Check the concrete package while parameters retain their remote roots."""
    checked = set()
    root = Path(remote_input)
    for key, value in parameters.items():
        if not isinstance(value, str) or value == 'none' or key == 'GLOBAL_RESTART_FILE':
            continue
        if not value.lower().endswith(('.fits', '.txt', '.dat')):
            continue
        path = Path(value)
        if path.is_absolute():
            path = path.relative_to(root)
        candidate = local_export / path
        if not candidate.is_file():
            raise FileNotFoundError(f'{key}: missing local counterpart {candidate}')
        checked.add(str(path))
    return checked


def generate_suite(manifest_path, output_root, cluster_root=DEFAULT_CLUSTER_ROOT, check_inputs=False):
    manifest_path, output_root = Path(manifest_path), Path(output_root)
    manifest = json.loads(manifest_path.read_text())
    if not manifest.get('primary_cmb_anisotropy') or not manifest.get('freefree_model', '').startswith('draine'):
        raise ValueError('Require the explicit primary-CMB/Draine matched export')
    profile = Path(generator.__file__).with_name('channel_selection_cmb_draine_64.json')
    profile_sha = sha(profile)
    selected = generator.select_channels(manifest, 64)
    inputs = str(Path(cluster_root) / 'fossil_inputs_cmb_draine')
    fullsky = manifest.get('analysis_mask') == 'fullsky'
    analysis_signature = manifest.get('analysis_signature')
    if fullsky and not analysis_signature:
        raise ValueError('A full-sky analysis requires explicit analysis_signature provenance')
    analysis_tag = '_fullsky' if fullsky else ''
    folder64, folder131 = output_root / 'cmb_draine64', output_root / 'cmb_draine131'
    folder64.mkdir(parents=True, exist_ok=True)
    folder131.mkdir(parents=True, exist_ok=True)
    if list(folder64.glob('*.par')) or list(folder131.glob('*.par')):
        raise ValueError('Versioned suites use .txt files; remove obsolete .par copies explicitly')
    records64, records131, all_inputs = [], [], set()

    def write(case, stage, nside, seed, count=64, label=None, iterations=None, first=1, restart='none', dependency=None):
        label = label or stage
        chain = str(Path(cluster_root) / f'tests_{count}channels_cmb_draine{analysis_tag}' / f'nside{nside}' /
                    f'{label}_{count}bands_seed{seed}')
        p = generator.build_parameters(manifest, nside, stage, inputs, chain, seed,
                                       iterations=iterations, band_count=count)
        p['GLOBAL_RESTART_FILE'] = restart
        generator.preflight_paths(p)
        if check_inputs:
            all_inputs.update(verify_inputs(p, manifest_path.parent, inputs))
        name = f'{case}_{label}_nside{nside:04d}_{count}bands_seed{seed}.txt'
        folder = folder64 if count == 64 else folder131
        header = (f"# Active model: {manifest['model']}\n"
                  f"# Export manifest signature: {manifest['signature']}\n"
                  + (f'# Analysis signature: {analysis_signature}\n' if analysis_signature else '')
                  + (f'# Fixed 64-channel profile SHA256: {profile_sha}\n' if count == 64 else '')
                  + f'# Suite case: {case}; seed: {seed}; first iteration: {first}\n')
        path = folder / name
        path.write_text(header + generator.render(p, label, nside))
        record = dict(file=name, case=case, stage=stage, label=label, nside=nside, seed=seed,
            operation=p['OPERATION'], data='noisy' if generator.STAGES[stage]['observed'] else 'noiseless',
            band_count=count, mpi_ranks=p['NUM_GROUPS'] * count * p['NUM_PROC_PER_BAND'],
            first_iteration=first, final_iteration=p['NUM_GIBBS_ITER'],
            expected_iterations=p['NUM_GIBBS_ITER']-first+1, input_directory=inputs,
            chain_directory=chain, global_restart_file=restart, depends_on=[dependency] if dependency else [],
            manifest_signature=manifest['signature'], sha256=sha(path),
            remote_file=str(Path(inputs) / name), verbosity=p['VERBOSITY'],
            file_format='sectioned_aligned_English_comments', pixel_window=p['PIXEL_WINDOW'],
            cmb_mode='known truth' if p['ENFORCE_ZERO_CL'] else 'harmonic realization with fixed D_l prior',
            cmb_min_l=p['FOSSIL_CMB_MIN_L'], cmb_si_units=p['FOSSIL_SI_CMB'], cmb_prior=p['IN_POWSPECFILE'],
            freefree_model=manifest['freefree_model'])
        if analysis_signature:
            record.update(analysis_signature=analysis_signature, analysis_mask=manifest.get('analysis_mask'),
                          postprocess_mask=manifest['resolutions'][str(nside)].get('postprocess_mask'))
        if count == 64:
            record['selection_profile_sha256'] = profile_sha
            records64.append(record)
        else:
            records131.append(record)
        return record

    for nside in (64, 128):
        for case, stage in CASES:
            write(case, stage, nside, SEEDS[0])
        for seed in SEEDS[1:]:
            write('09', 'forecast_sample', nside, seed)
        split = write('11', 'globals_sample', nside, SEEDS[0], label='globals_sample_split1000', iterations=1000)
        write('12', 'globals_sample', nside, SEEDS[0], label='globals_sample_resume2000', iterations=2000,
              first=1001, restart=str(Path(split['chain_directory']) / 'global_restart_c0001.txt'), dependency=split['file'])
        write('07', 'forecast_sample', nside, SEEDS[0], count=131)
    assert len(records64) == 26 and len(records131) == 2
    for folder, records in ((folder64, records64), (folder131, records131)):
        (folder / 'configurations.json').write_text(json.dumps(records, indent=2) + '\n')
    shutil.copyfile(profile, folder64 / 'channel_selection_64.json')
    selected_ids = {int(c['commander_band']) for c in selected}
    for filename, channels in (('selected_channels.csv', selected), ('removed_channels.csv',
            [c for c in manifest['channels'] if int(c['commander_band']) not in selected_ids])):
        with (folder64 / filename).open('w', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['active_commander_band', 'original_commander_band', 'frequency_GHz', 'sigma_MJy_sr'])
            for active, channel in enumerate(channels, 1):
                writer.writerow([active if filename.startswith('selected') else '', channel['commander_band'],
                                 channel['frequency_GHz'], channel['sigma_MJy_sr']])
    (folder64 / 'README.md').write_text(readme64(records64, cluster_root, manifest, profile_sha))
    (folder131 / 'README.md').write_text(readme131(records131, cluster_root, manifest))
    report = dict(status='passed', manifest_signature=manifest['signature'], manifest_sha256=sha(manifest_path),
        selection_profile_sha256=profile_sha, selected_static_files=26, full_static_files=2,
        checked_local_input_files=len(all_inputs), local_input_files_verified=check_inputs,
        primary_cmb_anisotropy=True, freefree_model=manifest['freefree_model'],
        selected_mpi_ranks=64, full_mpi_ranks=131, required_recompile=True,
        commander_or_mpi_executed=False, cluster_acceptance='pending',
        original_pdf_only_static_suite_modified=False)
    if analysis_signature:
        report.update(analysis_signature=analysis_signature, analysis_mask=manifest['analysis_mask'],
                      postprocess_masks={n: r.get('postprocess_mask') for n, r in manifest['resolutions'].items()})
    (folder64 / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def readme64(records, root, manifest, profile_sha):
    analysis_tag = '_fullsky' if manifest.get('analysis_mask') == 'fullsky' else ''
    text = f'''# Primary CMB and Draine free-free: 64-channel component separation

This versioned suite contains 26 complete `.txt` parameter files for the new
`{manifest['model']}` injection. The original PDF-only maps and parameter files
remain available separately. These files require the updated Commander source
to be rebuilt; full HEALPix/MPI acceptance of this extension is pending.

## Contents

- [Model and fixed quantities](#model-and-fixed-quantities)
- [Paths and manual execution](#paths-and-manual-execution)
- [Every configuration](#every-configuration)
- [Expected checks and limitations](#expected-checks-and-limitations)
- [Restart and posterior convergence](#restart-and-posterior-convergence)

## Model and fixed quantities

The six local entries are dust, synchrotron, free-free, AME1, AME2, and primary
CMB. Every intensity and RMS map is MJy/sr in Galactic HEALPix RING coordinates.
The Gaussian FWHM is 96 arcmin. Both NSIDE values retain the original per-pixel
TXT sensitivity and noise realizations. Electron temperature maps use K.

The likelihood, processing and overview masks now use the **full sky**. GAL070
is retained only for optional post-processing of recovered maps, spectra and
summary statistics; it does not remove pixels from the fit. The input
normalizations originally defined over GAL070 are unchanged, as are every
injected intensity map and per-pixel sensitivity. The scientific injection
signature is preserved; a separate analysis signature records the full-sky
mask and recomputed sufficient statistics. The active analysis entry point is
`fossil_fullsky.export_fullsky_commander`; it wraps the preserved physical
injection and configures this full-sky analysis.

Free-free uses the normalized Draine Gaunt-factor spectrum with one reference
intensity amplitude at 50 GHz and one electron temperature Te per pixel. Its
reference morphology comes from Planck's EM/Te maps, with one common amplitude
normalization at 30 GHz. EM is archived as input provenance; it is not an
additional fitted amplitude. Dust beta/temperature, synchrotron beta, Te, and
the two AME peak frequencies are the six free local spectral parameters in full
forecasts. Synchrotron curvature remains fixed at zero.

Primary CMB uses a shared Gaussian realization from the Planck TT spectrum.
Its harmonic monopole and dipole are exactly zero. In amplitude/full-forecast
cases Commander samples or optimizes this harmonic realization jointly with
the Galactic amplitudes, using a fixed, unwindowed TT covariance prior.
`CONDITION_ON_CL=true` fixes the spectrum; it does not fix the CMB realization.
`IN_POWSPECFILE` contains ell, D_ell_TT, and zero initialization uncertainty in
uK_CMB squared. The operator applies the Gaussian and native pixel window once.
`FOSSIL_CMB_MIN_L=2` excludes ell0/1; `FOSSIL_SI_CMB=true` matches the exact SI
derivative used by the injection and physical DeltaT monopole. The CMB pixel
entry is a synchronized reference-intensity representation, not a second free
pixel CMB field. Standard uK_CMB map diagnostics are disabled; its saved
component amplitude is MJy/sr at 50 GHz.

Cases 01, 02, 04, 08, 11, and 12 hold the exact injected CMB reference map fixed
and disable the harmonic solve. Cases 03, 05, 06, 07, and 09 fit the CMB
realization. The physical DeltaT/T0, mu, tSZ, CIB, CRR, and extragalactic CO
monopoles remain in the input. Instrumental offsets/dipoles stay zero, gains
one, bandpass shifts zero, and noise scaling fixed. The reference 2.725-K
blackbody is subtracted. The global DeltaT/T0 monopole is fitted alongside the
anisotropies: the anisotropy field has no harmonic monopole/dipole, so it cannot
absorb that independent constant sky amplitude. The CMB spectrum is fixed;
these files do not sample C_ell. Apply any GAL070 cut after fitting, retaining
the corresponding masked CMB mean rather than forcing it to zero.

The fixed 64-channel profile is in [channel_selection_64.json](channel_selection_64.json),
with active/original IDs and unchanged sensitivities in [selected_channels.csv](selected_channels.csv)
and the omitted bands in [removed_channels.csv](removed_channels.csv).
Its SHA256 is `{profile_sha}`. This fiducial spectral-design proxy does not
replace a spatial Fisher calculation or posterior convergence checks.

## Paths and manual execution

All files contain expanded absolute cluster paths:

```sh
FOSSIL_ROOT={root}
INPUT_DIRECTORY="$FOSSIL_ROOT/fossil_inputs_cmb_draine"
cd "$FOSSIL_ROOT"
PARFILE="$INPUT_DIRECTORY/01_globals_optimize_nside0064_64bands_seed20261102.txt"
./comm_run.sh --check "$PARFILE"
./comm_run.sh "$PARFILE"
```

Copy these parameter files into the new input directory beside the new matched
map/template package. Chains go into `tests_64channels_cmb_draine{analysis_tag}/nside64/` or
`nside128/`, each case and seed having its own directory. Use fresh directories;
scalar traces append. The updated launcher creates the chain parents and the
required native HEALPix pixel-window link before MPI. Its `--check` is read-only
and reports prospective links. Keep the original full-length
`pixel_window_n0064.fits` and `pixel_window_n0128.fits` in INPUT_DIRECTORY; they
contain 257/513 rows and both temperature/polarization columns. A custom launcher
must provide the standard filename in Commander's working directory too.

Each file requires 64 MPI ranks with one rank per selected frequency. On one
64-core owl node use one OpenMP/MKL/OpenBLAS thread per rank. Start with case 01,
then 03, then 06, then the full case 07; insert case 02 to test the analytic
global conditional before expensive sampling. Select cases individually rather
than automatically running every file. Four full seeds are supplied for each
resolution; iteration counts do not guarantee convergence.

## Every configuration

“Global shapes” means tSZ kTe and CIB beta/temperature. Their amplitudes are
spatially constant. “Harmonic CMB” always uses the fixed TT spectrum and ell>=2.
Optimize mode returns prior-dependent point estimates; sample mode draws from
the configured posterior.

| Case and purpose | Data; iterations | Free quantities | NSIDE 64 | NSIDE 128 |
| --- | --- | --- | --- | --- |
'''
    groupings = [(case, stage, SEEDS[0]) for case, stage in CASES]
    groupings += [('09', 'forecast_sample', seed) for seed in SEEDS[1:]]
    groupings += [('11', 'globals_sample_split1000', SEEDS[0]), ('12', 'globals_sample_resume2000', SEEDS[0])]
    for case, label, seed in groupings:
        purpose, data, free = DESCRIPTIONS[case]
        if case == '09':
            purpose += f' (seed {seed})'
        links = []
        for nside in (64, 128):
            record = next(r for r in records if r['case'] == case and r['label'] == label
                          and r['seed'] == seed and r['nside'] == nside)
            links.append(f"[{case} / NSIDE {nside}]({record['file']})")
        text += f'| {case}. {purpose} | {data} | {free} | {links[0]} | {links[1]} |\n'
    text += '''
## Expected checks and limitations

Case 01 should recover injected global amplitudes after the exact known
Galactic/CMB sky is subtracted. Case 02 is compared with the exact selected-band
Gaussian conditional. The exported global-only sufficient statistics remain
valid because that reference excludes the known local/CMB sky. Do not treat
this fixed-sky test as the uncertainty of full component separation.
For this analysis the sufficient statistics use every full-sky pixel; the
previous GAL070 Gaussian reference is not the correct acceptance target.

Case 03 must start and finish its CG solve at each accepted iteration. Inspect
CG convergence, finite outputs, residuals and physical harmonic ell0/1 zero.
It uses a CMB covariance prior even though the data are noiseless. Its optimum
is a Wiener/prior MAP solution, so exact zero residuals or exact recovery of
every CMB/Galactic/global amplitude are not valid acceptance requirements.
Case 06 also has nonlinear spectral priors and degeneracies. These cases
exercise the complete coupled model; investigate numerical failures and
systematic residual patterns rather than assuming noiseless means unregularized.

Cases 04 and 05 release nonlinear global shapes and uncertain local/CMB
amplitudes in stages. Inspect their correlations and mixing. Case 08 omits
physical global signals still present in the data; a substantial residual is
expected. It is an optional software-path control, not a foreground-only
forecast. No cluster run is certified by generation or local compiled fixtures.

## Restart and posterior convergence

Run case 11 before 12 at the same resolution. The resume file reads
`global_restart_c0001.txt` from the split1000 chain and writes elsewhere.
`NUM_GIBBS_ITER=2000` is the final iteration, generating rows 1001–2000.
Concatenate cases 11/12 and compare exactly with uninterrupted case 02 using
the same seed. This is a fixed-local global-state/RNG restart test. A full
forecast restart also needs matching Galactic maps, CMB alms, spectral
parameters and legacy RNG state; the global checkpoint alone is insufficient.

For final inference use case 07 and all three case 09 seeds for each resolution.
Check rank-normalized split R-hat <= 1.01, bulk ESS >= 1000 and tail ESS >= 400 for all
reported global amplitudes/shapes, plus representative Galactic parameters and
CMB modes. Inspect joint degeneracies, priors, boundaries and residuals. Extend
or improve blocking when these checks fail. Bias/coverage claims additionally
need independent sky/noise realizations. The current matched output-pixel
Te model does not assert that smoothing a nonlinear native SED commutes with
averaging its parameters.
'''
    return text


def readme131(records, root, manifest):
    text = '''# Full 131-channel primary-CMB/Draine forecasts

These two static files use all 131 original FOSSIL bands, the same versioned
input package as the 64-channel suite, and separate output chains. They require
131 MPI ranks with NUM_PROC_PER_BAND=1. Rebuild Commander with the opt-in CMB
fixes before running. Full HEALPix/MPI acceptance is pending.

See [the complete model, stage guide and diagnostics](../cmb_draine64/README.md).
The harmonic CMB realization is free with a fixed TT prior and ell>=2;
Galactic/global amplitudes and their allowed shapes are also free.
The likelihood uses the full sky; GAL070 is for optional post-processing only.
The CMB spectrum remains fixed and the global DeltaT monopole is fitted
separately from the zero-monopole/dipole anisotropy field. New `_fullsky` chain
roots preserve outputs from the earlier GAL070 analysis.

| Resolution | Static parameter file | Chain directory |
| --- | --- | --- |
'''
    for record in records:
        text += f"| NSIDE {record['nside']} | [{record['file']}]({record['file']}) | `{record['chain_directory']}` |\n"
    text += f'''
Copy the desired file into `{root}/fossil_inputs_cmb_draine/` and pass its full
path to your launcher. The launcher must allocate 131 ranks. Each file provides
10000 iterations with seed 20261102, not a convergence guarantee; create
separate files/seeds/chains for independent full-band inference.
'''
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, default=Path(__file__).parent)
    parser.add_argument('--cluster-root', default=DEFAULT_CLUSTER_ROOT)
    parser.add_argument('--verify-inputs', action='store_true')
    args = parser.parse_args()
    report = generate_suite(args.manifest, args.output_root, args.cluster_root, args.verify_inputs)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
