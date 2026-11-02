# FOSSIL component separation with 64 selected channels

This directory contains 26 complete, static Commander parameter files (`.txt`) for NSIDE 64 and
128. The full forecast cases fit Galactic pixel amplitudes and spectral
parameters together with spatially constant physical components. Each job uses
64 channels and 64 MPI ranks when `NUM_PROC_PER_BAND=1`.

## Contents

- [Reading the parameter files](#reading-the-parameter-files)
- [Channel selection and sky conventions](#channel-selection-and-sky-conventions)
- [Paths and manual runs](#paths-and-manual-runs)
- [Native HEALPix pixel windows](#native-healpix-pixel-windows)
- [Purpose of every configuration](#purpose-of-every-configuration)
- [What to inspect after each run](#what-to-inspect-after-each-run)
- [Restart and forecast interpretation](#restart-and-forecast-interpretation)

## Reading the parameter files

The `.txt` files use named sections, aligned parameter/value columns, and English
comments explaining the data paths, units, fixed quantities, and sampled
parameters. Channel inputs and individual sky components have their own blocks,
so their associated settings can be read together. They use the same Commander
parameter syntax as the earlier files.

| File section | Settings grouped there |
| --- | --- |
| 01. Run mode and logging | Optimize/sample operation and diagnostic verbosity |
| 02. MPI layout | Groups, active bands, pixel partitions, and required rank count |
| 03. Map resolution and spectral-parameter resolution | NSIDE, multipole limits, and 96-arcmin beam |
| 04. Chain length, random seed, and spectral sampling | Final iteration, seed, thinning, and spectral controls |
| 05. Amplitude solver and CG convergence | Linear solver, preconditioner, and convergence tolerance |
| 06. Polarization and CMB harmonic controls | Intensity-only setup and inactive primary-CMB harmonic model |
| 07. Output products | Maps, component outputs, and chi-square diagnostics |
| 08. Input/output paths and analysis masks | Absolute input/chain roots and relative mask filenames |
| 09. Noise, beam, units, and fixed instrument | Known RMS convention, beam, units, and fixed calibration |
| 10. Frequency map sets | One block per selected frequency, with its data/RMS files and fixed offsets |
| 11. Galactic pixel components | Local amplitude maps, spectral parameters, priors, and freeze settings |
| 12. Spatially constant physical sky components | Global amplitudes, spectral shapes, templates, and restart settings |

`VERBOSITY = 3` enables detailed diagnostic output, including the CG information
needed to inspect case 03. The readable suite preserves the scientific model,
data/noise maps, priors, seeds, free/fixed parameters, iteration counts, and
restart dependencies. The pixel-window correction changes only `PIXEL_WINDOW`
to the original full-length HEALPix filename and refreshes provenance comments.
A `.txt` filename can be passed directly to the existing compiled Commander
executable.

The source copies are kept in `examples/fossil/owl64/` in the Commander
repository. The cluster copies are stored directly in `fossil_inputs/`, beside
the existing input package; the full paths are shown below.

## Channel selection and sky conventions

The channel profile is fixed in the local `channel_selection_64.json`, with the
active channel mapping and unchanged per-pixel RMS values in
`selected_channels.csv`. Selection uses noise-weighted spectral information for
20 currently free amplitude/shape directions, retaining the 50, 545, and
2000-GHz anchors and coverage across the frequency range. It is a fiducial,
physics-informed channel choice; the 64-channel posterior has its own information
content and must be diagnosed independently. No new maps or recompilation are
needed for these parameter files.

| Channel-selection record | Contents |
| --- | --- |
| [Selected channels](selected_channels.csv) | All 64 active indices, original channel IDs, frequencies, and per-pixel RMS values |
| [Removed channels](removed_channels.csv) | The 67 original channels omitted from these configurations |
| [Fixed selection profile](channel_selection_64.json) | Exact selection, source hashes, method, and spectral-information diagnostics |

The four equal frequency intervals retain 28, 12, 12, and 12 channels; the largest
gap is 165 GHz. The selector uses a backward D-optimal heuristic for fiducial
spectral derivatives. Its information metrics are spectral-design proxies, not
the marginalized uncertainties of the spatial component-separation posterior.
The local profile is an exact copy of the version consumed by the repository's
parameter generator. This folder can be copied on its own for manual runs.

Sky-intensity and RMS maps use MJy/sr, Galactic coordinates, and HEALPix RING
ordering, with a 1.6-degree Gaussian FWHM. NSIDE 64 and 128 use the same TXT noise
standard deviation per pixel. The reference
2.725-K blackbody has been subtracted, while the physical temperature-offset
monopole, mu, tSZ, CIB, CRR, and aggregate extragalactic CO signals remain.
The current input sky has no primary CMB anisotropy realization. Instrumental
offsets/dipoles are fixed at zero, gains at one, and bandpass shifts at zero.

## Paths and manual runs

The supplied files use these exact cluster locations:

```sh
FOSSIL_ROOT=/mn/stornext/d23/cmbco/AST9240/2024/gabrieah/fossil_component_separation
INPUT_DIRECTORY=$FOSSIL_ROOT/fossil_inputs
```

The filtered parameter deployment can replace the 26 old parameter filenames
with the 26 `.txt` files and copy the suite guide, index, selection records,
and the exact native pixel-window/provenance files.
[`rsync_parameters.rules`](rsync_parameters.rules) lists the exact suite
filenames; its final exclusion protects the other contents of `fossil_inputs`.
The dedicated pixel-window fix bundle copies its contents into FOSSIL_ROOT
without deleting other files. It includes the updated launcher at the root and
the corrected inputs under `fossil_inputs/`.

Each parameter file contains the expanded absolute paths. Its chain directory is recorded
in `configurations.json`, below
`$FOSSIL_ROOT/tests_64channels/nside64/` or `nside128/`. Prepare the parents once:

```sh
mkdir -p "$FOSSIL_ROOT/tests_64channels/nside64" "$FOSSIL_ROOT/tests_64channels/nside128"
```

Choose a case from the table and pass that `.txt` file to your 64-rank launcher.
For example, use the supplied launcher and parameter copy on the cluster:

```sh
PARFILE="$INPUT_DIRECTORY/01_globals_optimize_nside0064_64bands_seed20261102.txt"
cd "$FOSSIL_ROOT"
./comm_run.sh --check "$PARFILE"
./comm_run.sh "$PARFILE"
```

The supplied launcher creates the selected file's full `CHAIN_DIRECTORY` and
its required subdirectories before a run. With a different launcher, create
them before invoking MPI. Start fresh cases in empty chain directories: scalar
traces append to existing files. The files
below are independent experiments, except for the explicitly linked split/resume
pair. Selecting an optimize file does not initialize another forecast from that
result automatically. Choose the desired cases individually; running every parameter file
includes long forecast chains.

For one shared owl31–35 node, use 64 MPI ranks and set `OMP_NUM_THREADS=1`,
`MKL_NUM_THREADS=1`, and `OPENBLAS_NUM_THREADS=1` in the launcher environment.
Check node availability before starting. Begin with NSIDE 64 cases 01, 02, and
03; continue to the more expensive cases after inspecting their outputs.

## Native HEALPix pixel windows

The input root must contain the unchanged original files
`pixel_window_n0064.fits` and `pixel_window_n0128.fits`. Each native binary table
has the `TEMPERATURE` and `POLARIZATION` columns through ell = 4 * NSIDE: 257 rows
for NSIDE 64 and 513 for NSIDE 128. These windows are dimensionless multipole
transfer functions. File hashes and original source URLs are recorded in
`pixel_window_provenance.json`, and `manifest.json` identifies the table for
each resolution. The channel profile's updated source-manifest hash records this
metadata addition; the selected channels and scientific simulation signature
are unchanged.

`PIXEL_WINDOW` now names the appropriate native table relative to
INPUT_DIRECTORY. Commander's legacy `read_pixwin` also opens the standard
filename in its process working directory, independently of that parameter.
The supplied launcher validates the full native table, then creates the
working-directory link before starting MPI. An existing working-directory
table must match the input file byte for byte. `--check` validates the inputs
and reports the prospective link; it creates no directories or links and does
not run MPI. A custom launcher must provide this same standard file/link in
the directory where Commander starts.

Do not use the shorter `nside64/pixel_window.fits` or
`nside128/pixel_window.fits` files for this implicit lookup: the legacy reader
requests ell = 0 through 4 * NSIDE. `NSIDE_CORR=16` does not trigger a separate
pixel-window lookup in these runs. The correction requires no Fortran rebuild,
new sky/noise maps, or changes to the per-pixel sensitivity.

## Purpose of every configuration

“Local” means the five Galactic amplitude maps: dust, synchrotron, free-free, and
two AME populations. Free local shapes are dust beta/temperature, synchrotron
beta, free-free beta, and both AME peaks. Synchrotron/free-free curvature stays
fixed at zero. “Global shapes” means tSZ kTe and CIB beta/temperature; the other
global templates have fixed spectral shapes.

Optimize mode produces point estimates with the configured priors. Sample mode
produces posterior draws. Iteration counts are initial run lengths, not a
convergence guarantee.

| Case and purpose | Data; iterations | Free parameters | NSIDE 64 TXT | NSIDE 128 TXT |
| --- | --- | --- | --- | --- |
| 01. Recover physical monopole amplitudes with a known local sky | Noiseless; 12 optimize | Global amplitudes; all local quantities and global shapes fixed at truth | [01 / NSIDE 64](01_globals_optimize_nside0064_64bands_seed20261102.txt) | [01 / NSIDE 128](01_globals_optimize_nside0128_64bands_seed20261102.txt) |
| 02. Check the exact Gaussian global-amplitude posterior | Noisy; 2000 sample | Global amplitudes; local quantities and global shapes fixed at truth | [02 / NSIDE 64](02_globals_sample_nside0064_64bands_seed20261102.txt) | [02 / NSIDE 128](02_globals_sample_nside0128_64bands_seed20261102.txt) |
| 03. Exercise the local CG solve and its coupling to global amplitudes | Noiseless; 12 optimize | Local and global amplitudes; all shapes fixed at truth | [03 / NSIDE 64](03_amplitudes_optimize_nside0064_64bands_seed20261102.txt) | [03 / NSIDE 128](03_amplitudes_optimize_nside0128_64bands_seed20261102.txt) |
| 04. Explore global spectral parameters with a known local sky | Noisy; 5000 sample | Global amplitudes and shapes; local quantities fixed at truth | [04 / NSIDE 64](04_global_shapes_sample_nside0064_64bands_seed20261102.txt) | [04 / NSIDE 128](04_global_shapes_sample_nside0128_64bands_seed20261102.txt) |
| 05. Explore amplitude degeneracies before releasing local shapes | Noisy; 5000 sample | Local and global amplitudes; all shapes fixed at truth | [05 / NSIDE 64](05_amplitudes_sample_nside0064_64bands_seed20261102.txt) | [05 / NSIDE 128](05_amplitudes_sample_nside0128_64bands_seed20261102.txt) |
| 06. Fit the complete matched model as a point-estimate test | Noiseless; 50 optimize | Local amplitudes/shapes and global amplitudes/shapes | [06 / NSIDE 64](06_forecast_optimize_nside0064_64bands_seed20261102.txt) | [06 / NSIDE 128](06_forecast_optimize_nside0128_64bands_seed20261102.txt) |
| 07. Full component separation, forecast chain 1 | Noisy; 10000 sample | Local amplitudes/shapes and global amplitudes/shapes | [07 / NSIDE 64](07_forecast_sample_nside0064_64bands_seed20261102.txt) | [07 / NSIDE 128](07_forecast_sample_nside0128_64bands_seed20261102.txt) |
| 08. Optional disabled-global-components regression control | Noiseless; 2 optimize | No fitted global components; local quantities fixed | [08 / NSIDE 64](08_no_globals_regression_nside0064_64bands_seed20261102.txt) | [08 / NSIDE 128](08_no_globals_regression_nside0128_64bands_seed20261102.txt) |
| 09. Full forecast chain 2, independent seed 20261103 | Noisy; 10000 sample | Same free quantities as case 07 | [09 / NSIDE 64, seed 3](09_forecast_sample_nside0064_64bands_seed20261103.txt) | [09 / NSIDE 128, seed 3](09_forecast_sample_nside0128_64bands_seed20261103.txt) |
| 09. Full forecast chain 3, independent seed 20261104 | Noisy; 10000 sample | Same free quantities as case 07 | [09 / NSIDE 64, seed 4](09_forecast_sample_nside0064_64bands_seed20261104.txt) | [09 / NSIDE 128, seed 4](09_forecast_sample_nside0128_64bands_seed20261104.txt) |
| 09. Full forecast chain 4, independent seed 20261105 | Noisy; 10000 sample | Same free quantities as case 07 | [09 / NSIDE 64, seed 5](09_forecast_sample_nside0064_64bands_seed20261105.txt) | [09 / NSIDE 128, seed 5](09_forecast_sample_nside0128_64bands_seed20261105.txt) |
| 11. First half of the global RNG restart check | Noisy; 1000 sample | Same fixed-local conditional model as case 02 | [11 / NSIDE 64](11_globals_sample_split1000_nside0064_64bands_seed20261102.txt) | [11 / NSIDE 128](11_globals_sample_split1000_nside0128_64bands_seed20261102.txt) |
| 12. Resume case 11 to iteration 2000 | Noisy; iterations 1001–2000 | Same model; restore global parameters and global RNG | [12 / NSIDE 64](12_globals_sample_resume2000_nside0064_64bands_seed20261102.txt) | [12 / NSIDE 128](12_globals_sample_resume2000_nside0128_64bands_seed20261102.txt) |

Case 08 intentionally omits physical global signals present in the data. A large
residual is expected and is not a component-recovery failure. It tests the
disabled-global path with this 64-channel configuration; comparison to a separate
old-binary 99-band control is a different test.

## What to inspect after each run

| Case | Evidence and interpretation |
| --- | --- |
| 01 | The final six global amplitudes recover their injected fiducials with the local sky and all shapes fixed. The optional acceptance tool requires maximum amplitude error divided by its configured `GLOBAL_SCALE` <= 1e-5. This isolates template units, normalization, band mapping, and the global linear solve. |
| 02 | The sampled mean and covariance agree with the exact six-dimensional Gaussian conditional computed from the selected channels, mask, RMS, and priors. Individual draws need not equal the injection. The acceptance tool checks every whitened mean and covariance entry at six standard errors and rejects constant draws. |
| 03 | Local CG searches start and finish at each accepted iteration, all global trace entries remain finite, and residuals/CG convergence are inspected. This exposes integration errors when the pixel amplitudes and global block are both active. It does not require a unique amplitude decomposition in degenerate directions. |
| 04 | kTe and CIB beta/temperature actually vary, stay within their bounds, and have sensible joint distributions with the global amplitudes. This exercises nonlinear spectral sampling in isolation; one exploratory chain does not establish convergence. |
| 05 | Shapes remain fixed, both amplitude blocks run, and local/global correlations and residuals are inspected. Compare its broader conditional distributions with case 02 to see the effect of uncertain foreground amplitudes. |
| 06 | The full nonlinear fit completes without invalid parameters or numerical failures. Inspect the fit residuals, CG convergence, priors, and boundaries. Degeneracies and priors can move point estimates away from individual injected parameters. |
| 07 + 09 | Use all four seeds for posterior convergence and component-separation uncertainties. Inspect global and local parameters, joint degeneracies, residuals, R-hat, and effective sample sizes. The settings below define diagnostic targets; the initial chain length is not a guarantee. |
| 08 | The run completes with the global block disabled. Nonzero residuals are expected because the input still contains the omitted physical signals. This is an optional software-path control, not a valid foreground-only forecast. |
| 11 + 12 | Concatenating the split and resumed global traces reproduces case 02 exactly with the same seed and fixed local model. This tests restoration of the new global parameters and RNG state. |

Save the terminal output with each chain, especially the CG messages for case
03. The parameter files alone configure the experiment; inspecting actual
Commander outputs is what establishes a pass. No cluster run has been certified
by the local parameter-file checks.

## Restart and forecast interpretation

For each resolution, run case 11 before case 12. The resume file reads
`global_restart_c0001.txt` from that resolution's
`globals_sample_split1000_64bands_seed20261102` chain directory and writes into a
separate `globals_sample_resume2000_64bands_seed20261102` directory. Compare the
concatenated scalar rows from cases 11 and 12 to the uninterrupted case 02.
They should match exactly for this fixed-local model with the same seed.
`NUM_GIBBS_ITER=2000` in the resume file is the final iteration number, so it
generates the remaining 1000 draws.

This is a global-state/RNG restart check. Resuming a full forecast also requires
the corresponding local maps, spectral parameters, and legacy RNG state; the
global checkpoint alone does not provide that complete restart.

For final uncertainties use case 07 plus all three case-09 seeds, separately for
each resolution. Check rank-normalized split R-hat <= 1.01, bulk ESS >= 1000,
and tail ESS >= 400 for every reported global amplitude and free global shape.
Inspect representative local parameters, joint contours, residuals, and prior
boundaries too. Slow mixing or degeneracy can require longer runs or improved
blocking; four seeds alone do not establish convergence. A bias/coverage claim
additionally needs multiple independent noise realizations.

The Galactic maps in this package form the matched, constant-fiducial-SED export.
Recovery of this dataset validates the specified model and sampler. Original PySM
maps with additional spectral structure require separate model-mismatch tests.
The configured local priors and global shape bounds remain explicit analysis
choices, documented in `docs/fossil.md` in the Commander repository.
