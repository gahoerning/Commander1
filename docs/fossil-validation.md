# FOSSIL validation and remote acceptance

## Operative native-unit validation

Use the `examples/fossil/native_units64/` or `native_units131/` parameter
files with the derived `fossil_inputs_native_units` package. The physical
master remains in MJy/sr; band data/RMS use `uK_cmb`, Galactic amplitudes use
`uK_ant`, and the CMB reference amplitude uses `uK_cmb`. Original internal RJ
units and the complete foreground/free-free source are unchanged from
commit8ffea73. SI-matched CMB band calibration is an opt-in boundary conversion.

The global interface converts cached and trial-shape MJy/sr SEDs to the band
unit. Priors/scalar parameters remain physical, and global spectral snapshots
remain explicitly MJy/sr. Isolated production tests cover Gaussian conditional
means/covariances, native/physical likelihood invariance, spectral-shape moves,
exact checkpoint continuation, and physical spectral-file units. Standard
free-free spline/derivative and CMB low-mode tests exercise the unchanged path.

```sh
python3 tests/run_fossil_cmb_tests.py
python3 tests/run_fossil_native_unit_tests.py /data/fossil_physical_master
```

The native test's positional package is the archived physical MJy/sr master,
used as an independent reference. It does not require running MPI. The
diagnostic notebook uses the native package's verified conversion metadata and
physical truth maps for actual chain inspection. Conditional references include
the tiny deterministic difference between the analytic injection and the
unchanged 100-point free-free temperature table. `NUMGRID_SPLINE` does not
control that hardcoded foreground table.

A source-hash check and compiled feature marker reject stale installations;
they do not prove executable build identity or full MPI acceptance. Recompile
after updating source, run case01 in the new chain directory, and inspect its
outputs before progressing.

## What local checks establish

`tests/fossil_parameter_tests.py` uses the standard library to audit generated
131-band configs, every direct required legacy key in the relevant production
readers, units, fixed instrumental settings, local freeze semantics, the release
sequence, and rejection of mismatched beam/CMB assumptions.

`tests/fossil_analysis_tests.py` needs NumPy. It checks analytic posterior
acceptance calculations and scalar-chain diagnostics using synthetic chains.
It verifies that repeated ML values, incorrect covariance, sticky chains, and
chains in different modes fail the appropriate criteria. These are tests of the
acceptance tooling, not runs of Commander.

```sh
python3 tests/fossil_parameter_tests.py
python3 tests/fossil_analysis_tests.py
python3 tests/fossil_index_tests.py --static-only
```

With a configured Fortran compiler, run `tests/fossil_index_tests.py` without
`--static-only`. It compiles the actual production formatting procedures and
checks identifiers 01 through 1000, reallocation, and legacy spelling. The static
audit verifies deferred-length tag callers and rejects direct fixed-width band-ID
formatting. Neither check replaces a full 131-band MPI run.

The production-kernel and extracted integration tests use the configured Fortran
compiler:

```sh
python3 tests/run_fossil_sed_tests.py /data/fossil/spectra --require-compiler
python3 tests/run_fossil_global_tests.py /data/fossil/spectra
python3 tests/run_fossil_local_sed_tests.py /data/fossil --require-compiler
python3 tests/run_fossil_integration_tests.py
python3 tests/fossil_index_tests.py
```

The integration runner exercises actual production residual/CG procedures with
controlled dependencies. A full Commander build still needs MPI, HEALPix, FITS,
and numerical-library toolchain. Record the compiler and source revision used
for remote acceptance. Do not mark the following tests passed until their actual
executable outputs exist.

## Opt-in primary CMB and Draine free-free extension

The versioned `outputs/commander_cmb_draine` physical injection adds primary
temperature anisotropy and the normalized Draine
free-free amplitude/Te model. The original PDF-only injection and its static
parameter suite are preserved. Recompile Commander before using the extension;
full HEALPix/MPI acceptance of this new path remains pending until actual
cluster outputs are inspected.

`FOSSIL_SI_CMB=true` makes the physical CMB intensity derivative use exact SI
constants, consistently with the injected CMB and global DeltaT monopole. Its
default is false. The internal legacy RJ definition and existing local
foreground normalization remain unchanged. `FOSSIL_CMB_MIN_L=2` fixes the
physical harmonic monopole and dipole exactly to zero through the covariance
and CG transformations; its default zero preserves the earlier behavior. Keep
`CMB_MONO_DIPOLE_MASK=none` so no removed CMB modes are transferred to fixed
instrumental offsets. The unwindowed TT input is D_ell in uK_CMB squared;
`CONDITION_ON_CL=true` fixes this covariance while the harmonic realization is
free in the coupled amplitude/full-forecast stages.

The existing free-free component has one native RJ reference amplitude and one electron
temperature in K. It uses `T_E_PRIOR_*` keys and its original spectral table;
EM is input provenance, not another independently fitted amplitude. Primary CMB
is fixed to its known truth in the analytic global-only controls. Their global
Gaussian sufficient statistics remain valid after the known local/CMB sky is
subtracted. In case03 the noiseless joint solve is still a prior-dependent CMB
Wiener/MAP solve, so exact recovery of every amplitude or identically zero
residuals is not an appropriate acceptance requirement.

```sh
python3 tests/test_fossil_cmb_parameters.py
python3 tests/run_fossil_cmb_tests.py
```

The compiled CMB fixtures exercise the production SI conversions, equivalence
to the physical DeltaT response, normalized free-free Te derivative, covariance
and D_ell initialization, exact low-mode CG output, and harmonic/pixel CMB
preconditioner separation. They use explicit small dependencies and do not
certify a complete MPI/HEALPix run. Repeat remote CG, residual, prior and
posterior-mixing checks with the new input package.

## Progressive MPI runs

The [`examples/fossil/owl64/`](../examples/fossil/owl64/README.md) directory
contains all static parameter files for a single owl31–35 node, using 64 MPI
ranks and a resource-limited 64-channel separation/forecast. Use your own
launcher with these files, or `--band-count 64` in the optional generator and
acceptance tools. The frequency set is selected by a noise-weighted fiducial
spectral-information heuristic, including global and Galactic amplitude/shape
directions, fixed reference channels and coverage constraints. Both the
map/reference-band IDs and analytic posterior statistics are mapped to the
selected original frequencies. Per-pixel sensitivity and the original maps do
not change. Conclusions for the reduced analysis refer to these 64 channels.

The supplied files progressively exercise the global block, local CG amplitudes,
global spectral shapes, the complete nonlinear model, four independent seeds,
and global RNG checkpoint continuation. The directory README explains the free
and fixed quantities and expected evidence for every case. The 131-channel
configuration remains the default when no band-count option is specified.

Run in a fresh directory on the configured machine. FOSSIL path buffers support
511 bytes; the generator/remote runner rejects longer resolved input paths.
The older executable used for regression still has 127-character limits, so
choose short input and output roots for that control. Never reuse a chain directory for a fresh run: the scalar trace
is append-only. The runner refuses to append to an existing nonempty directory.

```sh
python3 tests/fossil_remote_acceptance.py run \
  --commander /path/to/commander \
  --manifest /data/fossil/manifest.json \
  --run-root /scratch/fossil_acceptance \
  --nside 64 --stage globals_optimize globals_sample \
  --launcher 'srun --ntasks={ranks}'
```

The optimize fixture holds local amplitudes and local/global shapes at exact
truth and starts the global amplitudes away from truth. It must recover the
noiseless amplitudes. The sample fixture uses the noisy data and the same known
local sky. `global_response.csv` and `global_reference.json` contain independent
exported sufficient statistics under the actual binary mask. The acceptance
tool computes the exact Gaussian conditional, whitens the actual trace, and
checks means and every covariance entry at a six-standard-error threshold.
At least 500 retained draws are required. Identical ML values cannot pass.

Repeat both resolutions and then release parameters progressively:

```sh
python3 tests/fossil_remote_acceptance.py run --commander /path/to/commander \
  --manifest /data/fossil/manifest.json --run-root /scratch/fossil_shapes \
  --nside 128 64 --stage global_shapes_sample amplitudes_sample \
  --launcher 'srun --ntasks={ranks}'
```

For these exploratory stages the runner checks finite output, fixed parameters,
actual variation of free shapes, and complete scalar traces. It does not claim
nonlinear convergence from one chain. Run four independent full forecast chains
and evaluate the scalar traces, for example:

```sh
python3 tests/fossil_remote_acceptance.py diagnose \
  --chain-directory /scratch/fossil_seed1 /scratch/fossil_seed2 \
                    /scratch/fossil_seed3 /scratch/fossil_seed4 \
  --warmup 1000 --output forecast_diagnostics.json
```

Each directory must contain `global_components_c0001.txt` from a different seed.
All reported free scalars must have split rank-normalized R-hat <= 1.01, bulk
ESS >= 1000, and tail ESS >= 400. Also assess representative local amplitudes
and spectral parameters from the map outputs. The generator's prior widths and
global bounds must be reported and varied in prior-sensitivity tests.
For a conditional control with fixed shapes, declare those columns explicitly
with `--fixed-parameter y_p01 cib_p01 cib_p02`. A constant free parameter is a
failure, not an automatically accepted fixed quantity.

## MPI weighting, regression, restart, and noise

1. Repeat `globals_sample` with `--proc-per-band 2` in a separate run root. It
   then uses 262 ranks. Both decompositions must match the same analytic
   covariance. This detects duplicated or missing pixel/band contributions;
   different reduction order need not produce bit-identical random sequences.
2. Keep a known legacy executable and run `no_globals_regression` with
   `--legacy-commander /path/to/legacy_commander`. The tool compares actual FITS
   data arrays exactly. This explicit old-behavior control uses 99 channels and
   disables `FOSSIL_EXACT_LOCAL_SEDS` in both runs: the old executable cannot
   address band 100 or reproduce the new direct-SED option. Generated matched
   FOSSIL configs retain 131 channels and the exact-SED flag; the global stages
   exercise those capabilities separately. The sky contains monopoles omitted by this
   deliberately incomplete regression model, so a large chi-square is not a
   recovery failure. This check tests unchanged behavior with globals disabled.
3. For exact globals-only restart, compare a 2000-iteration `globals_sample`
   run with a 1000-iteration run resumed using its absolute
   `GLOBAL_RESTART_FILE`, with total `NUM_GIBBS_ITER=2000`. Keep all local maps,
   masks, responses, priors, and seed settings the same. Concatenated scalar
   rows should match uninterrupted rows exactly. The global restart contains
   its RNG state. A full pixel-chain restart needs additional legacy state and
   is not covered by this claim.
4. Confirm every exported RMS value is `sensitivity_TXT * 1e-6` at both
   resolutions. Global precision matrices should follow the exported binary
   pixel counts with unchanged per-pixel noise; they should not be identical
   between resolutions. Do not smooth RMS maps or noisy realizations.
5. Inspect residual maps, the total residual spectrum, and pulls over many
   independent noise realizations before publishing bias/coverage statements.
   A single injection's truth need not lie within one posterior standard
   deviation. Matched-model recovery and model-mismatch bias are separate runs.

## Interpreting slow mixing

The separated global amplitude block is an exact conditional transition. It
does not remove local/global degeneracies. Record ESS/hour and compare the
frozen-local, free-local-amplitude, and free-local-shape stages. If releasing local
quantities produces integrated autocorrelation times above 100 or prevents the
four-chain ESS targets within the intended runtime, inspect coupled local means
and global amplitudes and prioritize a joint amplitude/Schur-complement or exact
interweaving update. Shape correlations can additionally require blocked shape
updates. Check rank/conditioning and prior boundaries before assuming a longer
chain or a different amplitude sampler can resolve non-identifiability.
