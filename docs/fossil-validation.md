# FOSSIL validation and remote acceptance

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

## Progressive MPI runs

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
