# FOSSIL configurations

These configs target the matched simulation export described in
[`docs/fossil.md`](../../docs/fossil.md). They contain no primary CMB anisotropy.
All map and RMS inputs are MJy/sr, Gaussian FWHM is 96 arcmin, and instrumental
offsets, gains, shifts, and noise scaling are fixed.

Copy the self-contained `outputs/commander/` export to the machine running
Commander. Then, from the repository root, generate complete configs with paths
appropriate to that machine:

```sh
python3 examples/fossil/generate_parameters.py \
  --manifest /data/fossil/manifest.json \
  --input-directory /data/fossil \
  --output examples/fossil/generated --verify-inputs
```

For a portable layout at `outputs/commander`, omit `--input-directory`. Input
paths resolve from INPUT_DIRECTORY; CHAIN_DIRECTORY resolves from the process
working directory. The generator does not expand environment variables inside
the Fortran parameter file. Absolute input roots, as above, work from any process
directory. Relative roots require running from the directory containing that
layout. Regenerate on relocation or set INPUT_DIRECTORY to the new path. The checked-in `fossil_nside128_forecast_sample.par` and
`fossil_nside64_forecast_sample.par` use that portable layout.

| Stage | Data | Free quantities | Purpose |
| --- | --- | --- | --- |
| `no_globals_regression` | Noiseless | No new global components; frozen locals | Compare patched and legacy outputs; not a recovery test of this intentionally incomplete sky model |
| `globals_optimize` | Noiseless | Global amplitudes | Exact fixed-shape amplitude recovery; 12 optimize iterations |
| `globals_sample` | Noisy | Global amplitudes | Analytic Gaussian posterior mean/covariance; 2000 actual draws |
| `global_shapes_sample` | Noisy | Global amplitudes and shapes | Validate global shape exploration with a fixed local sky; 5000 iterations |
| `amplitudes_sample` | Noisy | Global and local amplitudes | Measure mixing introduced by local/global constant-mode correlations; 5000 iterations |
| `forecast_sample` | Noisy | Global amplitudes/shapes and local amplitudes/shapes | Full matched-model forecast; 10000 iterations per independent chain |

A single 131-band chain requires exactly 131 MPI ranks when
`NUM_PROC_PER_BAND=1`. The legacy decomposition requires
`NUM_GROUPS * NUMBAND * NUM_PROC_PER_BAND` ranks. A second pixel partition uses
262 ranks. An arbitrary one-rank invocation cannot run all 131 bands.

For one forecast chain:

```sh
mpirun -np 131 /path/to/commander examples/fossil/generated/fossil_nside128_forecast_sample.par
```

For final inference generate four independent seeds, each in a separate output
folder. Keep the corresponding scalar traces and timing information:

```sh
python3 examples/fossil/generate_parameters.py --manifest /data/fossil/manifest.json \
  --input-directory /data/fossil --output configs/seed20261103 \
  --stage forecast_sample --seed 20261103
```

The other three seeds can be 20261102, 20261104, and 20261105. NSIDE 64 uses the
same per-pixel TXT noise, not a degradation of the NSIDE 128 noisy maps. Compare
resolution-dependent forecast uncertainties with this explicit information-content
change in mind.

Use `tests/fossil_remote_acceptance.py` for staged remote runs and
[`docs/fossil-validation.md`](../../docs/fossil-validation.md) for acceptance
criteria and restart/MPI checks. Sample mode has no permitted silent ML fallback.
