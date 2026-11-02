# FOSSIL configurations

## Active native-unit CMB/Draine forecast

For new runs use `native_units64/` (64 channels, 64 MPI ranks) or
`native_units131/` (131 channels). Transfer the derived simulation package as
`fossil_inputs_native_units/`. Band maps and RMS use `uK_cmb`, Galactic
reference amplitudes use native `uK_ant`, and the CMB reference amplitude uses
`uK_cmb`. The physical master remains unchanged in MJy/sr.

The existing free-free equation, lookup table and foreground module are
unchanged. The new global block converts physical MJy/sr SED predictions into
the declared band unit; its physical parameters and saved spectral snapshots
retain their original units. Both maps and RMS receive the same conversion.
Use the updated diagnostic notebook to display native maps in MJy/sr.

This interface requires updating and recompiling Commander. The new chain
directories contain `native_units_fullsky`, preserving earlier runs. Start
with case01 at NSIDE64, inspect diagnostics, then case03 and case02. Only the
restart case12 requires a previous case11 execution.

## Archived PDF-only baseline configurations

These configs target the matched simulation export described in
[`docs/fossil.md`](../../docs/fossil.md). They contain no primary CMB anisotropy.
All map and RMS inputs are MJy/sr, Gaussian FWHM is 96 arcmin, and instrumental
offsets, gains, shifts, and noise scaling are fixed.

For one owl31–35 node with 64 physical cores, use the ready-made
[`owl64/`](owl64/README.md) parameter files. They include every validation stage,
optimization, sampling, four independent forecast chains, and global restart
controls for both NSIDE 64 and 128. Each file uses the same 64-channel selection
and requires 64 MPI ranks. They can be executed directly with your own launcher.
No Commander rebuild or new map generation is needed. The `.txt` inputs use
12 named sections, aligned values, English explanations, and `VERBOSITY=3`.

The selection in [`channel_selection_64.json`](channel_selection_64.json) removes
67 channels using the noise-weighted fiducial spectral Jacobian, preserving
amplitude references and frequency coverage. It includes Galactic and global
spectral directions. This is a spectral-information design heuristic; a full
spatial likelihood and priors can change the relative value of each channel.
The existing 131-channel forecast configurations remain available separately.

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
layout. Regenerate on relocation or set INPUT_DIRECTORY to the new path. The checked-in `fossil_nside128_forecast_sample.txt` and
`fossil_nside64_forecast_sample.txt` use that portable layout.

Keep the original full-length HEALPix tables `pixel_window_n0064.fits` and
`pixel_window_n0128.fits` in INPUT_DIRECTORY. They contain temperature and
polarization windows through ell = 4 * NSIDE (257 and 513 rows). The
`PIXEL_WINDOW` parameter points to the corresponding native table. Commander
also calls its legacy `read_pixwin` routine using that standard filename in the
process working directory. The supplied `comm_run.sh` validates the native
input table and creates a working-directory link before MPI starts; its
`--check` mode reports the planned link without creating it. A different launcher
must provide the same working-directory file or link. The shorter
`nside64/pixel_window.fits` and `nside128/pixel_window.fits` export tables do not
have enough rows for that legacy lookup. This deployment correction requires no
Commander rebuild or sky/noise-map regeneration. Original-file hashes and sources
are recorded in `pixel_window_provenance.json` in the input package.

| Stage | Data | Free quantities | Purpose |
| --- | --- | --- | --- |
| `no_globals_regression` | Noiseless | No new global components; frozen locals | Compare patched and legacy outputs; not a recovery test of this intentionally incomplete sky model |
| `globals_optimize` | Noiseless | Global amplitudes | Exact fixed-shape amplitude recovery; 12 optimize iterations |
| `globals_sample` | Noisy | Global amplitudes | Analytic Gaussian posterior mean/covariance; 2000 actual draws |
| `amplitudes_optimize` | Noiseless | Global and local amplitudes, fixed shapes | Exercise the local CG/global conditional integration; 12 optimize iterations |
| `global_shapes_sample` | Noisy | Global amplitudes and shapes | Validate global shape exploration with a fixed local sky; 5000 iterations |
| `amplitudes_sample` | Noisy | Global and local amplitudes | Measure mixing introduced by local/global constant-mode correlations; 5000 iterations |
| `forecast_optimize` | Noiseless | Global and local amplitudes and shapes | Check deterministic optimization of the complete model; 50 iterations |
| `forecast_sample` | Noisy | Global amplitudes/shapes and local amplitudes/shapes | Full matched-model forecast; 10000 iterations per independent chain |

A single 131-band chain requires exactly 131 MPI ranks when
`NUM_PROC_PER_BAND=1`. The legacy decomposition requires
`NUM_GROUPS * NUMBAND * NUM_PROC_PER_BAND` ranks. A second pixel partition uses
262 ranks. An arbitrary one-rank invocation cannot run all 131 bands.

For one forecast chain:

```sh
mpirun -np 131 /path/to/commander examples/fossil/generated/fossil_nside128_forecast_sample.txt
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
