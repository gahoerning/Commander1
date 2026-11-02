# Full 131-channel primary-CMB/Draine forecasts

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
| NSIDE 64 | [07_forecast_sample_nside0064_131bands_seed20261102.txt](07_forecast_sample_nside0064_131bands_seed20261102.txt) | `/mn/stornext/d23/cmbco/AST9240/2024/gabrieah/fossil_component_separation/tests_131channels_cmb_draine_fullsky/nside64/forecast_sample_131bands_seed20261102` |
| NSIDE 128 | [07_forecast_sample_nside0128_131bands_seed20261102.txt](07_forecast_sample_nside0128_131bands_seed20261102.txt) | `/mn/stornext/d23/cmbco/AST9240/2024/gabrieah/fossil_component_separation/tests_131channels_cmb_draine_fullsky/nside128/forecast_sample_131bands_seed20261102` |

Copy the desired file into `/mn/stornext/d23/cmbco/AST9240/2024/gabrieah/fossil_component_separation/fossil_inputs_cmb_draine/` and pass its full
path to your launcher. The launcher must allocate 131 ranks. Each file provides
10000 iterations with seed 20261102, not a convergence guarantee; create
separate files/seeds/chains for independent full-band inference.
