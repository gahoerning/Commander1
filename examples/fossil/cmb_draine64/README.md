# Primary CMB and Draine free-free: 64-channel component separation

This versioned suite contains 26 complete `.txt` parameter files for the new
`planck_draine_gaunt_plus_primary_cmb` injection. The original PDF-only maps and parameter files
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
Its SHA256 is `8a8f651e259d5c964a93864d0dc7f9f43e70b69d88c42001e6044bb720b91ad1`. This fiducial spectral-design proxy does not
replace a spatial Fisher calculation or posterior convergence checks.

## Paths and manual execution

All files contain expanded absolute cluster paths:

```sh
FOSSIL_ROOT=/mn/stornext/d23/cmbco/AST9240/2024/gabrieah/fossil_component_separation
INPUT_DIRECTORY="$FOSSIL_ROOT/fossil_inputs_cmb_draine"
cd "$FOSSIL_ROOT"
PARFILE="$INPUT_DIRECTORY/01_globals_optimize_nside0064_64bands_seed20261102.txt"
./comm_run.sh --check "$PARFILE"
./comm_run.sh "$PARFILE"
```

Copy these parameter files into the new input directory beside the new matched
map/template package. Chains go into `tests_64channels_cmb_draine_fullsky/nside64/` or
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
| 01. Recover the six physical monopole amplitudes with a completely known local sky | Noiseless; 12 optimize | Global amplitudes; Galactic maps/shapes and primary CMB fixed at truth | [01 / NSIDE 64](01_globals_optimize_nside0064_64bands_seed20261102.txt) | [01 / NSIDE 128](01_globals_optimize_nside0128_64bands_seed20261102.txt) |
| 02. Validate the exact Gaussian conditional for physical monopole amplitudes | Noisy; 2000 sample | Global amplitudes; all local quantities including primary CMB fixed at truth | [02 / NSIDE 64](02_globals_sample_nside0064_64bands_seed20261102.txt) | [02 / NSIDE 128](02_globals_sample_nside0128_64bands_seed20261102.txt) |
| 03. Exercise joint local-amplitude/CMB CG coupled to the global block | Noiseless; 12 optimize | Galactic amplitudes, harmonic CMB realization, global amplitudes; all spectral shapes fixed | [03 / NSIDE 64](03_amplitudes_optimize_nside0064_64bands_seed20261102.txt) | [03 / NSIDE 128](03_amplitudes_optimize_nside0128_64bands_seed20261102.txt) |
| 04. Explore global spectral shapes with a completely known local sky | Noisy; 5000 sample | Global amplitudes and tSZ/CIB shapes; Galactic quantities and primary CMB fixed at truth | [04 / NSIDE 64](04_global_shapes_sample_nside0064_64bands_seed20261102.txt) | [04 / NSIDE 128](04_global_shapes_sample_nside0128_64bands_seed20261102.txt) |
| 05. Study amplitude degeneracies with a free primary-CMB realization | Noisy; 5000 sample | Galactic amplitudes, harmonic CMB, global amplitudes; all spectral shapes fixed | [05 / NSIDE 64](05_amplitudes_sample_nside0064_64bands_seed20261102.txt) | [05 / NSIDE 128](05_amplitudes_sample_nside0128_64bands_seed20261102.txt) |
| 06. Optimize the complete matched component-separation model | Noiseless; 50 optimize | Galactic amplitudes/shapes, harmonic CMB realization, global amplitudes/shapes | [06 / NSIDE 64](06_forecast_optimize_nside0064_64bands_seed20261102.txt) | [06 / NSIDE 128](06_forecast_optimize_nside0128_64bands_seed20261102.txt) |
| 07. Run the first full posterior forecast chain | Noisy; 10000 sample | Galactic amplitudes/shapes, harmonic CMB realization, global amplitudes/shapes | [07 / NSIDE 64](07_forecast_sample_nside0064_64bands_seed20261102.txt) | [07 / NSIDE 128](07_forecast_sample_nside0128_64bands_seed20261102.txt) |
| 08. Optional disabled-global software-path control | Noiseless; 2 optimize | No globals fitted; all local quantities including primary CMB fixed | [08 / NSIDE 64](08_no_globals_regression_nside0064_64bands_seed20261102.txt) | [08 / NSIDE 128](08_no_globals_regression_nside0128_64bands_seed20261102.txt) |
| 09. Run an additional independent full forecast chain (seed 20261103) | Noisy; 10000 sample | Same complete model as case 07; independent seed | [09 / NSIDE 64](09_forecast_sample_nside0064_64bands_seed20261103.txt) | [09 / NSIDE 128](09_forecast_sample_nside0128_64bands_seed20261103.txt) |
| 09. Run an additional independent full forecast chain (seed 20261104) | Noisy; 10000 sample | Same complete model as case 07; independent seed | [09 / NSIDE 64](09_forecast_sample_nside0064_64bands_seed20261104.txt) | [09 / NSIDE 128](09_forecast_sample_nside0128_64bands_seed20261104.txt) |
| 09. Run an additional independent full forecast chain (seed 20261105) | Noisy; 10000 sample | Same complete model as case 07; independent seed | [09 / NSIDE 64](09_forecast_sample_nside0064_64bands_seed20261105.txt) | [09 / NSIDE 128](09_forecast_sample_nside0128_64bands_seed20261105.txt) |
| 11. Run the first half of the global-state/RNG restart check | Noisy; 1000 sample | Same fixed-local conditional as case 02 | [11 / NSIDE 64](11_globals_sample_split1000_nside0064_64bands_seed20261102.txt) | [11 / NSIDE 128](11_globals_sample_split1000_nside0128_64bands_seed20261102.txt) |
| 12. Resume the global-state/RNG check to iteration 2000 | Noisy; iterations 1001–2000 | Restore the global state and RNG from case 11; fixed local sky | [12 / NSIDE 64](12_globals_sample_resume2000_nside0064_64bands_seed20261102.txt) | [12 / NSIDE 128](12_globals_sample_resume2000_nside0128_64bands_seed20261102.txt) |

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
