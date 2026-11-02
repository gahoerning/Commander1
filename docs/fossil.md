# FOSSIL physical sky monopoles in Commander

This implementation adds spatially constant **physical sky components** to the
intensity likelihood. Each global component has one amplitude for the sky and,
where applicable, a small number of global spectral parameters. The Galactic
foreground amplitudes remain pixel dependent.

The current supported configuration has diagonal known RMS noise, monochromatic
channel responses, Stokes I, and all inputs in MJy/sr. Instrumental monopoles and
dipoles are fixed at zero, gains at one, and bandpass shifts at zero. Their prior
widths are zero. Empirical colour corrections and noise-amplitude fitting are
disabled. The code rejects unsupported instrumental settings when globals are
enabled; it does not silently reinterpret them.

## Scientific model and conventions

The reference blackbody at T0 = 2.725 K has already been subtracted from the
simulation. Its physical temperature-offset monopole and all distortion signals
remain in the maps and are fitted. Subtracting this known reference is an exact
change of data origin; removing a fitted per-band mean would destroy the signals
of interest and is not performed.

The present exported baseline has **no primary CMB anisotropy realization**.
`ENFORCE_ZERO_CL = .true.` describes that baseline. Adding primary CMB anisotropy
requires an explicit simulated realization, its input C_l, a consistent harmonic
CMB configuration, and new acceptance tests. It must not be claimed as present
merely because Commander contains a harmonic CMB sampler.

The matched Commander export uses separable Galactic morphologies and constant
fiducial SEDs: one modified-blackbody dust component, synchrotron and free-free
power laws, and two shifted AME templates. This is an exact-model validation
dataset. The original smoothed PySM maps are preserved separately; they have
additional spectral structure and require separate model-mismatch forecasts.
Successful recovery of the matched export does not demonstrate unbiased recovery
of the original maps.

| Global type | Amplitude and fiducial | Global spectral parameters |
| --- | --- | --- |
| `delta_t` | DeltaT/T0 = 1.2e-4 | None; linearized temperature derivative at T0 |
| `mu` | Chemical-potential amplitude = 2e-8 | None; fixed distortion shape |
| `tsz` | Compton y = 1.77e-6 | Parameter 01: kTe = 1.245 keV |
| `cib` | Modified-blackbody coefficient A = 0.346 MJy/sr | Parameter 01: beta = 0.86; parameter 02: T = 18.8 K |
| `crr` | A_CRR = 1 | Fixed recombination-radiation template |
| `co` | A_CO = 1 | Fixed aggregate extragalactic CO template |

The tSZ component contains the nonrelativistic y shape plus the Itoh corrections
through fourth order in theta_e. These terms share one y amplitude; they are not
independent components. The CIB coefficient multiplies
`x**(3+beta)/(exp(x)-1)`, with `x = h*nu/(k*T)`, and is not the intensity at a
reference frequency. CRR and aggregate CO amplitudes scale their complete
tabulated shapes. Template files contain frequency in GHz and intensity in
MJy/sr per unit amplitude. The export records its CO high-frequency continuation;
that continuation is part of the simulation assumption.

The full likelihood is

```
d[b,p] = local_sky[b,p] + sum_i F[b,i](theta) * a[i] + noise[b,p]
```

For fixed global shapes, all enabled global amplitudes are drawn jointly from
their Gaussian conditional posterior. With residual `r = d - local_sky`, the
sufficient statistics are

```
w[b] = sum_p mask[p] / sigma[b,p]**2
t[b] = sum_p mask[p] * r[b,p] / sigma[b,p]**2
Q = F.T @ diag(w) @ F + prior_precision
q = F.T @ t + prior_precision @ prior_mean
```

The posterior mean solves `Q @ mean = q`; its covariance is `Q**(-1)`. The code
uses a Cholesky solve and Gaussian fluctuations without explicitly inverting Q.
Global spectral shapes use conditional slice updates in sample mode. Optimize
mode returns conditional optima and must not be interpreted as a posterior
chain. A numerical sampling failure is fatal in the FOSSIL path; replacing failed
samples with ML values is not acceptable.

`FOSSIL_EXACT_LOCAL_SEDS = .true.` evaluates the supported local dust, power-law,
and shifted-AME spectra directly at the delta channel frequencies. This avoids
the ordinary foreground interpolation grid's errors imitating faint spectral
distortions. The generator enables it in every FOSSIL config, including the
no-globals regression fixture. It is independent of whether a local shape is
held fixed or sampled; existing non-FOSSIL configs keep their prior default.

Each sweep updates local/harmonic amplitudes conditional on the current globals,
then the global amplitudes jointly, then local spectral parameters, then global
spectral parameters. Residuals use the latest state at every block. The harmonic
block is inactive in the current no-primary-CMB baseline; the instrumental
sampling blocks are inactive throughout this experiment.

The instrument-fixed baseline is still susceptible to physical degeneracies.
A constant part of a free Galactic amplitude map can resemble a global SED.
Priors, frequency coverage, and spatial spectral variation determine whether the
data identify a component. A separate exact Gibbs block preserves the posterior
but can mix slowly across such degeneracies.

## Parameters

The two-digit component suffixes remain `01`, `02`, and so on. Band suffixes are
`01` through `99`, then `100` through `131`; they must never be truncated.

| Key | Meaning/default |
| --- | --- |
| `NUM_GLOBAL_COMPONENTS` | Number of physical global components; absent or zero preserves the legacy path |
| `GLOBAL_TYPE01` | One of `delta_t`, `mu`, `tsz`, `cib`, `crr`, `co` |
| `GLOBAL_LABEL01` | Output column label; use unique whitespace-free names |
| `GLOBAL_INIT01` | Initial physical amplitude |
| `GLOBAL_SAMPLE_AMPLITUDE01` | Enable amplitude fitting/sampling; default true; false fixes the initialized value |
| `GLOBAL_PRIOR_MEAN01` | Gaussian prior mean; defaults to the initial amplitude |
| `GLOBAL_PRIOR_SIGMA01` | Positive: Gaussian prior width; zero or negative: no Gaussian prior |
| `GLOBAL_SCALE01` | Positive numerical amplitude scale, with no change in the physical prior |
| `GLOBAL_TEMPLATE01` | Required template path for `crr` and `co`, relative to INPUT_DIRECTORY |
| `GLOBAL_INIT_PARAM03_01` | Initial spectral parameter 01 of global component 03 |
| `GLOBAL_SAMPLE_PARAM03_01` | Enable spectral-parameter fitting/sampling |
| `GLOBAL_PARAM_LOW03_01`, `GLOBAL_PARAM_HIGH03_01` | Required finite bounds when this shape parameter is enabled |
| `GLOBAL_PARAM_PRIOR_MEAN03_01` | Gaussian shape-prior mean |
| `GLOBAL_PARAM_PRIOR_SIGMA03_01` | Positive: Gaussian shape prior; zero or negative: flat within bounds |
| `GLOBAL_SLICE_STEPS` | Conditional slice transitions per shape update; default 8 |
| `GLOBAL_FIXED_INSTRUMENT` | Must be true for this release |
| `GLOBAL_RESTART_FILE` | Global checkpoint path, or `none` |
| `FOSSIL_EXACT_LOCAL_SEDS` | Required true for matched local spectra, avoiding coarse SED-grid interpolation errors |

**Legacy local spectral parameters use a different convention:** their Gaussian
prior standard deviation of zero freezes the parameter. Do not copy that
convention to `GLOBAL_PRIOR_SIGMA` or `GLOBAL_PARAM_PRIOR_SIGMA`; use the explicit
global sample flags to freeze global values.

The generated forecast uses flat global amplitude priors, kTe bounds [0,10] keV,
CIB beta bounds [0,3], and CIB temperature bounds [5,50] K. Gaussian local shape
priors are dust beta 1.53 +/- 0.2, dust temperature 21.7 +/- 3 K, synchrotron
RJ beta -3 +/- 0.3, free-free RJ beta -2.14 +/- 0.1, and AME peak frequencies
30 and 33.35 +/- 5 GHz. Synchrotron/free-free curvature stays fixed at zero.
These widths and bounds are explicit **analysis choices**, not extra fiducials
from the FOSSIL PDF. Change them for prior-sensitivity studies and report them.
The initial Galactic maps are perturbed from truth; fixed-local acceptance stages
instead select the exact truth maps from the export.

## Input data and noise

The matched export is self-contained under `outputs/commander/` in the simulation
project. Its `manifest.json` identifies maps, templates, truth/initial maps,
beams, masks, and channels. Portable parameter files use a single
`INPUT_DIRECTORY`, with relative paths from this root. Maps have one FITS field
in MJy/sr, Galactic coordinates, and HEALPix RING ordering. A filename example is

```
map_FOSSIL_0050GHz_fwhm096arcmin_nside0128_MJy-sr_noisy.fits
map_FOSSIL_0050GHz_fwhm096arcmin_nside0128_MJy-sr_noiseless.fits
rms_FOSSIL_0050GHz_fwhm096arcmin_nside0128_MJy-sr.fits
```

Both resolutions have a 1.6-degree Gaussian FWHM. Each RMS map contains the TXT
value converted from Jy/sr to MJy/sr, constant over its pixels. Noise is added
independently per pixel after smoothing the signal, and is not smoothed or
rescaled at NSIDE 64. Consequently NSIDE 128 has four times as many pixels and
approximately four times the Fisher information for a global signal at equal
mask coverage. Its ideal standard error is about half the NSIDE 64 value.

Commander uses a binary GAL070 mask, obtained by thresholding the fractional
mask at 0.5. The export records the resulting pixel count. This likelihood mask
is distinct from the fractional weights used to display a 70%-sky mean spectrum.

The channels are monochromatic centers. `BANDPASS_TYPE = 'delta'` exactly matches
that assumption; no finite instrumental response has been invented. All
regularization and fitted noise-amplitude widths are zero, so they cannot change
the supplied per-pixel RMS levels.

## Generated configurations and outputs

Run `examples/fossil/generate_parameters.py --help` and follow the example
workflow in `examples/fossil/README.md`. The generator supplies complete legacy
keys as well as the global keys; a short list of global settings alone is not a
runnable Commander configuration. It generates both definitive resolution
configs and progressive optimize/sample checks.

`global_components_c0001.txt` contains one row per accepted iteration, independent
of expensive map-output thinning. Its named header describes component amplitudes
and their shape parameters. Read the header rather than assuming a column order.
`global_spectra_c0001.txt` contains the latest total and component monopole spectra
in MJy/sr. These scalars define spatially constant full-sky maps without repeated
copies of the same value for every pixel. Local map products follow the existing
Commander output path, with map thinning of 100 in long configurations.

`global_restart_c0001.txt` preserves the global parameters, iteration, and global
RNG state. It supports exact continuation of the globals-only fixtures with a
fixed local sky. It is **not a complete checkpoint of a full pixel Gibbs chain**:
that additionally requires matching local amplitudes, local spectral maps,
harmonic state when used, and the legacy samplers' RNG states. A full forecast
must not be resumed by loading only this global file.

## Forecast acceptance and mixing

Use at least four independent seeds for final diagnostics. A production result
requires rank-normalized split R-hat <= 1.01, bulk ESS >= 1000, and tail ESS >= 400
for every reported global amplitude and free shape parameter. Also inspect traces,
joint contours, selected local spectral pixels/regions, prior boundaries, and
residual spectra/maps. The supplied diagnostics check scalar global chains;
they do not certify convergence of every local pixel map.

Record ESS per wall-clock hour as local amplitudes and then local shapes are
released. If a 10000-iteration post-warmup chain has integrated autocorrelation
time above 100, or the four-chain run fails the ESS thresholds while frozen-local
tests pass, treat this as a mixing failure requiring investigation. Extend a run
only when the measured ESS rate makes the target feasible. Persistent coupled
global/local mean modes justify a joint amplitude CG/Schur-complement update or
an exact interweaving move; global shape degeneracies may need blocked nonlinear
updates as well. Such an update improves exploration; it cannot identify an
intrinsically degenerate model or replace justified priors.

The remote validation checklist is in `docs/fossil-validation.md`. Generating
configs, passing Python checks, or compiling isolated kernels does not constitute
a successful end-to-end Commander/MPI run.
