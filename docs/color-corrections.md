# Colour corrections from the reconstructed foreground sky

Set `APPLY_COLOR_CORRECTIONS = T` and enable each desired band with
`USE_COLOR_CORRECTIONxx = T`, where `xx` is its two-digit band index.
When the global switch is absent or false, corrections are disabled.

## Update interval

```text
COLOR_CORR_UPDATE_INTERVAL = 10
```

This positive integer counts **outer Commander iterations**, in both `sample`
and `optimize` modes, not individual parameter proposals or internal substeps.
The default is `1`. Corrections are always calculated from the initial model
before the first fit. With an interval of 10, they are refreshed before iterations
11, 21, 31, and so on, using the current amplitudes and smoothed parameter maps.
With an interval larger than the run length, only the initial calculation occurs.
No extra refresh is performed when writing outputs or changing a single parameter.
A retried iteration uses the same iteration number and restores the corresponding
model state before its scheduled refresh.

The correction is held fixed between these updates to retain the existing linear
amplitude solvers and conditional spectral samplers. **This is an iterative
plug-in approximation, not exact posterior sampling for a response that depends
on the sampled amplitudes.** In optimize mode it is an alternating update scheme;
convergence of the correction and fitted parameters still needs checking. In
sample mode, evaluate sensitivity to the interval and the approximation before
using chains for scientific inference. Interval 1 does not remove the approximation.

## Reconstructed spectrum

For each enabled band and full-sky pixel, sum the physical continuum foreground
components at two frequencies in the internal `uK_RJ` units:

```text
I_fg(nu) = sum_c I_c(nu)
Q_fg(nu) = sum_c Q_c(nu)
U_fg(nu) = sum_c U_c(nu)
P_fg(nu) = sqrt(Q_fg(nu)^2 + U_fg(nu)^2)

alpha_I = 2 + log(I_fg(nu_high)/I_fg(nu_low)) / log(nu_high/nu_low)
alpha_P = 2 + log(P_fg(nu_high)/P_fg(nu_low)) / log(nu_high/nu_low)
```

The CMB is excluded. Component masks apply to each Stokes parameter. Faraday
components are reconstructed as rotated Q and U, using their actual Q0, U0,
beta and RM, before summing the components. In particular, the algorithm does
not sum the individual polarized intensities. Internal foreground amplitudes
have already been converted to RJ units when the initial maps are loaded.

Band-specific nuisance templates (including monopoles/dipoles) and spectral-line
components have no smooth continuum SED at the two endpoints. They are excluded
from the index estimate and retain their original response. The CMB response
also remains unchanged.

The extra `+2` converts the RJ temperature index to the flux-density index
**alpha** used by the polynomial. This follows section 2.2.1 of the Roke
Cepeda-Arroita thesis supplied with this change, especially equations 2.5-2.6
and Table 2.2 (pages 66-69).

## Per-band parameters

For band `01`:

```text
USE_COLOR_CORRECTION01 = T

# Optional paired endpoint overrides, in GHz.
COLOR_CORR_FREQ_MIN01 = 4.5
COLOR_CORR_FREQ_MAX01 = 5.5

# Polynomial K_X(alpha) = C0_X + C1_X*alpha + C2_X*alpha^2.
COLOR_CORR_C0_I01 = 0.99970
COLOR_CORR_C1_I01 = -0.00052
COLOR_CORR_C2_I01 = 0.00127
COLOR_CORR_C0_Q01 = 0.99776
COLOR_CORR_C1_Q01 = -0.00716
COLOR_CORR_C2_Q01 = 0.00090
COLOR_CORR_C0_U01 = 1.00465
COLOR_CORR_C1_U01 = 0.01650
COLOR_CORR_C2_U01 = 0.00240

# Optional minimum endpoint I or P, in uK_RJ; default shown.
COLOR_CORR_MIN_SIGNAL01 = 1.0e-12
```

The example coefficients are specifically the **C-BASS North average I/Q/U
coefficients at 4.76 GHz** in that thesis. They are not universal coefficients
for other instruments, reference frequencies or bandpasses. Set `FREQ_C01` to
the effective frequency associated with your coefficients. The two illustrative
endpoint values above should be replaced with the limits appropriate to your data.

If both endpoint parameters are omitted, the code uses the minimum and maximum
positive-response frequencies of the loaded input bandpass after its existing
loading threshold. A delta band therefore needs two explicit endpoint frequencies.
Providing only one endpoint, reversed/nonpositive limits, non-finite settings,
a negative signal floor, or an update interval below 1 is an error.

Q and U use the **same alpha_P** but their own polynomial coefficients. For a run
with `POLARIZATION = F`, only the I coefficients are required. For polarization
runs, all three coefficient sets are required.

## Forward-model convention

The polynomial K is the **divisive data correction** described in the thesis:
`corrected_map = observed_map / K`. Consequently, the forward model compared
with the uncorrected observed map is:

```text
model_X = unchanged_CMB_and_line_terms
        + gain * RJ_to_data_units * K_X * foreground_X(nu_c)
        + unchanged_band_templates
```

For enabled bands, the polynomial replaces continuum bandpass integration; it
is not multiplied onto an already integrated continuum response. Q and U retain
their signs. The Faraday block sampler uses separate corrected response rows for
Q and U. This scalar-polynomial prescription does not itself model Faraday
bandwidth depolarization; use the original bandpass-integrated response where
that effect needs to be represented explicitly.

Polynomial coefficients describe a fixed bandpass. Bandpass-shift sampling
(`BP_RMS > 0`) is rejected for a colour-corrected band; other bands may still
sample their bandpasses. The nominal input endpoints stay fixed between updates.

## Undefined indices and diagnostics

If either endpoint I is nonpositive or below the configured floor, K_I is set
to 1. If either endpoint P is below the floor, K_Q and K_U are set to 1.
Non-finite inputs and nonpositive/non-finite polynomial outputs also use unity
for the affected Stokes parameter. A negative Q or U is not an invalid input.

Every refresh prints the largest absolute change in K and the number of unity
fallbacks for I, Q and U in each band. Masked/empty pixels contribute to these
counts. Disabled temperature sampling does not contribute to the I count.
The in-memory alpha maps mark invalid estimates with the HEALPix missing value.

## Migration from the OWLS implementation

`COLOR_CORR_COMPONENTxx` is no longer used. Replace the old shared
`COLOR_CORR_C0xx/C1xx/C2xx` keys with the explicit I/Q/U keys above. The old
implementation evaluated the polynomial at a single power-law component's beta;
the new implementation evaluates it at the reconstructed foreground alpha.

## Tests

Run `python3 tests/run_color_tests.py` with a Fortran 2008 compiler available as
`gfortran`, or set `FC` and optional `FFLAGS` for your compiler environment.
The tests compile the numerical kernel and actual production routines in small
integration fixtures. They cover mixed-component slopes, RJ-to-alpha conversion,
polarization cancellation and rotation, CMB exclusion, independent Q/U factors,
forward-model units/gain, frozen factors and refresh, update cadence, global pixel indexing,
Faraday block sampling, and invalid configurations. Cluster-specific I/O,
smoothing, spline routines and single-rank MPI are substituted in these fixtures;
a full Commander build and multi-rank cluster validation are still required.
