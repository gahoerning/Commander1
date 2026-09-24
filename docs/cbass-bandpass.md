# C-BASS North Stokes I bandpass

`BANDPASS_TYPExx = 'CBASS'` integrates the Stokes I model through a measured
bandpass using the flux-density calibration convention verified against
Cepeda-Arroita's thesis, section 2.2.1 and Table 2.2. It is restricted to
`POLARIZATION = F` and `FREQ_UNITxx = 'uK_ant'` (Rayleigh-Jeans microkelvin).
It does not introduce or enable any sky component.

Use `examples/bandpass-cbass-I.par` as a fragment in a complete parameter file.
Resolve `BANDPASSxx` relative to the run's configured input directory as for the
other bandpass types. The required `FREQ_Cxx` is in GHz, and the required
`BANDPASS_CAL_INDEXxx` is a **flux-density** spectral index, not an RJ index.
The example uses 4.76 GHz and -0.299, matching the confirmed map convention;
neither number is hard-coded in the implementation. Changing these values
changes the map convention represented by the model; it does not recalibrate an
existing data map. Supply values consistent with the input map.

Set `A2Txx`, `F2Txx`, and `A2SZxx` to finite nonpositive values for automatic
conversions. Positive overrides are rejected for CBASS. These factors keep
existing internal unit/component paths consistent; they do not select CMB or SZ
as fitted components. Other required parameters of a normal run are unchanged.

## Input file and provenance

`examples/bandpasses/cbass_I_20130625.dat` contains all 201 samples from 4 to
6 GHz, spaced by 0.01 GHz, derived from the user-supplied
`20130625_v02_Passband[49].csv`. The CSV header orders the gains as
I1, Q1, U1, Q2, U2, I2, followed by six uncertainty columns. The output header
records the SHA256 of the original input.

Roke Cepeda-Arroita's correction in the 11 April 2025 email confirms that the
CSV is already linear, despite the earlier email saying dB. Following the
supplied plotting recipe, the converter computes `(abs(I1)+abs(I2))/2`, then
normalizes this average to its maximum. It does not apply a power of frequency
or a calibration normalization. `10*log10(G)` is only for plotting.

The conversion is reproducible without third-party Python packages:

```sh
python3 tools/prepare_cbass_intensity_bandpass.py input.csv output.dat
```

The magnitudes reproduce that specific supplied recipe, including small signed
fluctuations in low-response regions. Uncertainty columns are checked for finite
values but are not propagated into an instrumental uncertainty model. This is
not a general prescription to rectify arbitrary signed bandpass measurements.
The original input CSV is not modified. Commander reads two whitespace-separated
columns: frequency in GHz and nonnegative linear response. Overall response
scale cancels in the normalization. No tails are thresholded for CBASS.

## Response and normalization

Let `nu0 = FREQ_Cxx`, `a = BANDPASS_CAL_INDEXxx`, and `G(nu)` be the input
response. The band predicts the calibrated RJ map signal

```text
D      = integral G(nu) * (nu/nu0)^a dnu
w(nu)  = G(nu) * (nu/nu0)^2 / D
T_map  = integral w(nu) * T_RJ(nu) dnu
```

All integrals use the existing trapezoidal quadrature, with frequencies in Hz
internally. Each continuum component uses the existing band-integrated SED
machinery, so curved spectra and component mixtures need no effective sky
spectral index. The integral is linear in component amplitudes.

For `T_RJ(nu) = T0*(nu/nu0)^beta`, the response is `T0*K(alpha)` with
`alpha = beta+2` and

```text
K(alpha) = integral G(nu)*(nu/nu0)^alpha dnu / D
```

Thus the calibration spectrum (`beta=a-2`) gives exactly `T0`. For the supplied
bandpass at 4.76 GHz, a=-0.299 and alpha=-0.7, K=1.00065665502205. The thesis
Table 2.2 intensity polynomial gives 1.0006863, differing by about 0.003%.
K multiplies a monochromatic model to predict the observed map; 1/K multiplies
the observed data to correct them. This follows the thesis table and the
multiplicative correction implemented by the public fastcc C-BASS calculation;
the ratio printed as C in thesis equation 2.4 appears inverted relative to its
prose and table. The supplied raw-stream average also differs slightly from
fastcc's choice to normalize each stream before averaging.

The new branch does not borrow the LFI/HFI CMB or IRAS calibration normalizations.
It recalculates its own weights and unit factors whenever `update_tau` changes
the bandpass. Invalid frequencies, signed/nonfinite responses, zero-area bands,
and invalid calibration integrals fail explicitly. For initial comparisons use
`BP_RMSxx=0` and no nonzero initial bandpass shift.

## Choose integration OR the polynomial, per band

`CBASS` requires `USE_COLOR_CORRECTIONxx=F` whenever color corrections are enabled
globally. Combining them is rejected before polynomial coefficients are read.
Other bands may still use polynomial corrections in the same run.

- Full bandpass: `BANDPASS_TYPExx='CBASS'`, polynomial off.
- Monochromatic approximation: `BANDPASS_TYPExx='delta'`, polynomial optionally on.

`COLOR_CORR_UPDATE_INTERVAL` has no effect on CBASS integration. Existing SED
lookup tables are used in the usual manner; no new amplitude-dependent correction
update or sampler approximation is added by this branch.

## Validation

Run `python3 tests/run_bandpass_tests.py` with a Fortran compiler specified by
`FC` if needed, and the existing `python3 tests/run_color_tests.py` regression
suite. The bandpass suite compiles the complete production `comm_bp_mod.f90`
with single-process dependency fixtures and the production trapezoidal routine.
It exercises initialization, real file loading, measured-band numerical
responses, arbitrary reference parameters, mixtures, unit conversions, narrow
bands, updated bandpasses, legacy conventions and rejected configurations.
It is not a full HEALPix/MPI Commander build or a cluster data fit.
