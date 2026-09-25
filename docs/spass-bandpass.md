# Approximate S-PASS DR1 Stokes I bandpass

`BANDPASS_TYPExx = 'SPASS_RJ'` assumes that the delivered map is a weighted
average of Rayleigh-Jeans temperature. This is an explicit approximation to
the effective calibrated response, not a measured instrument bandpass or a
claim that the DR1 channel weighting and temperature conversion are fully known.
Only Stokes I (`POLARIZATION = F`) with `FREQ_UNITxx = 'uK_ant'` is supported.
Here `uK_ant` is Commander's name for RJ microkelvin.

## Configuration

See [the parameter fragment](../examples/bandpass-spass-I.par). Replace `01`
with the S-PASS band number and resolve the bandpass path using the run's usual
base path. The standard band/map/likelihood parameters are still required.
`FREQ_Cxx` is supplied in GHz; use `2.303` for the nominal DR1 reference.
Set `A2Txx`, `F2Txx` and `A2SZxx` to finite nonpositive values to enable automatic
conversions; positive overrides are rejected. No `BANDPASS_CAL_INDEXxx` is used.
The input map is not rescaled: convert a Kelvin map to microkelvin before use.

Set `USE_COLOR_CORRECTIONxx = F` when polynomial corrections are enabled
globally. Applying a polynomial as well as this integration would count the
bandwidth twice. Other bands can still use polynomial corrections. The color
correction update interval has no effect on this integration, and no sampler
changes are needed. Start with `BP_RMSxx = 0` if bandpass corrections are enabled.

## Response file and normalization

The supplied [response table](../examples/bandpasses/spass_dr1_two_windows.dat)
has two columns: frequency in GHz and nonnegative, dimensionless linear weight.
Comment lines start with `#`. There is no response threshold or tail trimming.

The DR1 windows are 2.176--2.216 GHz and 2.272--2.400 GHz, with a combined width
of 168 MHz. The weights have equal height in both windows and are zero in the
56 MHz gap. Their integrated contributions are 40/168 and 128/168, not one half
each. The uniform-weight mean frequency is 2.3026666667 GHz.

These intervals follow Carretti et al. (2019), MNRAS 489, 2330, Table 1 and
Section 2 ([paper](https://doi.org/10.1093/mnras/stz806)), and the
[official DR1 description](https://sites.google.com/inaf.it/spass).
Earlier products using 184 MHz and a quoted 2307 MHz are not represented here.
The published channel-by-channel calibration motivates a flat approximation;
it does not establish the exact map weights or prove a uniform RJ average.

The file uses a 0.25 MHz integration grid inside each window and 0.1 Hz ramps
at the edges. These are numerical sampling choices, not measured filter edges
or instrumental channel widths. Explicit zero-weight points prevent linear
interpolation from bridging the gap. Duplicate frequencies must not be used
to encode steps, because the shared reader removes them.

For an RJ spectrum `T_RJ(nu)`, Commander evaluates

```
w(nu) = W(nu) / integral[W(nu) dnu]
T_map = integral[w(nu) * T_RJ(nu) dnu]
```

All continuum components are integrated through the same response and then
combined. A constant RJ spectrum is preserved. There is no additional `nu^2`
weight and no single-spectrum color correction. An overall scale on `W` cancels.
Frequencies must be finite, positive and strictly increasing after the shared
reader's duplicate removal; responses must be finite and nonnegative, with
positive integrated support and at least two samples.

The existing bandpass update machinery can shift or tilt the shape; weights and
unit conversion factors are recomputed together. The factors keep existing
special-component and reference-unit paths consistent with the integrated RJ
response. They do not add CMB/SZ components or convert the input map to CMB units.
Changing `FREQ_Cxx` alone does not change the average of a fixed physical RJ
spectrum, although it changes reference-frequency brightness conversions.

## Verification and use

Run `python3 tests/run_bandpass_tests.py` with a supported Fortran compiler
(`FC` and `FFLAGS` can be set). The tests compile the complete production bandpass
module with lightweight external-dependency fixtures and production quadrature.
They check analytic two-window power laws and mixtures, a constant sky, the gap,
reference units, scale invariance, repeated updates, shifts/tilts, the narrow-band
limit and invalid inputs, alongside the existing C-BASS and legacy checks.
This does not replace a full cluster build or an end-to-end Commander run.

Compare a run using this approximation with a delta at the same reference
frequency. For example, a uniform RJ average of `(nu/2.303 GHz)^(-3)` is
1.0058340452, relative to unity for that delta. This is a prediction of the
assumed response, not an independently measured S-PASS color correction.
