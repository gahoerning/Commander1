#!/usr/bin/env python3
"""Generate complete Commander parameter files from a matched FOSSIL export.

This program uses only the Python standard library. It does not run Commander.
Paths in the export manifest are relative to the export's INPUT_DIRECTORY.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


STAGES = {
    "no_globals_regression": dict(operation="optimize", observed=False, local_amplitudes=False,
                                  local_shapes=False, global_shapes=False, iterations=2),
    "globals_optimize": dict(operation="optimize", observed=False, local_amplitudes=False,
                              local_shapes=False, global_shapes=False, iterations=12),
    "globals_sample": dict(operation="sample", observed=True, local_amplitudes=False,
                            local_shapes=False, global_shapes=False, iterations=2000),
    "amplitudes_optimize": dict(operation="optimize", observed=False, local_amplitudes=True,
                                 local_shapes=False, global_shapes=False, iterations=12),
    "global_shapes_sample": dict(operation="sample", observed=True, local_amplitudes=False,
                                  local_shapes=False, global_shapes=True, iterations=5000),
    "amplitudes_sample": dict(operation="sample", observed=True, local_amplitudes=True,
                               local_shapes=False, global_shapes=False, iterations=5000),
    "forecast_optimize": dict(operation="optimize", observed=False, local_amplitudes=True,
                               local_shapes=True, global_shapes=True, iterations=50),
    "forecast_sample": dict(operation="sample", observed=True, local_amplitudes=True,
                             local_shapes=True, global_shapes=True, iterations=10000),
}

# Bounds and priors are analysis choices, not additional FOSSIL fiducials.
LOCAL_PRIORS = {
    "dust": [("EMISSIVITY", 0.5, 2.8, 1.53, 0.2), ("DUST_TEMP", 8.0, 45.0, 21.7, 3.0)],
    "synch": [("INDEX", -5.0, -1.5, -3.0, 0.3), ("CURVATURE", -0.5, 0.5, 0.0, 0.0)],
    "freefree": [("INDEX", -3.0, -1.5, -2.14, 0.1), ("CURVATURE", -0.5, 0.5, 0.0, 0.0)],
    "ame1": [("FREQUENCY", 10.0, 70.0, 30.0, 5.0)],
    "ame2": [("FREQUENCY", 10.0, 70.0, 33.35, 5.0)],
}


def local_priors(manifest: dict) -> dict:
    """Return priors for the declared injection, preserving the PDF-only profile."""
    priors = dict(LOCAL_PRIORS)
    if manifest.get("freefree_model", "").startswith("draine"):
        priors["freefree"] = [("T_E", 3000.0, 15000.0, 7000.0, 500.0)]
    if manifest.get("primary_cmb_anisotropy", False):
        priors["cmb"] = []
    return priors


def band_suffix(number: int) -> str:
    """Preserve 01..99 legacy spelling; use 100..131 without truncation."""
    return f"{number:02d}"


def encode(value) -> str:
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, str):
        if "'" in value or "\n" in value:
            raise ValueError("Commander parameter strings cannot contain apostrophes or newlines")
        return f"'{value}'"
    if isinstance(value, float):
        return f"{value:.16g}".replace("e", "d")
    return str(value)


def _map_pattern(pattern: str, channel: dict) -> str:
    band = int(channel["band"])
    frequency = float(channel.get("frequency_GHz", channel.get("frequency_Hz", 0) / 1e9))
    return pattern.format(band=band, index=band, band03=f"{band:03d}", frequency_GHz=frequency)


def validate_manifest(manifest: dict, nside: int) -> dict:
    resolution = manifest["resolutions"][str(nside)]
    channels = manifest["channels"]
    if len(channels) != 131:
        raise ValueError(f"Expected 131 FOSSIL channels; found {len(channels)}")
    if [int(c["commander_band"]) for c in channels] != list(range(1, 132)):
        raise ValueError("commander_band must run from 1 to 131 in frequency order")
    metadata = manifest.get("metadata", manifest)
    if abs(float(metadata["beam_fwhm_deg"]) - 1.6) > 1e-12:
        raise ValueError("This configuration requires a 1.6-degree Gaussian FWHM")
    if metadata.get("primary_cmb_anisotropy", False) and not manifest.get("cmb"):
        raise ValueError("CMB anisotropy requires explicit realization/prior metadata")
    if metadata.get("units", metadata.get("unit")) != "MJy/sr":
        raise ValueError("All input intensity and RMS maps must use MJy/sr")
    if [c["name"] for c in resolution["components"]] != list(local_priors(manifest)):
        raise ValueError("Component names/order do not match the declared foreground/CMB model")
    return resolution


def select_channels(manifest: dict, band_count: int = 131) -> list[dict]:
    """Select a deterministic frequency-spanning subset without changing maps.

    Retain endpoints and all amplitude-reference channels. The Commander channel
    IDs are assigned later; original manifest IDs remain the file/statistic keys.
    The 64-band profile is a resource-limited analysis using a documented
    fiducial spectral-information selection. The full profile retains 131 bands.
    """
    if band_count not in (64, 131):
        raise ValueError("--band-count must be 64 or 131")
    channels = manifest["channels"]
    if len(channels) != 131:
        raise ValueError("Band selection requires the complete 131-channel export manifest")
    if band_count == 131:
        return list(channels)
    required = {0, len(channels) - 1}
    for resolution in manifest["resolutions"].values():
        for component in resolution["components"]:
            index = int(component["reference_band"]) - 1
            if index < 0 or index >= len(channels):
                raise ValueError("Local amplitude reference band is outside the export")
            frequency = float(channels[index].get("frequency_GHz", channels[index].get("frequency_Hz", 0) / 1e9))
            if abs(frequency - float(component["reference_frequency_GHz"])) > 1e-9:
                raise ValueError("Local amplitude reference frequency disagrees with its export band")
            required.add(index)
    if len(required) > band_count:
        raise ValueError("Too many reference channels for the requested subset")
    profile_path = Path(__file__).with_name("channel_selection_cmb_draine_64.json" if manifest.get("primary_cmb_anisotropy", False) else "channel_selection_64.json")
    profile = json.loads(profile_path.read_text())
    numbers = profile["original_commander_bands"]
    if len(numbers) != band_count or len(set(numbers)) != band_count or numbers != sorted(numbers):
        raise ValueError("Invalid or unordered 64-channel selection profile")
    if any(not isinstance(number, int) or number < 1 or number > len(channels) for number in numbers):
        raise ValueError("Selected original Commander band lies outside the export")
    indices = {number - 1 for number in numbers}
    if not required.issubset(indices):
        raise ValueError("Selection profile omits an endpoint or local amplitude reference channel")
    selected = [channels[index] for index in sorted(indices)]
    frequencies = [float(c.get("frequency_GHz", c.get("frequency_Hz", 0) / 1e9)) for c in selected]
    expected = profile["frequencies_GHz"]
    if len(expected) != band_count or any(abs(a - float(b)) > 1e-9 for a, b in zip(frequencies, expected)):
        raise ValueError("Selection profile frequencies do not match this export manifest")
    return selected


def preflight_paths(parameters: dict, verify_inputs: bool = False) -> None:
    """Reject path truncation before Fortran I/O; optionally check input files."""
    root = Path(parameters["INPUT_DIRECTORY"])
    for key, value in parameters.items():
        if not isinstance(value, str) or value == "none":
            continue
        if key in ("INPUT_DIRECTORY", "CHAIN_DIRECTORY"):
            resolved = value
        elif value.lower().endswith((".fits", ".txt", ".dat")):
            path = Path(value)
            cwd_keys = {"IN_POWSPECFILE", "BINFILE", "CMB_MONO_DIPOLE_MASK"}
            if not path.is_absolute() and key not in cwd_keys:
                path = root / path
            resolved = str(path)
            if verify_inputs and not path.is_file():
                raise FileNotFoundError(f"{key}: missing input {path}")
        else:
            continue
        if len(resolved.encode("utf-8")) >= 512:
            raise ValueError(f"{key}: path exceeds Commander's 511-byte limit")


def build_parameters(manifest: dict, nside: int, stage: str, input_directory: str,
                     chain_directory: str, seed: int = 20261102, groups: int = 1,
                     proc_per_band: int = 1, iterations: int | None = None,
                     band_count: int = 131) -> dict:
    resolution = validate_manifest(manifest, nside)
    channels = select_channels(manifest, band_count)
    reference_bands = {int(channel["commander_band"]): i for i, channel in enumerate(channels, 1)}
    settings = STAGES[stage]
    lmax = 3 * nside - 1
    cmb_present = bool(manifest.get("primary_cmb_anisotropy", False))
    harmonic_cmb = cmb_present and settings["local_amplitudes"]
    priors = local_priors(manifest)
    mask = resolution["mask"]
    fullsky = resolution["fullsky_mask"]
    p = dict(
        OPERATION=settings["operation"], VERBOSITY=3, POLARIZATION=False,
        SAMPLE_TT_MODES=True, SAMPLE_TE_MODES=False, SAMPLE_TB_MODES=False,
        SAMPLE_EE_MODES=False, SAMPLE_EB_MODES=False, SAMPLE_BB_MODES=False,
        NUM_REALIZATIONS=1, NUM_CHAIN_PER_REALIZATION=groups, NUM_GROUPS=groups,
        NUMBAND=len(channels), NUM_PROC_PER_BAND=proc_per_band, NSIDE=nside, LMAX=lmax,
        LMAX_LOWRES=lmax, FWHM_LOWRES=96.0,
        NUM_GIBBS_ITER=iterations if iterations is not None else settings["iterations"],
        NUM_STEP_WITH_ML_SEARCH=0, NUM_FAST_SUBSTEPS=1,
        THINNING_FACTOR=1, OUTPUT_THINNING_FACTOR=100 if settings["iterations"] >= 2000 else 1,
        SPECIND_SAMPLER_BLOCKSIZE=10000000,
        RMS_SPEC_INDEX_RANDOM_WALK=0.0, NUMGRID_SPLINE=1000,
        NUM_JOINT_ALMS_CL_STEPS=0, JOINT_ALMS_CL_CHISQ_THRESHOLD=20.0,
        BASE_SEED=seed, PRECOMPUTE_PLMS=0, SAMPLING_METHOD="joint_sampler",
        SOLUTION_METHOD="CG", LMAX_PRECOND=min(50, lmax), PRECOND_MATRIX_METHOD="diagonal",
        CONDITION_ON_CL=True, ENFORCE_ZERO_CL=True, IN_POWSPECFILE="none", DISPERSION_FACTOR=1.0,
        CL_BINNING=False, BINFILE="none", SAMPLE_INSIDE_MASK=False,
        CMB_MONO_DIPOLE_MASK="none", OUTPUT_ALMS=False, OUTPUT_MAPS=False,
        OUTPUT_SPECTRA=False, OUTPUT_FOREGROUNDS=True, OUTPUT_CHISQ=True,
        OUTPUT_MIXING_MATRIX=False, OUTPUT_BAND_CHISQ=True, OUTPUT_GAIN_RESIDUALS=False,
        OUTPUT_ML_MAP_AND_COVMAT=False, OUTPUT_CMB_FREQUENCY_MAPS=False,
        OUTPUT_CROSS_CORRELATION_STATS=False, NUM_OBJECTS=0,
        INPUT_DIRECTORY=input_directory, CHAIN_DIRECTORY=chain_directory,
        OVERVIEW_MASK_HIGHFG=mask, OVERVIEW_MASK_LOWFG=mask,
        OVERVIEW_FREQ_MIN=50.0, OVERVIEW_FREQ_MAX=2000.0,
        CGD_MAXITER=10000, CGD_TOLERANCE=1e-9,
        REGULARIZATION_NOISE=0.0, MASKFILE=mask, MASKFILE_CORR=mask,
        PROCESSING_MASK=mask, GAIN_APODIZATION_MASK=fullsky, FWHM_APOD_GAINMASK=0.0,
        CORR_CHISQ_THRESHOLD=30.0, LMAX_CORR=lmax, NSIDE_CORR=min(16, nside),
        TEMPLATE_AMP_INPUT="none", GAIN_INIT="none", BANDPASS_INIT="none",
        GAIN_INIT_RMS=0.0, BP_INIT_RMS=0.0, NOISE_FORMAT="rms",
        BEAM_CONVOLUTION_METHOD="transfer_function", PIXEL_WINDOW=resolution["pixel_window"],
        T_CMB=2.725, MJYSR_CONVENTION="PSM", APPLY_BP_CORRECTIONS=False,
        APPLY_GAIN_CORRECTIONS=False, APPLY_COLOR_CORRECTIONS=False,
        NUM_BP_MCMC_SUBSTEPS=0, SOLAR_AMP=0.0, SOLAR_LON=0.0, SOLAR_LAT=0.0,
        NUM_FG_TEMP=0, NUM_FG_SIGNAL_COMPONENTS=len(resolution["components"]),
        APPLY_PAR_RMS_SMOOTHING=False, OPTIMIZE_PRIORS=False, FOSSIL_EXACT_LOCAL_SEDS=True,
        GLOBAL_FIXED_INSTRUMENT=True, GLOBAL_SLICE_STEPS=8, GLOBAL_RESTART_FILE="none",
    )
    if cmb_present:
        p.update(FOSSIL_CMB_MIN_L=2, FOSSIL_SI_CMB=True,
                 ENFORCE_ZERO_CL=not harmonic_cmb,
                 IN_POWSPECFILE=str(Path(input_directory) / resolution["cmb_prior"]),
                 CMB_MONO_DIPOLE_MASK="none")
    maps = resolution["maps"]
    for number, channel in enumerate(channels, 1):
        suffix = band_suffix(number)
        values = dict(
            FREQ_LABEL=f"FOSSIL_{int(channel['commander_band']):03d}", FREQ_UNIT="MJy/sr",
            **{f"MAP{suffix}_0001": _map_pattern(maps["observed" if settings["observed"] else "signal"], channel)},
            NOISE_RMS=_map_pattern(maps["sigma"], channel), MASKFILE_CALIB=fullsky,
            BEAM=resolution["beam"], BANDPASS_TYPE="delta", BP_RMS=0.0,
            A2T=0.0, F2T=0.0, CO2T=0.0, A2SZ=0.0,
            FREQ_C=float(channel.get("frequency_GHz", channel.get("frequency_Hz", 0) / 1e9)),
            FIX_MONOPOLE=True, FIX_DIPOLE=True, REG_NOISE_SCALE=0.0,
            NOISEAMP_PRIOR_MEAN=1.0, NOISEAMP_PRIOR_RMS=0.0,
            GAIN_PRIOR_MEAN=1.0, GAIN_PRIOR_RMS=0.0, GAIN_CALIB_COMPONENT=1,
            GAIN_LMIN=-1, GAIN_LMAX=-1,
        )
        for key, value in values.items():
            p[key if key.startswith("MAP") else key + suffix] = value
    for number, component in enumerate(resolution["components"], 1):
        suffix = band_suffix(number)
        name = component["name"]
        fixed_truth = not settings["local_amplitudes"] and not settings["local_shapes"]
        # Shape-fixed stages must use truth, even when local amplitudes are free.
        amplitude = component["truth_amplitude"] if fixed_truth else component["amplitude"]
        values = dict(COMP_TYPE=component["type"], FG_LABEL=name, FG_UNIT="MJy/sr",
                      REFERENCE_FREQUENCY=float(component["reference_frequency_GHz"]),
                      REFERENCE_BAND=reference_bands[int(component["reference_band"])],
                      COMPONENT_MASK="fullsky", INDEX_SAMPLING_MASK="fullsky",
                      PRIOR_OPTIMIZATION_MASK="fullsky", INITIALIZATION_MODE="input_map",
                      INITIAL_AMPLITUDE_MAP=amplitude, SAMPLE_AMPLITUDES=settings["local_amplitudes"],
                      ENFORCE_POSITIVE_AMPLITUDE=False, IMPOSE_FG_ORTHOGONALITY=False,
                      APPLY_JEFFREYS_PRIOR=False, OUTPUT_FREQUENCY_COMPONENT_MAPS=False)
        for key, value in values.items():
            p[key + suffix] = value
        if name == "dust":
            p.update({"NU_FLATTENING" + suffix: 100.0, "FRAC_FLATTENING" + suffix: 0.0,
                      "DEFAULT_EMISSIVITY" + suffix: 1.53, "DEFAULT_DUST_TEMP" + suffix: 21.7})
        elif component["type"] == "freefree":
            p["DEFAULT_T_E" + suffix] = 7000.0
        elif name in ("synch", "freefree"):
            p["DEFAULT_SPECTRAL_INDEX" + suffix] = -3.0 if name == "synch" else -2.14
            p["DEFAULT_CURVATURE" + suffix] = 0.0
        elif name.startswith("ame"):
            p["DEFAULT_FREQUENCY" + suffix] = 30.0 if name == "ame1" else 33.35
            p["SPECTRUM_FILENAME" + suffix] = component.get("sed_template", manifest.get("ame_template", "spectra/ame_emissivity_GHz.txt"))
        for index, (parameter, prior) in enumerate(zip(component["parameters"], priors[name]), 1):
            prefix, low, high, mean, sigma = prior
            psuffix = suffix + "_" + band_suffix(index)
            p["INIT_INDEX_MAP" + psuffix] = parameter["map"] if settings["local_shapes"] else parameter["truth_map"]
            p["REGION_DEFINITION" + psuffix] = f"single_pix TQU {nside}"
            p["FWHM_PAR" + psuffix] = 0.0
            p[prefix + "_PRIOR_UNIFORM_LOW" + suffix] = low
            p[prefix + "_PRIOR_UNIFORM_HIGH" + suffix] = high
            p[prefix + "_PRIOR_GAUSSIAN_MEAN" + suffix] = mean
            p[prefix + "_PRIOR_GAUSSIAN_STDDEV" + suffix] = sigma if settings["local_shapes"] else 0.0
            p["OPTIMIZE_PRIOR_MEAN" + psuffix] = False
    add_globals(p, manifest, settings, stage)
    preflight_paths(p)
    return p


def add_globals(p: dict, manifest: dict, settings: dict, stage: str) -> None:
    globals_ = [
        ("delta_t", "delta_t", 1.2e-4, 1e-4, []),
        ("mu", "mu", 2e-8, 1e-8, []),
        ("tsz", "y", 1.77e-6, 1e-6, [(1.245, 0.0, 10.0)]),
        ("cib", "cib", 0.346, 1.0, [(0.86, 0.0, 3.0), (18.8, 5.0, 50.0)]),
        ("crr", "crr", 1.0, 1.0, []),
        ("co", "co", 1.0, 1.0, []),
    ]
    p["NUM_GLOBAL_COMPONENTS"] = 0 if stage == "no_globals_regression" else len(globals_)
    if p["NUM_GLOBAL_COMPONENTS"] == 0:
        return
    templates = manifest.get("global_templates", {})
    for number, (kind, label, truth, scale, parameters) in enumerate(globals_, 1):
        suffix = band_suffix(number)
        for key, value in dict(GLOBAL_TYPE=kind, GLOBAL_LABEL=label, GLOBAL_INIT=0.95 * truth,
                               GLOBAL_SAMPLE_AMPLITUDE=True, GLOBAL_PRIOR_MEAN=truth,
                               GLOBAL_PRIOR_SIGMA=-1.0, GLOBAL_SCALE=scale).items():
            p[key + suffix] = value
        if kind in ("crr", "co"):
            p["GLOBAL_TEMPLATE" + suffix] = templates.get(kind, f"spectra/{kind}_MJy_sr.txt")
        for index, (truth_parameter, low, high) in enumerate(parameters, 1):
            psuffix = suffix + "_" + band_suffix(index)
            for key, value in dict(GLOBAL_INIT_PARAM=truth_parameter,
                                   GLOBAL_SAMPLE_PARAM=settings["global_shapes"],
                                   GLOBAL_PARAM_LOW=low, GLOBAL_PARAM_HIGH=high,
                                   GLOBAL_PARAM_PRIOR_MEAN=truth_parameter,
                                   GLOBAL_PARAM_PRIOR_SIGMA=-1.0).items():
                p[key + psuffix] = value


STAGE_PURPOSE = {
    "globals_optimize": "Recover physical monopole amplitudes with a known Galactic sky.",
    "globals_sample": "Check the global amplitude posterior against its exact Gaussian reference.",
    "globals_sample_split1000": "Run the first half of the global state/RNG restart control.",
    "globals_sample_resume2000": "Resume the global state/RNG control and compare with its uninterrupted chain.",
    "amplitudes_optimize": "Exercise the local CG solve coupled to physical monopole amplitudes.",
    "global_shapes_sample": "Sample global spectral shapes while holding the Galactic sky at truth.",
    "amplitudes_sample": "Sample Galactic and global amplitudes with all spectral shapes fixed.",
    "forecast_optimize": "Optimize the complete matched model, including Galactic spectral parameters.",
    "forecast_sample": "Sample the complete component-separation posterior for the matched forecast.",
    "no_globals_regression": "Exercise the disabled-global path; omitted physical signals leave residuals.",
}


def render(parameters: dict, stage: str, nside: int) -> str:
    """Write a sectioned, commented input file without altering parameter values."""
    p = parameters
    used = set()
    base_stage = "globals_sample" if stage.startswith("globals_sample_") else stage
    fullsky_fit = p["MASKFILE"] == p.get("GAIN_APODIZATION_MASK")
    width = max(38, max(map(len, p)))
    lines = ["# Parameter file for Commander - FOSSIL component separation",
             f"# Stage: {stage}; NSIDE: {nside}; active channels: {p['NUMBAND']}",
             "# Purpose: " + STAGE_PURPOSE[stage],
             "# Data: " + ("signal plus Gaussian noise." if STAGES[base_stage]['observed'] else "noiseless signal."),
             "# Sky-intensity maps and RMS maps: MJy/sr; Galactic coordinates; HEALPix RING.",
             "# Gaussian beam: 1.6 degrees FWHM (96 arcmin).",
             "# The TXT specifies sigma per pixel, unchanged between NSIDE 64 and 128.",
             "# The reference CMB blackbody is subtracted; physical DeltaT/distortions remain.",
             "# Primary CMB anisotropy: " + ("included; harmonic monopole/dipole fixed to zero." if p.get("FOSSIL_SI_CMB", False) else "absent in this PDF-only baseline."),
             "# Instrumental offsets/dipoles = 0, gains = 1, and bandpass shifts = 0.",
             "# Contents: run/MPI; resolution; chain; solver; CMB; outputs; paths/masks;",
             "#           instrument; frequency maps; Galactic components; physical monopoles."]
    if fullsky_fit:
        lines.append("# Likelihood: FULL SKY. Apply GAL070 only to output maps/statistics after fitting.")

    def note(text):
        lines.append("# " + text)

    def section(title, description=None):
        while lines and not lines[-1]:
            lines.pop()
        lines.extend(["", "# " + "-" * 76, "# " + title, "# " + "-" * 76])
        if description:
            for text in description:
                note(text)
        lines.append("")

    def emit(key, comment=None):
        if key not in p:
            return
        if key in used:
            raise ValueError(f"Duplicate rendered parameter: {key}")
        used.add(key)
        line = f"{key:<{width}} = {encode(p[key])}"
        if comment:
            line = f"{line:<{width + 21}} # {comment}"
        if len(line.encode()) >= 2048:
            raise ValueError(f"Rendered line exceeds Commander input buffer: {key}")
        lines.append(line)

    def group(keys, comments=None):
        before = len(used)
        for key in keys:
            emit(key, (comments or {}).get(key))
        if len(used) > before:
            lines.append("")

    section("01. Run mode and logging")
    emit("OPERATION", "optimize = point estimates; sample = posterior draws")
    emit("VERBOSITY", "0 = silent, 1 = chatty, 2 = lots, 3 = debug")

    section("02. MPI layout", ["Required MPI ranks = NUM_GROUPS * NUMBAND * NUM_PROC_PER_BAND.",
                              "Use one OpenMP/BLAS thread per MPI rank."])
    group(["NUM_REALIZATIONS", "NUM_CHAIN_PER_REALIZATION", "NUM_GROUPS", "NUMBAND", "NUM_PROC_PER_BAND"],
          {"NUM_CHAIN_PER_REALIZATION": "Equal to NUM_GROUPS in this setup", "NUMBAND": "Active frequency maps in the likelihood",
           "NUM_PROC_PER_BAND": "Pixel partitions per frequency band"})

    section("03. Map resolution and spectral-parameter resolution")
    group(["NSIDE", "LMAX", "LMAX_LOWRES", "FWHM_LOWRES", "PRECOMPUTE_PLMS"],
          {"NSIDE": "HEALPix resolution of input maps", "LMAX": "3 * NSIDE - 1",
           "FWHM_LOWRES": "arcmin; 96 arcmin = 1.6 degrees FWHM"})

    section("04. Chain length, random seed, and spectral sampling")
    group(["NUM_GIBBS_ITER", "NUM_STEP_WITH_ML_SEARCH", "NUM_FAST_SUBSTEPS", "THINNING_FACTOR",
           "OUTPUT_THINNING_FACTOR", "BASE_SEED"],
          {"NUM_GIBBS_ITER": "Final iteration number, including when resuming",
           "NUM_STEP_WITH_ML_SEARCH": "No initial ML-only segment in sample mode",
           "OUTPUT_THINNING_FACTOR": "Map output cadence; global scalar trace is written every iteration",
           "BASE_SEED": "Independent chain seed; keep unchanged for the restart control"})
    group(["SPECIND_SAMPLER_BLOCKSIZE", "RMS_SPEC_INDEX_RANDOM_WALK", "NUMGRID_SPLINE"],
          {"RMS_SPEC_INDEX_RANDOM_WALK": "Zero random-walk perturbation in this setup"})

    section("05. Amplitude solver and CG convergence")
    group(["SAMPLING_METHOD", "SOLUTION_METHOD", "LMAX_PRECOND", "PRECOND_MATRIX_METHOD", "CGD_MAXITER", "CGD_TOLERANCE"],
          {"CGD_MAXITER": "Maximum iterations in each conjugate-gradient solve",
           "CGD_TOLERANCE": "Fractional CG convergence criterion"})

    section("06. Polarization and CMB harmonic controls", ["Intensity-only forecast; CMB covariance is fixed when its harmonic realization is free."])
    group(["POLARIZATION", "SAMPLE_TT_MODES", "SAMPLE_TE_MODES", "SAMPLE_TB_MODES", "SAMPLE_EE_MODES", "SAMPLE_EB_MODES", "SAMPLE_BB_MODES"])
    group(["FOSSIL_CMB_MIN_L", "FOSSIL_SI_CMB", "CONDITION_ON_CL", "ENFORCE_ZERO_CL", "IN_POWSPECFILE", "DISPERSION_FACTOR", "CL_BINNING", "BINFILE",
           "NUM_JOINT_ALMS_CL_STEPS", "JOINT_ALMS_CL_CHISQ_THRESHOLD", "CMB_MONO_DIPOLE_MASK"],
          {"ENFORCE_ZERO_CL": "True = no harmonic CMB solve; a fixed pixel CMB map may still be present",
           "FOSSIL_CMB_MIN_L": "Lowest allowed harmonic multipole; exclude CMB monopole/dipole",
           "FOSSIL_SI_CMB": "Use the SI CMB derivative, consistent with the physical DeltaT monopole",
           "NUM_JOINT_ALMS_CL_STEPS": "Joint harmonic-spectrum moves disabled"})

    section("07. Output products")
    group(["OUTPUT_ALMS", "OUTPUT_MAPS", "OUTPUT_SPECTRA", "OUTPUT_FOREGROUNDS", "OUTPUT_CHISQ",
           "OUTPUT_MIXING_MATRIX", "OUTPUT_BAND_CHISQ", "OUTPUT_GAIN_RESIDUALS", "OUTPUT_ML_MAP_AND_COVMAT",
           "OUTPUT_CMB_FREQUENCY_MAPS", "OUTPUT_CROSS_CORRELATION_STATS", "NUM_OBJECTS"],
          {"OUTPUT_FOREGROUNDS": "Save component amplitude/spectral outputs",
           "OUTPUT_CHISQ": "Save chi-square diagnostics"})

    section("08. Input/output paths and analysis masks", ["Input map/template paths below are relative to INPUT_DIRECTORY.",
                                                          "Use a fresh CHAIN_DIRECTORY for an independent run."])
    group(["INPUT_DIRECTORY", "CHAIN_DIRECTORY"])
    group(["MASKFILE", "MASKFILE_CORR", "PROCESSING_MASK", "SAMPLE_INSIDE_MASK", "OVERVIEW_MASK_HIGHFG",
           "OVERVIEW_MASK_LOWFG", "OVERVIEW_FREQ_MIN", "OVERVIEW_FREQ_MAX", "GAIN_APODIZATION_MASK", "FWHM_APOD_GAINMASK",
           "CORR_CHISQ_THRESHOLD", "LMAX_CORR", "NSIDE_CORR"],
          {"MASKFILE": "Full-sky likelihood mask" if fullsky_fit else "Binary GAL070 likelihood mask",
           "SAMPLE_INSIDE_MASK": "No additional inside-mask sampling"})

    section("09. Noise, beam, units, and fixed instrument", ["Known per-pixel RMS; no extra regularization noise or fitted noise scaling.",
                                                            "Instrumental offsets are separate from the physical monopoles in section 12.",
                                                            "PIXEL_WINDOW uses the original full-length HEALPix table from INPUT_DIRECTORY.",
                                                            f"Legacy read_pixwin also needs pixel_window_n{nside:04d}.fits in the process working directory.",
                                                            "The supplied launcher provides this working-directory link before MPI."])
    group(["NOISE_FORMAT", "REGULARIZATION_NOISE", "BEAM_CONVOLUTION_METHOD", "PIXEL_WINDOW", "T_CMB", "MJYSR_CONVENTION"],
          {"REGULARIZATION_NOISE": "No extra noise added to the maps",
           "PIXEL_WINDOW": "Native HEALPix temperature/polarization window; ell = 0 through 4 * NSIDE",
           "T_CMB": "K; reference blackbody temperature",
           "MJYSR_CONVENTION": "Matched intensity convention of the input export"})
    group(["APPLY_BP_CORRECTIONS", "APPLY_GAIN_CORRECTIONS", "APPLY_COLOR_CORRECTIONS", "NUM_BP_MCMC_SUBSTEPS",
           "TEMPLATE_AMP_INPUT", "GAIN_INIT", "BANDPASS_INIT", "GAIN_INIT_RMS", "BP_INIT_RMS", "SOLAR_AMP", "SOLAR_LON", "SOLAR_LAT"])

    section("10. Frequency map sets", ["Each block is one active Commander band; the label retains its original FOSSIL ID.",
                                         "RMS maps retain the TXT standard deviation per pixel for every resolution.",
                                         "FIX_MONOPOLE/FIX_DIPOLE fix instrumental offsets, not physical sky monopoles."])
    for band in range(1, p["NUMBAND"] + 1):
        suffix = band_suffix(band)
        note(f"------------ Band {suffix}: {p['FREQ_C' + suffix]:g} GHz | {p['FREQ_LABEL' + suffix]}")
        group([key + suffix for key in ("FREQ_LABEL", "FREQ_UNIT", "FREQ_C")]
              + [f"MAP{suffix}_0001"] + [key + suffix for key in ("NOISE_RMS", "MASKFILE_CALIB")])
        group([key + suffix for key in ("BEAM", "BANDPASS_TYPE", "BP_RMS", "A2T", "F2T", "CO2T", "A2SZ")])
        group([key + suffix for key in ("FIX_MONOPOLE", "FIX_DIPOLE", "REG_NOISE_SCALE", "NOISEAMP_PRIOR_MEAN",
                                             "NOISEAMP_PRIOR_RMS", "GAIN_PRIOR_MEAN", "GAIN_PRIOR_RMS", "GAIN_CALIB_COMPONENT", "GAIN_LMIN", "GAIN_LMAX")])

    section("11. Galactic pixel components", ["Local amplitude maps use MJy/sr at each component's reference frequency.",
                                               "A zero spectral Gaussian STDDEV freezes that parameter at its input map.",
                                               "Positive STDDEV activates spectral fitting/sampling with the listed prior."])
    group(["NUM_FG_TEMP", "NUM_FG_SIGNAL_COMPONENTS", "APPLY_PAR_RMS_SMOOTHING", "OPTIMIZE_PRIORS", "FOSSIL_EXACT_LOCAL_SEDS"],
          {"NUM_FG_TEMP": "No fitted instrumental templates", "FOSSIL_EXACT_LOCAL_SEDS": "Use the exact matched local spectral equations"})
    local_titles = {"dust": "Thermal dust: modified blackbody", "synch": "Synchrotron: curved power law",
                    "freefree": "Free-free: curved power law", "ame1": "AME population 1: shifted spinning-dust template",
                    "ame2": "AME population 2: shifted spinning-dust template", "cmb": "Primary CMB temperature anisotropy"}
    local_shapes = {"dust": [("EMISSIVITY", "DEFAULT_EMISSIVITY", "Dust emissivity index beta"), ("DUST_TEMP", "DEFAULT_DUST_TEMP", "Dust temperature [K]")],
                    "synch": [("INDEX", "DEFAULT_SPECTRAL_INDEX", "RJ spectral index beta"), ("CURVATURE", "DEFAULT_CURVATURE", "Spectral curvature (fixed)")],
                    "freefree": [("INDEX", "DEFAULT_SPECTRAL_INDEX", "RJ spectral index beta"), ("CURVATURE", "DEFAULT_CURVATURE", "Spectral curvature (fixed)")],
                    "ame1": [("FREQUENCY", "DEFAULT_FREQUENCY", "Spinning-dust peak frequency [GHz]")],
                    "ame2": [("FREQUENCY", "DEFAULT_FREQUENCY", "Spinning-dust peak frequency [GHz]")]}
    local_shapes["cmb"] = []
    for component in range(1, p["NUM_FG_SIGNAL_COMPONENTS"] + 1):
        suffix = band_suffix(component)
        name = p["FG_LABEL" + suffix]
        title = local_titles[name]
        shapes = local_shapes[name]
        if p["COMP_TYPE" + suffix] == "freefree":
            title = "Draine free-free: amplitude plus electron temperature"
            shapes = [("T_E", "DEFAULT_T_E", "Electron temperature Te [K]")]
        note(f"------------ Sky component {suffix}: {title}")
        group([key + suffix for key in ("COMP_TYPE", "FG_LABEL", "FG_UNIT", "REFERENCE_FREQUENCY", "REFERENCE_BAND",
                                             "INITIALIZATION_MODE", "INITIAL_AMPLITUDE_MAP", "SAMPLE_AMPLITUDES", "ENFORCE_POSITIVE_AMPLITUDE")],
              {"REFERENCE_FREQUENCY" + suffix: "GHz; normalization of the amplitude map",
               "REFERENCE_BAND" + suffix: "Active band index at the reference frequency",
               "SAMPLE_AMPLITUDES" + suffix: "False = hold input amplitude map fixed"})
        group([key + suffix for key in ("COMPONENT_MASK", "INDEX_SAMPLING_MASK", "PRIOR_OPTIMIZATION_MASK",
                                             "IMPOSE_FG_ORTHOGONALITY", "APPLY_JEFFREYS_PRIOR", "OUTPUT_FREQUENCY_COMPONENT_MAPS")])
        emit("SPECTRUM_FILENAME" + suffix, "Two-column spectral template; frequency in GHz")
        if "SPECTRUM_FILENAME" + suffix in p:
            lines.append("")
        for index, (prefix, default, title) in enumerate(shapes, 1):
            psuffix = suffix + "_" + band_suffix(index)
            note(f"Spectral parameter {index:02d}: {title}")
            group([default + suffix, "INIT_INDEX_MAP" + psuffix, "REGION_DEFINITION" + psuffix,
                   "FWHM_PAR" + psuffix, "OPTIMIZE_PRIOR_MEAN" + psuffix])
            group([prefix + ending + suffix for ending in ("_PRIOR_UNIFORM_LOW", "_PRIOR_UNIFORM_HIGH", "_PRIOR_GAUSSIAN_MEAN", "_PRIOR_GAUSSIAN_STDDEV")],
                  {prefix + "_PRIOR_GAUSSIAN_STDDEV" + suffix: "0 = fixed; >0 = Gaussian prior width"})
        group(["NU_FLATTENING" + suffix, "FRAC_FLATTENING" + suffix])

    section("12. Spatially constant physical sky components", ["One amplitude per physical component, shared by every sky pixel.",
                                                               "These are sky signals with frequency dependence, not per-band offsets.",
                                                               "GLOBAL_PRIOR_SIGMA <= 0 disables the Gaussian amplitude prior.",
                                                               "GLOBAL_SCALE rescales the numerical solve; it does not change the physics."])
    group(["NUM_GLOBAL_COMPONENTS", "GLOBAL_FIXED_INSTRUMENT", "GLOBAL_SLICE_STEPS", "GLOBAL_RESTART_FILE"],
          {"GLOBAL_FIXED_INSTRUMENT": "Keep instrumental offsets/gains/bandpass shifts fixed",
           "GLOBAL_RESTART_FILE": "none = fresh global state; otherwise restore global state/RNG"})
    global_titles = {"delta_t": "CMB temperature-offset monopole; amplitude = DeltaT/T0",
                     "mu": "mu spectral distortion; dimensionless amplitude",
                     "tsz": "Thermal SZ distortion; amplitude = Compton y",
                     "cib": "Cosmic infrared background; amplitude coefficient in MJy/sr",
                     "crr": "Cosmological recombination radiation; multiplier of fiducial template",
                     "co": "Aggregate extragalactic CO; multiplier of fiducial template"}
    for component in range(1, p["NUM_GLOBAL_COMPONENTS"] + 1):
        suffix = band_suffix(component)
        kind = p["GLOBAL_TYPE" + suffix]
        note(f"------------ Global component {suffix}: {global_titles[kind]}")
        group([key + suffix for key in ("GLOBAL_TYPE", "GLOBAL_LABEL", "GLOBAL_INIT", "GLOBAL_SAMPLE_AMPLITUDE",
                                             "GLOBAL_PRIOR_MEAN", "GLOBAL_PRIOR_SIGMA", "GLOBAL_SCALE", "GLOBAL_TEMPLATE")],
              {"GLOBAL_SAMPLE_AMPLITUDE" + suffix: "Enable amplitude updates in the selected OPERATION",
               "GLOBAL_TEMPLATE" + suffix: "Frequency [GHz], intensity [MJy/sr]"})
        shapes = [(1, "Electron temperature kTe [keV]")] if kind == "tsz" else [(1, "CIB emissivity index beta"), (2, "CIB temperature [K]")] if kind == "cib" else []
        for index, title in shapes:
            psuffix = suffix + "_" + band_suffix(index)
            note(f"Global spectral parameter {index:02d}: {title}")
            group([key + psuffix for key in ("GLOBAL_INIT_PARAM", "GLOBAL_SAMPLE_PARAM", "GLOBAL_PARAM_LOW", "GLOBAL_PARAM_HIGH",
                                                  "GLOBAL_PARAM_PRIOR_MEAN", "GLOBAL_PARAM_PRIOR_SIGMA")],
                  {"GLOBAL_SAMPLE_PARAM" + psuffix: "False = hold spectral parameter fixed",
                   "GLOBAL_PARAM_PRIOR_SIGMA" + psuffix: "<=0 = no Gaussian prior; finite bounds still apply"})
    missing = set(p) - used
    if missing:
        raise ValueError(f"Parameters missing from the readable layout: {sorted(missing)}")
    return "\n".join(lines).rstrip() + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--input-directory", default="outputs/commander")
    parser.add_argument("--chain-root", default="chains/fossil")
    parser.add_argument("--nside", type=int, nargs="+", default=[128, 64], choices=[128, 64])
    parser.add_argument("--stage", nargs="+", default=list(STAGES), choices=list(STAGES))
    parser.add_argument("--seed", type=int, default=20261102)
    parser.add_argument("--groups", type=int, default=1)
    parser.add_argument("--proc-per-band", type=int, default=1)
    parser.add_argument("--iterations", type=int)
    parser.add_argument("--band-count", type=int, choices=(64, 131), default=131,
                        help="Use all 131 channels or the fiducial-information 64-channel analysis profile")
    parser.add_argument("--verify-inputs", action="store_true",
                        help="Check referenced files under INPUT_DIRECTORY before writing configs")
    args = parser.parse_args()
    if args.groups < 1 or args.proc_per_band < 1:
        parser.error("--groups and --proc-per-band must be positive")
    if args.iterations is not None and args.iterations < 1:
        parser.error("--iterations must be positive")
    manifest = json.loads(args.manifest.read_text())
    args.output.mkdir(parents=True, exist_ok=True)
    index = []
    for nside in args.nside:
        for stage in args.stage:
            tag = "" if args.band_count == 131 else f"_{args.band_count}bands"
            chain = f"{args.chain_root}/nside{nside}/{stage}{tag}_seed{args.seed}"
            parameters = build_parameters(manifest, nside, stage, args.input_directory,
                                          chain, args.seed, args.groups, args.proc_per_band, args.iterations,
                                          band_count=args.band_count)
            preflight_paths(parameters, args.verify_inputs)
            filename = f"fossil_nside{nside}_{stage}{tag}.txt"
            provenance = f"# Export manifest signature: {manifest.get('signature', 'unspecified')}\n"
            if manifest.get("analysis_signature"):
                provenance += f"# Analysis signature: {manifest['analysis_signature']}\n"
            (args.output / filename).write_text(provenance + render(parameters, stage, nside))
            record = dict(file=filename, stage=stage, nside=nside,
                              operation=parameters["OPERATION"], chain_directory=chain,
                              input_directory=args.input_directory,
                              manifest_signature=manifest.get("signature"),
                              mpi_ranks=args.groups * parameters["NUMBAND"] * args.proc_per_band,
                              band_count=parameters["NUMBAND"],
                              frequencies_GHz=[parameters["FREQ_C" + band_suffix(i)]
                                               for i in range(1, parameters["NUMBAND"] + 1)],
                              iterations=parameters["NUM_GIBBS_ITER"])
            if manifest.get("analysis_signature"):
                record.update(analysis_signature=manifest["analysis_signature"],
                              analysis_mask=manifest.get("analysis_mask"),
                              postprocess_mask=manifest["resolutions"][str(nside)].get("postprocess_mask"))
            index.append(record)
    (args.output / "configurations.json").write_text(json.dumps(index, indent=2) + "\n")
    print(f"Generated {len(index)} configurations in {args.output}")


if __name__ == "__main__":
    main()
