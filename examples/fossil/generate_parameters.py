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
    "global_shapes_sample": dict(operation="sample", observed=True, local_amplitudes=False,
                                  local_shapes=False, global_shapes=True, iterations=5000),
    "amplitudes_sample": dict(operation="sample", observed=True, local_amplitudes=True,
                               local_shapes=False, global_shapes=False, iterations=5000),
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
    if metadata.get("primary_cmb_anisotropy", False):
        raise ValueError("These baseline configs exclude primary CMB anisotropy; use an explicit extension")
    if metadata.get("units", metadata.get("unit")) != "MJy/sr":
        raise ValueError("All input intensity and RMS maps must use MJy/sr")
    if [c["name"] for c in resolution["components"]] != list(LOCAL_PRIORS):
        raise ValueError("Expected local components dust, synch, freefree, ame1, ame2")
    return resolution


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
            path = path if path.is_absolute() else root / path
            resolved = str(path)
            if verify_inputs and not path.is_file():
                raise FileNotFoundError(f"{key}: missing input {path}")
        else:
            continue
        if len(resolved.encode("utf-8")) >= 512:
            raise ValueError(f"{key}: path exceeds Commander's 511-byte limit")


def build_parameters(manifest: dict, nside: int, stage: str, input_directory: str,
                     chain_directory: str, seed: int = 20261102, groups: int = 1,
                     proc_per_band: int = 1, iterations: int | None = None) -> dict:
    resolution = validate_manifest(manifest, nside)
    settings = STAGES[stage]
    lmax = 3 * nside - 1
    mask = resolution["mask"]
    fullsky = resolution["fullsky_mask"]
    p = dict(
        OPERATION=settings["operation"], VERBOSITY=2, POLARIZATION=False,
        SAMPLE_TT_MODES=True, SAMPLE_TE_MODES=False, SAMPLE_TB_MODES=False,
        SAMPLE_EE_MODES=False, SAMPLE_EB_MODES=False, SAMPLE_BB_MODES=False,
        NUM_REALIZATIONS=1, NUM_CHAIN_PER_REALIZATION=groups, NUM_GROUPS=groups,
        NUMBAND=131, NUM_PROC_PER_BAND=proc_per_band, NSIDE=nside, LMAX=lmax,
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
        NUM_FG_TEMP=0, NUM_FG_SIGNAL_COMPONENTS=5,
        APPLY_PAR_RMS_SMOOTHING=False, OPTIMIZE_PRIORS=False, FOSSIL_EXACT_LOCAL_SEDS=True,
        GLOBAL_FIXED_INSTRUMENT=True, GLOBAL_SLICE_STEPS=8, GLOBAL_RESTART_FILE="none",
    )
    maps = resolution["maps"]
    for channel in manifest["channels"]:
        number, band = int(channel["commander_band"]), int(channel["band"])
        suffix = band_suffix(number)
        values = dict(
            FREQ_LABEL=f"FOSSIL_{number:03d}", FREQ_UNIT="MJy/sr",
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
                      REFERENCE_BAND=int(component["reference_band"]),
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
        elif name in ("synch", "freefree"):
            p["DEFAULT_SPECTRAL_INDEX" + suffix] = -3.0 if name == "synch" else -2.14
            p["DEFAULT_CURVATURE" + suffix] = 0.0
        else:
            p["DEFAULT_FREQUENCY" + suffix] = 30.0 if name == "ame1" else 33.35
            p["SPECTRUM_FILENAME" + suffix] = component.get("sed_template", manifest.get("ame_template", "spectra/ame_emissivity_GHz.txt"))
        for index, (parameter, prior) in enumerate(zip(component["parameters"], LOCAL_PRIORS[name]), 1):
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


def render(parameters: dict, stage: str, nside: int) -> str:
    lines = ["# Generated FOSSIL Commander configuration; all maps and RMS inputs are MJy/sr.",
             "# 1.6-degree Gaussian FWHM; Galactic RING maps; fixed offsets=0, gains=1, shifts=0.",
             "# The CMB blackbody reference has been subtracted; physical DeltaT/distortions remain.",
             "# This baseline contains no primary CMB anisotropy; see docs/fossil.md.",
             f"# Stage: {stage}; NSIDE={nside}; do not reinterpret optimize output as samples."]
    lines.extend(f"{key:<42} = {encode(value)}" for key, value in parameters.items())
    return "\n".join(lines) + "\n"


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
            chain = f"{args.chain_root}/nside{nside}/{stage}_seed{args.seed}"
            parameters = build_parameters(manifest, nside, stage, args.input_directory,
                                          chain, args.seed, args.groups, args.proc_per_band, args.iterations)
            preflight_paths(parameters, args.verify_inputs)
            filename = f"fossil_nside{nside}_{stage}.par"
            provenance = f"# Export manifest signature: {manifest.get('signature', 'unspecified')}\n"
            (args.output / filename).write_text(provenance + render(parameters, stage, nside))
            index.append(dict(file=filename, stage=stage, nside=nside,
                              operation=parameters["OPERATION"], chain_directory=chain,
                              input_directory=args.input_directory,
                              manifest_signature=manifest.get("signature"),
                              mpi_ranks=args.groups * 131 * args.proc_per_band,
                              iterations=parameters["NUM_GIBBS_ITER"]))
    (args.output / "configurations.json").write_text(json.dumps(index, indent=2) + "\n")
    print(f"Generated {len(index)} configurations in {args.output}")


if __name__ == "__main__":
    main()
