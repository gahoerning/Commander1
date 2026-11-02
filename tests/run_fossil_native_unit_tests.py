#!/usr/bin/env python3
"""Compile the unchanged foreground code and the native-unit global adapter.

Global tests compile the complete production SED/global/RNG modules. Small
parameter, band-metadata and single-process MPI fixtures do not constitute a
full Commander/HEALPix/MPI run. Native likelihood conversions preserve physical
global parameters and the original MJy/sr snapshot-file convention.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import tempfile

import numpy as np

import run_fossil_global_tests as global_tests

ROOT = Path(__file__).resolve().parents[1]
FOREGROUND_HEAD_SHA256 = "8060051539020fbd88f97d7f9a85ba3041e13683401ebe3f6d04d9a9e8dfb83a"


def conversion(frequency_GHz):
    """Independent SI intensity-to-thermodynamic-microkelvin calibration."""
    nu = np.asarray(frequency_GHz, dtype=float) * 1e9
    h, k, c, temperature = 6.62607015e-34, 1.380649e-23, 299792458., 2.725
    x = h * nu / (k * temperature)
    derivative = 2 * h**2 * nu**4 / (c**2 * k * temperature**2) * np.exp(x) / np.expm1(x)**2
    return 1e-14 / derivative


class NativeRunner:
    """Transform physical likelihood sufficient statistics, never parameters."""
    def __init__(self, physical_runner):
        self.physical_runner = physical_runner
        self.work = physical_runner.work
        self.templates = physical_runner.templates

    def __call__(self, config, rows, *args, **kwargs):
        native = np.array(rows, dtype=float, copy=True)
        factors = conversion(native[:, 0])
        native[:, 1] /= factors**2
        native[:, 2] /= factors
        config = dict(config)
        config.update({f"FREQ_UNIT{band:02d}": "'uK_cmb'"
                       for band in range(1, len(native) + 1)})
        return self.physical_runner(config, native, *args, **kwargs)


def full_model_config(manifest, frequencies, templates):
    truth = manifest["global_parameters"]
    values = np.array([truth[name] for name in
                       ("delta_T", "mu", "y", "A_cib_MJy_sr", "A_crr", "A_co")])
    scales = np.array([1e-4, 1e-8, 1e-6, 1., 1., 1.])
    config = global_tests.common_config(len(frequencies), "optimize")
    config["NUM_GLOBAL_COMPONENTS"] = "6"
    for i, kind in enumerate(("delta_t", "mu", "tsz", "cib", "crr", "co"), 1):
        shape = ([truth["kTe_keV"]] if kind == "tsz" else
                 [truth["beta_cib"], truth["T_cib_K"]] if kind == "cib" else [])
        template = templates / (kind + "_MJy_sr.txt") if kind in ("crr", "co") else None
        global_tests.add_component(config, i, kind, .95 * values[i - 1],
                                   scale=scales[i - 1], shape=shape, template=template)
    return config, values, scales


def full_model_invariance(physical, native, package, count):
    manifest = json.loads((package / "manifest.json").read_text())
    source = np.genfromtxt(package / "global_response.csv", delimiter=",", names=True)
    references = json.loads((package / "global_reference.json").read_text())["64"]["bands"]
    if count == 64:
        profile = json.loads((ROOT / "examples/fossil/channel_selection_cmb_draine_64.json").read_text())
        ids = set(profile["original_commander_bands"])
        channels = [row for row in manifest["channels"] if row["commander_band"] in ids]
    else:
        channels = list(manifest["channels"])
    channels.sort(key=lambda row: row["frequency_GHz"])
    if len(channels) != count:
        raise AssertionError(f"Expected {count} actual source channels, found {len(channels)}")
    frequencies = np.array([row["frequency_GHz"] for row in channels])
    indices = [np.flatnonzero(np.isclose(source["frequency_GHz"], f, rtol=0, atol=1e-8))[0]
               for f in frequencies]
    design = np.column_stack([source[name][indices] for name in ("delta_t", "mu", "y", "cib", "crr", "co")])
    weights = np.array([next(row["weight"] for row in references
                             if np.isclose(row["frequency_GHz"], f, rtol=0, atol=1e-8))
                        for f in frequencies])
    config, truth, scales = full_model_config(manifest, frequencies, physical.templates)
    rows = np.column_stack([frequencies, weights, weights * (design @ truth)])
    physical_dir = physical(config, rows, "chain", 12)
    native_dir = native(config, rows, "chain", 12)
    assert "FOSSIL_NATIVE_UNIT_INTERFACE_V1" in (native_dir / "driver_stdout.txt").read_text()
    assert "FOSSIL_NATIVE_UNIT_INTERFACE_V1" not in (physical_dir / "driver_stdout.txt").read_text()
    physical_signal = np.loadtxt(physical_dir / "draws.txt")
    native_signal = np.loadtxt(native_dir / "draws.txt") / conversion(frequencies)
    np.testing.assert_allclose(native_signal, physical_signal, rtol=2e-10, atol=2e-13)
    np.testing.assert_allclose(native_signal, np.broadcast_to(design @ truth, native_signal.shape),
                               rtol=2e-10, atol=2e-13)
    states = []
    for directory in (physical_dir, native_dir):
        trace = np.atleast_2d(np.loadtxt(directory / "global_components_c0001.txt"))[-1]
        amplitudes = trace[[1, 2, 3, 5, 8, 9]]
        np.testing.assert_allclose((amplitudes - truth) / scales, 0., rtol=0, atol=1e-5)
        states.append(amplitudes)
        spectra = np.loadtxt(directory / "global_spectra_c0001.txt")
        assert spectra.shape == (count, 8), spectra.shape
        np.testing.assert_allclose(spectra[:, 0], frequencies, rtol=0, atol=0)
        np.testing.assert_allclose(spectra[:, 1], design @ amplitudes, rtol=2e-10, atol=2e-13)
        np.testing.assert_allclose(spectra[:, 2:], design * amplitudes, rtol=2e-10, atol=2e-13)
        np.testing.assert_allclose(spectra[:, 1], spectra[:, 2:].sum(axis=1), rtol=2e-13, atol=1e-14)
        lines = [line for line in (directory / "global_spectra_c0001.txt").read_text().splitlines()
                 if line and not line.startswith("#")]
        assert len(lines) == count and all(len(line.split()) == 8 for line in lines)
    np.testing.assert_allclose((states[1] - states[0]) / scales, 0., rtol=0, atol=2e-6)
    return dict(channels=count, iterations=12,
                max_native_signal_error_MJy_sr=float(np.max(abs(native_signal - design @ truth))),
                max_physical_amplitude_difference_in_scale_units=float(np.max(abs((states[1]-states[0])/scales))),
                snapshots_unit="MJy/sr", likelihood_unit="uK_cmb", native_unit_factor_max=float(conversion(frequencies).max()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("package", type=Path, help="Existing physical MJy/sr package with response/reference tables")
    parser.add_argument("--summary", type=Path, help="Write a machine-readable regression summary")
    args = parser.parse_args()
    package = args.package.resolve()
    foreground = ROOT / "src/commander/comm_fg_component_mod.f90"
    digest = hashlib.sha256(foreground.read_bytes()).hexdigest()
    assert digest == FOREGROUND_HEAD_SHA256, "Foreground source differs from the required unchanged HEAD8ffea73"
    summary = dict(status="passed", scope="isolated production procedures; not full MPI/HEALPix acceptance",
                   foreground_sha256=digest, foreground_matches_HEAD8ffea73=True,
                   source_sha256={name: hashlib.sha256((ROOT / "src/commander" / name).read_bytes()).hexdigest()
                                  for name in ("comm_fossil_global_mod.f90", "comm_bp_mod.f90")}, checks=[])
    with tempfile.TemporaryDirectory(prefix="fossil-native-global-") as temporary:
        work = Path(temporary)
        executable = global_tests.compile_driver(work)
        physical = global_tests.Runner(executable, work, package / "spectra")
        native = NativeRunner(physical)
        reference = np.loadtxt(package / "spectra/fiducial_spectra_MJy_sr.txt")[[0, 4, 10, 20, 40, 70, 97, 130]]
        global_tests.amplitude_cases(native, reference, quick=True)
        summary["checks"].append("Native global Gaussian mean/covariance, priors, frozen amplitude and arbitrary scales")
        global_tests.shape_cases(native, reference, quick=True)
        summary["checks"].append("Native global spectral slice distribution, optimization, fixed shape and boundary optimum")
        global_tests.restart_cases(native, reference)
        summary["checks"].append("Native exact global/RNG checkpoint continuation")
        summary["full_model"] = [full_model_invariance(physical, native, package, count) for count in (64, 131)]
        summary["checks"].append("Physical/native six-component optimize and MJy/sr snapshots at 64 and 131 actual channels")
    if args.summary is not None:
        args.summary.parent.mkdir(parents=True, exist_ok=True)
        args.summary.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
