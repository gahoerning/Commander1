#!/usr/bin/env python3
"""Validate generated FOSSIL configs without compiling or running Commander."""
import importlib.util
from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("fossil_parameters", ROOT / "examples/fossil/generate_parameters.py")
generator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(generator)


def fixture_manifest():
    components = []
    for name in generator.LOCAL_PRIORS:
        npar = 2 if name in ("dust", "synch", "freefree") else 1
        components.append(dict(name=name, type="one-component_dust" if name == "dust" else
                               "AME_freq_shift" if name.startswith("ame") else "power_law",
                               reference_frequency_GHz=545.0 if name == "dust" else 50.0,
                               reference_band=34 if name == "dust" else 1,
                               amplitude=f"initial/{name}.fits", truth_amplitude=f"truth/{name}.fits",
                               parameters=[dict(map=f"initial/{name}_{p}.fits",
                                                truth_map=f"truth/{name}_{p}.fits") for p in range(npar)]))
    resolution = dict(mask="mask70.fits", fullsky_mask="fullsky.fits", beam="beam.fits",
                      pixel_window="pixel.fits", components=components,
                      maps={kind: f"maps/{kind}_{{band:03d}}.fits" for kind in
                            ("signal", "observed", "sigma", "noise")})
    return dict(metadata=dict(beam_fwhm_deg=1.6, primary_cmb_anisotropy=False, units="MJy/sr"),
                channels=[dict(band=i, commander_band=i + 1, frequency_GHz=50 + i * 15,
                               sigma_MJy_sr=1e-6) for i in range(131)],
                resolutions={"128": resolution, "64": resolution})


class ParameterTests(unittest.TestCase):
    def parameters(self, stage="forecast_sample", nside=128):
        return generator.build_parameters(fixture_manifest(), nside, stage, "input", "chains/test")

    def test_all_stages_use_131_distinct_bands_and_fixed_instruments(self):
        for nside in (128, 64):
            for stage in generator.STAGES:
                p = self.parameters(stage, nside)
                self.assertEqual(p["NUMBAND"], 131)
                self.assertEqual(len([k for k in p if re.fullmatch(r"MAP\d+_0001", k)]), 131)
                self.assertEqual(p["FREQ_C99"], 1520)
                self.assertEqual(p["FREQ_C100"], 1535)
                self.assertEqual(p["FREQ_C131"], 2000)
                self.assertNotIn("FREQ_C*0", p)
                for band in range(1, 132):
                    b = generator.band_suffix(band)
                    self.assertEqual(p["FREQ_UNIT" + b], "MJy/sr")
                    self.assertTrue(p["FIX_MONOPOLE" + b])
                    self.assertTrue(p["FIX_DIPOLE" + b])
                    self.assertEqual(p["GAIN_PRIOR_RMS" + b], 0)
                    self.assertEqual(p["BP_RMS" + b], 0)
                    self.assertEqual(p["NOISEAMP_PRIOR_RMS" + b], 0)
                self.assertEqual(p["REGULARIZATION_NOISE"], 0)
                self.assertFalse(p["APPLY_COLOR_CORRECTIONS"])
                self.assertFalse(p["APPLY_GAIN_CORRECTIONS"])
                self.assertFalse(p["APPLY_BP_CORRECTIONS"])
                self.assertTrue(p["FOSSIL_EXACT_LOCAL_SEDS"])
                self.assertEqual(p["FWHM_LOWRES"], 96)

    def test_global_only_stages_really_freeze_truth_local_sky(self):
        for stage in ("globals_optimize", "globals_sample", "global_shapes_sample"):
            p = self.parameters(stage)
            for i in range(1, 6):
                suffix = generator.band_suffix(i)
                self.assertFalse(p["SAMPLE_AMPLITUDES" + suffix])
                self.assertTrue(p["INITIAL_AMPLITUDE_MAP" + suffix].startswith("truth/"))
            for key, value in p.items():
                if key.endswith("STDDEV01") or "_PRIOR_GAUSSIAN_STDDEV" in key:
                    self.assertEqual(value, 0)
                if key.startswith("INIT_INDEX_MAP"):
                    self.assertTrue(value.startswith("truth/"))

    def test_64_band_subset_keeps_reference_channels_and_original_maps(self):
        manifest = fixture_manifest()
        subset = generator.select_channels(manifest, 64)
        frequencies = [channel["frequency_GHz"] for channel in subset]
        self.assertEqual(len(subset), 64)
        self.assertEqual(len(set(frequencies)), 64)
        self.assertEqual(frequencies, sorted(frequencies))
        self.assertEqual(frequencies[0], 50)
        self.assertEqual(frequencies[-1], 2000)
        self.assertIn(545, frequencies)
        counts = [0] * 4
        for frequency in frequencies:
            counts[min(3, int((frequency - 50) / ((2000 - 50) / 4)))] += 1
        self.assertTrue(all(count >= 8 for count in counts))
        for nside in (64, 128):
            for stage in generator.STAGES:
                p = generator.build_parameters(manifest, nside, stage, "input", "chains/subset", band_count=64)
                self.assertEqual(p["NUMBAND"], 64)
                self.assertEqual(p["NUM_PROC_PER_BAND"], 1)
                self.assertEqual(p["NUM_GROUPS"], 1)
                self.assertEqual(len([key for key in p if re.fullmatch(r"MAP\d+_0001", key)]), 64)
                self.assertNotIn("FREQ_C65", p)
                for i, source in enumerate(subset, 1):
                    suffix = generator.band_suffix(i)
                    self.assertEqual(p["FREQ_C" + suffix], source["frequency_GHz"])
                    kind = "observed" if generator.STAGES[stage]["observed"] else "signal"
                    self.assertEqual(p[f"MAP{suffix}_0001"], f"maps/{kind}_{source['band']:03d}.fits")
                    self.assertEqual(p["NOISE_RMS" + suffix], f"maps/sigma_{source['band']:03d}.fits")
                for i in range(1, 6):
                    suffix = generator.band_suffix(i)
                    reference = generator.band_suffix(p["REFERENCE_BAND" + suffix])
                    self.assertEqual(p["FREQ_C" + reference], p["REFERENCE_FREQUENCY" + suffix])

    def test_invalid_band_selection_fails(self):
        for count in (0, 63, 65, 132):
            with self.assertRaises(ValueError):
                generator.select_channels(fixture_manifest(), count)
        manifest = fixture_manifest()
        manifest["resolutions"]["64"]["components"][0]["reference_band"] = 33
        with self.assertRaises(ValueError):
            generator.select_channels(manifest, 64)

    def test_native_pixel_window_paths_and_implicit_cwd_dependency(self):
        manifest = fixture_manifest()
        for nside in (64, 128):
            filename = f"pixel_window_n{nside:04d}.fits"
            manifest["resolutions"][str(nside)]["pixel_window"] = filename
            for band_count in (64, 131):
                p = generator.build_parameters(manifest, nside, "globals_optimize", "input", "chains/test",
                                               band_count=band_count)
                self.assertEqual(p["PIXEL_WINDOW"], filename)
                text = generator.render(p, "globals_optimize", nside)
                self.assertIn(f"Legacy read_pixwin also needs {filename} in the process working directory.", text)
                self.assertIn("ell = 0 through 4 * NSIDE", text)

    def test_release_sequence_and_sample_semantics(self):
        p = self.parameters("amplitudes_sample")
        self.assertTrue(p["SAMPLE_AMPLITUDES01"])
        self.assertEqual(p["EMISSIVITY_PRIOR_GAUSSIAN_STDDEV01"], 0)
        self.assertFalse(p["GLOBAL_SAMPLE_PARAM03_01"])
        p = self.parameters("amplitudes_optimize")
        self.assertEqual(p["OPERATION"], "optimize")
        self.assertTrue(p["SAMPLE_AMPLITUDES01"])
        self.assertFalse(p["GLOBAL_SAMPLE_PARAM03_01"])
        self.assertEqual(p["EMISSIVITY_PRIOR_GAUSSIAN_STDDEV01"], 0)
        p = self.parameters("forecast_sample")
        self.assertGreater(p["EMISSIVITY_PRIOR_GAUSSIAN_STDDEV01"], 0)
        self.assertGreater(p["FREQUENCY_PRIOR_GAUSSIAN_STDDEV04"], 0)
        self.assertEqual(p["CURVATURE_PRIOR_GAUSSIAN_STDDEV02"], 0)
        self.assertTrue(p["GLOBAL_SAMPLE_PARAM03_01"])
        self.assertTrue(p["GLOBAL_SAMPLE_PARAM04_01"])
        self.assertTrue(p["GLOBAL_SAMPLE_PARAM04_02"])
        for stage, options in generator.STAGES.items():
            p = self.parameters(stage)
            self.assertEqual(p["OPERATION"], options["operation"])
            self.assertEqual(p["NUM_STEP_WITH_ML_SEARCH"], 0)
            self.assertEqual(p["RMS_SPEC_INDEX_RANDOM_WALK"], 0)
        self.assertEqual(self.parameters("no_globals_regression")["NUM_GLOBAL_COMPONENTS"], 0)

    def test_required_literal_legacy_keys_are_complete(self):
        # Inspect production readers. Exclude disabled commented code and keys for
        # conditional formats not used here; dynamically constructed band/component
        # keys are checked separately above.
        files = ["commander", "comm_data_mod", "comm_bp_mod", "comm_N_mult_mod",
                 "comm_beam_mod", "comm_Cl_util_mod", "comm_fg_mod", "comm_cgd_mod",
                 "comm_cgd_precond_mod", "comm_genvec_mod", "comm_chisq_mod",
                 "comm_fg_component_mod", "comm_Cl_sampling_mod", "comm_signal_mod"]
        conditional = {"INV_N_DENSE_ORDERING", "COLOR_CORR_UPDATE_INTERVAL",
                       "FOSSIL_SI_CMB", "FOSSIL_CMB_MIN_L"}
        required = set()
        for name in files:
            source = (ROOT / f"src/commander/{name}.f90").read_text()
            source = "\n".join(line for line in source.splitlines() if not line.lstrip().startswith("!"))
            required.update(re.findall(r"call get_parameter\(paramfile,\s*'([A-Z_]+)'\s*,", source, re.I))
        missing = required - self.parameters().keys() - conditional
        self.assertEqual(missing, set(), f"Missing literal legacy keys: {sorted(missing)}")

    def test_prevents_mismatched_beam_units_and_unannounced_cmb(self):
        for key, value in [("beam_fwhm_deg", 2.0), ("units", "uK_cmb"),
                           ("primary_cmb_anisotropy", True)]:
            manifest = fixture_manifest()
            manifest["metadata"][key] = value
            with self.assertRaises(ValueError):
                generator.validate_manifest(manifest, 128)

    def test_path_preflight_rejects_truncation(self):
        p = self.parameters()
        p["INPUT_DIRECTORY"] = "/" + "x" * 500
        with self.assertRaises(ValueError):
            generator.preflight_paths(p)


if __name__ == "__main__":
    unittest.main()
