#!/usr/bin/env python3
"""Test 64-band smoke references without running Commander or creating maps."""
from contextlib import redirect_stdout
import csv
import io
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

import fossil_remote_acceptance as acceptance
from fossil_parameter_tests import fixture_manifest


class SubsetTests(unittest.TestCase):
    def make_export(self, root):
        manifest = fixture_manifest()
        rng = np.random.default_rng(813)
        scaled_design = rng.normal(size=(131, 6))
        scales = np.array([1e-4, 1e-8, 1e-6, 1.0, 1.0, 1.0])
        design = scaled_design / scales
        weights = np.linspace(0.5, 3.0, 131)
        truth = np.array([1.2e-4, 2e-8, 1.77e-6, 0.346, 1.0, 1.0])
        weighted_sum = weights * (design @ truth) + rng.normal(size=131) * np.sqrt(weights)
        with (root / "global_response.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["band", "frequency_GHz", *acceptance.AMPLITUDES])
            # Source files need not have the same row ordering as active bands.
            for i in rng.permutation(131):
                writer.writerow([i, 50 + i * 15, *design[i]])
        bands = [dict(band=i, frequency_GHz=50 + i * 15, weight=float(weights[i]),
                      weighted_sum=float(weighted_sum[i])) for i in reversed(range(131))]
        (root / "global_reference.json").write_text(json.dumps({"64": dict(bands=bands)}))
        return manifest, scaled_design, weights, weighted_sum, scales

    def parameters(self, manifest, band_count=131):
        return acceptance.generator.build_parameters(manifest, 64, "globals_sample", "input", "chains/test",
                                                     band_count=band_count)

    def manual_posterior(self, selected_ids, parameters, scaled_design, weights, weighted_sum, scales):
        # This independent calculation selects the known original source rows.
        # It never reads the production frequency matcher or its chosen rows.
        matrix = scaled_design[selected_ids]
        precision = matrix.T @ (weights[selected_ids, None] * matrix)
        rhs = matrix.T @ weighted_sum[selected_ids]
        for i in range(6):
            sigma = parameters[f"GLOBAL_PRIOR_SIGMA{i + 1:02d}"]
            if sigma > 0:
                prior_precision = (scales[i] / sigma) ** 2
                precision[i, i] += prior_precision
                rhs[i] += prior_precision * parameters[f"GLOBAL_PRIOR_MEAN{i + 1:02d}"] / scales[i]
        return np.linalg.solve(precision, rhs), np.linalg.inv(precision), precision

    def assert_posterior_matches(self, reference, expected_mean, expected_covariance):
        covariance = np.linalg.solve(reference["chol"].T, np.linalg.solve(reference["chol"], np.eye(6)))
        np.testing.assert_allclose(reference["mean"], expected_mean, rtol=1e-12, atol=1e-14)
        np.testing.assert_allclose(covariance, expected_covariance, rtol=1e-12, atol=1e-14)

    def test_64_band_mean_covariance_use_selected_original_sources(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest, matrix, weights, sums, scales = self.make_export(root)
            parameters = self.parameters(manifest, 64)
            # Include Gaussian priors to check selection and prior accumulation.
            parameters["GLOBAL_PRIOR_SIGMA02"] = 3e-8
            parameters["GLOBAL_PRIOR_SIGMA04"] = 0.5
            selected = acceptance.generator.select_channels(manifest, 64)
            ids = np.array([channel["band"] for channel in selected])
            self.assertEqual(len(ids), 64)
            self.assertIn(0, ids)
            self.assertIn(33, ids)  # Exact 545-GHz dust reference.
            self.assertIn(130, ids)
            self.assertFalse(np.array_equal(ids, np.arange(64)))
            expected_mean, expected_covariance, expected_precision = self.manual_posterior(
                ids, parameters, matrix, weights, sums, scales)
            reference = acceptance.gaussian_reference(root, 64, parameters)
            self.assertEqual(reference["source_bands"], ids.tolist())
            self.assertEqual(reference["frequencies_GHz"], [50 + 15 * int(i) for i in ids])
            self.assert_posterior_matches(reference, expected_mean, expected_covariance)
            _, _, wrong_precision = self.manual_posterior(
                np.arange(64), parameters, matrix, weights, sums, scales)
            self.assertGreater(np.linalg.norm(expected_precision - wrong_precision), 1.0)

    def test_131_band_default_is_unchanged_and_matches_full_reference(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest, matrix, weights, sums, scales = self.make_export(root)
            default = acceptance.generator.build_parameters(manifest, 64, "globals_sample", "input", "chains/test")
            explicit = self.parameters(manifest, 131)
            self.assertEqual(default, explicit)
            actual = acceptance.gaussian_reference(root, 64, default)
            repeated = acceptance.gaussian_reference(root, 64, explicit)
            np.testing.assert_array_equal(actual["mean"], repeated["mean"])
            np.testing.assert_array_equal(actual["chol"], repeated["chol"])
            self.assertEqual(actual["source_bands"], list(range(131)))
            mean, covariance, _ = self.manual_posterior(np.arange(131), default, matrix, weights, sums, scales)
            self.assert_posterior_matches(actual, mean, covariance)

    def test_invalid_active_frequency_and_duplicates_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest, *_ = self.make_export(root)
            parameters = self.parameters(manifest, 64)
            parameters["FREQ_C03"] = 9999.0
            with self.assertRaisesRegex(ValueError, "source response matches"):
                acceptance.gaussian_reference(root, 64, parameters)
            parameters = self.parameters(manifest, 64)
            parameters["FREQ_C02"] = parameters["FREQ_C01"]
            with self.assertRaisesRegex(ValueError, "Duplicate active frequency"):
                acceptance.gaussian_reference(root, 64, parameters)

    def test_missing_or_inconsistent_source_statistics_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest, *_ = self.make_export(root)
            parameters = self.parameters(manifest, 64)
            path = root / "global_reference.json"
            original = json.loads(path.read_text())
            modified = json.loads(path.read_text())
            modified["64"]["bands"] = [row for row in modified["64"]["bands"] if row["band"] != 130]
            path.write_text(json.dumps(modified))
            with self.assertRaisesRegex(ValueError, "Missing sufficient statistics"):
                acceptance.gaussian_reference(root, 64, parameters)
            for row in original["64"]["bands"]:
                if row["band"] == 130:
                    row["frequency_GHz"] += 0.1
            path.write_text(json.dumps(original))
            with self.assertRaisesRegex(ValueError, "Inconsistent response/statistics frequency"):
                acceptance.gaussian_reference(root, 64, parameters)

    def test_ambiguous_source_frequency_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest, *_ = self.make_export(root)
            path = root / "global_response.csv"
            with path.open() as handle:
                reader = csv.DictReader(handle)
                fields, rows = reader.fieldnames, list(reader)
            for row in rows:
                if int(row["band"]) == 1:
                    row["frequency_GHz"] = 50.0
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fields)
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaisesRegex(ValueError, "2 source response matches"):
                acceptance.gaussian_reference(root, 64, self.parameters(manifest, 64))

    def test_runner_reports_64_ranks_and_actual_frequencies(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest, *_ = self.make_export(root)
            manifest_path = root / "manifest.json"
            manifest_path.write_text(json.dumps(manifest))
            args = SimpleNamespace(manifest=manifest_path, run_root=root / "runs", band_count=64,
                                   legacy_commander=None, nside=[64], stage=["globals_optimize"], seed=20261102,
                                   proc_per_band=1, iterations=12, check_only=False,
                                   commander=Path("commander"), launcher="mpirun -np {ranks}", timeout=100)
            with patch.object(acceptance.generator, "preflight_paths"), \
                    patch.object(acceptance, "execute") as execute, \
                    patch.object(acceptance, "check_stage", return_value={"stage": "globals_optimize"}), \
                    redirect_stdout(io.StringIO()):
                acceptance.run(args)
            report = json.loads((args.run_root / "acceptance_report.json").read_text())[0]
            expected_frequencies = [c["frequency_GHz"] for c in acceptance.generator.select_channels(manifest, 64)]
            self.assertEqual(report["ranks"], 64)
            self.assertEqual(report["band_count"], 64)
            self.assertEqual(report["frequencies_GHz"], expected_frequencies)
            self.assertEqual(Path(report["config"]).name, "n64_globals_optimize_64bands.txt")
            self.assertEqual(execute.call_args.args[2], 64)
            self.assertEqual(execute.call_args.args[4].name, "n64_globals_optimize_64bands_s20261102")

    def test_64_band_runner_rejects_legacy_99_band_control(self):
        args = SimpleNamespace(band_count=64, legacy_commander=Path("old_commander"))
        with self.assertRaisesRegex(ValueError, "99-band legacy control"):
            acceptance.run(args)

    def test_amplitudes_optimize_acceptance_keeps_shapes_fixed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            values = np.column_stack((np.arange(1, 13), np.ones((12, 6)), np.tile([1.245, 0.86, 18.8], (12, 1))))
            path = root / "global_components_c0001.txt"
            header = "iteration " + " ".join(acceptance.AMPLITUDES + acceptance.SHAPES)
            np.savetxt(path, values, header=header)
            (root / "commander.log").write_text("Chain no. 1 -- Starting CG search\n"
                                                "Chain no. 1 -- number of CGD iterations = 6\n" * 12)
            result = acceptance.check_stage(root, 64, "amplitudes_optimize", {"NUM_GIBBS_ITER": 12}, root)
            self.assertTrue(result["finite"])
            self.assertEqual(result["cg_searches_started"], 12)
            self.assertEqual(result["cg_searches_completed"], 12)
            values[-1, -1] = 19.0
            np.savetxt(path, values, header=header)
            with self.assertRaisesRegex(AssertionError, "fixed global shape"):
                acceptance.check_stage(root, 64, "amplitudes_optimize", {"NUM_GIBBS_ITER": 12}, root)

    def test_amplitudes_optimize_rejects_missing_or_incomplete_cg_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            values = np.column_stack((np.arange(1, 13), np.ones((12, 6)), np.tile([1.245, 0.86, 18.8], (12, 1))))
            header = "iteration " + " ".join(acceptance.AMPLITUDES + acceptance.SHAPES)
            np.savetxt(root / "global_components_c0001.txt", values, header=header)
            parameters = {"NUM_GIBBS_ITER": 12}
            with self.assertRaisesRegex(AssertionError, "requires commander.log"):
                acceptance.check_stage(root, 64, "amplitudes_optimize", parameters, root)
            (root / "commander.log").write_text("Chain no. 1 -- Starting CG search\n" * 12
                                                + "Chain no. 1 -- number of CGD iterations = 6\n" * 11)
            with self.assertRaisesRegex(AssertionError, "finish for every accepted iteration"):
                acceptance.check_stage(root, 64, "amplitudes_optimize", parameters, root)


if __name__ == "__main__":
    unittest.main()
