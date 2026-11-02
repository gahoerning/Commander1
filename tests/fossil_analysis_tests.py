#!/usr/bin/env python3
"""Test acceptance calculations with synthetic chains, not Commander outputs."""
import csv
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace

import numpy as np

from fossil_diagnostics import diagnostics
import fossil_remote_acceptance as acceptance


class DiagnosticTests(unittest.TestCase):
    def test_independent_chains(self):
        values = np.random.default_rng(41).normal(size=(4, 4000))
        d = diagnostics(values)
        self.assertLess(d["rhat"], 1.01)
        self.assertGreater(d["bulk_ess"], 8000)
        self.assertGreater(d["tail_ess"], 8000)

    def test_sticky_chains_fail_ess(self):
        rng = np.random.default_rng(77)
        values = rng.normal(size=(4, 4000))
        for i in range(1, values.shape[1]):
            values[:, i] = 0.995 * values[:, i - 1] + np.sqrt(1 - 0.995**2) * values[:, i]
        d = diagnostics(values)
        self.assertLess(d["bulk_ess"], 400)

    def test_different_modes_fail_rhat(self):
        values = np.random.default_rng(13).normal(size=(4, 4000))
        values[:2] += 2
        self.assertGreater(diagnostics(values)["rhat"], 1.05)

    def test_constant_free_parameter_cannot_pass_forecast(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            directories = []
            header = "iteration " + " ".join(acceptance.AMPLITUDES + acceptance.SHAPES)
            for seed in range(4):
                directory = root / str(seed)
                directory.mkdir()
                values = np.random.default_rng(seed).normal(size=(2500, 9))
                values[:, 1] = 2e-8  # A failed free-mu sampler returning ML forever.
                np.savetxt(directory / "global_components_c0001.txt",
                           np.column_stack((np.arange(1, 2501), values)), header=header)
                directories.append(directory)
            args = SimpleNamespace(chain_directory=directories, chain=1, warmup=1000,
                                   fixed_parameter=[], output=root / "diagnostics.json")
            with redirect_stdout(io.StringIO()), self.assertRaises(SystemExit):
                acceptance.diagnose(args)
            result = json.loads(args.output.read_text())
            self.assertFalse(result["accepted"])
            self.assertFalse(result["parameters"]["mu"]["accepted"])

    def test_duplicate_chains_are_not_independent(self):
        args = SimpleNamespace(chain_directory=[Path("same_chain")] * 4)
        with self.assertRaises(ValueError):
            acceptance.diagnose(args)


class AnalyticAcceptanceTests(unittest.TestCase):
    def fixture(self, root):
        rng = np.random.default_rng(62)
        design = rng.normal(size=(18, 6))
        scales = np.array([1e-4, 1e-8, 1e-6, 1, 1, 1])
        design = design / scales
        weights = np.linspace(0.5, 3, 18)
        mean = np.array([1.2e-4, 2e-8, 1.77e-6, 0.346, 1, 1])
        weighted_sum = weights * (design @ mean)
        with (root / "global_response.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["band", "frequency_GHz", *acceptance.AMPLITUDES])
            for i, row in enumerate(design):
                writer.writerow([i, 50 + i * 15, *row])
        reference = {"64": dict(N_unmasked=123,
                                  bands=[dict(band=i, weight=float(weights[i]),
                                              weighted_sum=float(weighted_sum[i])) for i in range(18)])}
        (root / "global_reference.json").write_text(json.dumps(reference))
        parameters = {"NUM_GIBBS_ITER": 2500}
        for i in range(6):
            suffix = f"{i + 1:02d}"
            parameters["GLOBAL_SCALE" + suffix] = scales[i]
            parameters["GLOBAL_PRIOR_MEAN" + suffix] = mean[i]
            parameters["GLOBAL_PRIOR_SIGMA" + suffix] = -1
        return rng, parameters, mean, scales

    def test_exact_gaussian_chain_passes_and_ml_repeats_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rng, parameters, truth, scales = self.fixture(root)
            reference = acceptance.gaussian_reference(root, 64, parameters)
            gaussian = rng.normal(size=(2500, 6))
            amplitudes = (reference["mean"] + np.linalg.solve(reference["chol"].T, gaussian.T).T) * scales
            shapes = np.tile([1.245, 0.86, 18.8], (2500, 1))
            trace = np.column_stack((np.arange(1, 2501), amplitudes, shapes))
            target = root / "global_components_c0001.txt"
            header = "iteration " + " ".join(acceptance.AMPLITUDES + acceptance.SHAPES)
            np.savetxt(target, trace, header=header)
            result = acceptance.check_stage(root, 64, "globals_sample", parameters, root)
            self.assertLess(result["max_mean_z"], 6)
            trace[:, 1:7] = truth
            np.savetxt(target, trace, header=header)
            with self.assertRaises(AssertionError):
                acceptance.check_stage(root, 64, "globals_sample", parameters, root)

    def test_wrong_posterior_covariance_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rng, parameters, truth, scales = self.fixture(root)
            reference = acceptance.gaussian_reference(root, 64, parameters)
            gaussian = 2 * rng.normal(size=(2500, 6))
            amplitudes = (reference["mean"] + np.linalg.solve(reference["chol"].T, gaussian.T).T) * scales
            trace = np.column_stack((np.arange(1, 2501), amplitudes,
                                     np.tile([1.245, 0.86, 18.8], (2500, 1))))
            header = "iteration " + " ".join(acceptance.AMPLITUDES + acceptance.SHAPES)
            np.savetxt(root / "global_components_c0001.txt", trace, header=header)
            with self.assertRaises(AssertionError):
                acceptance.check_stage(root, 64, "globals_sample", parameters, root)


if __name__ == "__main__":
    unittest.main()
