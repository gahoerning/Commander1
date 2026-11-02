#!/usr/bin/env python3
"""Run full Commander FOSSIL acceptance tests on a configured MPI machine.

This script does not build Commander. The default 131-band run requires 131 MPI
ranks. Use --launcher 'srun --ntasks={ranks}' for a Slurm allocation. All process
arguments are passed without a shell. A nonzero exit or invalid output fails.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
from pathlib import Path
import shlex
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("fossil_generator", ROOT / "examples/fossil/generate_parameters.py")
generator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(generator)

AMPLITUDES = ("delta_t", "mu", "y", "cib", "crr", "co")
SHAPES = ("y_p01", "cib_p01", "cib_p02")


def load_trace(directory: Path, chain: int = 1):
    path = directory / f"global_components_c{chain:04d}.txt"
    array = np.loadtxt(path, ndmin=2)
    expected = 1 + len(AMPLITUDES) + len(SHAPES)
    if array.shape[1] != expected or not np.isfinite(array).all():
        raise AssertionError(f"Invalid global trace {path}: shape {array.shape}, expected {expected} columns")
    if not np.all(np.diff(array[:, 0]) == 1):
        raise AssertionError(f"Nonconsecutive global trace iterations: {path}")
    headers = [line.lstrip("# ").split() for line in path.read_text().splitlines()
               if line.startswith("#") and line.lstrip("# ").startswith("iteration")]
    if not headers:
        raise AssertionError(f"Missing named trace header: {path}")
    names = headers[0]
    canonical = ("iteration",) + AMPLITUDES + SHAPES
    if len(set(names)) != expected or set(names) != set(canonical):
        raise AssertionError(f"Unexpected trace parameter names: {names}")
    return array[:, [names.index(name) for name in canonical]]


def gaussian_reference(export_root: Path, nside: int, parameters: dict):
    with (export_root / "global_response.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    responses = {int(row["band"]): [float(row[key]) for key in AMPLITUDES] for row in rows}
    reference = json.loads((export_root / "global_reference.json").read_text())[str(nside)]
    bands = reference["bands"]
    design = np.array([responses[int(b["band"])] for b in bands])
    weights = np.array([b["weight"] for b in bands], dtype=float)
    weighted_sum = np.array([b["weighted_sum"] for b in bands], dtype=float)
    scales = np.array([parameters[f"GLOBAL_SCALE{i:02d}"] for i in range(1, 7)])
    scaled = design * scales
    precision = scaled.T @ (weights[:, None] * scaled)
    rhs = scaled.T @ weighted_sum
    for i in range(6):
        sigma = parameters[f"GLOBAL_PRIOR_SIGMA{i + 1:02d}"]
        if sigma > 0:
            prior_precision = (scales[i] / sigma) ** 2
            precision[i, i] += prior_precision
            rhs[i] += prior_precision * parameters[f"GLOBAL_PRIOR_MEAN{i + 1:02d}"] / scales[i]
    diagonal = np.sqrt(np.diag(precision))
    correlation = precision / diagonal[:, None] / diagonal[None, :]
    correlation_chol = np.linalg.cholesky(correlation)
    mean = np.linalg.solve(correlation_chol.T, np.linalg.solve(correlation_chol, rhs / diagonal)) / diagonal
    chol = diagonal[:, None] * correlation_chol
    return dict(scales=scales, chol=chol, mean=mean,
                scaled_condition_number=float(np.linalg.cond(precision)),
                normalized_condition_number=float(np.linalg.cond(correlation)))


def check_stage(export_root: Path, nside: int, stage: str, parameters: dict, directory: Path):
    if stage == "no_globals_regression":
        if not directory.is_dir() or not any(directory.glob("*.fits")):
            raise AssertionError(f"No no-globals FITS outputs: {directory}")
        return dict(stage=stage, status="completed; compare with legacy binary separately")
    trace = load_trace(directory)
    if len(trace) < parameters["NUM_GIBBS_ITER"]:
        raise AssertionError("Trace omitted accepted Gibbs iterations; scalar output must not use map thinning")
    amplitudes = trace[:, 1:7]
    result = dict(stage=stage, iterations=len(trace), finite=True)
    if stage == "globals_optimize":
        truth = np.array([parameters[f"GLOBAL_PRIOR_MEAN{i:02d}"] for i in range(1, 7)])
        scales = np.array([parameters[f"GLOBAL_SCALE{i:02d}"] for i in range(1, 7)])
        error = np.max(np.abs((amplitudes[-1] - truth) / scales))
        if error > 1e-5:
            raise AssertionError(f"Noiseless fixed-shape amplitude recovery failed: scaled error {error}")
        result["max_scaled_recovery_error"] = float(error)
    elif stage == "globals_sample":
        reference = gaussian_reference(export_root, nside, parameters)
        # Every draw in this fixture is an independent exact conditional draw.
        # Discarding 100 makes the criterion insensitive to initialization policy.
        draws = amplitudes[min(100, len(amplitudes) // 10):]
        if len(draws) < 500:
            raise AssertionError("Gaussian covariance acceptance requires at least 500 retained draws")
        whitened = (draws / reference["scales"] - reference["mean"]) @ reference["chol"]
        mean_z = np.sqrt(len(draws)) * whitened.mean(axis=0)
        covariance = np.cov(whitened, rowvar=False)
        covariance_z = (covariance - np.eye(6)) / np.sqrt((1 + np.eye(6)) / (len(draws) - 1))
        if np.max(np.abs(mean_z)) > 6 or np.max(np.abs(covariance_z)) > 6:
            raise AssertionError(f"Gaussian posterior mismatch: max mean z={np.max(np.abs(mean_z))}, "
                                 f"max covariance z={np.max(np.abs(covariance_z))}")
        if np.any(np.std(draws, axis=0) == 0):
            raise AssertionError("Sample mode returned fixed amplitudes instead of stochastic posterior draws")
        result.update(max_mean_z=float(np.max(np.abs(mean_z))),
                      max_covariance_z=float(np.max(np.abs(covariance_z))),
                      scaled_condition_number=reference["scaled_condition_number"],
                      normalized_condition_number=reference["normalized_condition_number"])
    if stage in ("globals_optimize", "globals_sample", "amplitudes_sample"):
        expected_shapes = np.array([1.245, 0.86, 18.8])
        if not np.allclose(trace[:, 7:], expected_shapes, rtol=0, atol=1e-12):
            raise AssertionError("A fixed global shape parameter changed")
    elif stage in ("global_shapes_sample", "forecast_sample"):
        if np.any(np.std(trace[:, 7:], axis=0) == 0):
            raise AssertionError("A requested global shape parameter was never sampled")
    return result


def compare_legacy(patched: Path, legacy: Path):
    from astropy.io import fits
    files = sorted(patched.glob("*.fits"))
    if not files:
        raise AssertionError("No FITS outputs available for no-globals regression comparison")
    for current in files:
        old = legacy / current.name
        if not old.is_file():
            raise AssertionError(f"Legacy output missing: {old}")
        with fits.open(current) as lhs, fits.open(old) as rhs:
            if len(lhs) != len(rhs):
                raise AssertionError(f"FITS extension count changed: {current.name}")
            for left, right in zip(lhs, rhs):
                if left.data is None:
                    continue
                fields = left.data.dtype.names or (None,)
                for field in fields:
                    a = left.data if field is None else left.data[field]
                    b = right.data if field is None else right.data[field]
                    nan = a.dtype.kind in "fc"
                    if not np.array_equal(a, b, equal_nan=nan):
                        raise AssertionError(f"No-globals output changed: {current.name}, column {field}")
    return dict(files=len(files), exact_data_match=True)


def execute(binary: Path, launcher: str, ranks: int, config: Path, directory: Path,
            cwd: Path, timeout: int):
    command = shlex.split(launcher.format(ranks=ranks)) + [str(binary), str(config)]
    directory.mkdir(parents=True, exist_ok=True)
    log = directory / "commander.log"
    with log.open("w") as output:
        completed = subprocess.run(command, cwd=cwd, stdout=output, stderr=subprocess.STDOUT, timeout=timeout)
    if completed.returncode:
        raise RuntimeError(f"Commander failed with status {completed.returncode}; inspect {log}")
    # Some legacy Fortran STOP paths return zero; missing outputs still fail checks.
    return command


def run(args):
    manifest = json.loads(args.manifest.read_text())
    export_root = args.manifest.resolve().parent
    args.run_root.mkdir(parents=True, exist_ok=True)
    reports = []
    for nside in args.nside:
        for stage in args.stage:
            directory = args.run_root.resolve() / f"n{nside}_{stage}_s{args.seed}"
            config = args.run_root.resolve() / f"n{nside}_{stage}.par"
            p = generator.build_parameters(manifest, nside, stage, str(export_root), str(directory),
                                           seed=args.seed, proc_per_band=args.proc_per_band,
                                           iterations=args.iterations)
            legacy_control = stage == "no_globals_regression" and args.legacy_commander is not None
            if legacy_control:
                # Legacy executables cannot express 100+ band identifiers or the
                # newly requested exact-SED path. Compare their ordinary behavior
                # under an explicit <=99-band control instead.
                p["NUMBAND"] = 99
                p["FOSSIL_EXACT_LOCAL_SEDS"] = False
            generator.preflight_paths(p, verify_inputs=True)
            if len(str(config).encode("utf-8")) >= 512:
                raise ValueError("Parameter filename exceeds Commander's 511-byte limit")
            config.write_text(generator.render(p, stage, nside))
            ranks = p["NUMBAND"] * args.proc_per_band
            if not args.check_only:
                if directory.exists() and any(directory.iterdir()):
                    raise FileExistsError(f"Refusing to append to previous chain: {directory}")
                execute(args.commander.resolve(), args.launcher, ranks, config, directory,
                        args.run_root.resolve(), args.timeout)
            report = check_stage(export_root, nside, stage, p, directory)
            report.update(nside=nside, ranks=ranks, config=str(config),
                          manifest_signature=manifest.get("signature"))
            if legacy_control:
                report["control"] = "99 bands; exact-local-SED option disabled in both legacy and patched runs"
            if stage == "no_globals_regression" and args.legacy_commander:
                old_directory = directory.with_name(directory.name + "_legacy")
                old_config = config.with_name(config.stem + "_legacy.par")
                p["CHAIN_DIRECTORY"] = str(old_directory)
                old_config.write_text(generator.render(p, stage, nside))
                if not args.check_only:
                    if old_directory.exists() and any(old_directory.iterdir()):
                        raise FileExistsError(f"Refusing to append to previous legacy chain: {old_directory}")
                    execute(args.legacy_commander.resolve(), args.launcher, ranks, old_config,
                            old_directory, args.run_root.resolve(), args.timeout)
                report["legacy_comparison"] = compare_legacy(directory, old_directory)
            reports.append(report)
            print(json.dumps(report), flush=True)
    (args.run_root / "acceptance_report.json").write_text(json.dumps(reports, indent=2) + "\n")


def diagnose(args):
    from fossil_diagnostics import diagnostics
    if len(args.chain_directory) < 4:
        raise ValueError("Forecast acceptance requires four independent chains with distinct seeds")
    if len({path.resolve() for path in args.chain_directory}) != len(args.chain_directory):
        raise ValueError("Each independent chain must come from a distinct directory")
    traces = [load_trace(path, args.chain) for path in args.chain_directory]
    for i, trace in enumerate(traces):
        if any(np.array_equal(trace, previous) for previous in traces[:i]):
            raise ValueError("Duplicated traces cannot be counted as independent chains")
    size = min(len(trace) for trace in traces) - args.warmup
    if size < 20:
        raise ValueError("Too few post-warmup draws")
    results = {}
    passed = True
    for column, label in enumerate(AMPLITUDES + SHAPES, 1):
        values = np.array([trace[args.warmup:args.warmup + size, column] for trace in traces])
        constant = bool(np.all(values == values[0, 0]))
        if label in args.fixed_parameter:
            if not constant:
                raise AssertionError(f"Parameter declared fixed varies in the trace: {label}")
            results[label] = dict(fixed=True)
            continue
        if constant:
            results[label] = dict(accepted=False, reason="A free parameter has no variation")
            passed = False
            continue
        diagnostic = diagnostics(values)
        valid = (np.isfinite(diagnostic["rhat"]) and diagnostic["rhat"] <= 1.01
                 and diagnostic["bulk_ess"] >= 1000 and diagnostic["tail_ess"] >= 400)
        diagnostic["accepted"] = bool(valid)
        results[label] = diagnostic
        passed &= valid
    result = dict(accepted=bool(passed), warmup=args.warmup, chains=len(traces), parameters=results)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    if not passed:
        raise SystemExit("Forecast diagnostics did not pass; do not report the chain as converged")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    runner = subparsers.add_parser("run")
    runner.add_argument("--commander", type=Path, required=True)
    runner.add_argument("--manifest", type=Path, required=True)
    runner.add_argument("--run-root", type=Path, required=True)
    runner.add_argument("--launcher", default="mpirun -np {ranks}")
    runner.add_argument("--nside", type=int, nargs="+", choices=[128, 64], default=[64])
    runner.add_argument("--stage", nargs="+", choices=list(generator.STAGES),
                        default=["globals_optimize", "globals_sample"])
    runner.add_argument("--proc-per-band", type=int, choices=[1, 2, 4], default=1)
    runner.add_argument("--iterations", type=int)
    runner.add_argument("--seed", type=int, default=20261102)
    runner.add_argument("--timeout", type=int, default=86400)
    runner.add_argument("--legacy-commander", type=Path)
    runner.add_argument("--check-only", action="store_true")
    runner.set_defaults(function=run)
    diagnostic = subparsers.add_parser("diagnose")
    diagnostic.add_argument("--chain-directory", type=Path, nargs="+", required=True)
    diagnostic.add_argument("--chain", type=int, default=1)
    diagnostic.add_argument("--warmup", type=int, default=1000)
    diagnostic.add_argument("--output", type=Path, default=Path("forecast_diagnostics.json"))
    diagnostic.add_argument("--fixed-parameter", nargs="*", choices=list(AMPLITUDES + SHAPES), default=[],
                            help="Explicitly declare fixed parameters in conditional control runs")
    diagnostic.set_defaults(function=diagnose)
    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
