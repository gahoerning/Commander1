#!/usr/bin/env python3
"""Compile and statistically verify the complete production FOSSIL global block.

The production SED, global sampler and Planck RNG are compiled unchanged.
Only parameter-file I/O, band metadata and single-process MPI are lightweight
fixtures. This does not validate the distributed map/residual integration.
"""
import argparse
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def compile_driver(work):
    compiler = shlex.split(os.environ.get("FC", "gfortran"))
    if not compiler or shutil.which(compiler[0]) is None:
        raise RuntimeError("A Fortran compiler is required; set FC and optional FFLAGS")
    flags = shlex.split(os.environ.get("FFLAGS", ""))
    support = (ROOT / "tests/fossil_global_support.f90").read_text()
    types, dependencies = support.split("\nmodule comm_utils", 1)
    (work / "healpix_types.f90").write_text(types)
    (work / "dependencies.f90").write_text("module comm_utils" + dependencies)
    executable = work / "test_fossil_global"
    sources = [work / "healpix_types.f90", ROOT / "src/include/rngmod.f90", work / "dependencies.f90",
               ROOT / "src/commander/comm_fossil_sed_mod.f90",
               ROOT / "src/commander/comm_fossil_global_mod.f90", ROOT / "tests/test_fossil_global.f90"]
    command = compiler + flags + ["-std=f2008", "-ffree-line-length-none", "-fcheck=all", "-O0",
                                  "-J", str(work), "-I", str(work)]
    subprocess.run(command + list(map(str, sources)) + ["-o", str(executable)], cwd=work, check=True)
    return executable


def common_config(nband, operation="sample"):
    config = dict(OPERATION=f"'{operation}'", POLARIZATION="F", NOISE_FORMAT="'rms'", SAMPLE_INSIDE_MASK="F",
                  GLOBAL_FIXED_INSTRUMENT="T", GLOBAL_SLICE_STEPS="4", NUM_FG_TEMP="0")
    for band in range(1, nband + 1):
        tag = f"{band:02d}"
        config.update({"GAIN_PRIOR_RMS"+tag: "0", "NOISEAMP_PRIOR_RMS"+tag: "0",
                       "FIX_MONOPOLE"+tag: "T", "FIX_DIPOLE"+tag: "T"})
    return config


def add_component(config, index, kind, amplitude, *, scale=1.0, sample=True,
                  mean=None, sigma=-1.0, shape=None, sample_shape=None, bounds=None,
                  shape_mean=None, shape_sigma=None, template=None):
    tag = f"{index:02d}"
    config.update({"GLOBAL_TYPE"+tag: f"'{kind}'", "GLOBAL_LABEL"+tag: f"'{kind}_{index}'",
                   "GLOBAL_INIT"+tag: repr(float(amplitude)), "GLOBAL_SCALE"+tag: repr(float(scale)),
                   "GLOBAL_SAMPLE_AMPLITUDE"+tag: "T" if sample else "F",
                   "GLOBAL_PRIOR_MEAN"+tag: repr(float(amplitude if mean is None else mean)),
                   "GLOBAL_PRIOR_SIGMA"+tag: repr(float(sigma))})
    if template is not None:
        config["GLOBAL_TEMPLATE"+tag] = f"'{template}'"
    for p, value in enumerate(shape or [], 1):
        key = f"{tag}_{p:02d}"
        config["GLOBAL_INIT_PARAM"+key] = repr(float(value))
        active = bool(sample_shape and sample_shape[p-1])
        config["GLOBAL_SAMPLE_PARAM"+key] = "T" if active else "F"
        if active:
            config["GLOBAL_PARAM_LOW"+key] = repr(float(bounds[p-1][0]))
            config["GLOBAL_PARAM_HIGH"+key] = repr(float(bounds[p-1][1]))
        if shape_mean is not None:
            config["GLOBAL_PARAM_PRIOR_MEAN"+key] = repr(float(shape_mean[p-1]))
        if shape_sigma is not None:
            config["GLOBAL_PARAM_PRIOR_SIGMA"+key] = repr(float(shape_sigma[p-1]))


class Runner:
    def __init__(self, executable, work, templates):
        self.executable, self.work, self.templates = executable, work, templates
        self.counter = 0

    def __call__(self, config, rows, mode, draws=1, seed=9123, error=None):
        self.counter += 1
        directory = self.work / f"case_{self.counter:02d}"
        directory.mkdir()
        parameters, likelihood = directory / "params.txt", directory / "likelihood.txt"
        parameters.write_text("\n".join(f"{key} = {value}" for key, value in config.items()) + "\n")
        np.savetxt(likelihood, rows, header=str(len(rows)), comments="", fmt="%.17e")
        result = subprocess.run([str(self.executable), str(parameters), str(likelihood), mode, str(draws),
                                 str(self.templates), str(seed)], cwd=directory, capture_output=True, text=True)
        (directory/"driver_stdout.txt").write_text(result.stdout)
        (directory/"driver_stderr.txt").write_text(result.stderr)
        if error is not None:
            if result.returncode == 0 or error not in result.stdout:
                raise AssertionError(f"Expected rejection {error!r}:\n{result.stdout}\n{result.stderr}")
            return directory
        if result.returncode:
            raise RuntimeError(result.stdout + result.stderr)
        return directory


def mbb(nu_GHz, beta, temperature=18.8):
    x = 6.62607015e-34 * np.asarray(nu_GHz) * 1e9 / (1.380649e-23 * temperature)
    return x**(3+beta) / np.expm1(x)


def check_gaussian(draws, mean, covariance):
    # Check in whitened coordinates so degeneracy and native unit differences
    # cannot hide a covariance error behind an absolute tolerance.
    cholesky = np.linalg.cholesky(covariance)
    whitened = np.linalg.solve(cholesky, (draws-mean).T).T
    n = len(draws)
    error = np.max(np.abs(whitened.mean(axis=0)))
    assert error < 6/np.sqrt(n), ("posterior mean", error, 6/np.sqrt(n))
    empirical = np.cov(whitened.T)
    error = np.max(np.abs(empirical-np.eye(len(mean))))
    assert error < 6*np.sqrt(2/n), ("posterior covariance", error, 6*np.sqrt(2/n))


def amplitude_cases(run, reference, quick):
    nu = reference[:, 0]
    design = reference[:, 1:3]
    frozen = 1.7 * reference[:, 5]
    mean_prior = np.array([0.9e-4, -0.3e-6])
    sigma_prior = np.array([0.6e-6, 1.1e-6])
    truth = np.array([1.2e-4, 2e-8])
    weight = np.linspace(3500, 14000, len(nu))
    data = weight * (design @ truth + frozen + np.linspace(-0.003, 0.002, len(nu)))

    def config(operation="sample", scale=(7e-5, 3e-7)):
        result = common_config(len(nu), operation)
        result["NUM_GLOBAL_COMPONENTS"] = "3"
        add_component(result, 1, "delta_t", truth[0], scale=scale[0], mean=mean_prior[0], sigma=sigma_prior[0])
        add_component(result, 2, "mu", truth[1], scale=scale[1], mean=mean_prior[1], sigma=sigma_prior[1])
        add_component(result, 3, "crr", 1.7, sample=False, template=run.templates/"crr_MJy_sr.txt")
        return result

    for multiplicity in (1, 4):
        weights, residual = multiplicity*weight, multiplicity*data
        precision = design.T @ (weights[:, None]*design) + np.diag(1/sigma_prior**2)
        rhs = design.T @ (residual-weights*frozen) + mean_prior/sigma_prior**2
        covariance = np.linalg.inv(precision)
        expected = np.linalg.solve(precision, rhs)
        rows = np.column_stack((nu, weights, residual))
        directory = run(config("optimize"), rows, "amplitudes")
        optimized = np.atleast_2d(np.loadtxt(directory/"draws.txt"))
        np.testing.assert_allclose(optimized[0, :2], expected, rtol=2e-11, atol=1e-17)
        assert optimized[0, 2] < 1e-12
        # The same physical prior is used once at both effective pixel counts.
        directory = run(config(), rows, "amplitudes", 10000 if quick else 50000)
        samples = np.loadtxt(directory/"draws.txt")
        assert np.max(samples[:, 2]) < 1e-12
        check_gaussian(samples[:, :2], expected, covariance)
        last = np.atleast_2d(np.loadtxt(directory/"global_components_c0001.txt"))[-1]
        assert last[3] == 1.7, ("frozen amplitude changed", last)
        print(f"PASS: joint Gaussian mean/covariance, frozen amplitude and single global prior (pixel factor {multiplicity})")

    # Numerical conditioning scales change neither the physical mean nor priors.
    directory = run(config("optimize", (2e-7, 5e-4)), np.column_stack((nu, weight, data)), "amplitudes")
    scaled = np.atleast_2d(np.loadtxt(directory/"draws.txt"))[0, :2]
    expected = np.linalg.solve(design.T@(weight[:, None]*design)+np.diag(1/sigma_prior**2),
                               design.T@(data-weight*frozen)+mean_prior/sigma_prior**2)
    np.testing.assert_allclose(scaled, expected, rtol=2e-11, atol=1e-17)
    print("PASS: arbitrary positive numerical scales preserve physical amplitudes")


def shape_config(nband, operation="sample", initial=0.55, active=True, upper=1.6):
    config = common_config(nband, operation)
    config["NUM_GLOBAL_COMPONENTS"] = "1"
    add_component(config, 1, "cib", 0.346, sample=False, shape=[initial, 18.8],
                  sample_shape=[active, False], bounds=[(0.2, upper), (1, 100)],
                  shape_mean=[1.02, 18.8], shape_sigma=[0.055, -1.0])
    return config


def shape_cases(run, reference, quick):
    nu = reference[:, 0]
    weights = np.linspace(40, 160, len(nu))
    observation = 0.346*mbb(nu, 0.86) + np.linspace(-0.01, 0.01, len(nu))
    rows = np.column_stack((nu, weights, weights*observation))
    grid = np.linspace(0.2, 1.6, 50001)
    model = 0.346*mbb(nu[:, None], grid[None, :])
    logp = -0.5*np.sum(weights[:, None]*(model-observation[:, None])**2, axis=0)
    logp -= 0.5*((grid-1.02)/0.055)**2
    probability = np.exp(logp-logp.max())
    probability /= np.trapezoid(probability, grid)
    cumulative = np.concatenate(([0], np.cumsum((probability[:-1]+probability[1:])*np.diff(grid)/2)))
    cumulative /= cumulative[-1]
    expected_mean = np.trapezoid(grid*probability, grid)
    expected_sigma = np.sqrt(np.trapezoid((grid-expected_mean)**2*probability, grid))
    directory = run(shape_config(len(nu)), rows, "shapes", 6000 if quick else 30000)
    draws = np.loadtxt(directory/"draws.txt")[100:]
    assert np.all((draws>=grid[0]) & (draws<=grid[-1]))
    assert abs(draws.mean()-expected_mean) < 7*expected_sigma/np.sqrt(len(draws))
    sorted_draws = np.sort(draws)
    empirical = np.arange(1, len(draws)+1)/len(draws)
    distance = np.max(np.abs(empirical-np.interp(sorted_draws, grid, cumulative)))
    assert distance < (0.04 if quick else 0.023), ("slice conditional CDF", distance)
    print(f"PASS: CIB beta conditional slice distribution matches numerical grid (CDF distance {distance:.4g})")
    directory = run(shape_config(len(nu), "optimize"), rows, "shapes")
    optimized = float(np.loadtxt(directory/"draws.txt"))
    assert abs(optimized-grid[np.argmax(logp)]) <= 2*(grid[1]-grid[0])
    directory = run(shape_config(len(nu), active=False), rows, "shapes", 2)
    np.testing.assert_allclose(np.loadtxt(directory/"draws.txt"), 0.55, rtol=0, atol=2e-14)
    # The constrained optimum must remain allowed to sit exactly on a boundary.
    boundary = shape_config(len(nu), "optimize", upper=0.7)
    directory = run(boundary, rows, "shapes")
    assert abs(float(np.loadtxt(directory/"draws.txt"))-0.7) < 2e-12
    print("PASS: global spectral optimization, frozen shape and boundary optimum")


def restart_cases(run, reference):
    nu = reference[:, 0]
    weights = np.linspace(80, 200, len(nu))
    observation = reference[:, 1]*1.2e-4 + 0.346*mbb(nu, 0.86)
    rows = np.column_stack((nu, weights, weights*observation))
    config = common_config(len(nu))
    config["NUM_GLOBAL_COMPONENTS"] = "2"
    # One Gaussian per iteration exercises the cached Box--Muller draw across
    # an odd checkpoint, mixed with uniform draws from the spectral sampler.
    add_component(config, 1, "delta_t", 1.2e-4, scale=1e-4, mean=1.2e-4, sigma=1e-5)
    add_component(config, 2, "cib", 0.346, sample=False, shape=[0.8, 18.8],
                  sample_shape=[True, False], bounds=[(0.2, 1.6), (1, 100)],
                  shape_mean=[1.02, 18.8], shape_sigma=[0.08, -1.0])
    run(config, rows, "restore")
    full = run(config, rows, "chain", 33)
    first = run(config, rows, "chain", 13)
    resumed = dict(config, GLOBAL_RESTART_FILE=f"'{first/'global_restart_c0001.txt'}'")
    second = run(resumed, rows, "chain", 20, seed=123456)
    expected, actual = np.loadtxt(full/"draws.txt"), np.loadtxt(second/"draws.txt")
    assert np.array_equal(np.loadtxt(first/"draws.txt"), expected[:13])
    assert np.array_equal(actual, expected[13:]), np.max(np.abs(actual-expected[13:]))
    assert int(np.atleast_2d(np.loadtxt(second/"global_components_c0001.txt"))[-1, 0]) == 33
    print("PASS: exact in-memory RNG restore and cross-process checkpoint continuation with cached Gaussian")


def error_cases(run, reference):
    nu = reference[:, 0]
    rows = np.column_stack((nu, np.ones(len(nu)), np.zeros(len(nu))))
    base = shape_config(len(nu))
    for key, value, message in (("NOISE_FORMAT", "'covariance'", "Use NOISE_FORMAT=rms"),
                                ("POLARIZATION", "T", "Stokes I only"),
                                ("GLOBAL_FIXED_INSTRUMENT", "F", "GLOBAL_FIXED_INSTRUMENT=true"),
                                ("GLOBAL_SCALE01", "0", "GLOBAL_SCALE must be positive"),
                                ("GLOBAL_SCALE01", "NaN", "Non-finite global amplitude"),
                                ("GLOBAL_INIT01", "NaN", "Non-finite global amplitude"),
                                ("GLOBAL_PRIOR_MEAN01", "NaN", "Non-finite global amplitude"),
                                ("GLOBAL_PRIOR_SIGMA01", "NaN", "Non-finite global amplitude"),
                                ("GLOBAL_INIT_PARAM01_01", "NaN", "Non-finite global spectral parameter"),
                                ("GLOBAL_PARAM_PRIOR_SIGMA01_01", "NaN", "Non-finite global spectral parameter"),
                                ("GLOBAL_PARAM_LOW01_01", "NaN", "Non-finite global spectral parameter bounds"),
                                ("GLOBAL_PARAM_LOW01_01", "1.7", "Invalid global spectral"),
                                ("GLOBAL_SLICE_STEPS", "0", "GLOBAL_SLICE_STEPS must be positive")):
        run(dict(base, **{key: value}), rows, "initialize", error=message)
    run(base, rows*np.array([1, 0, 1]), "shapes", error="No weighted pixels")
    bad = rows.copy(); bad[0, 1] = -1
    run(base, bad, "shapes", error="Invalid likelihood statistics")
    singular = common_config(len(nu)); singular["NUM_GLOBAL_COMPONENTS"] = "2"
    add_component(singular, 1, "delta_t", 1e-4)
    add_component(singular, 2, "delta_t", 1e-4)
    run(singular, rows, "chain", error="singular/ill-conditioned")
    # Checkpoint validation must fail during initialization, before a corrupted
    # cached Gaussian or out-of-bounds shape can enter a Markov transition.
    checkpoint = run(base, rows, "shapes")
    original = (checkpoint/"global_restart_c0001.txt").read_text().splitlines()
    corruptions = []
    lines = list(original); parts = lines[0].split(); parts[2] = "-1"; lines[0] = " ".join(parts)
    corruptions.append((lines, "Incompatible global restart header"))
    lines = list(original); parts = lines[1].split(); parts[1] = "NaN"; lines[1] = " ".join(parts)
    corruptions.append((lines, "Non-finite global amplitude"))
    lines = list(original); parts = lines[1].split(); parts[2] = "1.7"; lines[1] = " ".join(parts)
    corruptions.append((lines, "Invalid global spectral parameter bounds/start or restart"))
    lines = list(original); parts = lines[2].split(); parts[4] = "NaN"; lines[2] = " ".join(parts)
    corruptions.append((lines, "Invalid global restart RNG state"))
    lines = list(original); parts = lines[2].split(); parts[:4] = ["0"]*4; lines[2] = " ".join(parts)
    corruptions.append((lines, "Invalid global restart RNG state"))
    corruptions.append((["invalid header"], "Invalid global restart header"))
    for index, (lines, message) in enumerate(corruptions):
        path = run.work/f"invalid_restart_{index}.txt"
        path.write_text("\n".join(lines)+"\n")
        run(dict(base, GLOBAL_RESTART_FILE=f"'{path}'"), rows, "initialize", error=message)
    print("PASS: invalid configuration, nonfinite metadata, corrupted restart, likelihood and singular precision rejection")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("templates", type=Path, help="Path to templates/commander_spectral")
    parser.add_argument("--quick", action="store_true", help="Use fewer Monte Carlo samples for development")
    args = parser.parse_args()
    templates = args.templates.resolve()
    reference = np.loadtxt(templates/"fiducial_spectra_MJy_sr.txt")[[0, 4, 10, 20, 40, 70, 97, 130]]
    with tempfile.TemporaryDirectory(prefix="commander-fossil-global-") as temporary:
        work = Path(temporary)
        executable = compile_driver(work)
        run = Runner(executable, work, templates)
        amplitude_cases(run, reference, args.quick)
        shape_cases(run, reference, args.quick)
        restart_cases(run, reference)
        error_cases(run, reference)
    print("PASS: all isolated production FOSSIL global sampler tests")


if __name__ == "__main__":
    main()
