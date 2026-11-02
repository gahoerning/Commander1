#!/usr/bin/env python3
"""Validate physical FOSSIL spectra against independent simulation references.

The numerical source check evaluates the production Itoh expressions directly;
it is explicitly not a Fortran execution. If FC (default gfortran) is available,
also compile and execute the complete production module with a minimal dp fixture.
No HEALPix, MPI or complete Commander build is required for that isolated test.
"""
import argparse
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "src/commander/comm_fossil_sed_mod.f90"


def production_kernels(x):
    """Evaluate the actual kernel RHS expressions in the production source."""
    text = SOURCE.read_text()
    block = text.split("subroutine fossil_itoh_kernels(x, kernel)", 1)[1]
    block = block.split("end subroutine fossil_itoh_kernels", 1)[0]
    xbar = x / np.tanh(x / 2)
    s = x / np.sinh(x / 2)
    result = []
    for index in range(5):
        pattern = rf"kernel\({index}\)\s*=\s*(.*?)(?=\n\s*kernel\(|\Z)"
        expression = re.search(pattern, block, re.DOTALL).group(1)
        expression = expression.replace("&", " ").replace("_dp", "")
        expression = " ".join(expression.split())
        result.append(eval(expression, {"__builtins__": {}}, {"Xbar": xbar, "S": s}))
    return np.asarray(result)


def source_comparison(directory):
    reference = np.loadtxt(directory / "fiducial_spectra_MJy_sr.txt")
    assert reference.shape == (131, 7)
    assert np.allclose(reference[:, 0], np.arange(50, 2001, 15), rtol=0, atol=1e-12)
    # Constants are fixed by the simulation's explicit scientific convention.
    h, k, c, me, ev, t0 = 6.62607015e-34, 1.380649e-23, 299792458.0, 9.1093837139e-31, 1.602176634e-19, 2.725
    text = SOURCE.read_text()
    for constant, expected in (("h", h), ("k", k), ("c", c), ("me", me), ("ev", ev), ("T0", t0)):
        actual = float(re.search(rf"fossil_{constant}\s*=\s*([\d.eE+-]+)_dp", text).group(1))
        assert actual == expected, (constant, actual, expected)
    nu = reference[:, 0] * 1e9
    x = h * nu / (k * t0)
    g = x**4 * np.exp(-x) / (-np.expm1(-x))**2
    i0 = 2 * (k * t0)**3 / (h * c)**2 * 1e20
    kernel = production_kernels(x)
    theta = 1.245e3 * ev / (me * c * c)
    tsz = sum(kernel[index] * theta**index for index in range(5))
    cib_x = h * nu / (k * 18.8)
    actual = np.column_stack((i0*g, i0*g*(1/2.1923-1/x), i0*g*tsz,
                              cib_x**3.86/np.expm1(cib_x)))
    np.testing.assert_allclose(actual, reference[:, 1:5], rtol=5e-12, atol=1e-22)
    crr = np.loadtxt(directory / "crr_MJy_sr.txt")
    np.testing.assert_allclose(np.interp(reference[:, 0], crr[:, 0], crr[:, 1]), reference[:, 5],
                               rtol=5e-12, atol=1e-22)
    co = np.loadtxt(directory / "co_MJy_sr.txt")
    # Reproduce the production binary-search endpoint-slope policy, including extrapolation.
    upper = np.clip(np.searchsorted(co[:, 0], reference[:, 0], side="right"), 1, len(co)-1)
    lower = upper-1
    weight = np.log(reference[:, 0]/co[lower, 0]) / np.log(co[upper, 0]/co[lower, 0])
    values = np.exp(np.log(co[lower, 1])+weight*np.log(co[upper, 1]/co[lower, 1]))
    np.testing.assert_allclose(values, reference[:, 6], rtol=5e-12, atol=1e-22)
    ame = np.loadtxt(directory / "ame_emissivity_GHz.txt")
    assert np.all(ame > 0) and np.all(np.diff(ame[:, 0]) > 0)
    assert ame[np.argmax(ame[:, 1]), 0] == 29.989554
    print("PASS: production Itoh expressions and exported spectra match the 131-band simulation reference")


def compiled_comparison(directory, require_compiler):
    compiler = shlex.split(os.environ.get("FC", "gfortran"))
    if not compiler or shutil.which(compiler[0]) is None:
        message = "Fortran compiler unavailable; production-module execution was not performed"
        if require_compiler:
            raise RuntimeError(message)
        print("SKIP:", message)
        return
    flags = shlex.split(os.environ.get("FFLAGS", ""))
    with tempfile.TemporaryDirectory(prefix="commander-fossil-sed-") as tmp:
        work = Path(tmp)
        dependency = work / "healpix_types.f90"
        dependency.write_text("module healpix_types\n  implicit none\n"
                              "  integer, parameter :: dp=kind(1.0d0)\nend module healpix_types\n")
        executable = work / "test_fossil_sed"
        command = compiler + flags + ["-std=f2008", "-ffree-line-length-none", "-fcheck=all", "-O0",
                                      "-J", tmp, "-I", tmp, str(dependency), str(SOURCE),
                                      str(ROOT / "tests/test_fossil_sed.f90"), "-o", str(executable)]
        subprocess.run(command, cwd=work, check=True)
        subprocess.run([str(executable), str(directory)], cwd=work, check=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("templates", type=Path, help="Path to templates/commander_spectral")
    parser.add_argument("--require-compiler", action="store_true", help="Fail if the Fortran execution cannot run")
    args = parser.parse_args()
    directory = args.templates.resolve()
    source_comparison(directory)
    compiled_comparison(directory, args.require_compiler)
