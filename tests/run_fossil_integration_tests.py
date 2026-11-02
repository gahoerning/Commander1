#!/usr/bin/env python3
"""Compile production residual/CG procedures with small, explicit test doubles.

These fixtures validate arithmetic and dispatch dependencies without HEALPix or
MPI. Real multi-rank execution remains a separate remote acceptance gate.
"""
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src/commander"


def procedure(filename, name):
    text = (SRC / filename).read_text()
    match = re.search(
        rf"^  *(?:real\(dp\) +)?(subroutine|function) +{name}\b.*?^  *end +\1(?: +{name})?[^\n]*",
        text, re.M | re.S | re.I,
    )
    if not match:
        raise RuntimeError(f"Cannot extract {name} from {filename}")
    return match.group()


def main():
    procedures = [procedure("comm_data_mod.f90", name) for name in [
        "compute_residual", "compute_signal_map", "compute_lowres_residual", "sample_fossil_global_block",
    ]]
    procedures += [procedure("comm_cgd_mod.f90", "draw_constrained_realization_by_CG")]
    procedures += [procedure("comm_genvec_mod.f90", name) for name in [
        "allocate_genvec", "deallocate_genvec", "genvec_set_equal", "genvec_minus", "genvec_plus_times",
        "genvec_dot_product",
    ]]
    procedures += [procedure("comm_utils.f90", name) for name in [
        "int2string", "index2string", "parse_parameter", "get_token", "get_tokens", "tokenize",
    ]]
    fixture = (ROOT / "tests/fossil_integration_support.f90.in").read_text()
    fixture = fixture.replace("@PROCEDURES@", "\n\n".join(procedures))
    compiler = shlex.split(os.environ.get("FC", "gfortran"))
    flags = shlex.split(os.environ.get("FFLAGS", ""))
    with tempfile.TemporaryDirectory(prefix="commander-fossil-integration-") as tmp:
        path = Path(tmp)
        (path / "fixture.f90").write_text(fixture)
        command = compiler + flags + ["-std=f2008", "-ffree-line-length-none", "-fcheck=all", "-O0",
                                      "-J", tmp, "-I", tmp, str(path / "fixture.f90"), "-o", str(path / "test")]
        subprocess.run(command, check=True, cwd=tmp)
        subprocess.run([str(path / "test")], check=True, cwd=tmp)


if __name__ == "__main__":
    main()
