#!/usr/bin/env python3
"""Compile actual production index formatting and audit band-tag callers.

Set FC and FFLAGS for the target Fortran compiler. --static-only requires no
compiler and only audits use sites; it does not test the production formatter.
"""
import argparse
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src/commander"


def extract(source, name):
    match = re.search(rf"^  subroutine {name}\b.*?^  end subroutine {name}\b[^\n]*",
                      source, re.M | re.S | re.I)
    if match is None:
        raise RuntimeError(f"Production procedure missing: {name}")
    return match.group()


def static_audit():
    checked = 0
    for path in SRC.glob("*.f90"):
        source = "\n".join(line for line in path.read_text().splitlines()
                           if not line.lstrip().startswith("!"))
        # Direct channel IDs must never pass through a fixed-width tag helper.
        bad = re.search(r"call\s+int2string\s*\(\s*(?:map_id|map_id_in|band|band_reg)\s*,", source, re.I)
        if bad:
            raise AssertionError(f"Fixed-width band identifier remains in {path.name}: {bad.group()}")
        for call in re.finditer(r"call\s+index2string\s*\([^,]+,\s*(\w+)\s*\)", source, re.I):
            variable = call.group(1).lower()
            # Locate the surrounding routine, including a program's main scope.
            starts = list(re.finditer(r"^\s*(?:[\w()=,: ]+\s+)?(?:subroutine|function)\s+\w+\s*\(",
                                      source[:call.start()], re.M | re.I))
            start = starts[-1].start() if starts else 0
            ends = list(re.finditer(r"^\s*end\s+(?:subroutine|function)\b", source[start:call.start()], re.M | re.I))
            if ends:
                # Main program calls can follow earlier contained definitions in
                # some layouts; declaration audit still checks the main prefix.
                start = 0
            prefix = source[start:call.start()]
            declarations = re.findall(r"^\s*character\([^\n]*allocatable[^\n]*::([^\n]*)", prefix, re.M | re.I)
            names = {part.strip().lower() for declaration in declarations for part in declaration.split(",")}
            if variable not in names:
                raise AssertionError(f"index2string output must be deferred-length allocatable: {path.name}:{variable}")
            checked += 1
    print(f"PASS: {checked} index2string call sites and direct band identifiers audited")


def numeric_test():
    source = (SRC / "comm_utils.f90").read_text()
    procedures = extract(source, "int2string") + "\n" + extract(source, "index2string")
    harness = """module formatting
  implicit none
  integer, parameter :: i4b=kind(1)
contains
@PROCEDURES@
end module
program verify
  use formatting
  implicit none
  integer :: i
  character(len=:), allocatable :: tag
  character(len=32) :: expected
  do i=1,1000
    call index2string(i,tag)
    if (i<100) then
      write(expected,'(i2.2)') i
    else
      write(expected,'(i0)') i
    end if
    if (tag/=trim(expected)) stop 1
    if (len(tag)/=len_trim(expected)) stop 2
  end do
  call index2string(131,tag)
  if (tag/='131') stop 3
  call index2string(1,tag)
  if (tag/='01') stop 4
  print *, 'PASS: production formatter preserves 01..99 and formats 100..1000'
end program
""".replace("@PROCEDURES@", procedures)
    compiler = shlex.split(os.environ.get("FC", "gfortran"))
    flags = shlex.split(os.environ.get("FFLAGS", ""))
    with tempfile.TemporaryDirectory(prefix="commander-fossil-index-") as tmp:
        work = Path(tmp)
        (work / "formatting.f90").write_text(harness)
        exe = work / "verify"
        subprocess.run(compiler + flags + ["-std=f2008", "-ffree-line-length-none", "-fcheck=all",
                                          "-O0", "-J", tmp, str(work / "formatting.f90"), "-o", str(exe)],
                       check=True, cwd=tmp)
        subprocess.run([str(exe)], check=True, cwd=tmp)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--static-only", action="store_true")
    args = parser.parse_args()
    static_audit()
    if not args.static_only:
        numeric_test()
