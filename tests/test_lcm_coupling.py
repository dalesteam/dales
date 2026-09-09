"""Focused scheduling tests, independent of an installed LCM library.

The small Fortran harness compiles the actual post-dynamics routine with
stubbed physics and timing dependencies. Full DALES runs remain necessary
to validate the attached atmospheric fields and scientific coupling.
"""
from pathlib import Path
import re
import shutil
import subprocess

import pytest

SRC = Path(__file__).resolve().parents[1] / "src"


def routine(name, filename):
    text = (SRC / filename).read_text()
    return re.search(
        rf"^  subroutine {name}\b.*?^  end subroutine {name}\b",
        text, re.M | re.S,
    ).group()


def test_call_order():
    text = (SRC / "program.f90").read_text().lower()
    sequence = ["call tstep_integrate", "call boundary", "call thermodynamics",
                "call lcm_after_dynamics", "call leibniztend",
                "call writerestartfiles"]
    positions = [text.index(s) for s in sequence]
    assert positions == sorted(positions)
    assert text.count("call lcm_after_dynamics") == 1
    assert "call lcm_microphysics" not in routine("microphysics", "modmicrophysics.f90")
    adapter = routine("lcm_microphysics", "modlcm_adapter.f90")
    assert "rk3step" not in adapter
    assert "call lcm_advance(real(rdt, kind=real64))" in adapter


@pytest.mark.parametrize("enabled", [True, False])
def test_post_dynamics_dispatch(tmp_path, enabled):
    compiler = shutil.which("gfortran")
    if compiler is None:
        pytest.skip("gfortran is required for the isolated Fortran dispatch test")
    hook = routine("lcm_after_dynamics", "modmicrophysics.f90")
    source = """
module modglobal
  integer :: rk3step
end module
module harness
  use modglobal
  implicit none
  integer :: imicro, advances=0, failures=0, starts=0, stops=0
  integer, parameter :: imicro_lcm=12
  character(len=*), parameter :: modname='harness'
contains
""" + hook + """
  subroutine lcm_microphysics
    advances=advances+1
  end subroutine
  subroutine timer_tic(name, step)
    character(*) :: name
    integer :: step
    starts=starts+1
  end subroutine
  subroutine timer_toc(name)
    character(*) :: name
    stops=stops+1
  end subroutine
  subroutine finish(name, message)
    character(*) :: name, message
    failures=failures+1
  end subroutine
end module
program test
  use harness
  integer :: step
  ! All non-LCM selections must be no-ops at all stages.
  do imicro=0,11
    do rk3step=1,3
      call lcm_after_dynamics
    end do
  end do
  if (advances /= 0 .or. failures /= 0) stop 1
  imicro=imicro_lcm
  do step=1,2
    do rk3step=1,3
      call lcm_after_dynamics
      if (rk3step < 3) then
#ifdef USE_LCM
        if (advances /= step-1) stop 2
#else
        if (failures /= step-1) stop 2
#endif
      endif
    end do
  end do
#ifdef USE_LCM
  if (advances /= 2 .or. failures /= 0) stop 3
  if (starts /= 2 .or. stops /= 2) stop 4
#else
  if (advances /= 0 .or. failures /= 2) stop 3
#endif
end program
"""
    path = tmp_path / "dispatch.F90"
    path.write_text(source)
    flags = ["-DUSE_LCM"] if enabled else []
    subprocess.run([compiler, "-cpp", "-fcheck=all", *flags, str(path),
                    "-o", "dispatch"], cwd=tmp_path, check=True)
    subprocess.run([str(tmp_path / "dispatch")], cwd=tmp_path, check=True)
