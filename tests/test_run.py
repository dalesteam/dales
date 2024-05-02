import subprocess
import os

import pytest
import f90nml

DALES = os.environ["DALES"]

@pytest.mark.parametrize("kmax", [64, 61])
@pytest.mark.parametrize(
    "nprocx,nprocy,itot,jtot", 
    [
        (1, 3, 32, 48),
        (2, 2, 32, 32),
        (1, 1, 31, 31)
    ]
)
def test_domains(
    tmp_case,
    nprocx,
    nprocy,
    itot,
    jtot,
    kmax
):
    """Run DALES for interesting domains, check if it runs succesfully.
    The environment variable DALES has to point to the DALES executable to run.
    """

    # Patch namoptions
    patch = {}

    patch["DOMAIN"] = {
        "itot": itot,
        "jtot": jtot,
        "kmax": kmax
    }

    patch["RUN"] = {
        "runtime": 100,
        "nprocx": nprocx,
        "nprocy": nprocy
    }

    namopts = f90nml.read(tmp_case / "namoptions.001")
    namopts.patch(patch)
    namopts.write(tmp_case / "namoptions.001", force=True)

    nprocs = nprocx * nprocy

    # Run DALES as a subprocess, capturing its output
    result = subprocess.run(
        ["mpirun", "-oversubscribe", "-np", str(nprocs), DALES, "namoptions.001"],
        cwd=tmp_case,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True 
    )

    # Print output for debugging purposes
    print(result.stdout)

    # Check if run was successful (sort of)
    assert "Run successful!" in str(result.stdout)

@pytest.mark.parametrize(
    "mom,tke,thl,qt",
    [
        (2, 2, 2, 2),
        (5, 5, 5, 5),
        (52, 52, 52, 52),
        (6, 6, 6, 6),
        (62, 62, 62, 62),
        (2, 55, 55, 55), # Hybrid scheme
        (2, 555, 555, 555), # Alternative hybrid scheme
        (2, 7, 7, 7), # Kappa scheme
        (2, 2, 1, 1) # Upwinding for thl and qt
    ]
)
def test_advection(tmp_case, mom, tke, thl, qt):
    patch = {}
    
    patch["DYNAMICS"] = {
        "iadv_mom": mom,
        "iadv_tke": tke,
        "iadv_thl": thl,
        "iadv_qt": qt
    }

    patch["DOMAIN"] = {
        "itot": 32,
        "jtot": 32
    }

    patch["RUN"] = {
        "runtime": 100,
        "nprocx": 2,
        "nprocy": 2
    }

    namopts = f90nml.read(tmp_case / "namoptions.001")
    namopts.patch(patch)
    namopts.write(tmp_case / "namoptions.001", force=True)
   
    result = subprocess.run(
        ["mpirun", "-oversubscribe", "-np", "4", DALES, "namoptions.001"],
        cwd=tmp_case,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True 
    )

    print(result.stdout)

    assert "Run successful!" in str(result.stdout)
