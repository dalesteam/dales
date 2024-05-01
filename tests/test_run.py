import subprocess
import os

import pytest
import f90nml

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
    dales = os.environ["DALES"]

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
        ["mpirun", "-oversubscribe", "-np", str(nprocs), dales, "namoptions.001"],
        cwd=tmp_case,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True 
    )

    # Print output for debugging purposes
    print(result.stdout)

    # Check if run was successful (sort of)
    assert "Run successful!" in str(result.stdout)
