import subprocess
import os

import pytest
import f90nml

@pytest.mark.parametrize(
    "nprocs", [1, 2]
)
def test_run(nprocs, tmp_case):
    """Run DALES, check if it runs succesfully.
    The environment variable DALES has to point to the DALES executable to run.
    """
    dales = os.environ["DALES"]

    namopts = f90nml.read(tmp_case / "namoptions.001")
    namopts["RUN"]["runtime"] = 100
    namopts.write(tmp_case / "namoptions.001", force=True)

    # Run DALES as a subprocess, capturing its output
    result = subprocess.run(
        ["mpirun", "-np", str(nprocs), dales, "namoptions.001"],
        cwd=tmp_case,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True 
    )

    # Print output for debugging purposes
    print(result.stdout)

    # Check if run was successful (sort of)
    assert "Run successful!" in str(result.stdout)
