import subprocess
import os

import pytest

@pytest.mark.parametrize(
    "nprocs", [1, 2, 3, 4]
)
def test_run(nprocs):
    """Run DALES, check if it runs succesfully (exit code equal to zero).
    The environment variable DALES has to point to the DALES executable to run.
    """
    dales = os.environ["DALES"]

    # Run DALES as a subprocess, capturing its output
    result = subprocess.run(
        ["mpirun", "-np", str(nprocs), dales, "namoptions.001"],
        cwd="cases/bomex",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True 
    )

    # Print output for debugging purposes
    print(result.stdout)

    # Check if run was successful (sort of)
    assert "Run successful!" in str(result.stdout)
