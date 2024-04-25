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
    result = subprocess.run(["mpirun", "-np", str(nprocs), dales, "namoptions.001"], cwd="cases/bomex")
    assert result.returncode == 0
