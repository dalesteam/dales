import shutil
from pathlib import Path

import pytest

CASES = Path("cases")


def pytest_addoption(parser):
    parser.addoption("--case", action="store")


@pytest.fixture()
def tmp_case(tmp_path, request):
    """Setup a temporary directory for the selected case"""
    # Source dir for input files
    case_path = CASES / request.config.getoption("--case")
    
    # Copy input files
    shutil.copytree(
        case_path,
        tmp_path,
        ignore=shutil.ignore_patterns("results*"), # Don't copy existing results
        dirs_exist_ok=True
    )

    return tmp_path
