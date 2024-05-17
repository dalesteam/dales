# Running tests

DALES can be tested automatically using [Pytest](https://docs.pytest.org/en/latest/). The tests that are performed can be seen as [smoke tests](https://en.wikipedia.org/wiki/Smoke_testing_(software)) and only test for catastrophic failure (reading/writing out of array bounds, over/underflow, et cetera). The output of the model is not tested for correctness, and should therefore be done manually.

Tests are performed in CI, but can also be done locally. To do so:

1. Make and activate a Python virtual environment (optional, but recommended):

```{bash}
python -m venv .venv
source .venv/bin/activate
```

2. Install requirements:

```{bash}
python -m pip install -r requirements.txt
```

3. Use an environment variable to point to the DALES executable:

```{bash}
export DALES=$(pwd)/build/src/dales4.4
```

4. Run the tests:

```{bash}
pytest -rf --assert=plain --case=bomex
```

Any of the cases in the `cases/` directory can be tested by changing the `--case` argument. Keep in mind that some cases require a custom `moduser.f90`.
