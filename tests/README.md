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

## LCM call scheduling

Run `pytest -q tests/test_lcm_coupling.py` for the isolated call-order and
Fortran dispatch tests. The dispatch tests require `gfortran`, but do not
require a DALES executable or LCM library. They compile the production
post-dynamics hook with stubbed dependencies, with and without `USE_LCM`.

LCM particle transport is dispatched only at stage 3, after `tstep_integrate`,
boundary updates and `thermodynamics`, and before diagnostics/restarts.
The adapter receives the full `rdt`. Other microphysics schemes retain their
within-RK calls. These tests check scheduling, not trajectory accuracy.

This placement currently covers particle transport. Saturation adjustment is
unchanged. Future two-way condensation coupling must define conservative
prognostic-field updates and refresh thermodynamic diagnostics/halos after
LCM feedback, without repeating saturation adjustment.
