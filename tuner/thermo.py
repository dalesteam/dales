#!/usr/bin/env python
from kernel_tuner import tune_kernel
from kernel_tuner.utils.directives import (
        Code,
        OpenACC,
        Fortran,
        extract_directive_signature,
        extract_directive_code,
        extract_preprocessor,
        generate_directive_function,
        extract_directive_data,
        allocate_signature_memory,
        extract_initialization_code,
        extract_deinitialization_code
)

with open("thermo.f90", "r") as file:
    code = file.read()

app = Code(OpenACC(), Fortran())
preprocessor = extract_preprocessor(code)
signature = extract_directive_signature(code, app)
body = extract_directive_code(code, app)
data = extract_directive_data(code, app)
init = extract_initialization_code(code, app)
deinit = extract_deinitialization_code(code, app)

args = allocate_signature_memory(data["icethermo"], preprocessor)
kernel_string = generate_directive_function(
        preprocessor,
        signature["icethermo"],
        body["icethermo"],
        app,
        data=data["icethermo"],
        initialization=init,
        deinitialization=deinit,
)

tune_params = {}
tune_params["nthreads"] = [32 * i for i in range(1, 33)]
tune_params["ncoll"] = [1, 2, 3]

metrics = {}
tune_kernel(
    "icethermo",
    kernel_string,
    0,
    args,
    tune_params,
    metrics=metrics,
    iterations=20,
    compiler_options=[ "-Mr8", "-O3", "-acc=gpu", "-gpu=cc80"],
    compiler="nvfortran"
)
print(metrics)
