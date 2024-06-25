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

def main():
    with open("./kernels/advection_2nd.f90", "r") as file:
        source = file.read()

    app = Code(OpenACC(), Fortran())
    preprocessor = extract_preprocessor(source)
    signature = extract_directive_signature(source, app)
    body = extract_directive_code(source, app)
    data = extract_directive_data(source, app)
    init = extract_initialization_code(source, app)
    deinit = extract_deinitialization_code(source, app)

    for func in signature.keys():
        args = allocate_signature_memory(data[func], preprocessor)

        kernel_string = generate_directive_function(
            preprocessor,
            signature[func],
            body[func],
            app,
            initialization=init,
            deinitialization=deinit
        )

        collapse_factor = 2 if func == "advecu_2nd_vert_bot" else 3

        tune_params = {}
        tune_params["vlength"] = [32 * i for i in range(1, 33)]
        tune_params["ngangs"] = [2 ** i for i in range(1, 8)]
        tune_params["ncollapse"] = [i for i in range(1, collapse_factor + 1)]

        metrics = {}

        restrictions = []
        
        tune_kernel(
            func,
            kernel_string,
            0,
            args,
            tune_params,
            metrics=metrics,
            restrictions=restrictions,
            compiler="nvfortran",
            compiler_options=["-acc=gpu", "-gpu=cc80"]
        )

if __name__ == "__main__":
    main()
