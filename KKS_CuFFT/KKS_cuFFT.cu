#include <cuda_runtime.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <fenv.h>
#include <cuda.h>
#include <cuComplex.h>
#include <cufft.h>
#include <cublas_v2.h>
#include <math_constants.h>
#include "cub/cub3/cub.cuh"
#include <unistd.h>

#include <arpa/inet.h>
#include <endian.h>
#include <stdint.h>

#define  IS_BIG_ENDIAN     (1 == htons(1))
#define  IS_LITTLE_ENDIAN  (!IS_BIG_ENDIAN)

#include "functions/global_vars.h"
#include "functions/functions.h"
#include "functions/helper_string.h"
#include "functions/helper_cuda.h"
#include "solverloop/kernels.cuh"
#include "solverloop/kernels_typeA.cuh"
#include "solverloop/kernels_typeB.cuh"
#include "functions/initialize_variables.h"
#include "functions/utility.cuh"
#include "functions/cuda_filling.cuh"
#include "functions/reading_input_parameters.h"
#include "functions/fill_domain.h"
#include "functions/read_restart.h"
#include "functions/calc_bn.h"
#include "solverloop/file_writer.h"
#include "solverloop/evolve_typeA.h"
#include "solverloop/evolve_typeB.h"

int main (int argc, char *argv[])
{
    // Check to ensure the program is run with all the necessary info specified
    if (argc < 4)
    {
        printf("Insufficient number of arguments in execution command\n"
        "Using default arguments - Input.in, Filling.in and Output\n");

        const int argc_const = 4;
        argc = argc_const;

        char *argv_temp[argc_const];

        for (int i = 0; i < argc; i++)
            argv_temp[i] = (char*)malloc(20*sizeof(char));

        strcpy(argv_temp[0], "kks.out");
        strcpy(argv_temp[1], "Input.in");
        strcpy(argv_temp[2], "Filling.in");
        strcpy(argv_temp[3], "Output");

        for (int i = 0; i < argc; i++)
            argv[i] = argv_temp[i];
    }

    mkdir("DATA", 0777);

    // Reads all the input data from the specified input file
    reading_input_parameters(argv);

    printf("Periodic boundary conditions will be applied for all variables at all boundaries\n");
    // Uses the information read by reading_input_parameters to allocate memory
    // for variables that are necessary for the simulation
    initialize_variables();

    // Sets the initial phase-fields using the geometric data
    // specified in the filling file
    if (restart == 0)
        fill_domain(argv);
    else if (restart == 1)
        read_domain(argv, initcount);

    // Executes the solution procedure
    evolve(argv);

    cudaDeviceReset();
    return EXIT_SUCCESS;
}
