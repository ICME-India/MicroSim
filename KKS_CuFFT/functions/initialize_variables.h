#ifndef INITIALIAZE_VARIABLES_H_
#define INITIALIAZE_VARIABLES_H_

void initialize_variables()
{
    // Check if homogenous elastic constants are well-defined
    // Run calculations to obtain elastic free-energy contribution
    if (elast_int == 1)
    {
        for (int i = 0; i < NUMPHASES; i++)
        {
            if (Stiffness_c[i].C11 < 0)
            {
                printf("\nStiffness modulus C11 for phase %s can not be negative. Taking absolute value\n", PHASES[i]);
                Stiffness_c[i].C11 *= -1;
            }
            if (Stiffness_c[i].C12 < 0)
            {
                printf("\nStiffness modulus C12 for phase %s can not be negative. Taking absolute value\n", PHASES[i]);
                Stiffness_c[i].C12 *= -1;
            }
            if (Stiffness_c[i].C44 < 0)
            {
                printf("\nStiffness modulus C44 for phase %s can not be negative. Taking absolute value\n", PHASES[i]);
                Stiffness_c[i].C44 *= -1;
            }
        }

        if (Stiffness_c[0].C11 != Stiffness_c[1].C11 || Stiffness_c[0].C12 != Stiffness_c[1].C12 || Stiffness_c[0].C44 != Stiffness_c[1].C44)
        {
            printf("\nAssuming homogenous elastic moduli, taking average of given input.\n");
            Stiffness_c[1].C11 = (Stiffness_c[0].C11 + Stiffness_c[1].C11) / 2;
            Stiffness_c[1].C12 = (Stiffness_c[0].C12 + Stiffness_c[1].C12) / 2;
            Stiffness_c[1].C44 = (Stiffness_c[0].C44 + Stiffness_c[1].C44) / 2;
        }

        double *B_h = (double*)malloc(sizeof(double)*MESH_X*MESH_Y*MESH_Z);
        checkCudaErrors(cudaMalloc((void**)&B, sizeof(double)*MESH_X*MESH_Y*MESH_Z));

        printf("Calculating eigenstrain\n");
        calc_strn(eigen_strn1, eigen_strain[1].xx, eigen_strain[1].yy, eigen_strain[1].zz,
                  eigen_strain[1].xy, eigen_strain[1].xz, eigen_strain[1].yz);
        calculate_Bn(B_h, eigen_strn1, eigen_strn1, stress1, stress1);
        printf("Done calculating eigenstrain\n");
        checkCudaErrors(cudaMemcpy(B, B_h, double_size, cudaMemcpyHostToDevice));
        free(B_h);
    }

    long index_count;
    long index;

    start[X]   = 0;
    start[Y]   = 0;
    start[Z]   = 0;

    rows_x     = MESH_X + 0;
    rows_y     = MESH_Y + 0;
    rows_z     = MESH_Z + 0;
    end[X]     = rows_x - 1;
    end[Y]     = rows_y - 1;
    end[Z]     = rows_z - 1;

    layer_size = rows_y*rows_z;

    // Setting up for mGPU data organization

    if (DIMENSION == 2)
    {
        checkCudaErrors(cufftPlan2d(&plan, MESH_X, MESH_Y, CUFFT_Z2Z));

        // Flatten 3rd dimension as a precaution
        MESH_Z     = 1;
        DELTA_Z    = 1.0;
        rows_z     = 1;
        start[Z]   = 0;
        end[Z]     = 0;
        layer_size = rows_y;

        NUM_THREADS_X = 16;
        NUM_THREADS_Y = 16;
        NUM_THREADS_Z = 1;

        Blocks_X = ceil(MESH_X/NUM_THREADS_X);
        Blocks_Y = ceil(MESH_Y/NUM_THREADS_Y);
        Blocks_Z = 1;

        Gridsize  = dim3(Blocks_X, Blocks_Y, Blocks_Z);
        Blocksize = dim3(NUM_THREADS_X, NUM_THREADS_Y, NUM_THREADS_Z);
    }
    else
    {
        checkCudaErrors(cufftPlan3d(&plan, MESH_X ,MESH_Y, MESH_Z, CUFFT_Z2Z));

        NUM_THREADS_X = 2;
        NUM_THREADS_Y = 2;
        NUM_THREADS_Z = 64;

        Blocks_X = ceil(MESH_X/NUM_THREADS_X);
        Blocks_Y = ceil(MESH_Y/NUM_THREADS_Y);
        Blocks_Z = ceil(MESH_Z/NUM_THREADS_Z);

        Gridsize = dim3(Blocks_X, Blocks_Y, Blocks_Z);
        Blocksize = dim3(NUM_THREADS_X, NUM_THREADS_Y, NUM_THREADS_Z);
    }

    double_size = MESH_X * MESH_Y * MESH_Z * sizeof(double);
    complex_size = MESH_X * MESH_Y * MESH_Z * sizeof(cufftDoubleComplex);

    // Memory for Fourier differentiation matrices on host and device
    checkCudaErrors(cudaMalloc((void**)&kx, MESH_X*sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&ky, MESH_Y*sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&kz, MESH_Z*sizeof(double)));

    kx_h = (double*)malloc(MESH_X*sizeof(double));
    ky_h = (double*)malloc(MESH_Y*sizeof(double));
    kz_h = (double*)malloc(MESH_Z*sizeof(double));

    printf("\nGrid Dimensions:  (%d, %d, %d)\n", Blocks_X, Blocks_Y, Blocks_Z);
    printf("Block Dimensions: (%d, %d, %d)\n", NUM_THREADS_X, NUM_THREADS_Y, NUM_THREADS_Z);

    index_count = layer_size*rows_x;

    // gridinfo stores fields data on the CPU
    gridinfo = (struct fields*)malloc(index_count*sizeof(*gridinfo));

    for (index = 0; index < index_count; index++)
    {
        gridinfo[index].phia = (double*)malloc(NUMPHASES*sizeof(double));
        gridinfo[index].compi = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));
    }

    // Allocating memory on the GPU for all the fields
    checkCudaErrors(cudaMalloc((void**)&compDev, sizeof(cufftDoubleComplex*)*(NUMCOMPONENTS-1)));
    checkCudaErrors(cudaMalloc((void**)&dfdcDev, sizeof(cufftDoubleComplex*)*(NUMCOMPONENTS-1)));
    checkCudaErrors(cudaMalloc((void**)&phiDev, sizeof(cufftDoubleComplex*)*(NUMPHASES)));
    checkCudaErrors(cudaMalloc((void**)&dfdphiDev, sizeof(cufftDoubleComplex*)*(NUMPHASES)));

    compHost = (cufftDoubleComplex**)malloc((NUMCOMPONENTS-1)*sizeof(cufftDoubleComplex*));
    dfdcHost = (cufftDoubleComplex**)malloc((NUMCOMPONENTS-1)*sizeof(cufftDoubleComplex*));
    phiHost = (cufftDoubleComplex**)malloc(NUMPHASES*sizeof(cufftDoubleComplex*));
    dfdphiHost = (cufftDoubleComplex**)malloc(NUMPHASES*sizeof(cufftDoubleComplex*));

    for (int j = 0; j < NUMCOMPONENTS-1; j++)
    {
        checkCudaErrors(cudaMalloc((void**)&compHost[j], sizeof(cufftDoubleComplex)*MESH_X*MESH_Y*MESH_Z));
        checkCudaErrors(cudaMemcpy(&compDev[j], &compHost[j], sizeof(cufftDoubleComplex*), cudaMemcpyHostToDevice));
    }

    for (int j = 0; j < NUMCOMPONENTS-1; j++)
    {
        checkCudaErrors(cudaMalloc((void**)&dfdcHost[j], sizeof(cufftDoubleComplex)*MESH_X*MESH_Y*MESH_Z));
        checkCudaErrors(cudaMemcpy(&dfdcDev[j], &dfdcHost[j], sizeof(cufftDoubleComplex*), cudaMemcpyHostToDevice));
    }

    for (int j = 1; j < NUMPHASES; j++)
    {
        checkCudaErrors(cudaMalloc((void**)&phiHost[j], sizeof(cufftDoubleComplex)*MESH_X*MESH_Y*MESH_Z));
        checkCudaErrors(cudaMemcpy(&phiDev[j], &phiHost[j], sizeof(cufftDoubleComplex*), cudaMemcpyHostToDevice));
    }

    for (int j = 1; j < NUMPHASES; j++)
    {
        checkCudaErrors(cudaMalloc((void**)&dfdphiHost[j], sizeof(cufftDoubleComplex)*MESH_X*MESH_Y*MESH_Z));
        checkCudaErrors(cudaMemcpy(&dfdphiDev[j], &dfdphiHost[j], sizeof(cufftDoubleComplex*), cudaMemcpyHostToDevice));
    }

    // Setting the solver type as read from Input.in

    solvertype = 1;

    if (solvertype == 0)
        evolve = evolve_FFT;
    else if (solvertype == 1)
        evolve = evolve_FD;
}

#endif
