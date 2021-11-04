#ifndef INITIALIAZE_VARIABLES_H_
#define INITIALIAZE_VARIABLES_H_

void initialize_variables()
{
    if (DIMENSION == 3 && elast_int == 1)
    {
        printf("Elasticity not yet supported for 3D simulations\n"
        "Setting elasticity off\n");
        elast_int = 0;
    }

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
    }


    double_size = MESH_X * MESH_Y * MESH_Z * sizeof(double);
    complex_size = MESH_X * MESH_Y * MESH_Z * sizeof(cufftDoubleComplex);

    checkCudaErrors(cudaMallocManaged(&kx, MESH_X*sizeof(double)));
    checkCudaErrors(cudaMallocManaged(&ky, MESH_Y*sizeof(double)));
    checkCudaErrors(cudaMallocManaged(&kz, MESH_Z*sizeof(double)));

    checkCudaErrors(cudaMalloc((void**)&gradphix, complex_size));
    checkCudaErrors(cudaMalloc((void**)&gradphiy, complex_size));
    checkCudaErrors(cudaMalloc((void**)&gradphiz, complex_size));

    if (DIMENSION == 3)
        cufftPlan3d(&plan, MESH_X ,MESH_Y, MESH_Z, CUFFT_Z2Z);
    else if (DIMENSION == 2)
        cufftPlan2d(&plan, MESH_X, MESH_Y, CUFFT_Z2Z);

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

    if (DIMENSION == 2)
    {
        MESH_Z     = 1;
        DELTA_Z    = 1.0;
        rows_z     = 1;
        start[Z]   = 0;
        end[Z]     = 0;
        layer_size = rows_y;

        NUM_THREADS_X = 16;
        NUM_THREADS_Y = 16;

        Blocks_X = ceil(MESH_X/NUM_THREADS_X);
        Blocks_Y = ceil(MESH_Y/NUM_THREADS_Y);

        Gridsize  = dim3(Blocks_X, Blocks_Y);
        Blocksize = dim3(NUM_THREADS_X, NUM_THREADS_Y);
    }
    else
    {
        NUM_THREADS_X = 8;
        NUM_THREADS_Y = 8;
        NUM_THREADS_Z = 8;

        Blocks_X = ceil(MESH_X/NUM_THREADS_X);
        Blocks_Y = ceil(MESH_Y/NUM_THREADS_Y);
        Blocks_Z = ceil(MESH_Z/NUM_THREADS_Z);

        Gridsize = dim3(Blocks_X, Blocks_Y, Blocks_Z);
        Blocksize = dim3(NUM_THREADS_X, NUM_THREADS_Y, NUM_THREADS_Z);
    }

    printf("\nGrid Dimensions:  (%d, %d, %d)\n", Blocks_X, Blocks_Y, Blocks_Z);
    printf("Block Dimensions: (%d, %d, %d)\n", NUM_THREADS_X, NUM_THREADS_Y, NUM_THREADS_Z);

    index_count = layer_size*rows_x;

    gridinfo = (struct fields*)malloc(index_count*sizeof(*gridinfo));

    for (index = 0; index < index_count; index++)
    {
        gridinfo[index].phia = (cufftDoubleComplex*)malloc(NUMPHASES*sizeof(cufftDoubleComplex));
        gridinfo[index].comp = (cufftDoubleComplex*)malloc((NUMCOMPONENTS-1)*sizeof(cufftDoubleComplex));
    }

    gridinfo_size = index_count * sizeof(cufftDoubleComplex) * (NUMPHASES + NUMCOMPONENTS - 1);

    checkCudaErrors(cudaMallocManaged(&CompBuff, sizeof(*CompBuff)*(NUMCOMPONENTS-1)));
    checkCudaErrors(cudaMallocManaged(&PhiBuff, sizeof(*PhiBuff)*(NUMPHASES)));
    checkCudaErrors(cudaMallocManaged(&dfdc, sizeof(*dfdc)*(NUMCOMPONENTS-1)));
    checkCudaErrors(cudaMallocManaged(&dfdphi, sizeof(*dfdphi)*(NUMPHASES)));

    for (int j = 0; j < NUMCOMPONENTS-1; j++)
        checkCudaErrors(cudaMallocManaged(&CompBuff[j], complex_size));

    for (int j = 0; j < NUMPHASES; j++)
        checkCudaErrors(cudaMallocManaged(&PhiBuff[j], complex_size));

    for (int j = 0; j < NUMCOMPONENTS-1; j++)
        checkCudaErrors(cudaMallocManaged(&dfdc[j], complex_size));

    for (int j = 0; j < NUMPHASES; j++)
        checkCudaErrors(cudaMallocManaged(&dfdphi[j], complex_size));


    checkCudaErrors(cudaMallocManaged(&B, sizeof(double)*MESH_X*MESH_Y));
}

#endif
