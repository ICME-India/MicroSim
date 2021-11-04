void Evolve2D(char *argv[])
{
    long x, y, z, index;
    long b, a;

    int     loop_condition, count = initcount;
    double  *tempComp, *tempPhi;
    double  *maxerrComp, *maxerrPhi;
    double  *maxVal, *minVal;
    double  f0AVminv, f0BVminv;
    void    *t_storage = NULL;
    size_t  t_storage_bytes = 0;

    kappa_phi = 3.0/(2.0*alpha) * (Gamma[0][1]*Ln);
    w = 6.0 * alpha * Gamma[0][1] / Ln;

    calc_strn(eigen_strn1, eigen_strain[1].xx, eigen_strain[1].yy, eigen_strain[1].zz,
              eigen_strain[1].xy, eigen_strain[1].xz, eigen_strain[1].yz);
    calculate_Bn(B, eigen_strn1, eigen_strn1, stress1, stress1);

    f0AVminv = F0[0][0][0] * (1.0/Vm);
    f0BVminv = F0[1][0][0] * (1.0/Vm);

    checkCudaErrors(cudaMallocManaged(&tempComp, double_size));
    checkCudaErrors(cudaMallocManaged(&tempPhi, double_size));
    checkCudaErrors(cudaMallocManaged(&maxerrComp, sizeof(double)));
    checkCudaErrors(cudaMallocManaged(&maxerrPhi, sizeof(double)));
    checkCudaErrors(cudaMallocManaged(&maxVal, sizeof(double)));
    checkCudaErrors(cudaMallocManaged(&minVal, sizeof(double)));

    *maxerrComp = 1;
    *maxerrPhi  = 1;

    cub::DeviceReduce::Max(t_storage, t_storage_bytes, tempComp, maxerrComp, MESH_X * MESH_Y);
    cub::DeviceReduce::Sum(t_storage, t_storage_bytes, tempComp, maxerrComp, MESH_X * MESH_Y);

    checkCudaErrors(cudaMallocManaged(&t_storage, t_storage_bytes));

    if (UMFlag)
    {
        cudaMemAdvise(CompBuff[0], complex_size, cudaMemAdviseSetAccessedBy, 0);
        cudaMemPrefetchAsync(CompBuff[0], complex_size, 0);
        cudaMemPrefetchAsync(tempComp, double_size, 0);
    }

    //save real part
    SaveReal_2D<<< Gridsize, Blocksize >>>(tempComp, CompBuff[0], MESH_X, MESH_Y);

    if(cufftExecZ2Z(plan, PhiBuff[1], PhiBuff[1], CUFFT_FORWARD)!= CUFFT_SUCCESS)
        printf("forward fft of phi failed\n");
    cudaDeviceSynchronize();

    loop_condition = 1;

    dkx = 2.0 * CUDART_PI / ((double) MESH_X * DELTA_X);
    dky = 2.0 * CUDART_PI / ((double) MESH_Y * DELTA_Y);

    sizescale = 1.0/(double)(MESH_X * MESH_Y);

    for (int i = 0; i < MESH_X; i++)
    {
        if (i < MESH_X/2)
            kx[i] = (double)i * dkx;
        else
            kx[i] = (double)(i - MESH_X) * dkx;
    }

    for (int j = 0; j < MESH_Y; j++)
    {
        if (j < MESH_Y/2)
            ky[j] = (double)j * dky;
        else
            ky[j] = (double)(j - MESH_Y) * dky;
    }

    if (UMFlag)
    {
        cudaMemPrefetchAsync(dfdc, complex_size, 0, NULL);
        cudaMemPrefetchAsync(dfdphi, complex_size, 0, NULL);
        cudaMemPrefetchAsync(gradphix, complex_size, 0, NULL);
        cudaMemPrefetchAsync(gradphiy, complex_size, 0, NULL);
        cudaMemPrefetchAsync(kx, MESH_X*sizeof(double), 0, NULL);
        cudaMemPrefetchAsync(ky, MESH_Y*sizeof(double), 0, NULL);
        cudaMemPrefetchAsync(tempComp, double_size, 0, NULL);
        cudaMemPrefetchAsync(tempPhi, double_size, 0, NULL);
    }

    //Time loop
    while (count <= numsteps && loop_condition == 1)
    {
        if ((count % saveT)==0 || (count == numsteps) ||  (loop_condition == 0))
        {
            if(sim_time == 0)
            {
                printf("\nTime = %6lf\n", sim_time);
                printf("Writing configuration to file!\n");
            }

            if (UMFlag)
            {
                for (int i = 0; i < NUMCOMPONENTS-1; i++)
                    cudaMemPrefetchAsync(CompBuff[i], complex_size, cudaCpuDeviceId);
                cudaMemPrefetchAsync(dfdphi, complex_size, cudaCpuDeviceId);
            }

            for(x = 0; x < rows_x; x++)
            {
                for(y = 0; y < rows_y; y++)
                {
                    for (z = 0; z < rows_z; z++)
                    {
                        index = x*layer_size + y*rows_z + z;

                        for(a = 0; a < NUMCOMPONENTS-1; a++)
                            gridinfo[index].comp[a].x = CompBuff[a][index].x;

                        gridinfo[index].phia[0].x = 1.0;

                        for (b = 1; b < NUMPHASES; b++)
                        {
                            gridinfo[index].phia[b].x = dfdphi[b][index].x;
                            gridinfo[index].phia[0].x -= gridinfo[index].phia[b].x;
                        }
                    }
                }
            }
            if (ASCII == 0)
                writetofile_serial2D_binary(gridinfo, argv, count);
            else
                writetofile_serial2D(gridinfo, argv, count);

            if (UMFlag)
            {
                for (int i = 0; i < NUMCOMPONENTS-1; i++)
                    cudaMemPrefetchAsync(CompBuff[i], complex_size, 0);
                for (int i = 0; i < NUMPHASES; i++)
                    cudaMemPrefetchAsync(dfdphi[i], complex_size, 0);
            }
        }

        //gradient of phi
        ComputeGradphi_2D<<< Gridsize, Blocksize >>>(kx, ky, MESH_X, MESH_Y, PhiBuff[1], gradphix, gradphiy);

        if (cufftExecZ2Z(plan, gradphix, gradphix, CUFFT_INVERSE) != CUFFT_SUCCESS)
            printf("inverse fft of gradphix failed\n");
        if (cufftExecZ2Z(plan, gradphiy, gradphiy, CUFFT_INVERSE) != CUFFT_SUCCESS)
            printf("inverse fft of gradphiy failed\n");

        Normalize_2D<<< Gridsize, Blocksize >>>(gradphix, sizescale, MESH_X, MESH_Y);
        Normalize_2D<<< Gridsize, Blocksize >>>(gradphiy, sizescale, MESH_X, MESH_Y);

        //Driving forces: varmob (saved in gradphi) and dfdphi
        ComputeDrivForce_2D<<< Gridsize, Blocksize >>>(CompBuff[0], dfdphi[1], gradphix,
                                                       gradphiy,f0AVminv,f0BVminv,
                                                       ceq[1][1][0], ceq[0][0][0], Diffusivity[0][0][0], w, MESH_X, MESH_Y);

        if (cufftExecZ2Z(plan,gradphix,gradphix,CUFFT_FORWARD) != CUFFT_SUCCESS)
            printf("forward fft of gradphix failed\n");
        if (cufftExecZ2Z(plan,gradphiy,gradphiy,CUFFT_FORWARD) != CUFFT_SUCCESS)
            printf("forward fft of gradphiy failed\n");

        //dfdc
        ComputeDfdc_2D<<< Gridsize, Blocksize >>>(dfdc[0], gradphix, gradphiy,
                                                  MESH_X, MESH_Y, kx, ky);

        if (cufftExecZ2Z(plan, CompBuff[0], CompBuff[0], CUFFT_FORWARD) != CUFFT_SUCCESS)
            printf("forward fft of comp failed\n");;
        if (cufftExecZ2Z(plan, dfdphi[1], dfdphi[1], CUFFT_FORWARD) != CUFFT_SUCCESS)
            printf("forward fft of dfphi failed\n");

        //evolve comp and phi
        Update_comp_phi_2D<<< Gridsize, Blocksize >>>(CompBuff[0], dfdc[0], PhiBuff[1], dfdphi[1],
                                                      kx, ky, DELTA_t, Diffusivity[0][0][0], kappa_phi,
                                                      relax_coeff, MESH_X, MESH_Y, elast_int, B);

        cufftExecZ2Z(plan, CompBuff[0], CompBuff[0], CUFFT_INVERSE);
        cufftExecZ2Z(plan, dfdphi[1], dfdphi[1], CUFFT_INVERSE);

        Normalize_2D<<< Gridsize, Blocksize >>>(CompBuff[0], sizescale, MESH_X, MESH_Y);
        Normalize_2D<<< Gridsize, Blocksize >>>(dfdphi[1], sizescale, MESH_X, MESH_Y);

        Find_err_matrix_2D<<< Gridsize, Blocksize >>>(tempComp, CompBuff[0], MESH_X, MESH_Y);

        //cudaDeviceSynchronize();
        cudaGetLastError();

        if (UMFlag)
            cudaMemAdvise(maxerrComp, sizeof(double), cudaMemAdviseSetAccessedBy, 0);

        cub::DeviceReduce::Max(t_storage, t_storage_bytes, tempComp, maxerrComp, MESH_X*MESH_Y);

        if (*maxerrComp <= COMPERR)
        {
            printf("!!!! CONVERGENCE ACHIEVED !!!!.\n");
            printf("Maxerror = %le\n", *maxerrComp);
            loop_condition = 0;
        }

        sim_time = sim_time + DELTA_t;
        count++;

        SaveReal_2D<<< Gridsize, Blocksize >>>(tempComp, CompBuff[0], MESH_X, MESH_Y);

        if (count % time_output == 0)
        {
            printf("\nTime=%6lf\n", sim_time);

            Find_err_matrix_2D<<< Gridsize, Blocksize >>>(tempPhi, dfdphi[1], MESH_X, MESH_Y);
            cub::DeviceReduce::Max(t_storage, t_storage_bytes, tempPhi, maxerrPhi, MESH_X*MESH_Y);
            SaveReal_2D<<< Gridsize, Blocksize >>>(tempPhi, dfdphi[1], MESH_X, MESH_Y);
            cub::DeviceReduce::Max(t_storage, t_storage_bytes, tempPhi, maxVal, MESH_X*MESH_Y);
            cub::DeviceReduce::Min(t_storage, t_storage_bytes, tempPhi, minVal, MESH_X*MESH_Y);
            cudaDeviceSynchronize();

            printf("%*s, Max = %lf, Min = %lf, Max. error = %lf\n", 5, PHASES[1], *maxVal, *minVal, *maxerrPhi);

            cub::DeviceReduce::Max(t_storage, t_storage_bytes, tempComp, maxVal, MESH_X*MESH_Y);
            cub::DeviceReduce::Min(t_storage, t_storage_bytes, tempComp, minVal, MESH_X*MESH_Y);
            cudaDeviceSynchronize();

            printf("%*s, Max = %lf, Min = %lf, Max. error = %lf\n", 5, COMPONENTS[0], *maxVal, *minVal, *maxerrComp);

            if (count % saveT == 0)
                printf("Writing configuration to file!\n");
        }
        else if (count % time_output == time_output-1)
            SaveReal_2D<<< Gridsize, Blocksize >>>(tempPhi, dfdphi[1], MESH_X, MESH_Y);

    }//time loop ends

    cudaFree(tempComp);
    cudaFree(tempPhi);
    cudaFree(maxerrComp);
    cudaFree(maxerrPhi);
    cudaFree(maxVal);
    cudaFree(minVal);
    cudaFree(t_storage);

    for (int j = 0; j < NUMCOMPONENTS-1; j++)
        cudaFree(CompBuff[j]);
    for(int j = 0; j < NUMPHASES; j++)
        cudaFree(PhiBuff[j]);

    for (int j = 0; j < NUMCOMPONENTS-1; j++)
        cudaFree(dfdc[j]);
    for(int j = 0; j < NUMPHASES; j++)
        cudaFree(dfdphi[j]);

    cudaFree(CompBuff);
    cudaFree(PhiBuff);
    cudaFree(dfdc);
    cudaFree(dfdphi);

    cudaFree(gradphix);
    cudaFree(gradphiy);
    cudaFree(gradphiz);

    cudaFree(kx);
    cudaFree(ky);
    cudaFree(kz);

    cudaFree(B);

    for (int i = 0; i < rows_x*rows_y*rows_z; i++)
    {
        free(gridinfo[i].phia);
        free(gridinfo[i].comp);
    }
    free(gridinfo);
}

void Evolve3D(char *argv[])
{
    long x, y, z, index;
    long b, a;

    int     loop_condition, count = initcount;
    double  *tempComp, *tempPhi;
    double  *maxerrComp, *maxerrPhi;
    double  *maxVal, *minVal;
    double  f0AVminv, f0BVminv;
    void    *t_storage = NULL;
    size_t  t_storage_bytes = 0;

    kappa_phi = 3.0/(2.0*alpha) * (Gamma[0][1]*Ln);
    w = 6.0 * alpha * Gamma[0][1] / Ln;

    f0AVminv = F0[0][0][0] * (1.0/Vm);
    f0BVminv = F0[1][0][0] * (1.0/Vm);

    checkCudaErrors(cudaMallocManaged(&tempComp, double_size));
    checkCudaErrors(cudaMallocManaged(&tempPhi, double_size));
    checkCudaErrors(cudaMallocManaged(&maxerrComp, sizeof(double)));
    checkCudaErrors(cudaMallocManaged(&maxerrPhi, sizeof(double)));
    checkCudaErrors(cudaMallocManaged(&maxVal, sizeof(double)));
    checkCudaErrors(cudaMallocManaged(&minVal, sizeof(double)));

    *maxerrComp = 1;
    *maxerrPhi  = 1;

    cub::DeviceReduce::Max(t_storage, t_storage_bytes, tempComp, maxerrComp, MESH_X * MESH_Y * MESH_Z);
    cub::DeviceReduce::Sum(t_storage, t_storage_bytes, tempComp, maxerrComp, MESH_X * MESH_Y * MESH_Z);

    checkCudaErrors(cudaMallocManaged(&t_storage, t_storage_bytes));

    if (UMFlag)
    {
        cudaMemAdvise(CompBuff[0], complex_size, cudaMemAdviseSetAccessedBy, 0);
        cudaMemPrefetchAsync(CompBuff[0], complex_size, 0);
        cudaMemPrefetchAsync(tempComp, double_size, 0);
    }

    //save real part
    SaveReal_3D<<< Gridsize, Blocksize >>>(tempComp, CompBuff[0], MESH_X, MESH_Y, MESH_Z);

    if(cufftExecZ2Z(plan, PhiBuff[1], PhiBuff[1], CUFFT_FORWARD)!= CUFFT_SUCCESS)
        printf("forward fft of phi failed\n");
    cudaDeviceSynchronize();

    loop_condition = 1;

    dkx = 2.0 * CUDART_PI / ((double) MESH_X * DELTA_X);
    dky = 2.0 * CUDART_PI / ((double) MESH_Y * DELTA_Y);
    dkz = 2.0 * CUDART_PI / ((double) MESH_Z * DELTA_Z);

    sizescale = 1.0/(double)(MESH_X * MESH_Y * MESH_Z);

    for (int i = 0; i < MESH_X; i++)
    {
        if (i < MESH_X/2)
            kx[i] = (double)i * dkx;
        else
            kx[i] = (double)(i - MESH_X) * dkx;
    }

    for (int j = 0; j < MESH_Y; j++)
    {
        if (j < MESH_Y/2)
            ky[j] = (double)j * dky;
        else
            ky[j] = (double)(j - MESH_Y) * dky;
    }

    for (int k = 0; k < MESH_Z; k++)
    {
        if (k < MESH_Z/2)
            kz[k] = (double)k * dkz;
        else
            kz[k] = (double)(k - MESH_Z) * dkz;
    }

    if (UMFlag)
    {
        cudaMemPrefetchAsync(dfdc, complex_size, 0, NULL);
        cudaMemPrefetchAsync(dfdphi, complex_size, 0, NULL);
        cudaMemPrefetchAsync(gradphix, complex_size, 0, NULL);
        cudaMemPrefetchAsync(gradphiy, complex_size, 0, NULL);
        cudaMemPrefetchAsync(gradphiz, complex_size, 0, NULL);
        cudaMemPrefetchAsync(kx, MESH_X*sizeof(double), 0, NULL);
        cudaMemPrefetchAsync(ky, MESH_Y*sizeof(double), 0, NULL);
        cudaMemPrefetchAsync(kz, MESH_Z*sizeof(double), 0, NULL);
        cudaMemPrefetchAsync(tempComp, double_size, 0, NULL);
        cudaMemPrefetchAsync(tempPhi, double_size, 0, NULL);
    }

    //Time loop
    while (count <= numsteps && loop_condition == 1)
    {
        if ((count % saveT)==0 || (count == numsteps) ||  (loop_condition == 0))
        {
            if(sim_time == 0)
            {
                printf("\nTime = %6lf\n", sim_time);
                printf("Writing configuration to file!\n");
            }

            if (UMFlag)
            {
                for (int i = 0; i < NUMCOMPONENTS-1; i++)
                    cudaMemPrefetchAsync(CompBuff[i], complex_size, cudaCpuDeviceId);
                cudaMemPrefetchAsync(dfdphi, complex_size, cudaCpuDeviceId);
            }

            for(x = 0; x < rows_x; x++)
            {
                for(y = 0; y < rows_y; y++)
                {
                    for (z = 0; z < rows_z; z++)
                    {
                        index = x*layer_size + y*rows_z + z;

                        for(a = 0; a < NUMCOMPONENTS-1; a++)
                            gridinfo[index].comp[a].x = CompBuff[a][index].x;

                        gridinfo[index].phia[0].x = 1.0;

                        for (b = 1; b < NUMPHASES; b++)
                        {
                            gridinfo[index].phia[b].x = dfdphi[b][index].x;
                            gridinfo[index].phia[0].x -= gridinfo[index].phia[b].x;
                        }
                    }
                }
            }

            if (ASCII == 0)
                writetofile_serial2D_binary(gridinfo, argv, count);
            else
                writetofile_serial2D(gridinfo, argv, count);

            if (UMFlag)
            {
                for (int i = 0; i < NUMCOMPONENTS-1; i++)
                    cudaMemPrefetchAsync(CompBuff[i], complex_size, 0);
                cudaMemPrefetchAsync(dfdphi, complex_size, 0);
            }
        }

        //gradient of phi
        ComputeGradphi_3D<<< Gridsize, Blocksize >>>(kx, ky, kz, MESH_X, MESH_Y, MESH_Z, PhiBuff[1], gradphix, gradphiy, gradphiz);

        if (cufftExecZ2Z(plan, gradphix, gradphix, CUFFT_INVERSE) != CUFFT_SUCCESS)
            printf("inverse fft of gradphix failed\n");
        if (cufftExecZ2Z(plan, gradphiy, gradphiy, CUFFT_INVERSE) != CUFFT_SUCCESS)
            printf("inverse fft of gradphiy failed\n");
        if (cufftExecZ2Z(plan, gradphiz, gradphiz, CUFFT_INVERSE) != CUFFT_SUCCESS)
            printf("inverse fft of gradphiz failed\n");

        Normalize_3D<<< Gridsize, Blocksize >>>(gradphix, sizescale, MESH_X, MESH_Y, MESH_Z);
        Normalize_3D<<< Gridsize, Blocksize >>>(gradphiy, sizescale, MESH_X, MESH_Y, MESH_Z);
        Normalize_3D<<< Gridsize, Blocksize >>>(gradphiz, sizescale, MESH_X, MESH_Y, MESH_Z);

        //Driving forces: varmob (saved in gradphi) and dfdphi
        ComputeDrivForce_3D<<< Gridsize, Blocksize >>>(CompBuff[0], dfdphi[1], gradphix,
                                                       gradphiy, gradphiz, f0AVminv, f0BVminv,
                                                       ceq[1][1][0], ceq[0][0][0], Diffusivity[0][0][0], w, MESH_X, MESH_Y, MESH_Z);

        if (cufftExecZ2Z(plan,gradphix,gradphix,CUFFT_FORWARD) != CUFFT_SUCCESS)
            printf("forward fft of gradphix failed\n");
        if (cufftExecZ2Z(plan,gradphiy,gradphiy,CUFFT_FORWARD) != CUFFT_SUCCESS)
            printf("forward fft of gradphiy failed\n");
        if (cufftExecZ2Z(plan,gradphiz,gradphiz,CUFFT_FORWARD) != CUFFT_SUCCESS)
            printf("forward fft of gradphiz failed\n");

        //dfdc
        ComputeDfdc_3D<<< Gridsize, Blocksize >>>(dfdc[0], gradphix, gradphiy, gradphiz,
                                                  MESH_X, MESH_Y, MESH_Z, kx, ky, kz);

        if (cufftExecZ2Z(plan, CompBuff[0], CompBuff[0], CUFFT_FORWARD) != CUFFT_SUCCESS)
            printf("forward fft of comp failed\n");;
        if (cufftExecZ2Z(plan, dfdphi[1], dfdphi[1], CUFFT_FORWARD) != CUFFT_SUCCESS)
            printf("forward fft of dfphi failed\n");

        //evolve comp and phi
        Update_comp_phi_3D<<< Gridsize, Blocksize >>>(CompBuff[0], dfdc[0], PhiBuff[1], dfdphi[1],
                                                      kx, ky, kz, DELTA_t, Diffusivity[0][0][0], kappa_phi,
                                                      relax_coeff, MESH_X, MESH_Y, MESH_Z);

        cufftExecZ2Z(plan, CompBuff[0], CompBuff[0], CUFFT_INVERSE);
        cufftExecZ2Z(plan, dfdphi[1], dfdphi[1], CUFFT_INVERSE);

        Normalize_3D<<< Gridsize, Blocksize >>>(CompBuff[0], sizescale, MESH_X, MESH_Y, MESH_Z);
        Normalize_3D<<< Gridsize, Blocksize >>>(dfdphi[1], sizescale, MESH_X, MESH_Y, MESH_Z);

        Find_err_matrix_3D<<< Gridsize, Blocksize >>>(tempComp, CompBuff[0], MESH_X, MESH_Y, MESH_Z);

        cudaGetLastError();

        if (UMFlag)
            cudaMemAdvise(maxerrComp, sizeof(double), cudaMemAdviseSetAccessedBy, 0);

        cub::DeviceReduce::Max(t_storage, t_storage_bytes, tempComp, maxerrComp, MESH_X*MESH_Y*MESH_Z);

        if(*maxerrComp <= COMPERR)
        {
            printf("!!!! CONVERGENCE ACHIEVED !!!!.\n");
            printf("Maxerror = %le\n", *maxerrComp);
            loop_condition = 0;
        }

        sim_time = sim_time + DELTA_t;
        count++;

        SaveReal_3D<<< Gridsize, Blocksize >>> (tempComp, CompBuff[0], MESH_X, MESH_Y, MESH_Z);
        if (count % time_output == 0)
        {
            printf("\nTime=%6lf\n", sim_time);

            Find_err_matrix_3D<<< Gridsize, Blocksize >>>(tempPhi, dfdphi[1], MESH_X, MESH_Y, MESH_Z);
            cub::DeviceReduce::Max(t_storage, t_storage_bytes, tempPhi, maxerrPhi, MESH_X*MESH_Y*MESH_Z);
            SaveReal_3D<<< Gridsize, Blocksize >>>(tempPhi, dfdphi[1], MESH_X, MESH_Y, MESH_Z);
            cub::DeviceReduce::Max(t_storage, t_storage_bytes, tempPhi, maxVal, MESH_X*MESH_Y*MESH_Z);
            cub::DeviceReduce::Min(t_storage, t_storage_bytes, tempPhi, minVal, MESH_X*MESH_Y*MESH_Z);
            cudaDeviceSynchronize();

            printf("%*s, Max = %lf, Min = %lf, Max. error = %lf\n", 5, PHASES[0], *maxVal, *minVal, *maxerrPhi);

            cub::DeviceReduce::Max(t_storage, t_storage_bytes, tempComp, maxVal, MESH_X*MESH_Y*MESH_Z);
            cub::DeviceReduce::Min(t_storage, t_storage_bytes, tempComp, minVal, MESH_X*MESH_Y*MESH_Z);
            cudaDeviceSynchronize();

            printf("%*s, Max = %lf, Min = %lf, Max. error = %lf\n", 5, COMPONENTS[0], *maxVal, *minVal, *maxerrComp);

            if (count % saveT == 0)
                printf("Writing configuration to file!\n");
        }
        else if (count % time_output == time_output-1)
            SaveReal_3D<<< Gridsize, Blocksize >>>(tempPhi, dfdphi[1], MESH_X, MESH_Y, MESH_Z);

    }//time loop ends

    cudaFree(tempComp);
    cudaFree(tempPhi);
    cudaFree(maxerrComp);
    cudaFree(maxerrPhi);
    cudaFree(maxVal);
    cudaFree(minVal);
    cudaFree(t_storage);

    for (int j = 0; j < NUMCOMPONENTS-1; j++)
        cudaFree(CompBuff[j]);
    for(int j = 0; j < NUMPHASES; j++)
        cudaFree(PhiBuff[j]);

    for (int j = 0; j < NUMCOMPONENTS-1; j++)
        cudaFree(dfdc[j]);
    for(int j = 0; j < NUMPHASES; j++)
        cudaFree(dfdphi[j]);

    cudaFree(CompBuff);
    cudaFree(PhiBuff);
    cudaFree(dfdc);
    cudaFree(dfdphi);

    cudaFree(gradphix);
    cudaFree(gradphiy);
    cudaFree(gradphiz);

    cudaFree(kx);
    cudaFree(ky);
    cudaFree(kz);

    cudaFree(B);

    for (int i = 0; i < rows_x*rows_y*rows_z; i++)
    {
        free(gridinfo[i].phia);
        free(gridinfo[i].comp);
    }
    free(gridinfo);
}
