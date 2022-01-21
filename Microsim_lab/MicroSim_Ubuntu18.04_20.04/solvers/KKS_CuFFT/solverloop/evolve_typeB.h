void evolve_FD(char *argv[])
{
    int     loop_condition, count = initcount;
    double  *tempComp, *tempPhi;
    double  *maxerrComp, *maxerrPhi;
    double  *maxerrComp_h, *maxerrPhi_h;
    double  *maxVal, *minVal;
    double  *maxVal_h, *minVal_h;
    double  f0AVminv, f0BVminv;
    void    *t_storage = NULL;
    size_t  t_storage_bytes = 0;
    double  *gradPhi;

    kappa_phi = 3.0/(2.0*alpha) * (Gamma[0][1]*Ln);
    w = 6.0 * alpha * Gamma[0][1] / Ln;

    f0AVminv = F0[0][0][0] * (1.0/Vm);
    f0BVminv = F0[1][0][0] * (1.0/Vm);

    cudaMalloc((void**)&tempComp, sizeof(double)*MESH_X*MESH_Y*MESH_Z);
    cudaMalloc((void**)&tempPhi, sizeof(double)*MESH_X*MESH_Y*MESH_Z);
    cudaMalloc((void**)&maxerrComp, sizeof(double));
    cudaMalloc((void**)&maxerrPhi, sizeof(double));
    cudaMalloc((void**)&maxVal, sizeof(double));
    cudaMalloc((void**)&minVal, sizeof(double));

    cudaMalloc((void**)&gradPhi, sizeof(double)*MESH_X*MESH_Y*MESH_Z);

    maxerrComp_h = (double*)malloc(sizeof(double));
    maxerrPhi_h = (double*)malloc(sizeof(double));

    maxVal_h = (double*)malloc(sizeof(double));
    minVal_h = (double*)malloc(sizeof(double));

    cudaMemset(maxerrComp, 1, sizeof(double));
    cudaMemset(maxerrPhi, 1, sizeof(double));

    cub::DeviceReduce::Max(t_storage, t_storage_bytes, tempComp, maxerrComp, MESH_X * MESH_Y * MESH_Z);
    cub::DeviceReduce::Sum(t_storage, t_storage_bytes, tempComp, maxerrComp, MESH_X * MESH_Y * MESH_Z);

    cudaMalloc((void**)&t_storage, t_storage_bytes);

    //save real part
    SaveReal<<< Gridsize, Blocksize >>>(tempComp, compHost[0], MESH_X, MESH_Y, MESH_Z);
    SaveReal<<< Gridsize, Blocksize >>>(tempPhi, phiHost[1], MESH_X, MESH_Y, MESH_Z);

    loop_condition = 1;

    dkx = 2.0 * CUDART_PI / ((double) MESH_X * DELTA_X);
    dky = 2.0 * CUDART_PI / ((double) MESH_Y * DELTA_Y);
    dkz = 2.0 * CUDART_PI / ((double) MESH_Z * DELTA_Z);

    sizescale = 1.0/(double)(MESH_X * MESH_Y * MESH_Z);

    for (int i = 0; i < MESH_X; i++)
    {
        if (i < MESH_X/2)
            kx_h[i] = (double)i * dkx;
        else
            kx_h[i] = (double)(i - MESH_X) * dkx;
    }

    for (int j = 0; j < MESH_Y; j++)
    {
        if (j < MESH_Y/2)
            ky_h[j] = (double)j * dky;
        else
            ky_h[j] = (double)(j - MESH_Y) * dky;
    }

    if (DIMENSION == 3)
    {
        for (int k = 0; k < MESH_Z; k++)
        {
            if (k < MESH_Z/2)
                kz_h[k] = (double)k * dkz;
            else
                kz_h[k] = (double)(k - MESH_Z) * dkz;
        }
    }
    else
        kz_h[0] = 0.0;

    cudaMemcpy(kx, kx_h, MESH_X*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(ky, ky_h, MESH_Y*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(kz, kz_h, MESH_Z*sizeof(double), cudaMemcpyHostToDevice);

    printf("\nTime = %6lf\n", sim_time);

    //Time loop
    while (count <= initcount+numsteps && loop_condition == 1)
    {
        if ((count % saveT) == 0 || (count == numsteps) ||  (loop_condition == 0) && (count == 0 || count != initcount))
        {
            gpu_to_gridinfo(gridinfo, compHost, phiHost);

            printf("Writing configuration to file! Time = %06lf\n", sim_time);

            if (ASCII == 0)
                writetofile_serial2D_binary(gridinfo, argv, count);
            else
                writetofile_serial2D(gridinfo, argv, count);

            if (count == initcount+numsteps)
                break;
        }

        //Driving forces: varmob (saved in gradphi) and dfdphi
        ComputeDrivForce_FD<<< Gridsize, Blocksize >>>(compHost[0], phiHost[1],
                                                       dfdcHost[0], dfdphiHost[1],
                                                       f0AVminv, f0BVminv,
                                                       ceq[0][0][0], ceq[1][1][0],
                                                       Diffusivity[0][0][0], w,
                                                       MESH_X, MESH_Y, MESH_Z);

        ComputeGradphi_x<<< Gridsize, Blocksize >>>(tempPhi, gradPhi,
                                                    MESH_X, MESH_Y, MESH_Z, DELTA_X);

        ComputeDfdc_x<<<Gridsize, Blocksize>>>(dfdcHost[0], gradPhi, MESH_X, MESH_Y, MESH_Z, DELTA_X);

        ComputeGradphi_y<<< Gridsize, Blocksize >>>(tempPhi, gradPhi,
                                                    MESH_X, MESH_Y, MESH_Z, DELTA_Y);

        ComputeDfdc_y<<<Gridsize, Blocksize>>>(dfdcHost[0], gradPhi, MESH_X, MESH_Y, MESH_Z, DELTA_Y);

        if (DIMENSION == 3)
        {
            ComputeGradphi_z<<< Gridsize, Blocksize >>>(tempPhi, gradPhi,
                                                        MESH_X, MESH_Y, MESH_Z, DELTA_Z);

            ComputeDfdc_z<<<Gridsize, Blocksize>>>(dfdcHost[0], gradPhi, MESH_X, MESH_Y, MESH_Z, DELTA_Z);
        }

        reset_dfdc<<<Gridsize, Blocksize>>>(dfdcHost[0], MESH_X, MESH_Y, MESH_Z);

        if (cufftExecZ2Z(plan, compHost[0], compHost[0], CUFFT_FORWARD) != CUFFT_SUCCESS)
            printf("forward fft of comp failed\n");
        if (cufftExecZ2Z(plan, phiHost[1], phiHost[1], CUFFT_FORWARD) != CUFFT_SUCCESS)
            printf("forward fft of phi failed\n");

        if (cufftExecZ2Z(plan, dfdcHost[0], dfdcHost[0], CUFFT_FORWARD) != CUFFT_SUCCESS)
            printf("forward fft of dfdc failed\n");
        if (cufftExecZ2Z(plan, dfdphiHost[1], dfdphiHost[1], CUFFT_FORWARD) != CUFFT_SUCCESS)
            printf("forward fft of dfdphi failed\n");

        //evolve comp and phi
        Update_comp_phi<<< Gridsize, Blocksize >>>(compHost[0], phiHost[1],
                                                   dfdcHost[0], dfdphiHost[1],
                                                   B, elast_int,
                                                   DELTA_t, Diffusivity[0][0][0],
                                                   kappa_phi, relax_coeff,
                                                   kx, ky, kz,
                                                   MESH_X, MESH_Y, MESH_Z);

        cufftExecZ2Z(plan, compHost[0], compHost[0], CUFFT_INVERSE);
        cufftExecZ2Z(plan, phiHost[1], phiHost[1], CUFFT_INVERSE);

        Normalize<<< Gridsize, Blocksize >>>(compHost[0], sizescale, MESH_X, MESH_Y, MESH_Z);
        Normalize<<< Gridsize, Blocksize >>>(phiHost[1], sizescale, MESH_X, MESH_Y, MESH_Z);

        Find_err_matrix<<< Gridsize, Blocksize >>>(tempComp, compHost[0], MESH_X, MESH_Y, MESH_Z);

        cudaDeviceSynchronize();

        cub::DeviceReduce::Max(t_storage, t_storage_bytes, tempComp, maxerrComp, MESH_X*MESH_Y*MESH_Z);

        cudaMemcpy(maxerrComp_h, maxerrComp, sizeof(double), cudaMemcpyDeviceToHost);

        if(*maxerrComp_h <= COMPERR)
        {
            printf("!!!! CONVERGENCE ACHIEVED !!!!.\n");
            printf("Maxerror = %le\n", *maxerrComp_h);
            loop_condition = 0;
        }

        sim_time = sim_time + DELTA_t;
        count++;

        SaveReal<<< Gridsize, Blocksize >>> (tempComp, compHost[0], MESH_X, MESH_Y, MESH_Z);
        if (count % time_output == 0)
        {
            printf("\nTime=%6lf\n", sim_time);

            Find_err_matrix<<< Gridsize, Blocksize >>>(tempPhi, phiHost[1], MESH_X, MESH_Y, MESH_Z);
            cub::DeviceReduce::Max(t_storage, t_storage_bytes, tempPhi, maxerrPhi, MESH_X*MESH_Y*MESH_Z);
            SaveReal<<< Gridsize, Blocksize >>>(tempPhi, phiHost[1], MESH_X, MESH_Y, MESH_Z);
            cub::DeviceReduce::Max(t_storage, t_storage_bytes, tempPhi, maxVal, MESH_X*MESH_Y*MESH_Z);
            cub::DeviceReduce::Min(t_storage, t_storage_bytes, tempPhi, minVal, MESH_X*MESH_Y*MESH_Z);

            cudaMemcpy(maxVal_h, maxVal, sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(minVal_h, minVal, sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(maxerrPhi_h, maxerrPhi, sizeof(double), cudaMemcpyDeviceToHost);

            printf("%*s, Max = %lf, Min = %lf, Relative_Change = %lf\n", 5, PHASES[1], *maxVal_h, *minVal_h, *maxerrPhi_h);

            cub::DeviceReduce::Max(t_storage, t_storage_bytes, tempComp, maxVal, MESH_X*MESH_Y*MESH_Z);
            cub::DeviceReduce::Min(t_storage, t_storage_bytes, tempComp, minVal, MESH_X*MESH_Y*MESH_Z);

            cudaMemcpy(maxVal_h, maxVal, sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(minVal_h, minVal, sizeof(double), cudaMemcpyDeviceToHost);

            printf("%*s, Max = %lf, Min = %lf, Relative_Change = %lf\n", 5, COMPONENTS[0], *maxVal_h, *minVal_h, *maxerrComp_h);
        }
        else
            SaveReal<<< Gridsize, Blocksize >>>(tempPhi, phiHost[1], MESH_X, MESH_Y, MESH_Z);
    }//time loop ends

    cudaFree(tempComp);
    cudaFree(tempPhi);
    cudaFree(maxerrComp);
    cudaFree(maxerrPhi);
    free(maxerrComp_h);
    free(maxerrPhi_h);
    cudaFree(maxVal);
    cudaFree(minVal);
    free(maxVal_h);
    free(minVal_h);
    cudaFree(t_storage);

    for (int j = 0; j < NUMCOMPONENTS-1; j++)
    {
        cudaFree(compHost[j]);
        cudaFree(dfdcHost[j]);
    }

    cudaFree(compDev);
    cudaFree(dfdcDev);
    free(compHost);
    free(dfdcHost);

    for(int j = 1; j < NUMPHASES; j++)
    {
        cudaFree(phiHost[j]);
        cudaFree(dfdphiHost[j]);
    }

    cudaFree(phiDev);
    cudaFree(dfdphiDev);
    free(phiHost);
    free(dfdphiHost);

    cudaFree(kx);
    cudaFree(ky);
    cudaFree(kz);

    cudaFree(gradPhi);

    free(kx_h);
    free(ky_h);
    free(kz_h);

    if (elast_int)
        cudaFree(B);

    //     for (int i = 0; i < NUMCOMPONENTS; i++)
    //         free(dfdc_thermo[i]);
    //
    //     free(dfdc_thermo);

    for (int i = 0; i < rows_x*rows_y*rows_z; i++)
    {
        free(gridinfo[i].phia);
        free(gridinfo[i].compi);
    }
    free(gridinfo);
}
