#ifndef READ_RESTART_H_
#define READ_RESTART_H_

int read_vtk_ascii(FILE *fp, struct fields *gridinfo)
{
    char temp[1000];
    long a = 0, k = 0, x, y, z;

    while(fscanf(fp, "%s", temp))
    {
        if (strcmp(temp, "DIMENSIONS") == 0)
        {
            long temp_mx = 0, temp_my = 0, temp_mz = 0;
            fscanf(fp, "%ld", &temp_mx);
            fscanf(fp, "%ld", &temp_my);
            fscanf(fp, "%ld", &temp_mz);

            printf("Read dimensions: %ld, %ld, %ld\n", temp_mx, temp_my, temp_mz);
            printf("Input dimensions: %ld, %ld, %ld\n", MESH_X, MESH_Y, MESH_Z);

            if (temp_mx != MESH_X || temp_my != MESH_Y || temp_mz != MESH_Z)
            {
                printf("Dimensions do not match. Filling using specified filling file\n");
                return 0;
            }
        }
        else if (strcmp(temp, "SPACING") == 0)
        {
            double temp_dx = 0, temp_dy = 0, temp_dz = 0;

            fscanf(fp, "%le", &temp_dx);
            fscanf(fp, "%le", &temp_dy);
            fscanf(fp, "%le", &temp_dz);

            if (temp_dx != DELTA_X || temp_dy != DELTA_Y || temp_dz != DELTA_Z)
            {
                printf("Dimensions do not match. Filling using specified filling file\n");
                return 0;
            }
        }
        else if (strcmp(temp, "default") == 0 && a < NUMPHASES)
        {
            printf("Reading %s\n", PHASES[a]);
            for (x = start[X]; x <= end[X]; x++) {
                for (y = start[Y]; y <= end[Y]; y++) {
                    for (z = start[Z]; z <= end[Z]; z++) {
                        fscanf(fp, "%le", &gridinfo[x*layer_size + y*rows_z + z].phia[a]);
                    }
                }
            }
            a++;
        }
        else if (strcmp(temp, "default") == 0 && a >= NUMPHASES && k < NUMCOMPONENTS-1)
        {
            printf("Reading %s\n", COMPONENTS[k]);
            for (x = start[X]; x <= end[X]; x++) {
                for (y = start[Y]; y <= end[Y]; y++) {
                    for (z = start[Z]; z <= end[Z]; z++) {
                        fscanf(fp, "%le", &gridinfo[x*layer_size + y*rows_z + z].compi[k]);
                    }
                }
            }
            k++;
        }
        else if (k == NUMCOMPONENTS -  1)
            break;
    }
    return 1;
}

int read_vtk_binary(FILE *fp, struct fields *gridinfo)
{
    char temp[1000];
    long a = 1, k = 0, x, y, z;
    double value;

    while(fscanf(fp, "%s", temp))
    {
        if (strcmp(temp, "DIMENSIONS") == 0)
        {
            long temp_mx = 0, temp_my = 0, temp_mz = 0;
            fscanf(fp, "%ld", &temp_mx);
            fscanf(fp, "%ld", &temp_my);
            fscanf(fp, "%ld", &temp_mz);

            printf("Read dimensions: %ld, %ld, %ld\n", temp_mx, temp_my, temp_mz);
            printf("Input dimensions: %ld, %ld, %ld\n", MESH_X, MESH_Y, MESH_Z);

            if (temp_mx != MESH_X || temp_my != MESH_Y || temp_mz != MESH_Z)
            {
                printf("Dimensions do not match. Filling using specified filling file\n");
                return 0;
            }
        }
        else if (strcmp(temp, "SPACING") == 0)
        {
            double temp_dx = 0, temp_dy = 0, temp_dz = 0;

            fscanf(fp, "%le", &temp_dx);
            fscanf(fp, "%le", &temp_dy);
            fscanf(fp, "%le", &temp_dz);

            if (temp_dx != DELTA_X || temp_dy != DELTA_Y || temp_dz != DELTA_Z)
            {
                printf("Dimensions do not match. Filling using specified filling file\n");
                return 0;
            }
        }
        else if (strcmp(temp, "default") == 0 && a < NUMPHASES)
        {
            printf("Reading %s\n", PHASES[a]);
            fscanf(fp, "%lf", &value);

            for (x = start[X]; x <= end[X]; x++) {
                for (y = start[Y]; y <= end[Y]; y++) {
                    for (z = start[Z]; z <= end[Z]; z++) {
                        fread(&value, sizeof(double), 1, fp);
                        gridinfo[x*layer_size + y*rows_z + z].phia[a] = value;
                        if (IS_LITTLE_ENDIAN)
                            gridinfo[x*layer_size + y*rows_z + z].phia[a] = swap_bytes(value);
                        else
                            gridinfo[x*layer_size + y*rows_z + z].phia[a] = value;
                    }
                }
            }
            a++;
        }
        else if (strcmp(temp, "default") == 0 && a >= NUMPHASES && k < NUMCOMPONENTS-1)
        {
            printf("Reading %s\n", COMPONENTS[k]);
            fscanf(fp, "%lf", &value);
            for (x = start[X]; x <= end[X]; x++) {
                for (y = start[Y]; y <= end[Y]; y++) {
                    for (z = start[Z]; z <= end[Z]; z++) {
                        fread(&value, sizeof(double), 1, fp);
                        if (IS_LITTLE_ENDIAN)
                            gridinfo[x*layer_size + y*rows_z + z].compi[k] = swap_bytes(value);
                        else
                            gridinfo[x*layer_size + y*rows_z + z].compi[k] = value;
                    }
                }
            }
            k++;
        }
        else if (k == NUMCOMPONENTS -  1)
            break;
    }
    return 1;
}

/******************************************************************************
 * Driver function for simulation resumption functionality. Calls appropriate
 * read functions depending on the type of file that is in the folder.
 *
 * Arguments:
 *  char *argv[]    : arguments passed to program execution call. Uses argv[2] to
 *                    read filling file if no data exists or if it can not be
 *                    read.
 *****************************************************************************/
void read_domain(char *argv[], long t)
{
    FILE *fp;
    char name[1000], writeformat[10];

    sprintf(name,"DATA/%s_%07ld.vtk", argv[3], t);

    if (fp = fopen(name, "rb"))
    {
        printf("\nReading from %s\n", name);

        fscanf(fp,"%*[^\n]\n");
        fscanf(fp,"%*[^\n]\n");

        fscanf(fp, "%s", writeformat);

        if (strcmp(writeformat, "ASCII") == 0)
        {
            if (read_vtk_ascii(fp, gridinfo))
            {
                gridinfo_to_gpu(gridinfo, compHost, phiHost);
                for (int i = 1; i < NUMPHASES; i++)
                    cudaMemcpy(&dfdphiHost[i], phiHost[i], sizeof(cufftDoubleComplex)*MESH_X*MESH_Y*MESH_Z, cudaMemcpyDeviceToDevice);
            }
            else
                fill_domain(argv);
        }
        else if (strcmp(writeformat, "BINARY") == 0)
        {
            if (read_vtk_binary(fp, gridinfo))
            {
                gridinfo_to_gpu(gridinfo, compHost, phiHost);
                for (int i = 1; i < NUMPHASES; i++)
                    cudaMemcpy(dfdphiHost[i], phiHost[i], sizeof(cufftDoubleComplex)*MESH_X*MESH_Y*MESH_Z, cudaMemcpyDeviceToDevice);
            }
            else
                fill_domain(argv);
        }
        fclose(fp);
    }
    else
    {
        printf("\nCould not find %s.\nFilling domain using %s\n", name, argv[2]);
        fill_domain(argv);
    }
}

#endif
