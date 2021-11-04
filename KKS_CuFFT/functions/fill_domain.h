#ifndef FILL_DOMAIN_H_
#define FILL_DOMAIN_H_

void fill_domain(char *argv[])
{
    int i;
    char tempbuff[100];

    char tmpstr1[100];
    char tmpstr2[100];
    char **tmp;

    char *str1, *token;
    char *saveptr1;

    long phase;

    FILE *fr;
    if (fr = fopen(argv[2], "r"))
        printf("\nReading filling parameters from %s\n", argv[2]);
    else
        printf("\nFile %s not found\n", argv[2]);

    cudaMemPrefetchAsync(gridinfo, gridinfo_size, cudaCpuDeviceId);

    while (fgets(tempbuff, 100, fr))
    {
        sscanf(tempbuff, "%100s = %100[^;];", tmpstr1, tmpstr2);
        if (tmpstr1[0] != '#')
        {
            if ((strcmp(tmpstr1, "FILLCYLINDER") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) > 0))
            {
                printf("Filling cylinder\n");
                tmp = (char**)malloc(sizeof(char*)*6);
                for (i = 0; i < 6; i++)
                {
                    tmp[i] = (char*)malloc(sizeof(char)*10);
                }
                for (i = 0, str1 = tmpstr2; ; i++, str1 = NULL)
                {
                    token = strtok_r(str1, "{,}", &saveptr1);
                    if (token == NULL)
                        break;
                    strcpy(tmp[i],token);
                }
                phase = atol(tmp[0]);

                fill_cylinder_parameters.x_center = atol(tmp[1]) + start[X];
                fill_cylinder_parameters.y_center = atol(tmp[2]) + start[Y];
                fill_cylinder_parameters.z_start  = atol(tmp[3]) + start[Z];
                fill_cylinder_parameters.z_end    = atol(tmp[4]) + start[Z];
                fill_cylinder_parameters.radius   = atof(tmp[5]);

                fill_phase_cylinder_cuda<<<Gridsize, Blocksize>>>(fill_cylinder_parameters, PhiBuff, dfdphi, phase, MESH_X, MESH_Y, MESH_Z, NUMPHASES, NUMCOMPONENTS);

                for (i = 0; i < 6; i++)
                {
                    free(tmp[i]);
                }
                free(tmp);
                printf("End filling cylinder\n");
            }
            else if ((strcmp(tmpstr1, "FILLSPHERE") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) > 0))
            {
                printf("Filling sphere\n");
                tmp = (char**)malloc(sizeof(char*)*5);
                for (i = 0; i < 5; i++)
                {
                    tmp[i] = (char*)malloc(sizeof(char)*10);
                }
                for (i = 0, str1 = tmpstr2; ; i++, str1 = NULL)
                {
                    token = strtok_r(str1, "{,}", &saveptr1);
                    if (token == NULL)
                        break;
                    strcpy(tmp[i],token);
                }
                phase = atol(tmp[0]);

                fill_sphere_parameters.x_center = atol(tmp[1]) + start[X];
                fill_sphere_parameters.y_center = atol(tmp[2]) + start[Y];
                fill_sphere_parameters.z_center = atol(tmp[3]) + start[Z];
                fill_sphere_parameters.radius   = atof(tmp[4]);

                fill_phase_sphere_cuda<<<Gridsize, Blocksize>>>(fill_sphere_parameters, PhiBuff, dfdphi, phase, MESH_X, MESH_Y, MESH_Z, NUMPHASES, NUMCOMPONENTS);
                for (i = 0; i < 5; i++)
                {
                    free(tmp[i]);
                }
                free(tmp);
                printf("End filling sphere\n");
            }
            else if ((strcmp(tmpstr1, "FILLCYLINDERRANDOM") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) > 0))
            {
                printf("Filling cylinders at random\n");
                tmp = (char**)malloc(sizeof(char*)*4);
                for (i = 0; i < 4; i++)
                {
                    tmp[i] = (char*)malloc(sizeof(char)*10);
                }
                for (i = 0, str1 = tmpstr2; ; i++, str1 = NULL)
                {
                    token = strtok_r(str1, "{,}", &saveptr1);
                    if (token == NULL)
                        break;
                    strcpy(tmp[i],token);
                }

                phase           = atol(tmp[0]);
                ppt_radius      = atol(tmp[1]);
                volume_fraction = atof(tmp[2]);
                shield_dist     = atol(tmp[3]);

                if (shield_dist > 8)
                    shield_dist = 8;
                else if (shield_dist == 1)
                    shield_dist = 2;

                fill_phase_cylinder_random(phase);

                for (i = 0; i < 4; i++)
                {
                    free(tmp[i]);
                }
                free(tmp);
                printf("End filling cylinders at random\n");
            }
            else if ((strcmp(tmpstr1, "FILLSPHERERANDOM") == 0) && (NUMPHASES > 0) && ((NUMCOMPONENTS-1) > 0))
            {
                printf("Filling spheres at random\n");
                tmp = (char**)malloc(sizeof(char*)*4);
                for (i = 0; i < 4; i++)
                {
                    tmp[i] = (char*)malloc(sizeof(char)*10);
                }
                for (i = 0, str1 = tmpstr2; ; i++, str1 = NULL)
                {
                    token = strtok_r(str1, "{,}", &saveptr1);
                    if (token == NULL)
                        break;
                    strcpy(tmp[i],token);
                }

                phase           = atol(tmp[0]);
                ppt_radius      = atol(tmp[1]);
                volume_fraction = atof(tmp[2]);
                shield_dist     = atol(tmp[3]);

                if (shield_dist > 8)
                    shield_dist = 8;
                else if (shield_dist == 1)
                    shield_dist = 2;

                fill_phase_sphere_random(phase);

                for (i = 0; i < 4; i++)
                {
                    free(tmp[i]);
                }
                free(tmp);
                printf("End filling spheres at random\n");
            }
        }
    }

    fclose(fr);
    printf("Filling composition\n");
    fill_composition_cuda<<<Gridsize, Blocksize>>>(PhiBuff, CompBuff, MESH_X, MESH_Y, MESH_Z, NUMPHASES, NUMCOMPONENTS, c0, ceq[1][1][0]);
    printf("Filled composition\n");
}
#endif
