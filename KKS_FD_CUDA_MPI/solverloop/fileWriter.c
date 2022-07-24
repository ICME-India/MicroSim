#include "fileWriter.h"

double swap_bytes(double value)
{
    double  src_num = value;
    int64_t tmp_num = htobe64(le64toh(*(int64_t*)&src_num));
    double  dst_num = *(double*)&tmp_num;
    return  dst_num;
}

void writeVTK_ASCII(double *phi, double *comp,
                    domainInfo simDomain, subdomainInfo subdomain,
                    int t, int rank, MPI_Comm comm,
                    char *argv[])
{
    FILE *fp;
    char name[100];

    int layer_size = simDomain.MESH_X * simDomain.MESH_Y;

    sprintf(name, "DATA/Processor_%d/%s_%07d.vtk", rank, argv[3], t);
    fp = fopen(name, "w");

    /*
     * Metadata required for the reader to understand the data's topology and format
     */
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "Microsim_fields\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_POINTS\n");
    fprintf(fp, "DIMENSIONS %d %d %d\n", subdomain.xE-subdomain.xS+1, subdomain.yE-subdomain.yS+1, subdomain.zE-subdomain.zS+1);
    fprintf(fp, "ORIGIN 0 0 0\n");
    fprintf(fp, "SPACING %le %le %le\n", simDomain.DELTA_X, simDomain.DELTA_Y, simDomain.DELTA_Z);
    fprintf(fp, "POINT_DATA %d\n", subdomain.numCells);

    for (int b = 0; b < simDomain.numComponents-1; b++)
    {
        fprintf(fp, "SCALARS %s double 1\n", simDomain.componentNames[b]);
        fprintf(fp, "LOOKUP_TABLE default\n");

        for (int z = 0; z <= subdomain.zE - subdomain.zS; z++)
        {
            for (int y = 0; y <= subdomain.yE - subdomain.yS; y++)
            {
                for (int x = 0; x <= subdomain.xE - subdomain.xS; x++)
                {
                    fprintf(fp, "%le\n", comp[b*subdomain.numCells + z*layer_size + y*simDomain.MESH_X + x]);
                }
            }
        }
        fprintf(fp, "\n");
    }

    for (int a = 0; a < simDomain.numPhases; a++)
    {
        fprintf(fp, "SCALARS %s double 1\n", simDomain.phaseNames[a]);
        fprintf(fp, "LOOKUP_TABLE default\n");


        for (int z = 0; z <= subdomain.zE - subdomain.zS; z++)
        {
            for (int y = 0; y <= subdomain.yE - subdomain.yS; y++)
            {
                for (int x = 0; x <= subdomain.xE - subdomain.xS; x++)
                {
                    fprintf(fp, "%le\n", phi[a*subdomain.numCells + z*layer_size + y*simDomain.MESH_X + x]);
                }
            }
        }

        fprintf(fp, "\n");
    }

    fclose(fp);
}

void writeVTK_BINARY(double *phi, double *comp,
                     domainInfo simDomain, subdomainInfo subdomain,
                     int t, int rank, MPI_Comm comm,
                     char *argv[])
{
    FILE *fp;
    char name[100];

    int layer_size = simDomain.MESH_X * simDomain.MESH_Y;

    double value;

    sprintf(name, "DATA/Processor_%d/%s_%07d.vtk", rank, argv[3], t);
    fp = fopen(name, "w");

    /*
     * Metadata required for the reader to understand the data's topology and format
     */
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "Microsim_fields\n");
    fprintf(fp, "BINARY\n");
    fprintf(fp, "DATASET STRUCTURED_POINTS\n");
    fprintf(fp, "DIMENSIONS %d %d %d\n", subdomain.xE-subdomain.xS+1, subdomain.yE-subdomain.yS+1, subdomain.zE-subdomain.zS+1);
    fprintf(fp, "ORIGIN 0 0 0\n");
    fprintf(fp, "SPACING %le %le %le\n", simDomain.DELTA_X, simDomain.DELTA_Y, simDomain.DELTA_Z);
    fprintf(fp, "POINT_DATA %d\n", subdomain.numCells);

    for (int b = 0; b < simDomain.numComponents-1; b++)
    {
        fprintf(fp, "SCALARS %s double 1\n", simDomain.componentNames[b]);
        fprintf(fp, "LOOKUP_TABLE default\n");

        for (int z = 0; z <= subdomain.zE - subdomain.zS; z++)
        {
            for (int y = 0; y <= subdomain.yE - subdomain.yS; y++)
            {
                for (int x = 0; x <= subdomain.xE - subdomain.xS; x++)
                {
                    if (IS_LITTLE_ENDIAN)
                        value = swap_bytes(comp[b*subdomain.numCells + z*layer_size + y*simDomain.MESH_X + x]);
                    else
                        value = comp[b*subdomain.numCells + z*layer_size + y*simDomain.MESH_X + x];
                    fwrite(&value, sizeof(double), 1, fp);
                }
            }
        }
        fprintf(fp, "\n");
    }

    for (int a = 0; a < simDomain.numPhases; a++)
    {
        fprintf(fp, "SCALARS %s double 1\n", simDomain.phaseNames[a]);
        fprintf(fp, "LOOKUP_TABLE default\n");

        for (int z = 0; z <= subdomain.zE - subdomain.zS; z++)
        {
            for (int y = 0; y <= subdomain.yE - subdomain.yS; y++)
            {
                for (int x = 0; x <= subdomain.xE - subdomain.xS; x++)
                {
                    if (IS_LITTLE_ENDIAN)
                        value = swap_bytes(phi[a*subdomain.numCells + z*layer_size + y*simDomain.MESH_X + x]);
                    else
                        value = phi[a*subdomain.numCells + z*layer_size + y*simDomain.MESH_X + x];
                    fwrite(&value, sizeof(double), 1, fp);
                }
            }
        }

        fprintf(fp, "\n");
    }

    fclose(fp);
}

int readVTK_ASCII(FILE *fp, double *phi, double *comp,
                  domainInfo simDomain, subdomainInfo subdomain,
                  int t)
{
    char temp[1000];
    int a = 0, k = 0, x, y, z;
    double value;

    int temp_mx = 0, temp_my = 0, temp_mz = 0;

    int layer_size = simDomain.MESH_X*simDomain.MESH_Y;

    while(fscanf(fp, "%s", temp))
    {
        if (strcmp(temp, "DIMENSIONS") == 0)
        {
            if(fscanf(fp, "%d", &temp_mx));
            if(fscanf(fp, "%d", &temp_my));
            if(fscanf(fp, "%d", &temp_mz));

//             printf("Read dimensions: %d, %d, %d\n", temp_mx, temp_my, temp_mz);
//             printf("Required dimensions: %d, %d, %d\n", (subdomain.xE-subdomain.xS+1), (subdomain.yE-subdomain.yS+1), (subdomain.zE-subdomain.zS+1));

            if (temp_mx != (subdomain.xE-subdomain.xS+1) || temp_my != (subdomain.yE-subdomain.yS+1) || temp_mz != (subdomain.zE-subdomain.zS+1))
            {
                printf("Dimensions do not match. Filling using specified filling file\n");
                return 0;
            }
        }
        else if (strcmp(temp, "SPACING") == 0)
        {
            double temp_dx = 0, temp_dy = 0, temp_dz = 0;

            if(fscanf(fp, "%le", &temp_dx));
            if(fscanf(fp, "%le", &temp_dy));
            if(fscanf(fp, "%le", &temp_dz));

            if (temp_dx != simDomain.DELTA_X || temp_dy != simDomain.DELTA_Y || temp_dz != simDomain.DELTA_Z)
            {
                printf("Dimensions do not match. Filling using specified filling file\n");
                return 0;
            }
        }
        else if (strcmp(temp, "default") == 0 && k < simDomain.numComponents-1)
        {
           // printf("Reading %s\n", simDomain.componentNames[k]);

            for (z = 0; z <= subdomain.zE - subdomain.zS; z++)
            {
                for (y = 0; y <= subdomain.yE - subdomain.yS; y++)
                {
                    for (x = 0; x <= subdomain.xE - subdomain.xS; x++)
                    {
                        if(fscanf(fp, "%le", &comp[k*subdomain.numCells + z*layer_size + y*simDomain.MESH_X + x]));
                    }
                }
            }
            k++;
        }
        else if (strcmp(temp, "default") == 0 && a < simDomain.numPhases && k >= simDomain.numComponents-1)
        {
           // printf("Reading %s\n", simDomain.phaseNames[a]);

            for (z = 0; z <= subdomain.zE - subdomain.zS; z++)
            {
                for (y = 0; y <= subdomain.yE - subdomain.yS; y++)
                {
                    for (x = 0; x <= subdomain.xE - subdomain.xS; x++)
                    {
                        if(fscanf(fp, "%le", &phi[a*subdomain.numCells + z*layer_size + y*simDomain.MESH_X + x]));
                    }
                }
            }
            a++;
        }
        else if (a == simDomain.numPhases)
            break;
    }

    return 1;
}

int readVTK_BINARY(FILE *fp, double *phi, double *comp,
                   domainInfo simDomain, subdomainInfo subdomain,
                   int t)
{
    char temp[1000];
    int a = 0, k = 0, x, y, z;
    double value;

    int temp_mx = 0, temp_my = 0, temp_mz = 0;

    int layer_size = simDomain.MESH_X*simDomain.MESH_Y;

    while(fscanf(fp, "%s", temp))
    {
        if (strcmp(temp, "DIMENSIONS") == 0)
        {
            if(fscanf(fp, "%d", &temp_mx));
            if(fscanf(fp, "%d", &temp_my));
            if(fscanf(fp, "%d", &temp_mz));

//             printf("Read dimensions: %d, %d, %d\n", temp_mx, temp_my, temp_mz);
//             printf("Required dimensions: %d, %d, %d\n", (subdomain.xE-subdomain.xS+1), (subdomain.yE-subdomain.yS+1), (subdomain.zE-subdomain.zS+1));

            if (temp_mx != (subdomain.xE-subdomain.xS+1) || temp_my != (subdomain.yE-subdomain.yS+1) || temp_mz != (subdomain.zE-subdomain.zS+1))
            {
                printf("Dimensions do not match. Filling using specified filling file\n");
                return 0;
            }
        }
        else if (strcmp(temp, "SPACING") == 0)
        {
            double temp_dx = 0, temp_dy = 0, temp_dz = 0;

            if(fscanf(fp, "%le", &temp_dx));
            if(fscanf(fp, "%le", &temp_dy));
            if(fscanf(fp, "%le", &temp_dz));

            if (temp_dx != simDomain.DELTA_X || temp_dy != simDomain.DELTA_Y || temp_dz != simDomain.DELTA_Z)
            {
                printf("Dimensions do not match. Filling using specified filling file\n");
                return 0;
            }
        }
        else if (strcmp(temp, "default") == 0 && k < simDomain.numComponents-1)
        {
           // printf("Reading %s\n", simDomain.componentNames[k]);
            if(fscanf(fp, "%le", &value));
            for (z = 0; z <= subdomain.zE - subdomain.zS; z++)
            {
                for (y = 0; y <= subdomain.yE - subdomain.yS; y++)
                {
                    for (x = 0; x <= subdomain.xE - subdomain.xS; x++)
                    {
                        if(fread(&value, sizeof(double), 1, fp));
                        if (IS_LITTLE_ENDIAN)
                            comp[k*subdomain.numCells + z*layer_size + y*simDomain.MESH_X + x] = swap_bytes(value);
                        else
                            comp[k*subdomain.numCells + z*layer_size + y*simDomain.MESH_X + x] = value;
                    }
                }
            }
            k++;
        }
        else if (strcmp(temp, "default") == 0 && a < simDomain.numPhases && k >= simDomain.numComponents-1)
        {
           // printf("Reading %s\n", simDomain.phaseNames[a]);
            if(fscanf(fp, "%le", &value));

            for (z = 0; z <= subdomain.zE - subdomain.zS; z++)
            {
                for (y = 0; y <= subdomain.yE - subdomain.yS; y++)
                {
                    for (x = 0; x <= subdomain.xE - subdomain.xS; x++)
                    {
                        if(fread(&value, sizeof(double), 1, fp));
                        if (IS_LITTLE_ENDIAN)
                            phi[a*subdomain.numCells + z*layer_size + y*simDomain.MESH_X + x] = swap_bytes(value);
                        else
                            phi[a*subdomain.numCells + z*layer_size + y*simDomain.MESH_X + x] = value;
                    }
                }
            }
            a++;
        }
        else if (a == simDomain.numPhases)
            break;
    }

    return 1;
}

void read_domain(double *phi, double *comp,
                 domainInfo simDomain, subdomainInfo subdomain,
                 int t, int rank, MPI_Comm comm,
                 char *argv[])
{
    FILE *fp;
    char name[1000], writeformat[10];

    int size;
    MPI_Comm_size(comm, &size);

    // Ensure numprocs = numfiles
    if (rank == size-1)
    {
        sprintf(name,"DATA/Processor_%d/%s_%07d.vtk", rank+1, argv[3], t);
        if (fp = fopen(name, "rb"))
        {
            printf("Found files that are not being read. Either move those files or change the number of processors\n");
            printf("File name in violation: %s\n", name);
            exit(-1);
        }
    }
    MPI_Barrier(comm);

    sprintf(name,"DATA/Processor_%d/%s_%07d.vtk", rank, argv[3], t);

    if (fp = fopen(name, "rb"))
    {
        if (rank == 0 || rank == size-1)
            printf("\nReading from %s\n", name);

        if(fscanf(fp, "%*[^\n]\n"));
        if(fscanf(fp, "%*[^\n]\n"));

        if(fscanf(fp, "%s", writeformat));

        if (strcmp(writeformat, "ASCII") == 0)
        {
            if (readVTK_ASCII(fp, phi, comp, simDomain, subdomain, t))
            {
                printf("Read %s successfully\n", name);
            }
            else
            {
                printf("Verify input and saved files\n");
                exit(-1);
            }
        }
        else if (strcmp(writeformat, "BINARY") == 0)
        {
            if (readVTK_BINARY(fp, phi, comp, simDomain, subdomain, t))
            {
                printf("Read %s successfully\n", name);
            }
            else
            {
                printf("Verify input and saved files\n");
                exit(-1);
            }
        }
        fclose(fp);
    }
    else
    {
        printf("\nCould not find %s.", name);
        exit(-1);
    }
}

// void writeHDF5(double *phi, double *comp,
//                domainInfo simDomain, subdomainInfo subdomain,
//                int t, int rank, MPI_Comm comm,
//                char *argv[])
// {
//     hid_t file_id;
//     hid_t plist_id;
//
//     char filename_hdf5[1000];
//
//     herr_t status_h;
//
//     //every processor creates a file collectively
//     //   sprintf(filename_hdf5, "DATA/%s_%d_%ld.h5", argv[3], numtasks, t);
//     sprintf(filename_hdf5, "DATA/%s_%d.h5", argv[3], t);
//     /* Set up file access property list with parallel I/O access*/
//     plist_id = H5Pcreate(H5P_FILE_ACCESS);
//     H5Pset_fapl_mpio(plist_id, comm, MPI_INFO_NULL);
//
//     file_id = H5Fcreate(filename_hdf5, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
//     H5Pclose(plist_id); // we'll use this plist_id again later
//
//     hsize_t dims[3];
//     hsize_t count[3];
//     hsize_t offset_slab[3];
//     int rank_hdf5 = 3;
//
//     int i, a, b, k, index, index_to;
//     int x, y, z;
//
//     double **buffer;
//     double composition;
//
//     int index_count;
//
//     hid_t dset_id; //handles
//     hid_t dataspace_id, memspace_id;
//
//     dims[2] = simDomain.MESH_Z;
//     dims[1] = simDomain.MESH_Y;
//     dims[0] = simDomain.MESH_X;
//
//     if (subdomain.sizeZ > 1)
//     {
//         count[2] = subdomain.sizeZ - 2*subdomain.padding;
//         count[1] = subdomain.sizeY;
//         count[0] = subdomain.sizeX;
//
//         index_count = count[2]*count[1]*count[0];
//     }
//     else
//     {
//         count[2] = subdomain.sizeZ;
//         count[1] = subdomain.sizeY - 2*subdomain.padding;
//         count[0] = subdomain.sizeX;
//
//         index_count = count[1]*count[0];
//     }
//
//     int size_fields = simDomain.numPhases + simDomain.numComponents-1;
//
//     char **coordNames   = (char**)malloc(sizeof(char*)*(size_fields));
//     char name[100];
//     i = 0;
//     for (k = 0; k < simDomain.numComponents-1; k++)
//     {
//         sprintf(name, "/%s", simDomain.componentNames[k]);
//         coordNames[i] = (char*)malloc(sizeof(char)*(strlen(name)+1));
//         strcpy(coordNames[i], name);
//         i++;
//     }
//     for (k = 0; k < simDomain.numPhases; k++)
//     {
//         sprintf(name, "/%s", simDomain.phaseNames[k]);
//         coordNames[i] = (char*)malloc(sizeof(char)*(strlen(name)+1));
//         strcpy(coordNames[i], name);
//         i++;
//     }
//
//     buffer = (double**)malloc(size_fields*sizeof(double*));
//     for (i = 0; i < size_fields; i++)
//         buffer[i] = (double*)malloc(index_count*sizeof(double));
//
//     for (z = 0; z < count[2]; z++)
//     {
//         for (y = 0; y < count[1]; y++)
//         {
//             for (x = 0; x < count[0]; x++)
//             {
//                 index    = (y + z*simDomain.MESH_Y)*simDomain.MESH_Z + x;
//                 index_to =  (y + z*count[1])*count[0] + x;
//
//                 for (a = 0; a < simDomain.numComponents-1; a++)
//                     buffer[a][index_to] = comp[a*simDomain.numCells + index];
//
//                 for (k = 0; k < simDomain.numPhases; k++)
//                     buffer[simDomain.numComponents-1+k][index_to] = phi[a*simDomain.numCells + index];
//             }
//         }
//     }
//
//     for (i = 0; i < size_fields; i++)
//     {
//         dataspace_id = H5Screate_simple(rank_hdf5, dims, NULL);
//         /*
//          * Create the dataset with default properties and close dataspace_id.
//          */
//         dset_id = H5Dcreate(file_id, coordNames[i], H5T_NATIVE_DOUBLE, dataspace_id,
//                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//         H5Sclose(dataspace_id);
//         //
//         memspace_id = H5Screate_simple(rank_hdf5, count, NULL); //count is of shape rank_hdf5
//         //
//         offset_slab[2] = subdomain.zS;
//         offset_slab[1] = subdomain.yS;
//         offset_slab[0] = subdomain.xS;
//         //
//         //     /* Select hyperslab in the file.
//         //             */
//         dataspace_id = H5Dget_space(dset_id); // filespace is the handle given by dset_id
//         //     // 			H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
//         //     // 			printf("dims[0] = %ld, dims[1] = %ld\n", (long)dims[0], (long)dims[1]);
//         H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset_slab, NULL, count, NULL);
//         //
//         //     /* Create property list for collective dataset write.
//         //     */
//         plist_id = H5Pcreate(H5P_DATASET_XFER);
//         H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
//         //
//         status_h = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id,
//                             plist_id, buffer[i]);
//         //
//         status_h = H5Dclose(dset_id);
//         status_h = H5Sclose(dataspace_id);
//         status_h = H5Sclose(memspace_id);
//         status_h = H5Pclose(plist_id);
//     }
//
//     for (i = 0; i < size_fields; i++)
//     {
//         free(buffer[i]);
//         free(coordNames[i]);
//     }
//
//     free(buffer);
//     free(coordNames);
//
//     status_h = H5Fclose(file_id);
// }
//
// void readHDF5(double *phi, double *comp,
//               domainInfo simDomain, subdomainInfo subdomain,
//               int t, int rank, MPI_Comm comm,
//               char *argv[])
// {
//     hid_t file_id;
//     hid_t plist_id;
//
//     char filename_hdf5[1000];
//
//     herr_t status_h;
//
//     //every processor creates a file collectively
//     //   sprintf(filename_hdf5, "DATA/%s_%ld_%ld.h5", argv[3], numworkers, t);
//     sprintf(filename_hdf5, "DATA/%s_%d.h5", argv[3], t);
//     /* Set up file access property list with parallel I/O access*/
//     plist_id = H5Pcreate(H5P_FILE_ACCESS);
//     H5Pset_fapl_mpio(plist_id, comm, MPI_INFO_NULL);
//
//     //   file_id = H5Fcreate(filename_hdf5, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
//     file_id = H5Fopen(filename_hdf5, H5F_ACC_RDONLY, plist_id);
//     H5Pclose(plist_id); // we'll use this plist_id again later
//
//     hsize_t dims[3];
//     hsize_t count[3];
//     hsize_t offset_slab[3];
//     int rank_hdf5 = 3;
//
//     long i, a, b, k, index, index_to;
//     long x, y, z;
//
//     double **buffer;
//     double composition;
//
//     hid_t dset_id; //handles
//     hid_t dataspace_id, memspace_id;
//
//     dims[2] = simDomain.MESH_Z;
//     dims[1] = simDomain.MESH_Y;
//     dims[0] = simDomain.MESH_X;
//
//     int index_count;
//
//     if (subdomain.sizeZ > 1)
//     {
//         count[2] = subdomain.sizeZ - 2*subdomain.padding;
//         count[1] = subdomain.sizeY;
//         count[0] = subdomain.sizeX;
//
//         index_count = count[2]*count[1]*count[0];
//     }
//     else
//     {
//         count[2] = subdomain.sizeZ;
//         count[1] = subdomain.sizeY - 2*subdomain.padding;
//         count[0] = subdomain.sizeX;
//
//         index_count = count[1]*count[0];
//     }
//
//     int size_fields = simDomain.numPhases + simDomain.numComponents-1;
//
//     char **coordNames   = (char**)malloc(sizeof(char*)*(size_fields));
//     char name[100];
//     i = 0;
//     for (k = 0; k < simDomain.numComponents-1; k++)
//     {
//         sprintf(name, "/%s", simDomain.componentNames[k]);
//         coordNames[i] = (char*)malloc(sizeof(char)*(strlen(name)+1));
//         strcpy(coordNames[i], name);
//         i++;
//     }
//     for (k = 0; k < simDomain.numPhases; k++)
//     {
//         sprintf(name, "/%s", simDomain.phaseNames[k]);
//         coordNames[i] = (char*)malloc(sizeof(char)*(strlen(name)+1));
//         strcpy(coordNames[i], name);
//         i++;
//     }
//
//     buffer = (double **)malloc(size_fields*sizeof(double*));
//     for (i=0; i < size_fields; i++) {
//         buffer[i] = (double *)malloc(index_count*sizeof(double));
//     }
//
//     for (i = 0; i < size_fields; i++)
//     {
//         /* open the dset_id collectively */
//         dset_id = H5Dopen(file_id, coordNames[i], H5P_DEFAULT);
//         /* create a file dataspace independently */
//         dataspace_id = H5Dget_space(dset_id);
//
//         offset_slab[2] = subdomain.zS;
//         offset_slab[1] = subdomain.yS;
//         offset_slab[0] = subdomain.xS;
//
//         H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset_slab, NULL, count, NULL);
//         /* create a memory dataspace independently */
//         memspace_id = H5Screate_simple(rank_hdf5, count, NULL);
//
//         //     /* Create property list for collective dataset write.
//         //     */
//         plist_id = H5Pcreate(H5P_DATASET_XFER);
//         H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
//         //
//         status_h = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id,
//                            plist_id, buffer[i]);
//         //     status_h = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id,	plist_id, vec[i]);
//         status_h = H5Dclose(dset_id);
//         status_h = H5Sclose(dataspace_id);
//         status_h = H5Sclose(memspace_id);
//         status_h = H5Pclose(plist_id);
//     }
//
//     for (z = 0; z < count[2]; z++)
//     {
//         for (y = 0; y < count[1]; y++)
//         {
//             for (x = 0; x < count[0]; x++)
//             {
//                 index    = (y + z*simDomain.MESH_Y)*simDomain.MESH_Z + x;
//                 index_to =  (y + z*count[1])*count[0] + x;
//
//                 for (a = 0; a < simDomain.numComponents-1; a++)
//                     comp[a*simDomain.numCells + index] = buffer[a][index_to];
//
//                 for (k = 0; k < simDomain.numPhases; k++)
//                     phi[a*simDomain.numCells + index] = buffer[simDomain.numComponents-1+k][index_to];
//             }
//         }
//     }
//
//     for (i = 0; i < size_fields; i++)
//     {
//         free(buffer[i]);
//         free(coordNames[i]);
//     }
//
//     free(buffer);
//     free(coordNames);
//
//     status_h = H5Fclose(file_id);
// }
