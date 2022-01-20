#ifndef FILE_WRITER_H_
#define FILE_WRITER_H_

void writetofile_serial2D(struct fields* gridinfo, char *argv[], long t)
{
    FILE *fp;
    char name[1000];
    sprintf(name, "DATA/%s_%07ld.vtk", argv[3], t);
    fp = fopen(name, "w");
    write_cells_vtk_2D(fp, gridinfo);
    fclose(fp);
}

void writetofile_serial2D_binary(struct fields* gridinfo, char *argv[], long t)
{
    FILE *fp;
    char name[1000];
    sprintf(name,"DATA/%s_%07ld.vtk", argv[3], t);
    fp = fopen(name, "wb");
    write_cells_vtk_2D_binary(fp, gridinfo);
    fclose(fp);
}

void write_cells_vtk_2D(FILE *fp, struct fields *gridinfo)
{
    long index;
    long a, k;
    long x, y, z;

    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "Microsim_fields\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_POINTS\n");
    fprintf(fp, "DIMENSIONS %ld %ld %ld\n", MESH_Y, MESH_X, MESH_Z);
    fprintf(fp, "ORIGIN 0 0 0\n");
    fprintf(fp, "SPACING %le %le %le\n", DELTA_X, DELTA_Y, DELTA_Z);
    fprintf(fp, "POINT_DATA %ld\n", (long)MESH_X*(long)MESH_Y*(long)MESH_Z);

    for (a = 1; a < NUMPHASES; a++) {
        fprintf(fp,"SCALARS %s double 1\n", PHASES[a]);
        fprintf(fp,"LOOKUP_TABLE default\n");
        for (x = start[X]; x <= end[X]; x++) {
            for (y = start[Y]; y <= end[Y]; y++) {
                for (z = start[Z]; z <= end[Z]; z++) {
                    index = x*layer_size + y*rows_z + z;
                    fprintf(fp, "%le\n",gridinfo[index].phia[a]);
                }
            }
        }
        fprintf(fp,"\n");
    }

    for (k = 0; k < NUMCOMPONENTS-1; k++) {
        fprintf(fp, "SCALARS %s double 1\n", COMPONENTS[k]);
        fprintf(fp, "LOOKUP_TABLE default\n");
        for (x = start[X]; x <= end[X]; x++) {
            for (y = start[Y]; y <= end[Y]; y++) {
                for (z = start[Z]; z <= end[Z]; z++) {
                    index = x*layer_size + y*rows_z + z;
                    fprintf(fp, "%le\n", gridinfo[index].compi[k]);
                }
            }
        }
        fprintf(fp, "\n");
    }
}

void write_cells_vtk_2D_binary(FILE *fp, struct fields *gridinfo) {
    long x, y, z, index;
    long a;
    long k;
    double value;

    fprintf(fp,"# vtk DataFile Version 3.0\n");
    fprintf(fp,"Microsim_fields\n");
    fprintf(fp,"BINARY\n");
    fprintf(fp,"DATASET STRUCTURED_POINTS\n");
    fprintf(fp,"DIMENSIONS %ld %ld %ld\n", MESH_Y, MESH_X, MESH_Z);
    fprintf(fp,"ORIGIN 0 0 0\n");
    fprintf(fp,"SPACING %le %le %le\n", DELTA_X, DELTA_Y, DELTA_Z);
    fprintf(fp,"POINT_DATA %ld\n",(long)MESH_X*(long)MESH_Y*(long)MESH_Z);

    for (a = 1; a < NUMPHASES; a++) {
        fprintf(fp,"SCALARS %s double 1\n", PHASES[a]);
        fprintf(fp,"LOOKUP_TABLE default\n");
        for (x = start[X]; x <= end[X]; x++) {
            for (y = start[Y]; y <= end[Y]; y++) {
                for (z = start[Z]; z <= end[Z]; z++) {
                    index = x*layer_size + y*rows_z + z;
                    if (IS_LITTLE_ENDIAN) {
                        value = swap_bytes(gridinfo[index].phia[a]);
                    } else {
                        value = gridinfo[index].phia[a];
                    }
                    fwrite(&value, sizeof(double), 1, fp);
                }
            }
        }
        fprintf(fp,"\n");
    }
    for (k = 0; k < NUMCOMPONENTS-1; k++) {
        fprintf(fp,"SCALARS %s double 1\n",COMPONENTS[k]);
        fprintf(fp,"LOOKUP_TABLE default\n");
        for (x = start[X]; x <= end[X]; x++) {
            for (y = start[Y]; y <= end[Y]; y++) {
                for (z = start[Z]; z <= end[Z]; z++) {
                    index = x*layer_size + y*rows_z + z;
                    if (IS_LITTLE_ENDIAN) {
                        value = swap_bytes(gridinfo[index].compi[k]);
                    } else {
                        value = gridinfo[index].compi[k];
                    }
                    fwrite(&value, sizeof(double), 1, fp);
                }
            }
        }
        fprintf(fp,"\n");
    }
}

#endif
