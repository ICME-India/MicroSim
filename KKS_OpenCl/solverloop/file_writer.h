#ifndef FILE_WRITER_H_
#define FILE_WRITER_H_

#include <arpa/inet.h>
#include <endian.h>
#include <stdint.h>

#define  IS_BIG_ENDIAN     (1 == htons(1))
#define  IS_LITTLE_ENDIAN  (!IS_BIG_ENDIAN)

double swap_bytes(double value) {
  double  src_num = value;
  int64_t tmp_num = htobe64(le64toh(*(int64_t*)&src_num));
  double  dst_num = *(double*)&tmp_num;
  return dst_num;
}

void writetofile_serial2D(struct fields* gridinfo, char *argv[], long t);
void writetofile_serial2D_binary(struct fields* gridinfo, char *argv[], long t);
void write_cells_vtk_2D(FILE *fp, struct fields *gridinfo);
void write_cells_vtk_2D_binary(FILE *fp, struct fields *gridinfo);
void writetofile_mpi(struct fields* gridinfo, char *argv[], long t);
void write_cells_vtk_mpi(FILE *fp, struct fields* gridinfo);
void writetofile_mpi_binary(struct fields* gridinfo, char *argv[], long t);
void write_cells_vtk_mpi_binary(FILE *fp, struct fields* gridinfo);

void writetofile_serial2D(struct fields* gridinfo, char *argv[], long t) {
  long x,y,z;
  long gidy;
  FILE *fp;
  char name[1000];
  double composition;
  long b, k;
  sprintf(name,"DATA/%s_%ld.vtk",argv[3], t);
  fp=fopen(name,"w");
  write_cells_vtk_2D(fp, gridinfo);
  fclose(fp);
}
void writetofile_serial2D_binary(struct fields* gridinfo, char *argv[], long t) {
  long x,y,z;
  long gidy;
  FILE *fp;
  char name[1000];
  double composition;
  long b, k;
  sprintf(name,"DATA/%s_%ld.vtk",argv[3], t);
  fp=fopen(name,"wb");
  write_cells_vtk_2D_binary(fp, gridinfo);
  fclose(fp);
}
void readfromfile_serial2D(struct fields* gridinfo, char *argv[], long t) {
  long x,y,z;
  long gidy;
  FILE *fp;
  char name[1000];
  double composition;
  long b, k;
  sprintf(name,"DATA/%s_%ld.vtk",argv[3], t);
  fp=fopen(name,"r");
  read_cells_vtk_2D(fp, gridinfo);
  fclose(fp);
}
void readfromfile_serial2D_binary(struct fields* gridinfo, char *argv[], long t) {
  long x,y,z;
  long gidy;
  FILE *fp;
  char name[1000];
  double composition;
  long b, k;
  sprintf(name,"DATA/%s_%ld.vtk",argv[3], t);
  fp=fopen(name,"rb");
  read_cells_vtk_2D_binary(fp, gridinfo);
  fclose(fp);
}
void readfromfile_serialmpi(struct fields* gridinfo, char *argv[], long t) {
  long x,y,z;
  long gidy;
  FILE *fp;
  char name[1000];
  double composition;
  long b, k;
  sprintf(name,"DATA/Processor_%d/%s_%ld.vtk",rank, argv[3], t);
  fp=fopen(name,"rb");
  read_cells_vtk_mpi(fp, gridinfo);
  fclose(fp);
}
void readfromfile_serialmpi_binary(struct fields* gridinfo, char *argv[], long t) {
  long x,y,z;
  long gidy;
  FILE *fp;
  char name[1000];
  double composition;
  long b, k;
  sprintf(name,"DATA/Processor_%d/%s_%ld.vtk",rank, argv[3], t);
  fp=fopen(name,"rb");
  read_cells_vtk_mpi_binary(fp, gridinfo);
  fclose(fp);
}
void write_cells_vtk_2D(FILE *fp, struct fields *gridinfo) {
  long x, y, z, index;
  long a, b;
  long k;
  double composition;
  fprintf(fp,"# vtk DataFile Version 3.0\n");
  fprintf(fp,"Microsim_fields\n");
  fprintf(fp,"ASCII\n");
  fprintf(fp,"DATASET STRUCTURED_POINTS\n");
  fprintf(fp,"DIMENSIONS %ld %ld %ld\n",MESH_Y, MESH_X, MESH_Z);
  fprintf(fp,"ORIGIN 0 0 0\n");
  fprintf(fp,"SPACING %le %le %le\n",deltax, deltay, deltaz);
  fprintf(fp,"POINT_DATA %ld\n",(long)MESH_X*(long)MESH_Y*(long)MESH_Z);
//   fprintf(fp,"FIELD PHASEFIELD 2\n");
//   fprintf(fp,"PHI[0] 1 %ld double\n", (long)Nx*(long)Ny);
  for (a=0; a < NUMPHASES; a++) {
    fprintf(fp,"SCALARS %s double 1\n",Phases[a]);
    fprintf(fp,"LOOKUP_TABLE default\n");
    for (x=start[X]; x<=end[X]; x++) {
      for (z=start[Z]; z <= end[Z]; z++) {
        for (y=start[Y]; y <= end[Y]; y++) {
          index = x*layer_size + z*rows_y + y;
          fprintf(fp, "%le\n",gridinfo[index].phia[a]);
        }
      }
    }
    fprintf(fp,"\n");
  }

  for (k=0; k < NUMCOMPONENTS-1; k++) {
    fprintf(fp,"SCALARS Mu_%s double 1\n",Components[k]);
    fprintf(fp,"LOOKUP_TABLE default\n");
    for (x=start[X]; x<=end[X]; x++) {
      for (z=start[Z]; z <= end[Z]; z++) {
        for (y=start[Y]; y <= end[Y]; y++) {
          index = x*layer_size + z*rows_y + y;
          fprintf(fp, "%le\n",gridinfo[index].compi[k]);
        }
      }
    }
    fprintf(fp,"\n");
  }

  for (k=0; k < NUMCOMPONENTS-1; k++) {
    fprintf(fp,"SCALARS comp_%s double 1\n",Components[k]);
    fprintf(fp,"LOOKUP_TABLE default\n");
    for (x=start[X]; x<=end[X]; x++) {
      for (z=start[Z]; z <= end[Z]; z++) {
        for (y=start[Y]; y <= end[Y]; y++) {
          index = x*layer_size + z*rows_y + y;
          fprintf(fp, "%le\n",gridinfo[index].composition[k]);
        }
      }
    }
    fprintf(fp,"\n");
  }
  if (!ISOTHERMAL) {
    fprintf(fp,"SCALARS Temperature double 1\n");
    fprintf(fp,"LOOKUP_TABLE default\n");
    for (x=start[X]; x<=end[X]; x++) {
      for (z=start[Z]; z <= end[Z]; z++) {
        for (y=start[Y]; y <= end[Y]; y++) {
          index = x*layer_size + z*rows_y + y;
          fprintf(fp, "%le \n",gridinfo[index].temperature);
        }
      }
    }
  }
}
void write_cells_vtk_2D_binary(FILE *fp, struct fields *gridinfo) {
  long x, y, z, index;
  long a, b;
  long k;
  double composition;
  char first_lines[50];
  int length;
  double value;
  
  fprintf(fp,"# vtk DataFile Version 3.0\n");
  fprintf(fp,"Microsim_fields\n");
  fprintf(fp,"BINARY\n");
  fprintf(fp,"DATASET STRUCTURED_POINTS\n");
  fprintf(fp,"DIMENSIONS %ld %ld %ld\n",MESH_Y, MESH_X, (long)1);
  fprintf(fp,"ORIGIN 0 0 0\n");
  fprintf(fp,"SPACING %le %le %le\n",deltax, deltay, deltaz);
  fprintf(fp,"POINT_DATA %ld\n",(long)MESH_X*(long)MESH_Y);
  
  for (a=0; a < NUMPHASES; a++) {
    fprintf(fp,"SCALARS %s double 1\n",Phases[a]);
    fprintf(fp,"LOOKUP_TABLE default\n");
    for (x=start[X]; x<=end[X]; x++) {
      for (z=start[Z]; z <= end[Z]; z++) {
        for (y=start[Y]; y <= end[Y]; y++) {
          index = x*layer_size + z*rows_y + y;
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
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    fprintf(fp,"SCALARS Mu_%s double 1\n",Components[k]);
    fprintf(fp,"LOOKUP_TABLE default\n");
    for (x=start[X]; x<=end[X]; x++) {
      for (z=start[Z]; z <= end[Z]; z++) {
        for (y=start[Y]; y <= end[Y]; y++) {
          index = x*layer_size + z*rows_y + y;
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
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    fprintf(fp,"SCALARS Composition_%s double 1\n",Components[k]);
    fprintf(fp,"LOOKUP_TABLE default\n");
    for (x=start[X]; x<=end[X]; x++) {
      for (z=start[Z]; z <= end[Z]; z++) {
        for (y=start[Y]; y <= end[Y]; y++) {
          index = x*layer_size + z*rows_y + y;
          if (IS_LITTLE_ENDIAN) {
            value = swap_bytes(gridinfo[index].compi[k]);
          } else {
            value = gridinfo[index].composition[k];
          }
          fwrite(&value, sizeof(double), 1, fp);
        }
      }
    }
   fprintf(fp,"\n");
  }
  if (!ISOTHERMAL) {
    fprintf(fp,"SCALARS Temperature double 1\n");
    fprintf(fp,"LOOKUP_TABLE default\n");
    for (x=start[X]; x<=end[X]; x++) {
      for (z=start[Z]; z <= end[Z]; z++) {
        for (y=start[Y]; y <= end[Y]; y++) {
          index = x*layer_size + z*rows_y + y;
          if(IS_LITTLE_ENDIAN) {
            value = swap_bytes(gridinfo[index].temperature);
          } else {
            value = gridinfo[index].temperature;
          }
          fwrite(&value, sizeof(double), 1, fp);
        }
      }
    }
  }
}
void read_cells_vtk_2D(FILE *fp, struct fields *gridinfo) {
  long x, y, z, index;
  long a, b;
  long k;
  char name[1000];
  long mesh_x, mesh_y, mesh_z;
  long ox, oy, oz;
  long total_points;
  double dx, dy, dz;
  double composition;
  long size;
    
  fscanf(fp,"%*[^\n]\n");
  fscanf(fp,"%*[^\n]\n");
  fscanf(fp,"%*[^\n]\n");
  fscanf(fp,"%*[^\n]\n");
  fscanf(fp,"%*[^\n]\n");
  fscanf(fp,"%*[^\n]\n");
  fscanf(fp,"%*[^\n]\n");
  fscanf(fp,"%*[^\n]\n");

  for (a=0; a < NUMPHASES; a++) {
    fscanf(fp,"%*[^\n]\n");
    fscanf(fp,"%*[^\n]\n");
    for (x=start[X]; x<=end[X]; x++) {
      for (z=start[Z]; z <= end[Z]; z++) {
        for (y=start[Y]; y <= end[Y]; y++) {
          index = x*layer_size + z*rows_y + y;
          fscanf(fp, "%le \n",&gridinfo[index].phia[a]);
        }
      }
    }
  }
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    fscanf(fp,"%*[^\n]\n");
    fscanf(fp,"%*[^\n]\n");
    for (x=start[X]; x<=end[X]; x++) {
      for (z=start[Z]; z <= end[Z]; z++) {
        for (y=start[Y]; y <= end[Y]; y++) {
          index = x*layer_size + z*rows_y + y;
          fscanf(fp, "%le \n",&gridinfo[index].compi[k]);
        }
      }
    }
  }
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    fscanf(fp,"%*[^\n]\n");
    fscanf(fp,"%*[^\n]\n");
    for (x=start[X]; x<=end[X]; x++) {
      for (z=start[Z]; z <= end[Z]; z++) {
        for (y=start[Y]; y <= end[Y]; y++) {
          index = x*layer_size + z*rows_y + y;
          fscanf(fp, "%le \n",&gridinfo[index].composition[k]);
        }
      }
    }
  }
  if (!ISOTHERMAL) {
    fscanf(fp,"%*[^\n]\n");
    fscanf(fp,"%*[^\n]\n");
    for (x=start[X]; x<=end[X]; x++) {
      for (z=start[Z]; z <= end[Z]; z++) {
        for (y=start[Y]; y <= end[Y]; y++) {
          index = x*layer_size + z*rows_y + y;
          fscanf(fp, "%le \n",&gridinfo[index].temperature);
        }
      }
    }
  }
}
void read_cells_vtk_2D_binary(FILE *fp, struct fields *gridinfo) {
  long x, y, z, index;
  long a, b;
  long k;
  char name[1000];
  long mesh_x, mesh_y, mesh_z;
  long ox, oy, oz;
  long total_points;
  double dx, dy, dz;
  double composition;
  long size;
  double value;
  
  fscanf(fp,"%*[^\n]\n");
  fscanf(fp,"%*[^\n]\n");
  fscanf(fp,"%*[^\n]\n");
  fscanf(fp,"%*[^\n]\n");
  fscanf(fp,"%*[^\n]\n");
  fscanf(fp,"%*[^\n]\n");
  fscanf(fp,"%*[^\n]\n");
  fscanf(fp,"%*[^\n]\n");
  
  
  for (a=0; a < NUMPHASES; a++) {
    fscanf(fp,"%*[^\n]\n");
    fscanf(fp,"%*[^\n]\n");
    for (x=start[X]; x<=end[X]; x++) {
      for (z=start[Z]; z <= end[Z]; z++) {
        for (y=start[Y]; y <= end[Y]; y++) {
          index = x*layer_size + z*rows_y + y;
          fread(&value, sizeof(double), 1, fp);
          if (IS_LITTLE_ENDIAN) {
            gridinfo[index].phia[a] = swap_bytes(value);
          } else {
            gridinfo[index].phia[a] = value;
          }
        }
      }
    }
  }
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    fscanf(fp,"%*[^\n]\n");
    fscanf(fp,"%*[^\n]\n");
    for (x=start[X]; x<=end[X]; x++) {
      for (z=start[Z]; z <= end[Z]; z++) {
        for (y=start[Y]; y <= end[Y]; y++) {
          index = x*layer_size + z*rows_y + y;
          fread(&value, sizeof(double), 1, fp);
          if (IS_LITTLE_ENDIAN) {
            gridinfo[index].compi[k] = swap_bytes(value);
          } else {
            gridinfo[index].compi[k] = value;
          }
        }
      }
    }
  }
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    fscanf(fp,"%*[^\n]\n");
    fscanf(fp,"%*[^\n]\n");
    for (x=start[X]; x<=end[X]; x++) {
      for (z=start[Z]; z <= end[Z]; z++) {
        for (y=start[Y]; y <= end[Y]; y++) {
          index = x*layer_size + z*rows_y + y;
          fread(&value, sizeof(double), 1, fp);
          if (IS_LITTLE_ENDIAN) {
            gridinfo[index].composition[k] = swap_bytes(value);
          } else {
            gridinfo[index].composition[k] = value;
          }
        }
      }
    }
  }
  if (!ISOTHERMAL) {
    fscanf(fp,"%*[^\n]\n");
    fscanf(fp,"%*[^\n]\n");
    for (x=start[X]; x<=end[X]; x++) {
      for (z=start[Z]; z <= end[Z]; z++) {
        for (y=start[Y]; y <= end[Y]; y++) {
          index = x*layer_size + z*rows_y + y;
          fread(&value, sizeof(double), 1, fp);
          if (IS_LITTLE_ENDIAN) {
            gridinfo[index].temperature = swap_bytes(value);
          } else {
            gridinfo[index].temperature = value;
          }
        }
      }
    }
  }
}
void read_cells_vtk_mpi(FILE *fp, struct fields *gridinfo) {
  long x, y, z, index;
  long a, b;
  long k;
  char name[1000];
  long mesh_x, mesh_y, mesh_z;
  long ox, oy, oz;
  long total_points;
  double dx, dy, dz;
  double composition;
  long size,dim;

  fscanf(fp,"%*[^\n]\n");
  fscanf(fp,"%*[^\n]\n");
  fscanf(fp,"%*[^\n]\n");

  for (a=0; a < NUMPHASES; a++) {
    for (x=1; x < mpiparam.rows_x-1; x++) {
      for (z=1; z < mpiparam.rows_z-1; z++) {
        for (y=1; y < mpiparam.rows_y-1; y++) {
          index        = y + mpiparam.rows_y * ( z + mpiparam.rows_z * x);
          fscanf(fp, "%le \n",&gridinfo[index].phia[a]);
        }
      }
    }
    fprintf(fp,"\n");
  }
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    for (x=1; x < mpiparam.rows_x-1; x++) {
      for (z=1; z < mpiparam.rows_z-1; z++) {
        for (y=1; y < mpiparam.rows_y-1; y++) {
          index        = y + mpiparam.rows_y * ( z + mpiparam.rows_z * x);
          fscanf(fp, "%le \n",&gridinfo[index].compi[k]);
        }
      }
    }
    fprintf(fp,"\n");
  }
//   if(WRITECOMPOSITION) {
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      for (x=1; x < mpiparam.rows_x-1; x++) {
        for (z=1; z < mpiparam.rows_z-1; z++) {
          for (y=1; y < mpiparam.rows_y-1; y++) {
            index        = y + mpiparam.rows_y * ( z + mpiparam.rows_z * x);
            fscanf(fp, "%le \n",&gridinfo[index].composition[k]);
          }
        }
      }
      fprintf(fp,"\n");
    }
//   }
  if (ELASTICITY) {
    for (dim=0; dim < DIMENSION; dim++) {
      for (x=1; x < mpiparam.rows_x-1; x++) {
        for (z=1; z < mpiparam.rows_z-1; z++) {
          for (y=1; y < mpiparam.rows_y-1; y++) {
            index        = y + mpiparam.rows_y * ( z + mpiparam.rows_z * x);
            fscanf(fp, "%le \n",&iter_gridinfom[index].disp[dim][2]);
          }
        }
      }
      fprintf(fp,"\n");
    }
  }

  if (!ISOTHERMAL) {
    for (x=1; x < mpiparam.rows_x-1; x++) {
      for (z=1; z < mpiparam.rows_z-1; z++) {
        for (y=1; y < mpiparam.rows_y-1; y++) {
          index        = y + mpiparam.rows_y * ( z + mpiparam.rows_z * x);
          fscanf(fp, "%le \n",&gridinfo[index].temperature);
        }
      }
    }
  }

}
void read_cells_vtk_mpi_binary(FILE *fp, struct fields *gridinfo) {
  long x, y, z, index;
  long a, b;
  long k;
  char name[1000];
  long mesh_x, mesh_y, mesh_z;
  long ox, oy, oz;
  long total_points;
  double dx, dy, dz;
  double composition;
  long size, dim;
  double value;

  fscanf(fp,"%*[^\n]\n");
  fscanf(fp,"%*[^\n]\n");
  fscanf(fp,"%*[^\n]\n");

  for (a=0; a < NUMPHASES; a++) {
    for (x=1; x < mpiparam.rows_x-1; x++) {
      for (z=1; z < mpiparam.rows_z-1; z++) {
        for (y=1; y < mpiparam.rows_y-1; y++) {
          index        = y + mpiparam.rows_y * ( z + mpiparam.rows_z * x);
          fread(&value, sizeof(double), 1, fp);
          if (IS_LITTLE_ENDIAN) {
            gridinfo[index].phia[a] = swap_bytes(value);
          } else {
            gridinfo[index].phia[a] = value;
          }
        }
      }
    }
    fprintf(fp,"\n");
  }
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    for (x=1; x < mpiparam.rows_x-1; x++) {
      for (z=1; z < mpiparam.rows_z-1; z++) {
        for (y=1; y < mpiparam.rows_y-1; y++) {
          index        = y + mpiparam.rows_y * ( z + mpiparam.rows_z * x);
          fread(&value, sizeof(double), 1, fp);
          if (IS_LITTLE_ENDIAN) {
            gridinfo[index].composition[k] = swap_bytes(value);
          } else {
            gridinfo[index].composition[k] = value;
          }
        }
      }
    }
    fprintf(fp,"\n");
  }
//   if(WRITECOMPOSITION) {
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      for (x=1; x < mpiparam.rows_x-1; x++) {
        for (z=1; z < mpiparam.rows_z-1; z++) {
          for (y=1; y < mpiparam.rows_y-1; y++) {
            index        = y + mpiparam.rows_y * ( z + mpiparam.rows_z * x);
          if (IS_LITTLE_ENDIAN) {
            value = swap_bytes(gridinfo[index].compi[k]);
          } else {
            value = gridinfo[index].composition[k];
          }
          fwrite(&value, sizeof(double), 1, fp);
          }
        }
      }
      fprintf(fp,"\n");
    }
//   }
  if (ELASTICITY) {
    for (dim=0; dim < DIMENSION; dim++) {
      for (x=1; x < mpiparam.rows_x-1; x++) {
        for (z=1; z < mpiparam.rows_z-1; z++) {
          for (y=1; y < mpiparam.rows_y-1; y++) {
            index        = y + mpiparam.rows_y * ( z + mpiparam.rows_z * x);

          fread(&value, sizeof(double), 1, fp);
          if (IS_LITTLE_ENDIAN) {
            iter_gridinfom[index].disp[dim][2] = swap_bytes(value);
          } else {
            iter_gridinfom[index].disp[dim][2] = value;
          }

          }
        }
      }
      fprintf(fp,"\n");
    }
  }

  if (!ISOTHERMAL) {
    for (x=1; x < mpiparam.rows_x-1; x++) {
      for (z=1; z < mpiparam.rows_z-1; z++) {
        for (y=1; y < mpiparam.rows_y-1; y++) {
          index        = y + mpiparam.rows_y * ( z + mpiparam.rows_z * x);
          fread(&value, sizeof(double), 1, fp);
          if (IS_LITTLE_ENDIAN) {
            gridinfo[index].temperature = swap_bytes(value);
          } else {
            gridinfo[index].temperature = value;
          }
        }
      }
    }
  }
}

void populate_table_names(){
  long i, a, b, k;
  char chempot_name[100];
  char composition_name[100]; 
  char phase_name[100];
  
  size_fields = NUMPHASES + (NUMCOMPONENTS-1);
  
//   if (WRITECOMPOSITION) {
    size_fields += (NUMCOMPONENTS-1);
//   }
  if(!ISOTHERMAL) {
    size_fields += 1;
  }
  
  coordNames   = (char**)malloc(sizeof(char*)*(size_fields));
  i=0;
  for (a = 0; a < NUMPHASES; a++) {
    sprintf(phase_name, "/%s",Phases[a]);
    coordNames[i] = (char*)malloc(sizeof(char)*strlen(phase_name));
    strcpy(coordNames[i], phase_name);
    i++;
  }
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    sprintf(chempot_name, "/Mu_%s",Components[k]);
    coordNames[i] = (char*)malloc(sizeof(char)*strlen(chempot_name));
    strcpy(coordNames[i], chempot_name);
    i++;
  }
//   if (WRITECOMPOSITION) {
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      sprintf(composition_name, "/Composition_%s",Components[k]);
      coordNames[i] = (char*)malloc(sizeof(char)*strlen(composition_name));
      strcpy(coordNames[i], composition_name);
      i++;
    }
//   }
  if (!ISOTHERMAL) {
    coordNames[i] = (char*)malloc(sizeof(char)*(strlen("/T")+1));
    strcpy(coordNames[i], "/T");
  }
}

void writetofile_mpi(struct fields* gridinfo, char *argv[], long t) {
  long x,y,z;
  long gidy;
  FILE *fp;
  char name[1000];
  double composition;
  long b, k;
  sprintf(name,"DATA/Processor_%d/%s_%ld.vtk",rank, argv[3], t);
  fp=fopen(name,"w");
  if(fp == NULL) {
    printf("file not found inside  DATA directory: %s", name);
  }
  write_cells_vtk_mpi(fp, gridinfo);
  fclose(fp);

}
void write_cells_vtk_mpi(FILE *fp, struct fields* gridinfo) {
  long x, y, z, index;
  long a, b;
  long k;
  long global_index;
  long dim;
  FILE *fpr;
  char name[1000];

  //printf("%ld\n", mpiparam.rows_x);
  //printf("%ld\n", mpiparam.rows_y);
  //printf("%ld\n", mpiparam.rows_z);

  fprintf(fp,"%ld\n", mpiparam.rows_x-2);
  fprintf(fp,"%ld\n", mpiparam.rows_y-2);
  fprintf(fp,"%ld\n", mpiparam.rows_z-2);

  for (a=0; a < NUMPHASES; a++) {
    for (x=1; x < mpiparam.rows_x-1; x++) {
      for (z=1; z < mpiparam.rows_z-1; z++) {
        for (y=1; y < mpiparam.rows_y-1; y++) {
          index        = y + mpiparam.rows_y * ( z + mpiparam.rows_z * x);
          fprintf(fp, "%le\n", gridinfo[index].phia[a]);
        }
      }
    }
    fprintf(fp,"\n");
  }
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    for (x=1; x < mpiparam.rows_x-1; x++) {
      for (z=1; z < mpiparam.rows_z-1; z++) {
        for (y=1; y < mpiparam.rows_y-1; y++) {
          index        = y + mpiparam.rows_y * ( z + mpiparam.rows_z * x);
          fprintf(fp, "%le\n", gridinfo[index].compi[k]);
        }
      }
    }
    fprintf(fp,"\n");
  }
//   if(WRITECOMPOSITION) {
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      for (x=1; x < mpiparam.rows_x-1; x++) {
        for (z=1; z < mpiparam.rows_z-1; z++) {
          for (y=1; y < mpiparam.rows_y-1; y++) {
            index        = y + mpiparam.rows_y * ( z + mpiparam.rows_z * x);
            fprintf(fp, "%le\n", gridinfo[index].composition[k]);
          }
        }
      }
      fprintf(fp,"\n");
    }
//   }
  if (ELASTICITY) {
    for (dim=0; dim < DIMENSION; dim++) {
      for (x=1; x < mpiparam.rows_x-1; x++) {
        for (z=1; z < mpiparam.rows_z-1; z++) {
          for (y=1; y < mpiparam.rows_y-1; y++) {
            index        = y + mpiparam.rows_y * ( z + mpiparam.rows_z * x);
            fprintf(fp, "%le\n", iter_gridinfom[index].disp[dim][2]);
          }
        }
      }
      fprintf(fp,"\n");
    }
  }

  if (!ISOTHERMAL) {
    for (x=1; x < mpiparam.rows_x-1; x++) {
      for (z=1; z < mpiparam.rows_z-1; z++) {
        for (y=1; y < mpiparam.rows_y-1; y++) {
          index        = y + mpiparam.rows_y * ( z + mpiparam.rows_z * x);
          fprintf(fp, "%le\n", gridinfo[index].temperature);
        }
      }
    }
  }

//   sprintf(name,"DATA/Processor_%d/u2_%ld.dat",rank, t);
//   fpr=fopen(name,"w");
//   if(fpr == NULL) {
//     printf("file not found inside  DATA directory: %s", name);
//   }
//
//     for (k=0; k < 1; k++) {
//       for (x=1; x < mpiparam.rows_x-1; x++) {
//         //for (z=1; z < 2; z++) {
//         z = (int) mpiparam.rows_y/2;
//           for (y=1; y < mpiparam.rows_y-1; y++) {
//             index        = y + mpiparam.rows_y * ( z + mpiparam.rows_z * x);
//             fprintf(fpr, "%le\t", iter_gridinfom[index].disp[2][2]);
//           }
//         //}
//         fprintf(fpr,"\n");
//       }
//
//     }
}

void writetofile_mpi_binary(struct fields* gridinfo, char *argv[], long t) {
  long x,y,z;
  long gidy;
  FILE *fp;
  char name[1000];
  double composition;
  long b, k;
  sprintf(name,"DATA/Processor_%d/%s_%ld.vtk",rank, argv[3], t);
  fp=fopen(name,"wb");
  write_cells_vtk_mpi_binary(fp, gridinfo);
  fclose(fp);

}
void write_cells_vtk_mpi_binary(FILE *fp, struct fields* gridinfo) {
  long x, y, z, index;
  long a, b;
  long k;
  long global_index;
  long dim;
  FILE *fpr;
  char name[1000];
  double value;

  fprintf(fp,"%ld\n", mpiparam.rows_x-2);
  fprintf(fp,"%ld\n", mpiparam.rows_y-2);
  fprintf(fp,"%ld\n", mpiparam.rows_z-2);

  for (a=0; a < NUMPHASES; a++) {
    for (x=1; x < mpiparam.rows_x-1; x++) {
      for (z=1; z < mpiparam.rows_z-1; z++) {
        for (y=1; y < mpiparam.rows_y-1; y++) {
          index        = y + mpiparam.rows_y * ( z + mpiparam.rows_z * x);
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
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    for (x=1; x < mpiparam.rows_x-1; x++) {
      for (z=1; z < mpiparam.rows_z-1; z++) {
        for (y=1; y < mpiparam.rows_y-1; y++) {
          index        = y + mpiparam.rows_y * ( z + mpiparam.rows_z * x);
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
//   if(WRITECOMPOSITION) {
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      for (x=1; x < mpiparam.rows_x-1; x++) {
        for (z=1; z < mpiparam.rows_z-1; z++) {
          for (y=1; y < mpiparam.rows_y-1; y++) {
            index        = y + mpiparam.rows_y * ( z + mpiparam.rows_z * x);
          if (IS_LITTLE_ENDIAN) {
            value = swap_bytes(gridinfo[index].compi[k]);
          } else {
            value = gridinfo[index].composition[k];
          }
          fwrite(&value, sizeof(double), 1, fp);
          }
        }
      }
      fprintf(fp,"\n");
    }
//   }
  if (ELASTICITY) {
    for (dim=0; dim < DIMENSION; dim++) {
      for (x=1; x < mpiparam.rows_x-1; x++) {
        for (z=1; z < mpiparam.rows_z-1; z++) {
          for (y=1; y < mpiparam.rows_y-1; y++) {
            index        = y + mpiparam.rows_y * ( z + mpiparam.rows_z * x);
          if(IS_LITTLE_ENDIAN) {
            value = swap_bytes(iter_gridinfom[index].disp[dim][2]);
          } else {
            value = iter_gridinfom[index].disp[dim][2];
          }
          fwrite(&value, sizeof(double), 1, fp);
          }
        }
      }
      fprintf(fp,"\n");
    }
  }

  if (!ISOTHERMAL) {
    for (x=1; x < mpiparam.rows_x-1; x++) {
      for (z=1; z < mpiparam.rows_z-1; z++) {
        for (y=1; y < mpiparam.rows_y-1; y++) {
          index        = y + mpiparam.rows_y * ( z + mpiparam.rows_z * x);
          if(IS_LITTLE_ENDIAN) {
            value = swap_bytes(gridinfo[index].temperature);
          } else {
            value = gridinfo[index].temperature;
          }
          fwrite(&value, sizeof(double), 1, fp);
        }
      }
    }
  }
}

#endif
