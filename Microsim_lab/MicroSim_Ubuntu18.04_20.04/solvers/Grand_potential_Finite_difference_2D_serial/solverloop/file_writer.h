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
void write_cells_vtk_2D(FILE *fp, struct fields *gridinfo) {
  long x, y, z, index;
  long a, b;
  long k;
  double composition;
  fprintf(fp,"# vtk DataFile Version 3.0\n");
  fprintf(fp,"Microsim_fields\n");
  fprintf(fp,"ASCII\n");
  fprintf(fp,"DATASET STRUCTURED_POINTS\n");
  fprintf(fp,"DIMENSIONS %ld %ld %ld\n",MESH_Y, MESH_X, (long)1);
  fprintf(fp,"ORIGIN 0 0 0\n");
  fprintf(fp,"SPACING %le %le %le\n",deltax, deltay, deltaz);
  fprintf(fp,"POINT_DATA %ld\n",(long)MESH_X*(long)MESH_Y);
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
  if(WRITECOMPOSITION) {
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      fprintf(fp,"SCALARS Composition_%s double 1\n",Components[k]);
      fprintf(fp,"LOOKUP_TABLE default\n");
      for (x=start[X]; x<=end[X]; x++) {
        for (z=start[Z]; z <= end[Z]; z++) {
          for (y=start[Y]; y <= end[Y]; y++) {
            index = x*layer_size + z*rows_y + y;
            composition=0.0;
            if(ISOTHERMAL) {
              for (b=0; b < NUMPHASES; b++) {
                composition += c_mu(gridinfo[index].compi, T, b, k)*hphi(gridinfo[index].phia, b);
              }
              fprintf(fp,"%le \n",composition);
            } else {
              for (b=0; b < NUMPHASES; b++) {
                composition += c_mu(gridinfo[index].compi, gridinfo[index].temperature, b, k)*hphi(gridinfo[index].phia, b);
              }
              fprintf(fp,"%le \n",composition);
            }
          }
        }
      }
      fprintf(fp,"\n");
    }
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
  }
  if(WRITECOMPOSITION) {
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      fprintf(fp,"SCALARS Composition_%s double 1\n",Components[k]);
      fprintf(fp,"LOOKUP_TABLE default\n");
      for (x=start[X]; x<=end[X]; x++) {
        for (z=start[Z]; z <= end[Z]; z++) {
          for (y=start[Y]; y <= end[Y]; y++) {
            index = x*layer_size + z*rows_y + y;
            composition=0.0;
            if(ISOTHERMAL) {
              for (b=0; b < NUMPHASES; b++) {
                composition += c_mu(gridinfo[index].compi, T, b, k)*hphi(gridinfo[index].phia, b);
              }
              if (IS_LITTLE_ENDIAN) {
                value = swap_bytes(composition);
              } else {
                value = composition;
              }
              fwrite(&value, sizeof(double), 1, fp);
            } else {
              for (b=0; b < NUMPHASES; b++) {
                composition += c_mu(gridinfo[index].compi, gridinfo[index].temperature, b, k)*hphi(gridinfo[index].phia, b);
              }
              if (IS_LITTLE_ENDIAN) {
                value = swap_bytes(composition);
              } else {
                value = composition;
              }
              fwrite(&value, sizeof(double), 1, fp);
            }
          }
        }
      }
    }
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
  if(WRITECOMPOSITION) {
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      fscanf(fp,"%*[^\n]\n");
      fscanf(fp,"%*[^\n]\n");
      for (x=start[X]; x<=end[X]; x++) {
        for (z=start[Z]; z <= end[Z]; z++) {
          for (y=start[Y]; y <= end[Y]; y++) {
            index = x*layer_size + z*rows_y + y;
            fscanf(fp,"%le \n",&composition);
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
  if(WRITECOMPOSITION) {
    for (k=0; k < NUMCOMPONENTS-1; k++) {
     fscanf(fp,"%*[^\n]\n");
     fscanf(fp,"%*[^\n]\n");
     for (x=start[X]; x<=end[X]; x++) {
       for (z=start[Z]; z <= end[Z]; z++) {
          for (y=start[Y]; y <= end[Y]; y++) {
            index = x*layer_size + z*rows_y + y;
            fread(&value, sizeof(double), 1, fp);
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
#endif