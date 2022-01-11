// #include <mpi.h>
#include <hdf5.h>
#include<assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h> 
#include <stdbool.h>
#include <sys/stat.h>
#include "functions/global_vars.h"
#include "functions/functions.h"
#include "functions/matrix.h"
#include "functions/utility_functions.h"
#include "functions/reading_input_parameters.h"

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

void write_to_file_2D_ascii(FILE *fp, double **buffer);
void write_to_file_2D_binary(FILE *fp, double **buffer);

long get_global_index(long index, long *workers_rows, long* workers_offset) {
  long x;
  long y;
  x = index/workers_rows[Y];
  y = index%workers_rows[Y];
  return (x+workers_offset[X])*(MESH_Y+6) + (y+workers_offset[Y]);
}


int main (int argc, char * argv[]) {
  FILE *fp;
  char filename[1000];
  long x, y, z;
  long index;
  long global_index;
  double *values;
  int numworkers;
  long n, num_lines;
//   long rows[2];
  long workers_rows[2];
  long workers_offset[2];
  double value;
  double value_bin;
  long numx;
  long numy;
  long t;
  
  long t_start;
  long t_end;
//   long size_fields;
  
  
  reading_input_parameters(argv);
//   serialinfo_xy();

  start[X]   = 3;
  start[Y]   = 3;
  start[Z]   = 3;
  
  rows_x     = MESH_X + 6;
  rows_y     = MESH_Y + 6;
  rows_z     = MESH_Z + 6;
  end[X]     = rows_x - 4;
  end[Y]     = rows_y - 4;
  end[Z]     = rows_z - 4;
  
  layer_size = rows_y*rows_z;
  
  if (DIMENSION == 2) {
    rows_z     = 1;
    start[Z]   = 0; 
    end[Z]     = 0;
    layer_size = rows_y;
  }
  
  long index_count = (MESH_X+6)*(MESH_Y+6);
  

  double **buffer;
  
  buffer = (double **)malloc((MESH_X+6)*(MESH_Y+6)*sizeof(double*));
  
  size_fields = NUMPHASES + (NUMCOMPONENTS-1);
  
  if(WRITECOMPOSITION) {
    size_fields += (NUMCOMPONENTS-1);
  }
  if(!ISOTHERMAL) {
    size_fields += 1;
  }
  
  for (index=0; index < index_count; index++) {
    buffer[index] = (double *)malloc(size_fields*sizeof(double));
  }
  
  
  numworkers =  atoi(argv[3]);
  t_start    =  atol(argv[4]);
  t_end      =  atol(argv[5]);
  
  
  for (t=t_start; t<=t_end; t+=saveT) {
     for (n=0; n < numworkers; n++) {
      sprintf(filename, "DATA/Processor_%ld/%s_%ld.vtk",n, argv[2],t);
      if (ASCII) {
        fp = fopen(filename, "r");
        fscanf(fp,"%ld", &workers_rows[X]);
        fscanf(fp,"%ld", &workers_rows[Y]);
        fscanf(fp,"%ld", &workers_offset[X]);
        fscanf(fp,"%ld", &workers_offset[Y]);  
      } else {
        fp = fopen(filename, "rb");
        fread(&workers_rows[X],   sizeof(long), 1, fp);
        fread(&workers_rows[Y],   sizeof(long), 1, fp);
        fread(&workers_offset[X], sizeof(long), 1, fp);
        fread(&workers_offset[Y], sizeof(long), 1, fp);
      }
      
      for (a=0; a<NUMPHASES; a++) {
        if (ASCII) {
          for (num_lines=0; num_lines < (workers_rows[X]*workers_rows[Y]); num_lines++) {
            fscanf(fp,"%ld %le", &global_index, &value);
            global_index = get_global_index(num_lines, workers_rows, workers_offset);
            buffer[global_index][a] = value;
          }
        } else {
          for (num_lines=0; num_lines < (workers_rows[X]*workers_rows[Y]); num_lines++) {
            fread(&value_bin, sizeof(double), 1, fp);
            global_index = get_global_index(num_lines, workers_rows, workers_offset);
            buffer[global_index][a] = value_bin;
          }
        }
      }
      for (k=0; k<(NUMCOMPONENTS-1); k++) {
        if(ASCII) {
          for (num_lines=0; num_lines < (workers_rows[X]*workers_rows[Y]); num_lines++) {
            fscanf(fp,"%ld %le", &global_index, &value);
            buffer[global_index][NUMPHASES+k] = value;
          }
        } else {
          for (num_lines=0; num_lines < (workers_rows[X]*workers_rows[Y]); num_lines++) {
            fread(&value_bin, sizeof(double), 1, fp);
            global_index = get_global_index(num_lines, workers_rows, workers_offset);
            buffer[global_index][NUMPHASES+k] = value_bin;
          }
        }
      }
      if(WRITECOMPOSITION) {
        for (k=0; k<(NUMCOMPONENTS-1); k++) {
          if (ASCII) {
            for (num_lines=0; num_lines < (workers_rows[X]*workers_rows[Y]); num_lines++) {
              fscanf(fp,"%ld %le", &global_index, &value);
              buffer[global_index][NUMPHASES+(NUMCOMPONENTS-1)+k] = value;
            }
          } else {
            for (num_lines=0; num_lines < (workers_rows[X]*workers_rows[Y]); num_lines++) {
              fread(&value_bin, sizeof(double), 1, fp);
              global_index = get_global_index(num_lines, workers_rows, workers_offset);
              buffer[global_index][NUMPHASES+(NUMCOMPONENTS-1)+k] = value_bin;
            }
          }
        }
      }
      if (!ISOTHERMAL) {
        if (ASCII) {
          for (num_lines=0; num_lines < (workers_rows[X]*workers_rows[Y]); num_lines++) {
            fscanf(fp,"%ld %le", &global_index, &value);
            buffer[global_index][size_fields-1] = value;
          }
        } else {
          for (num_lines=0; num_lines < (workers_rows[X]*workers_rows[Y]); num_lines++) {
            fread(&value_bin, sizeof(double), 1, fp);
            global_index = get_global_index(num_lines, workers_rows, workers_offset);
            buffer[global_index][size_fields-1] = value_bin;
          }
        }
      }
      fclose(fp);
    }
    
    if (ASCII) {
      sprintf(filename, "DATA/%s_%ld.vtk", argv[2],t);
      fp = fopen(filename, "w");
      write_to_file_2D_ascii(fp, buffer);
      fclose(fp);
    } else {
      sprintf(filename, "DATA/%s_%ld.vtk", argv[2],t);
      fp = fopen(filename, "wb");
      write_to_file_2D_binary(fp, buffer);
      fclose(fp);
    }
  }
  for (index=0; index < index_count; index++) {
   free(buffer[index]);
  }
  free(buffer);
  Free3M(Diffusivity, NUMPHASES, NUMCOMPONENTS-1);
  Free3M(ceq,         NUMPHASES, NUMPHASES);
  Free3M(cfill,       NUMPHASES, NUMPHASES);
  Free3M(ceq_coeffs,  NUMPHASES, NUMCOMPONENTS-1);
  Free3M(slopes,      NUMPHASES, NUMPHASES);
  Free3M(A,           NUMPHASES, NUMCOMPONENTS-1);
  Free3M(dcbdT,       NUMPHASES, NUMPHASES);

  FreeM(DELTA_T,      NUMPHASES);
  FreeM(DELTA_C,      NUMPHASES);
  FreeM(dcbdT_phase,  NUMPHASES);
  FreeM(B,            NUMPHASES);
  FreeM(Beq,          NUMPHASES);
  FreeM(dBbdT,        NUMPHASES);

  free(C);
  free(eigen_strain);
  free(Stiffness_c);
  free(Stiffness_t);
  Free3M(cmu,NUMCOMPONENTS-1,NUMCOMPONENTS-1);
  Free3M(muc,NUMCOMPONENTS-1,NUMCOMPONENTS-1);
  Free4M(Rotation_matrix,NUMPHASES, NUMPHASES, DIMENSION);
  Free4M(Inv_Rotation_matrix,NUMPHASES, NUMPHASES, DIMENSION);
  free(Rotated_qab);
}
void write_to_file_2D_ascii(FILE *fp, double **buffer) {
  long x, y, z, index;
  long a, k;
//   long size_fields;
  
  
  fprintf(fp,"# vtk DataFile Version 3.0\n");
  fprintf(fp,"Microsim_fields\n");
  fprintf(fp,"ASCII\n");
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
          fprintf(fp, "%le\n",buffer[index][a]);
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
          fprintf(fp, "%le\n",buffer[index][NUMPHASES+k]);
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
            fprintf(fp, "%le\n",buffer[index][NUMPHASES+(NUMCOMPONENTS-1)+k]);
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
          fprintf(fp, "%le \n",buffer[index][size_fields-1]);
        }
      }
    }
  }
}
void write_to_file_2D_binary(FILE *fp, double **buffer) {
  long x, y, z, index;
  long a, k;
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
            value = swap_bytes(buffer[index][a]);
          } else {
            value = buffer[index][a];
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
            value = swap_bytes(buffer[index][NUMPHASES+k]);
          } else {
            value = buffer[index][NUMPHASES+k];
          }
          fwrite(&value, sizeof(double), 1, fp);
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
            if (IS_LITTLE_ENDIAN) {
              value     = swap_bytes(buffer[index][NUMPHASES+(NUMCOMPONENTS-1)+k]);
            } else {
              value = buffer[index][NUMPHASES+(NUMCOMPONENTS-1)+k];
            }
            fwrite(&value, sizeof(double), 1, fp);
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
          if (IS_LITTLE_ENDIAN) {
            value = swap_bytes(buffer[index][size_fields-1]);
          } else {
            value = buffer[index][size_fields-1];
          }
          fwrite(&value, sizeof(double), 1, fp);
        }
      }
    }
  }
}
