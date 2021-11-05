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
}

void writetofile_serial2D(struct fields* gridinfo, char *argv[], long t);
void writetofile_serial2D_binary(struct fields* gridinfo, char *argv[], long t);
void write_cells_vtk_2D(FILE *fp, struct fields *gridinfo);
void write_cells_vtk_2D_binary(FILE *fp, struct fields *gridinfo);

void writetofile_serial2D(struct fields* gridinfo, char *argv[], long t) {
  long x,y,z;
  long gidy;
  long index;
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
  for (a=0; a < NUMPHASES-1; a++) {
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
    fprintf(fp,"SCALARS comp_%s double 1\n",Components[k]);
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
  
  for (a=0; a < NUMPHASES-1; a++) {
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
    fprintf(fp,"SCALARS comp_%s double 1\n",Components[k]);
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
}
#endif
