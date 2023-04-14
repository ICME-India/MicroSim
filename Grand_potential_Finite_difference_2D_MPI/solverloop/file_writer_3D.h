#ifndef FILE_WRITER_3D_H_
#define FILE_WRITER_3D_H_

// #include <arpa/inet.h>
// #include <endian.h>
// #include <stdint.h>
// 
// #define  IS_BIG_ENDIAN     (1 == htons(1))
// #define  IS_LITTLE_ENDIAN  (!IS_BIG_ENDIAN)

// double swap_bytes(double value) {
//   double  src_num = value;
//   int64_t tmp_num = htobe64(le64toh(*(int64_t*)&src_num));
//   double  dst_num = *(double*)&tmp_num;
// }

// void writetofile_serial2D(struct fields* gridinfo, char *argv[], long t);
// void writetofile_serial2D_binary(struct fields* gridinfo, char *argv[], long t);
// void write_cells_vtk_3D(FILE *fp, struct fields *gridinfo);
// void write_cells_vtk_3D_binary(FILE *fp, struct fields *gridinfo);
void writetofile_mpi3D(struct fields* gridinfo, char *argv[], long t);
void writetofile_mpi3D_binary(struct fields* gridinfo, char *argv[], long t);
void write_cells_vtk_3D_mpi(FILE *fp, struct fields* gridinfo);
void write_cells_vtk_3D_mpibinary(FILE *fp, struct fields* gridinfo);
void writetofile_mpi3D_hdf5(struct fields* gridinfo, char *argv[], long t);
void write_cells_hdf5_3D_mpi(hid_t file_id, struct fields* gridinfo);

// void writetofile_serial2D(struct fields* gridinfo, char *argv[], long t) {
//   long x,y,z;
//   long gidy;
//   FILE *fp;
//   char name[1000];
//   double composition;
//   long b, k;
//   sprintf(name,"DATA/%s_%ld.vtk",argv[3], t);
//   fp=fopen(name,"w");
//   write_cells_vtk_2D(fp, gridinfo);
//   fclose(fp);
// }
// void writetofile_serial2D_binary(struct fields* gridinfo, char *argv[], long t) {
//   long x,y,z;
//   long gidy;
//   FILE *fp;
//   char name[1000];
//   double composition;
//   long b, k;
//   sprintf(name,"DATA/%s_%ld.vtk",argv[3], t);
//   fp=fopen(name,"wb");
//   write_cells_vtk_2D_binary(fp, gridinfo);
//   fclose(fp);
// }
void writetofile_mpi3D(struct fields* gridinfo, char *argv[], long t) {
  long x,y,z;
  long gidy;
  FILE *fp;
  char name[1000];
  double composition;
  long b, k;
  sprintf(name,"DATA/Processor_%d/%s_%ld.vtk",taskid, argv[3], t);
  fp=fopen(name,"w");
  write_cells_vtk_3D_mpi(fp, gridinfo);
  fclose(fp);
}
void writetofile_mpi3D_binary(struct fields* gridinfo, char *argv[], long t) {
  long x,y,z;
  long gidy;
  FILE *fp;
  char name[1000];
  double composition;
  long b, k;
  sprintf(name,"DATA/Processor_%d/%s_%ld.vtk",taskid, argv[3], t);
  fp=fopen(name,"wb");
  write_cells_vtk_3D_mpibinary(fp, gridinfo);
  fclose(fp);
}
void writetofile_mpi3D_hdf5(struct fields* gridinfo, char *argv[], long t) {
  hid_t file_id;
  hid_t plist_id;
  
  char filename_hdf5[1000];
  
  //every processor creates a file collectively
//   sprintf(filename_hdf5, "DATA/%s_%d_%ld.h5", argv[3], numtasks, t);
  sprintf(filename_hdf5, "DATA/%s_%ld.h5", argv[3], t);
  /* Set up file access property list with parallel I/O access*/
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  
  file_id = H5Fcreate(filename_hdf5, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id); // we'll use this plist_id again later 
  
  write_cells_hdf5_3D_mpi(file_id, gridinfo);
  
  status_h = H5Fclose(file_id);
}

void readfromfile_mpi3D_hdf5(struct fields* gridinfo, char *argv[], long numworkers, long t) {
  hid_t file_id;
  hid_t plist_id;
  
  char filename_hdf5[1000];
  
  //every processor creates a file collectively
//   sprintf(filename_hdf5, "DATA/%s_%ld_%ld.h5", argv[3], numworkers, t);
  sprintf(filename_hdf5, "DATA/%s_%ld.h5", argv[3], t);
  /* Set up file access property list with parallel I/O access*/
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  
//   file_id = H5Fcreate(filename_hdf5, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  file_id = H5Fopen(filename_hdf5, H5F_ACC_RDONLY, plist_id);
  H5Pclose(plist_id); // we'll use this plist_id again later 
  
  read_cells_hdf5_3D_mpi(file_id, gridinfo);
  
  status_h = H5Fclose(file_id);
}

void readfromfile_mpi3D(struct fields* gridinfo, char *argv[], long t) {
  long x,y,z;
  long gidy;
  FILE *fp;
  char name[1000];
  double composition;
  long b, k;
  sprintf(name,"DATA/Processor_%d/%s_%ld.vtk",taskid, argv[3], t);
  fp=fopen(name,"r");
  read_cells_vtk_3D_mpi(fp, gridinfo);
  fclose(fp);
}
void readfromfile_mpi3D_binary(struct fields* gridinfo, char *argv[], long t) {
  long x,y,z;
  long gidy;
  FILE *fp;
  char name[1000];
  double composition;
  long b, k;
  sprintf(name,"DATA/Processor_%d/%s_%ld.vtk",taskid, argv[3], t);
  fp=fopen(name,"rb");
  read_cells_vtk_3D_mpibinary(fp, gridinfo);
  fclose(fp);
}

// void write_cells_vtk_3D(FILE *fp, struct fields *gridinfo) {
//   long x, y, z, index;
//   long a, b;
//   long k;
//   double composition;
//   fprintf(fp,"# vtk DataFile Version 3.0\n");
//   fprintf(fp,"Microsim_fields\n");
//   fprintf(fp,"ASCII\n");
//   fprintf(fp,"DATASET STRUCTURED_POINTS\n");
//   fprintf(fp,"DIMENSIONS %ld %ld %ld\n",MESH_Y, MESH_X, (long)1);
//   fprintf(fp,"ORIGIN 0 0 0\n");
//   fprintf(fp,"SPACING %le %le %le\n",deltax, deltay, deltaz);
//   fprintf(fp,"POINT_DATA %ld\n",(long)MESH_X*(long)MESH_Y);
// //   fprintf(fp,"FIELD PHASEFIELD 2\n");
// //   fprintf(fp,"PHI[0] 1 %ld double\n", (long)Nx*(long)Ny);
//   for (a=0; a < NUMPHASES; a++) {
//     fprintf(fp,"SCALARS %s double 1\n",Phases[a]);
//     fprintf(fp,"LOOKUP_TABLE default\n");
//     for (x=start[X]; x<=end[X]; x++) {
//       for (z=start[Z]; z <= end[Z]; z++) {
//         for (y=start[Y]; y <= end[Y]; y++) {
//           index = x*layer_size + z*rows_y + y;
//           fprintf(fp, "%le\n",gridinfo[index].phia[a]);
//         }
//       }
//     }
//     fprintf(fp,"\n");
//   }
//   for (k=0; k < NUMCOMPONENTS-1; k++) {
//     fprintf(fp,"SCALARS Mu_%s double 1\n",Components[k]);
//     fprintf(fp,"LOOKUP_TABLE default\n");
//     for (x=start[X]; x<=end[X]; x++) {
//       for (z=start[Z]; z <= end[Z]; z++) {
//         for (y=start[Y]; y <= end[Y]; y++) {
//           index = x*layer_size + z*rows_y + y;
//           fprintf(fp, "%le\n",gridinfo[index].compi[k]);
//         }
//       }
//     }
//     fprintf(fp,"\n");
//   }
// //   if(WRITECOMPOSITION) {
//     for (k=0; k < NUMCOMPONENTS-1; k++) {
//       fprintf(fp,"SCALARS Composition_%s double 1\n",Components[k]);
//       fprintf(fp,"LOOKUP_TABLE default\n");
//       for (x=start[X]; x<=end[X]; x++) {
//         for (z=start[Z]; z <= end[Z]; z++) {
//           for (y=start[Y]; y <= end[Y]; y++) {
//             index = x*layer_size + z*rows_y + y;
//             fprintf(fp, "%le\n",gridinfo[index].composition[k]);
// //             composition=0.0;
// //             if(ISOTHERMAL) {
// //               for (b=0; b < NUMPHASES; b++) {
// //                 composition += c_mu(gridinfo[index].compi, T, b, k)*hphi(gridinfo[index].phia, b);
// //               }
// //               fprintf(fp,"%le \n",composition);
// //             } else {
// //               for (b=0; b < NUMPHASES; b++) {
// //                 composition += c_mu(gridinfo[index].compi, gridinfo[index].temperature, b, k)*hphi(gridinfo[index].phia, b);
// //               }
// //               fprintf(fp,"%le \n",composition);
// //             }
//           }
//         }
//       }
//       fprintf(fp,"\n");
//     }
// //   }
//   if (!ISOTHERMAL) {
//     fprintf(fp,"SCALARS Temperature double 1\n");
//     fprintf(fp,"LOOKUP_TABLE default\n");
//     for (x=start[X]; x<=end[X]; x++) {
//       for (z=start[Z]; z <= end[Z]; z++) {
//         for (y=start[Y]; y <= end[Y]; y++) {
//           index = x*layer_size + z*rows_y + y;
//           fprintf(fp, "%le \n",gridinfo[index].temperature);
//         }
//       }
//     }
//   }
// }
// void write_cells_vtk_2D_binary(FILE *fp, struct fields *gridinfo) {
//   long x, y, z, index;
//   long a, b;
//   long k;
//   double composition;
//   char first_lines[50];
//   int length;
//   double value;
//   
//   fprintf(fp,"# vtk DataFile Version 3.0\n");
//   fprintf(fp,"Microsim_fields\n");
//   fprintf(fp,"BINARY\n");
//   fprintf(fp,"DATASET STRUCTURED_POINTS\n");
//   fprintf(fp,"DIMENSIONS %ld %ld %ld\n",MESH_Y, MESH_X, (long)1);
//   fprintf(fp,"ORIGIN 0 0 0\n");
//   fprintf(fp,"SPACING %le %le %le\n",deltax, deltay, deltaz);
//   fprintf(fp,"POINT_DATA %ld\n",(long)MESH_X*(long)MESH_Y);
//   
//   for (a=0; a < NUMPHASES; a++) {
//     fprintf(fp,"SCALARS %s double 1\n",Phases[a]);
//     fprintf(fp,"LOOKUP_TABLE default\n");
//     for (x=start[X]; x<=end[X]; x++) {
//       for (z=start[Z]; z <= end[Z]; z++) {
//         for (y=start[Y]; y <= end[Y]; y++) {
//           index = x*layer_size + z*rows_y + y;
//           if (IS_LITTLE_ENDIAN) {
//             value = swap_bytes(gridinfo[index].phia[a]);
//           } else {
//             value = gridinfo[index].phia[a];
//           }
//           fwrite(&value, sizeof(double), 1, fp);
//         }
//       }
//     }
//     fprintf(fp,"\n");
//   }
//   for (k=0; k < NUMCOMPONENTS-1; k++) {
//     fprintf(fp,"SCALARS Mu_%s double 1\n",Components[k]);
//     fprintf(fp,"LOOKUP_TABLE default\n");
//     for (x=start[X]; x<=end[X]; x++) {
//       for (z=start[Z]; z <= end[Z]; z++) {
//         for (y=start[Y]; y <= end[Y]; y++) {
//           index = x*layer_size + z*rows_y + y;
//           if (IS_LITTLE_ENDIAN) {
//             value = swap_bytes(gridinfo[index].compi[k]);
//           } else {
//             value = gridinfo[index].compi[k];
//           }
//           fwrite(&value, sizeof(double), 1, fp);
//         }
//       }
//     }
//    fprintf(fp,"\n");
//   }
// //   if(WRITECOMPOSITION) {
//     for (k=0; k < NUMCOMPONENTS-1; k++) {
//       fprintf(fp,"SCALARS Composition_%s double 1\n",Components[k]);
//       fprintf(fp,"LOOKUP_TABLE default\n");
//       for (x=start[X]; x<=end[X]; x++) {
//         for (z=start[Z]; z <= end[Z]; z++) {
//           for (y=start[Y]; y <= end[Y]; y++) {
//             index = x*layer_size + z*rows_y + y;
// //             composition=0.0;
// //             if(ISOTHERMAL) {
// //               for (b=0; b < NUMPHASES; b++) {
// //                 composition += c_mu(gridinfo[index].compi, T, b, k)*hphi(gridinfo[index].phia, b);
// //               }
// //               if (IS_LITTLE_ENDIAN) {
// //                 value = swap_bytes(composition);
// //               } else {
// //                 value = composition;
// //               }
// //               fwrite(&value, sizeof(double), 1, fp);
// //             } else {
// //               for (b=0; b < NUMPHASES; b++) {
// //                 composition += c_mu(gridinfo[index].compi, gridinfo[index].temperature, b, k)*hphi(gridinfo[index].phia, b);
// //               }
// //               if (IS_LITTLE_ENDIAN) {
// //                 value = swap_bytes(composition);
// //               } else {
// //                 value = composition;
// //               }
// //               fwrite(&value, sizeof(double), 1, fp);
// //             }
//             if (IS_LITTLE_ENDIAN) {
//               value = swap_bytes(gridinfo[index].composition[k]);
//             } else {
//               value = gridinfo[index].composition[k];
//             }
//             fwrite(&value, sizeof(double), 1, fp);
//           }
//         }
//       }
//       fprintf(fp,"\n");
//     }
// //   }
//   if (!ISOTHERMAL) {
//     fprintf(fp,"SCALARS Temperature double 1\n");
//     fprintf(fp,"LOOKUP_TABLE default\n");
//     for (x=start[X]; x<=end[X]; x++) {
//       for (z=start[Z]; z <= end[Z]; z++) {
//         for (y=start[Y]; y <= end[Y]; y++) {
//           index = x*layer_size + z*rows_y + y;
//           if(IS_LITTLE_ENDIAN) {
//             value = swap_bytes(gridinfo[index].temperature);
//           } else {
//             value = gridinfo[index].temperature;
//           }
//           fwrite(&value, sizeof(double), 1, fp);
//         }
//       }
//     }
//   }
// }
void write_cells_vtk_3D_mpi(FILE *fp, struct fields* gridinfo) {
  long x, y, z, index;
  long a, b;
  long k;
  double composition;
  long global_index;
  long dim;
  
  fprintf(fp,"%ld\n",workers_mpi.rows[X]);
  fprintf(fp,"%ld\n",workers_mpi.rows[Y]);
  fprintf(fp,"%ld\n",workers_mpi.rows[Z]);
  fprintf(fp,"%ld\n",workers_mpi.offset[X]);
  fprintf(fp,"%ld\n",workers_mpi.offset[Y]);
  fprintf(fp,"%ld\n",workers_mpi.offset[Z]);
  
  for (a=0; a < NUMPHASES; a++) {
    for (x=0; x < workers_mpi.rows[X]; x++) {
      for (z=0; z < workers_mpi.rows[Z]; z++) {
        for (y=0; y < workers_mpi.rows[Y]; y++) {
          index        = (x + workers_mpi.offset_x)*workers_mpi.layer_size   + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y + workers_mpi.offset_y);
          global_index = (x + workers_mpi.offset[X])*(MESH_Y+6)*(MESH_Z+6) + (z + workers_mpi.offset[Z])*(MESH_Y+6) + (y + workers_mpi.offset[Y]); 
          fprintf(fp, "%ld %le\n",global_index, gridinfo[index].phia[a]);
        }
      }
    }
    fprintf(fp,"\n");
  }
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    for (x=0; x < workers_mpi.rows[X]; x++) {
      for (z=0; z < workers_mpi.rows[Z]; z++) {
        for (y=0; y < workers_mpi.rows[Y]; y++) {
          index        = (x + workers_mpi.offset_x)*workers_mpi.layer_size   + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y + workers_mpi.offset_y);
          global_index = (x + workers_mpi.offset[X])*(MESH_Y+6)*(MESH_Z+6) + (z + workers_mpi.offset[Z])*(MESH_Y+6) + (y + workers_mpi.offset[Y]); 
          fprintf(fp, "%ld %le\n",global_index, gridinfo[index].compi[k]);
        }
      }
    }
    fprintf(fp,"\n");
  }
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    for (x=0; x < workers_mpi.rows[X]; x++) {
      for (z=0; z < workers_mpi.rows[Z]; z++) {
        for (y=0; y < workers_mpi.rows[Y]; y++) {
          index        = (x + workers_mpi.offset_x)*workers_mpi.layer_size   + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y + workers_mpi.offset_y);
          global_index = (x + workers_mpi.offset[X])*(MESH_Y+6)*(MESH_Z+6) + (z + workers_mpi.offset[Z])*(MESH_Y+6) + (y + workers_mpi.offset[Y]); 
          fprintf(fp,"%ld %le \n",global_index, gridinfo[index].composition[k]);
        }
      }
    }
    fprintf(fp,"\n");
  }
  if (ELASTICITY) {
    for (dim=0; dim < DIMENSION; dim++) {
      for (x=0; x < workers_mpi.rows[X]; x++) {
        for (z=0; z < workers_mpi.rows[Z]; z++) {
          for (y=0; y < workers_mpi.rows[Y]; y++) {
            index        = (x + workers_mpi.offset_x)*workers_mpi.layer_size   + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y + workers_mpi.offset_y);
            global_index = (x + workers_mpi.offset[X])*(MESH_Y+6)*(MESH_Z+6) + (z + workers_mpi.offset[Z])*(MESH_Y+6) + (y + workers_mpi.offset[Y]); 
            fprintf(fp,"%ld %le \n",global_index, iter_gridinfo_w[index].disp[dim][2]);
          }
        }
      }
      fprintf(fp,"\n");
    }
  }
  
  if (!ISOTHERMAL) {
    for (x=0; x < workers_mpi.rows[X]; x++) {
      for (z=0; z < workers_mpi.rows[Z]; z++) {
        for (y=0; y < workers_mpi.rows[Y]; y++) {
          index        = (x + workers_mpi.offset_x)*workers_mpi.layer_size + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y + workers_mpi.offset_y);
          global_index = (x + workers_mpi.offset[X])*(MESH_Y+6)*(MESH_Z+6) + (z + workers_mpi.offset[Z])*(MESH_Y+6)+ (y + workers_mpi.offset[Y]); 
          fprintf(fp, "%ld %le \n",global_index, gridinfo[index].temperature);
        }
      }
    }
  }
}
void write_cells_vtk_3D_mpibinary(FILE *fp, struct fields* gridinfo) {
  long x, y, z, index;
  long a, b;
  long k;
  double composition;
  long global_index;
  double value;
  long dim;
  
  fwrite(&workers_mpi.rows[X],   sizeof(long), 1, fp);
  fwrite(&workers_mpi.rows[Y],   sizeof(long), 1, fp);
  fwrite(&workers_mpi.rows[Z],   sizeof(long), 1, fp);
  fwrite(&workers_mpi.offset[X], sizeof(long), 1, fp);
  fwrite(&workers_mpi.offset[Y], sizeof(long), 1, fp);
  fwrite(&workers_mpi.offset[Z], sizeof(long), 1, fp);
  
  for (a=0; a < NUMPHASES; a++) {
    for (x=0; x < workers_mpi.rows[X]; x++) {
      for (z=0; z < workers_mpi.rows[Z]; z++) {
        for (y=0; y < workers_mpi.rows[Y]; y++) {
          index        = (x + workers_mpi.offset_x)*workers_mpi.layer_size   + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y + workers_mpi.offset_y);
          value        = gridinfo[index].phia[a];
          fwrite(&value, sizeof(double), 1, fp);
        }
      }
    }
  }
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    for (x=0; x < workers_mpi.rows[X]; x++) {
      for (z=0; z < workers_mpi.rows[Z]; z++) {
        for (y=0; y < workers_mpi.rows[Y]; y++) {
          index        = (x + workers_mpi.offset_x)*workers_mpi.layer_size   + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y + workers_mpi.offset_y);
          value        = gridinfo[index].compi[k];
          fwrite(&value, sizeof(double), 1, fp);
        }
      }
    }
  }
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    for (x=0; x < workers_mpi.rows[X]; x++) {
      for (z=0; z < workers_mpi.rows[Z]; z++) {
        for (y=0; y < workers_mpi.rows[Y]; y++) {
          index        = (x + workers_mpi.offset_x)*workers_mpi.layer_size   + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y + workers_mpi.offset_y);
          value        = gridinfo[index].composition[k];
          fwrite(&value, sizeof(double), 1, fp);
        }
      }
    }
  }
  if (ELASTICITY) {
    for (dim=0; dim < DIMENSION; dim++) {
      for (x=0; x < workers_mpi.rows[X]; x++) {
        for (z=0; z < workers_mpi.rows[Z]; z++) {
          for (y=0; y < workers_mpi.rows[Y]; y++) {
            index        = (x + workers_mpi.offset_x)*workers_mpi.layer_size   + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y + workers_mpi.offset_y);
            value        = iter_gridinfo_w[index].disp[dim][2];
            fwrite(&value, sizeof(double), 1, fp);
          }
        }
      }
    }
  }
  if (!ISOTHERMAL) {
    for (x=0; x < workers_mpi.rows[X]; x++) {
      for (z=0; z < workers_mpi.rows[Z]; z++) {
        for (y=0; y < workers_mpi.rows[Y]; y++) {
          index        = (x + workers_mpi.offset_x)*workers_mpi.layer_size   + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y + workers_mpi.offset_y);
          value        = gridinfo[index].temperature;
          fwrite(&value, sizeof(double), 1, fp);
        }
      }
    }
  }
}
void write_cells_hdf5_3D_mpi(hid_t file_id, struct fields* gridinfo_w) {
  hsize_t dims[DIMENSION];
  hsize_t count[DIMENSION];
  hsize_t offset_slab[DIMENSION];
  int rank_hdf5 = DIMENSION;
  
  long i, a, b, k, index, index_to;
  long x, y, z;
  
  double **buffer;
  double composition;
  
  hid_t dset_id; //handles
  hid_t dataspace_id, memspace_id; 
  hid_t plist_id;
  
  long dim;
  
  dims[2] = MESH_Y+6;
  dims[1] = MESH_Z+6;
  dims[0] = MESH_X+6;
  
  count[2] = workers_mpi.rows[Y];
  count[1] = workers_mpi.rows[Z];
  count[0] = workers_mpi.rows[X];
  
  offset_slab[2] = workers_mpi.offset[Y];
  offset_slab[1] = workers_mpi.offset[Z];
  offset_slab[0] = workers_mpi.offset[X];
  
  long index_count = workers_mpi.rows[X]*workers_mpi.rows[Y]*workers_mpi.rows[Z];
  
  buffer = (double **)malloc(size_fields*sizeof(double*));
  for (i=0; i < size_fields; i++) {
    buffer[i] = (double *)malloc(index_count*sizeof(double));
  }
  for (x=0; x < workers_mpi.rows[X]; x++) {
    for (z=0; z < workers_mpi.rows[Z]; z++) {
      for (y=0; y < workers_mpi.rows[Y]; y++) {
        index        = (x + workers_mpi.offset_x)*workers_mpi.layer_size   + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y + workers_mpi.offset_y);
        index_to     =  x*(workers_mpi.rows[Y]*workers_mpi.rows[Z])        + z*workers_mpi.rows[Y] + y;
        for (a=0; a<NUMPHASES; a++) {
          buffer[a][index_to] = gridinfo_w[index].phia[a];
        }
        for (k=0; k<(NUMCOMPONENTS-1); k++) {
          buffer[NUMPHASES+k][index_to] = gridinfo_w[index].compi[k];
        }
        for (k=0; k<(NUMCOMPONENTS-1); k++) {
          buffer[NUMPHASES+(NUMCOMPONENTS-1)+k][index_to] = gridinfo_w[index].composition[k];
        }
        if (ELASTICITY) {
          for (dim=0; dim < DIMENSION; dim ++) {
            buffer[NUMPHASES+2*(NUMCOMPONENTS-1)+dim][index_to] = iter_gridinfo_w[index].disp[dim][2];
          }
        }
        if (!ISOTHERMAL) {
          buffer[size_fields-1][index_to] = gridinfo_w[index].temperature;
        }
      }
    }
  }
  for (i = 0; i < size_fields; i++) {
    dataspace_id = H5Screate_simple(rank_hdf5, dims, NULL);
    /*
    * Create the dataset with default properties and close dataspace_id.
    */
    dset_id = H5Dcreate(file_id, coordNames[i], H5T_NATIVE_DOUBLE, dataspace_id,
                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(dataspace_id);
//     
    memspace_id = H5Screate_simple(rank_hdf5, count, NULL); //count is of shape rank_hdf5
//         
//     /* Select hyperslab in the file.
//             */
    dataspace_id = H5Dget_space(dset_id); // filespace is the handle given by dset_id
//     // 			H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
//     // 			printf("dims[0] = %ld, dims[1] = %ld\n", (long)dims[0], (long)dims[1]);
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset_slab, NULL, count, NULL);
//     
//     /* Create property list for collective dataset write.
//     */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
//     
    status_h = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id,
				plist_id, buffer[i]);
//     
    status_h = H5Dclose(dset_id);
    status_h = H5Sclose(dataspace_id);
    status_h = H5Sclose(memspace_id);
    status_h = H5Pclose(plist_id);
  }
  for(i=0; i < size_fields; i++) {
    free(buffer[i]);
  }
  free(buffer);
}
void read_cells_hdf5_3D_mpi(hid_t file_id, struct fields* gridinfo_w) {
  hsize_t dims[DIMENSION];
  hsize_t count[DIMENSION];
  hsize_t offset_slab[DIMENSION];
  int rank_hdf5 = DIMENSION;
  
  long i, a, b, k, index, index_to;
  long x, y, z;
  
  double **buffer;
  double composition;
  
  hid_t dset_id; //handles
  hid_t dataspace_id, memspace_id; 
  hid_t plist_id;

  long dim;
  
  dims[2] = MESH_Y + 6;
  dims[1] = MESH_Z + 6;
  dims[0] = MESH_X + 6;
  
  count[2] = workers_mpi.rows[Y];
  count[1] = workers_mpi.rows[Z];
  count[0] = workers_mpi.rows[X];
  
  offset_slab[2] = workers_mpi.offset[Y];
  offset_slab[1] = workers_mpi.offset[Z];
  offset_slab[0] = workers_mpi.offset[X];
  
  long index_count = workers_mpi.rows[X]*workers_mpi.rows[Y]*workers_mpi.rows[Z];
  
  buffer = (double **)malloc(size_fields*sizeof(double*));
  for (i=0; i < size_fields; i++) {
    buffer[i] = (double *)malloc(index_count*sizeof(double));
  }
  
  for (i = 0; i < size_fields; i++) {
    /* open the dset_id collectively */
    dset_id = H5Dopen(file_id, coordNames[i], H5P_DEFAULT);
    /* create a file dataspace independently */
    dataspace_id = H5Dget_space(dset_id);
    
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset_slab, NULL, count, NULL);				
    /* create a memory dataspace independently */
    memspace_id = H5Screate_simple(rank_hdf5, count, NULL);
    
//     /* Create property list for collective dataset write.
//     */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
//     
    status_h = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id,
				plist_id, buffer[i]);
//     status_h = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id,	plist_id, vec[i]);
    status_h = H5Dclose(dset_id);
    status_h = H5Sclose(dataspace_id);
    status_h = H5Sclose(memspace_id);
    status_h = H5Pclose(plist_id);
  }
  
  for (x=0; x < workers_mpi.rows[X]; x++) {
    for (z=0; z < workers_mpi.rows[Z]; z++) {
      for (y=0; y < workers_mpi.rows[Y]; y++) {
        index_to     = (x + workers_mpi.offset_x)*workers_mpi.layer_size   + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y + workers_mpi.offset_y);
        index        =  x*(workers_mpi.rows[Y]*workers_mpi.rows[Z])        + z*workers_mpi.rows[Y] + y;
        for (a=0; a<NUMPHASES; a++) {
          gridinfo_w[index_to].phia[a]    = buffer[a][index];
        }
        for (k=0; k<(NUMCOMPONENTS-1); k++) {
          gridinfo_w[index_to].compi[k]    = buffer[NUMPHASES+k][index];
        }
        for (k=0; k<(NUMCOMPONENTS-1); k++) {
          gridinfo_w[index_to].composition[k] = buffer[NUMPHASES+(NUMCOMPONENTS-1)+k][index];
        }
        if (ELASTICITY) {
          for (dim=0; dim < DIMENSION; dim++) {
            iter_gridinfo_w[index_to].disp[dim][2] = buffer[NUMPHASES+2*(NUMCOMPONENTS-1)+dim][index]; 
          }
        }
        if (!ISOTHERMAL) {
          gridinfo_w[index_to].temperature = buffer[size_fields-1][index];
        }
      }
    }
  }
  
  for(i=0; i < size_fields; i++) {
    free(buffer[i]);
  }
  free(buffer);
}
void read_cells_vtk_3D_mpi(FILE *fp, struct fields* gridinfo_w) {
  long x, y, z, index;
  long a, b;
  long k;
  double composition;
  long global_index;
  long dim;
  
  fscanf(fp,"%ld\n",&workers_mpi.rows[X]);
  fscanf(fp,"%ld\n",&workers_mpi.rows[Y]);
  fscanf(fp,"%ld\n",&workers_mpi.rows[Z]);
  fscanf(fp,"%ld\n",&workers_mpi.offset[X]);
  fscanf(fp,"%ld\n",&workers_mpi.offset[Y]);
  fscanf(fp,"%ld\n",&workers_mpi.offset[Z]);
  
  for (a=0; a < NUMPHASES; a++) {
    for (x=0; x < workers_mpi.rows[X]; x++) {
      for (z=0; z < workers_mpi.rows[Z]; z++) {
        for (y=0; y < workers_mpi.rows[Y]; y++) {
          index        = (x + workers_mpi.offset_x)*workers_mpi.layer_size   + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y + workers_mpi.offset_y);
          fscanf(fp, "%ld %le\n",&global_index, &gridinfo_w[index].phia[a]);
        }
      }
    }
//     fprintf(fp,"\n");
  }
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    for (x=0; x < workers_mpi.rows[X]; x++) {
      for (z=0; z < workers_mpi.rows[Z]; z++) {
        for (y=0; y < workers_mpi.rows[Y]; y++) {
          index        = (x + workers_mpi.offset_x)*workers_mpi.layer_size   + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y + workers_mpi.offset_y);
          fscanf(fp, "%ld %le\n",&global_index, &gridinfo_w[index].compi[k]);
        }
      }
    }
//     fprintf(fp,"\n");
  }
//   if(WRITECOMPOSITION) {
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      for (x=0; x < workers_mpi.rows[X]; x++) {
        for (z=0; z < workers_mpi.rows[Z]; z++) {
          for (y=0; y < workers_mpi.rows[Y]; y++) {
            index        = (x + workers_mpi.offset_x)*workers_mpi.layer_size   + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y + workers_mpi.offset_y);
            fscanf(fp,"%ld %le \n",&global_index, &gridinfo_w[index].composition[k]);
          }
        }
      }
//       fprintf(fp,"\n");
    }
    if (ELASTICITY) {
      for (dim=0; dim < DIMENSION; dim++) {
        for (x=0; x < workers_mpi.rows[X]; x++) {
          for (z=0; z < workers_mpi.rows[Z]; z++) {
            for (y=0; y < workers_mpi.rows[Y]; y++) {
              index        = (x + workers_mpi.offset_x)*workers_mpi.layer_size   + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y + workers_mpi.offset_y);
              fscanf(fp,"%ld %le \n",&global_index, &iter_gridinfo_w[index].disp[dim][2]);
            }
          }
        }
  //       fprintf(fp,"\n");
      }
    }
//   }
  if (!ISOTHERMAL) {
    for (x=0; x < workers_mpi.rows[X]; x++) {
      for (z=0; z < workers_mpi.rows[Z]; z++) {
        for (y=0; y < workers_mpi.rows[Y]; y++) {
          index        = (x + workers_mpi.offset_x)*workers_mpi.layer_size   + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y + workers_mpi.offset_y);
          fscanf(fp, "%ld %le \n",&global_index, &gridinfo_w[index].temperature);
        }
      }
    }
  }
}
void read_cells_vtk_3D_mpibinary(FILE *fp, struct fields* gridinfo_w) {
  long x, y, z, index;
  long a, b;
  long k;
  double composition;
  long global_index;
  double value;
  long dim;
  

  fread(&workers_mpi.rows[X],   sizeof(long), 1, fp);
  fread(&workers_mpi.rows[Y],   sizeof(long), 1, fp);
  fread(&workers_mpi.rows[Z],   sizeof(long), 1, fp);
  fread(&workers_mpi.offset[X], sizeof(long), 1, fp);
  fread(&workers_mpi.offset[Y], sizeof(long), 1, fp);
  fread(&workers_mpi.offset[Z], sizeof(long), 1, fp);
  
  for (a=0; a < NUMPHASES; a++) {
    for (x=0; x < workers_mpi.rows[X]; x++) {
      for (z=0; z < workers_mpi.rows[Z]; z++) {
        for (y=0; y < workers_mpi.rows[Y]; y++) {
          index        = (x + workers_mpi.offset_x)*workers_mpi.layer_size   + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y + workers_mpi.offset_y);
          fread(&value, sizeof(double), 1, fp);
          gridinfo_w[index].phia[a] = value;
        }
      }
    }
  }
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    for (x=0; x < workers_mpi.rows[X]; x++) {
      for (z=0; z < workers_mpi.rows[Z]; z++) {
        for (y=0; y < workers_mpi.rows[Y]; y++) {
          index        = (x + workers_mpi.offset_x)*workers_mpi.layer_size   + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y + workers_mpi.offset_y);
          fread(&value, sizeof(double), 1, fp);
          gridinfo_w[index].compi[k] = value;
        }
      }
    }
  }
//   if(WRITECOMPOSITION) {
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    for (x=0; x < workers_mpi.rows[X]; x++) {
      for (z=0; z < workers_mpi.rows[Z]; z++) {
        for (y=0; y < workers_mpi.rows[Y]; y++) {
          index        = (x + workers_mpi.offset_x)*workers_mpi.layer_size   + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y + workers_mpi.offset_y);
          fread(&value, sizeof(double), 1, fp);
          gridinfo_w[index].composition[k] = value;
        }
      }
    }
  }
  if (ELASTICITY) {
    for (dim=0; dim < DIMENSION; dim++) {
      for (x=0; x < workers_mpi.rows[X]; x++) {
        for (z=0; z < workers_mpi.rows[Z]; z++) {
          for (y=0; y < workers_mpi.rows[Y]; y++) {
            index        = (x + workers_mpi.offset_x)*workers_mpi.layer_size   + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y + workers_mpi.offset_y);
            fread(&value, sizeof(double), 1, fp);
            iter_gridinfo_w[index].disp[dim][2] = value;
          }
        }
      }
    }
  }
//   }
  if (!ISOTHERMAL) {
    for (x=0; x < workers_mpi.rows[X]; x++) {
      for (z=0; z < workers_mpi.rows[Z]; z++) {
        for (y=0; y < workers_mpi.rows[Y]; y++) {
          index        = (x + workers_mpi.offset_x)*workers_mpi.layer_size   + (z + workers_mpi.offset_z)*workers_mpi.rows_y + (y + workers_mpi.offset_y);
          fread(&value, sizeof(double), 1, fp);
          gridinfo_w[index].temperature = value;
        }
      }
    }
  }
}
// void populate_table_names(){
//   long i, a, b, k;
//   char chempot_name[100];
//   char composition_name[100]; 
//   char phase_name[100];
//   
//   size_fields = NUMPHASES + (NUMCOMPONENTS-1);
//   
// //   if (WRITECOMPOSITION) {
//     size_fields += (NUMCOMPONENTS-1);
// //   }
//   if(!ISOTHERMAL) {
//     size_fields += 1;
//   }
//   
//   coordNames   = (char**)malloc(sizeof(char*)*(size_fields));
//   i=0;
//   for (a = 0; a < NUMPHASES; a++) {
//     sprintf(phase_name, "/%s",Phases[a]);
//     coordNames[i] = (char*)malloc(sizeof(char)*strlen(phase_name));
//     strcpy(coordNames[i], phase_name);
//     i++;
//   }
//   for (k=0; k < NUMCOMPONENTS-1; k++) {
//     sprintf(chempot_name, "/Mu_%s",Components[k]);
//     coordNames[i] = (char*)malloc(sizeof(char)*strlen(chempot_name));
//     strcpy(coordNames[i], chempot_name);
//     i++;
//   }
// //   if (WRITECOMPOSITION) {
//     for (k=0; k < NUMCOMPONENTS-1; k++) {
//       sprintf(composition_name, "/Composition_%s",Components[k]);
//       coordNames[i] = (char*)malloc(sizeof(char)*strlen(composition_name));
//       strcpy(coordNames[i], composition_name);
//       i++;
//     }
// //   }
//   if (!ISOTHERMAL) {
//     coordNames[i] = (char*)malloc(sizeof(char)*(strlen("/T")+1));
//     strcpy(coordNames[i], "/T");
//   }
// }
#endif
