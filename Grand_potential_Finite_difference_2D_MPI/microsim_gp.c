// #include <mpi.h>
#include <hdf5.h>
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
#include "functions/functionH.h"
#include "functions/functionF_01.h"
#include "functions/functionQ.h"
#include "functions/anisotropy_01.h"
#include "functions/functionW_01.h"
#include "functions/functionW_02.h"
#include "functions/function_A_00.h"
#include "functions/function_A_01.h"
// #include "functions/anisotropy_01.h"
#include "functions/functionTau.h"
#include "functions/functionD.h"
#include "functions/filling.h"
#include "functions/reading_input_parameters.h"
#include "functions/read_boundary_conditions.h"
#include "functions/print_input_parameters.h"
#include "functions/print_boundary_conditions.h"
#include "functions/initialize_variables.h"
#include "functions/free_variables.h"
#include "functions/fill_domain.h"
#include "solverloop/serialinfo_xy.h"
#include "solverloop/gradients.h"
#include "solverloop/simplex_projection.h"
#include "solverloop/calculate_gradients.h"
#include "solverloop/calculate_fluxes_concentration.h"
#include "solverloop/calculate_divergence_phasefield.h"
#include "solverloop/calculate_divergence_concentration.h"
#include "solverloop/initialize_functions_solverloop.h"
#include "solverloop/solverloop.h"
#include "solverloop/mpiinfo_xy.h"
#include "solverloop/boundary_mpi.h"
#include "solverloop/file_writer.h"

// void writetofile_worker();

int main(int argc, char * argv[]) {
  
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  
  reading_input_parameters(argv);
  initialize_variables();
  initialize_functions_solverloop();
  read_boundary_conditions(argv);
  init_propertymatrices(T);
  
  
  Build_derived_type(gridinfo_instance, &MPI_gridinfo);
  
  if(taskid==MASTER) {
    print_input_parameters(argv);
    print_boundary_conditions(argv);
    serialinfo_xy();
    fill_domain(argv);
    mkdir("DATA",0777);
  }
  populate_table_names();
  
  sprintf(dirname,"DATA/Processor_%d",taskid);
  
  mkdir(dirname, 0777);
  
  numworkers_x = atol(argv[4]);
  numworkers_y = atol(argv[5]);
  
  if(numtasks != numworkers_x*numworkers_y) {
    fprintf(stdin,"The domain decomposition does not correspond to the number of spawned tasks!!\n: number of processes=numworkers_x*numworkers_y");
    exit(0);
  }
  
  Mpiinfo(taskid);
    
  mpiexchange_left_right(taskid);
  mpiexchange_top_bottom(taskid);
 
 if(boundary_worker) {
   apply_boundary_conditions(taskid);
 }
  if (!WRITEHDF5) {
      if ((ASCII == 0)) {
        writetofile_mpi2D_binary(gridinfo_w, argv, 0);
      } else {
        writetofile_mpi2D(gridinfo_w, argv, 0);
      }
  } else {
    writetofile_mpi2D_hdf5(gridinfo_w, argv, 0);
  }

//   writetofile_worker();
  
  //Preconditioning
  for(t=1; t<nsmooth; t++) {
    smooth(workers_mpi.start, workers_mpi.end);
    mpiexchange_left_right(taskid);
    mpiexchange_top_bottom(taskid);
    if(boundary_worker) {
      apply_boundary_conditions(taskid);
    }
  }
  printf("Finished Smoothing:%d\n", taskid);
 
  for(t=1;t<=ntimesteps;t++) {
    mpiexchange_left_right(taskid);
    mpiexchange_top_bottom(taskid);
    
    solverloop_phasefield(workers_mpi.start, workers_mpi.end);
    
    if(boundary_worker) {
      apply_boundary_conditions(taskid);
    }
    mpiexchange_left_right(taskid);
    mpiexchange_top_bottom(taskid);
    
    solverloop_concentration(workers_mpi.start,workers_mpi.end);
    
    if(boundary_worker) {
      apply_boundary_conditions(taskid);
    }
    
    if(t%saveT == 0) {
      if (!WRITEHDF5) {
        if ((ASCII == 0)) {
          writetofile_mpi2D_binary(gridinfo_w, argv, t);
        } else {
          writetofile_mpi2D(gridinfo_w, argv, t);
        }
      } else {
        writetofile_mpi2D_hdf5(gridinfo_w, argv, t);
      }
    }
    //printf("Iteration=%d\n",t);
    if (t%time_output == 0) {
      for (b=0; b<NUMPHASES; b++) {
        MPI_Reduce(&workers_max_min.phi_max[b],        &global_max_min.phi_max[b],         1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&workers_max_min.phi_min[b],        &global_max_min.phi_min[b],         1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&workers_max_min.rel_change_phi[b], &global_max_min.rel_change_phi[b],  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      }
      for (k=0; k<(NUMCOMPONENTS-1); k++) {
        MPI_Reduce(&workers_max_min.mu_max[k],         &global_max_min.mu_max[k],          1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&workers_max_min.mu_min[k],         &global_max_min.mu_min[k],          1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&workers_max_min.rel_change_mu[k],  &global_max_min.rel_change_mu[k],   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      }
      if (taskid == MASTER) {
        fprintf(stdout, "Time=%le\n", t*deltat);
        for (b=0; b<NUMPHASES; b++) {
          fprintf(stdout, "%*s, Max = %le, Min = %le, Relative_Change=%le\n", max_length, Phases[b], global_max_min.phi_max[b], global_max_min.phi_min[b], sqrt(global_max_min.rel_change_phi[b]));
        }
        for (k=0; k<NUMCOMPONENTS-1; k++) {
          fprintf(stdout, "%*s, Max = %le, Min = %le, Relative_Change=%le\n", max_length, Components[k], global_max_min.mu_max[k], global_max_min.mu_min[k], sqrt(global_max_min.rel_change_mu[k]));
        }
        fprintf(stdout, "\n");
      }
    }
  }
  free_variables();
  
  
  if (taskid == MASTER) {
    index_count = layer_size*rows_x;
    for (index_=0; index_ < index_count; index_++) {
      if ((&gridinfo[index_]) !=NULL) {
        free_memory_fields(&gridinfo[index_]);
      }
    }
    free(gridinfo);
  }
  for(i=0; i<size_fields; i++) {
    free(coordNames[i]);
  }
  MPI_Type_free(&MPI_gridinfo_vector_b);
  
  MPI_Type_free(&MPI_gridinfo);
  MPI_Finalize();
}

// void writetofile_worker() {
//   FILE *fp;
//   long index;
//   long x,y;
//   char name[100];
//   long start_x;
//   long start_y;
//   sprintf(name,"Worker_%d.dat",taskid);
//   fp = fopen(name,"w");
//   
//   fprintf(fp,"%ld\n",workers_mpi.rows[X]);
//   fprintf(fp,"%ld\n",workers_mpi.rows[Y]);
//   for(x=0;x < workers_mpi.rows[X]; x++) {
//     for(y=0;y < workers_mpi.rows[Y]; y++) {
//       index = (x + workers_mpi.offset_x)*workers_mpi.layer_size + (y + workers_mpi.offset_y);
//       fprintf(fp,"%ld %ld %le\n",x + workers_mpi.offset[X],y + workers_mpi.offset[Y], gridinfo_w[index].phia[0]);
//     }
//     fprintf(fp,"\n");
//   }
//   fprintf(fp,"\n");
//   fclose(fp);
// }



