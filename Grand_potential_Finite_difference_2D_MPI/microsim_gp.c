// #include <mpi.h>
#include <hdf5.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h> 
#include <stdbool.h>
#include <sys/stat.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_roots.h>
#include "tdbs/Thermo.h"
#include "tdbs/Thermo.c"
#include "functions/global_vars.h"
#include "functions/functions.h"
#include "functions/matrix.h"
#include "functions/utility_functions.h"
#include "functions/functionH.h"
#include "functions/functionF_01.h"
#include "functions/functionF_02.h"
#include "functions/functionF_03.h"
#include "functions/functionF_04.h"
#include "functions/functionF_05.h"
#include "functions/functionF_elast.h"
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
#include "functions/shift.h"
#include "functions/Temperature_gradient.h"
#include "solverloop/serialinfo_xy.h"
#include "solverloop/gradients.h"
#include "solverloop/simplex_projection.h"
#include "solverloop/calculate_gradients.h"
#include "solverloop/calculate_fluxes_concentration.h"
#include "solverloop/calculate_divergence_phasefield.h"
#include "solverloop/calculate_divergence_concentration.h"
#include "solverloop/calculate_divergence_stress.h"
#include "solverloop/solverloop.h"
#include "solverloop/file_writer.h"
#include "solverloop/file_writer_3D.h"
// #include "solverloop/mpiinfo_xy.h"
#include "solverloop/mpiinfo_xyz.h"
#include "solverloop/boundary_mpi.h"
#include "solverloop/initialize_functions_solverloop.h"

// void writetofile_worker();

MPI_Request request;
int main(int argc, char * argv[]) {
  
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  
  reading_input_parameters(argv);
  initialize_variables();
  initialize_functions_solverloop();
  read_boundary_conditions(argv);
  if (!((FUNCTION_F == 2) || (FUNCTION_F==5) || (GRAIN_GROWTH))) {
    init_propertymatrices(T);
  }
  
  for (a=0; a<NUMPHASES-1; a++) {
    for(k=0; k<NUMCOMPONENTS-1; k++) {
      printf("slopes[a][NUMPHASES-1][k]=%le\n", slopes[a][NUMPHASES-1][k]);
      printf("slopes[a][a][k]=%le\n", slopes[a][a][k]);
    }
  }

  Build_derived_type(gridinfo_instance, &MPI_gridinfo);
  if (ELASTICITY) {
    Build_derived_type_stress(iter_gridinfo_w_instance, &MPI_iter_gridinfo);
  }
  
  
  if(DIMENSION == 2) {
    numworkers_x = atol(argv[4]);
    numworkers_y = atol(argv[5]);
    numworkers_z = 1;
    if(numtasks != numworkers_x*numworkers_y) {
      fprintf(stdin,"The domain decomposition does not correspond to the number of spawned tasks!!\n: number of processes=numworkers_x*numworkers_y");
      exit(0);
    }
  } else {
    numworkers_x = atol(argv[4]);
    numworkers_y = atol(argv[5]);
    numworkers_z = atol(argv[6]);
    printf("numworkers_x=%d, numworkers_y=%d, numworkers_z=%d\n",numworkers_x, numworkers_y,  numworkers_z);
    if(numtasks != numworkers_x*numworkers_y*numworkers_z) {
      fprintf(stdin,"The domain decomposition does not correspond to the number of spawned tasks!!\n: number of processes=numworkers_x*numworkers_y");
      exit(0);
    }
  }
  
  if(taskid==MASTER) {
    mkdir("DATA",0777);
    if (!WRITEHDF5){
     for (n=0; n < numtasks; n++) {
       sprintf(dirname,"DATA/Processor_%d",n);
       mkdir(dirname, 0777);
     }
    }
    print_input_parameters(argv);
    print_boundary_conditions(argv);
    serialinfo_xy();
    if ((STARTTIME == 0) && (RESTART ==0)) {
      fill_domain(argv);
    }
  }
 
  populate_table_names();
  
  if(ELASTICITY) {
    for (b=0; b<NUMPHASES; b++) {
      stiffness_phase_n[b].C11 = stiffness_phase[b].C11/stiffness_phase[NUMPHASES-1].C44;
      stiffness_phase_n[b].C12 = stiffness_phase[b].C12/stiffness_phase[NUMPHASES-1].C44;
      stiffness_phase_n[b].C44 = stiffness_phase[b].C44/stiffness_phase[NUMPHASES-1].C44;
    }
  }
  
//   if (DIMENSION == 2) {
//     Mpiinfo(taskid);
//   } else {
  Mpiinfo(taskid);
  
 
//   }
//   exit(0);
//   if (FUNCTION_F == 2) {
  
  
  if ((FUNCTION_F != 5) && (!GRAIN_GROWTH)) {
    Calculate_Tau();
  } else {
    if ((FUNCTION_F ==5) && (!GRAIN_GROWTH)) {
      for (a=0; a < NUMPHASES-1; a++) {
        tau_ab[a][NUMPHASES-1] = (Lf*Lf)*epsilon*0.2222/(V*V*therm_cond*Teq);
        tau_ab[NUMPHASES-1][a] = tau_ab[a][NUMPHASES-1];
      }
      for (a=0; a < NUMPHASES-1; a++) {
        for (b=0; b < NUMPHASES-1; b++) {
          tau_ab[a][b] = (Lf*Lf)*epsilon*0.2222/(V*V*therm_cond*Teq);
        }
      }
    }
//     if (GRAIN_GROWTH) {
//       
//     }
    printf("tau[0][NUMPHASES-1]=%le\n",tau_ab[0][NUMPHASES-1]);
//     exit(0);
  }

//   }
  
  //Checking tdb functions
  
//   if (taskid == MASTER) {
//     double c_x;
//     double c[NUMCOMPONENTS-1];
//     double c_calc[NUMCOMPONENTS-1];
//     double mu[NUMCOMPONENTS-1];
//     double dpsi;    
//     char filename[1000];
//     double fe;
//     FILE *fp_check;
//     for (a=0; a<NUMPHASES; a++) {
//       sprintf(filename, "Thermodynamic_functions_%ld.dat", a);
//       fp_check = fopen(filename, "w");
//       for(c_x=0.01; c_x < 0.99;) {
//         c[0] = c_x;
//         Mu(c, T, a, mu);
//         dc_dmu(mu, c, T, a, dcdmu);
//         fe = free_energy(c, T, a);
//         dpsi = fe - mu[0]*c[0];
//         c_mu(mu, c, T, a, ceq[a][a]);
//         fprintf(fp_check, "%le %le %le %le %le %le\n", c_x, mu[0], dcdmu[0][0], fe,  dpsi, c_calc[0]);
//         c_x += 0.05;
//       }
//       fclose(fp_check);
//     }
//   }
  
  
  if ((STARTTIME !=0) || (RESTART !=0)) {
    if (WRITEHDF5) {
      readfromfile_mpi_hdf5(gridinfo_w, argv, numworkers, STARTTIME);
    } else {
      if (ASCII) {
        readfromfile_mpi(gridinfo_w, argv, STARTTIME);
      } else {
        readfromfile_mpi_binary(gridinfo_w, argv, STARTTIME);
      }
    }
    if (SHIFT) {
      FILE *fp;
      fp = fopen("DATA/shift.dat","r");
      
      for(file_iter=0; file_iter <= STARTTIME/saveT; file_iter++) {
        fscanf(fp,"%ld %ld\n",&time_file, &position);
      }
      fclose(fp);
      shift_position = position;
    }
    if (!ISOTHERMAL) {
      if (SHIFT) {
         temperature_gradientY.gradient_OFFSET = (temperature_gradientY.gradient_OFFSET) + floor((temperature_gradientY.velocity)*(STARTTIME*deltat)) - shift_position*deltay;
      } else {
        temperature_gradientY.gradient_OFFSET  = (temperature_gradientY.gradient_OFFSET) + floor((temperature_gradientY.velocity)*(STARTTIME*deltat));
      }
    }
  }
  
  if(boundary_worker) {
   apply_boundary_conditions(taskid);
  }
    
  mpiexchange_left_right(taskid);
  mpiexchange_top_bottom(taskid);
  if (DIMENSION==3) {
    mpiexchange_front_back(taskid);
  }
  
  if (TEMPGRADY) {
    BASE_POS    = (temperature_gradientY.gradient_OFFSET/deltay) - shift_OFFSET;
    GRADIENT    = (temperature_gradientY.DeltaT)*deltay/(temperature_gradientY.Distance);
    temp_bottom = temperature_gradientY.base_temp - BASE_POS*GRADIENT + (workers_mpi.offset[Y]-workers_mpi.offset_y)*GRADIENT;
    apply_temperature_gradientY(gridinfo_w, shift_OFFSET, 0);
  }
   
  
  if (!WRITEHDF5) {
    if ((ASCII == 0)) {
      writetofile_mpi_binary(gridinfo_w, argv, 0 + STARTTIME);
    } else {
      writetofile_mpi(gridinfo_w, argv, 0 + STARTTIME);
    }
  } else {
    writetofile_mpi_hdf5(gridinfo_w, argv, 0 + STARTTIME);
  }
  
//   printf("I am coming here\n");
//   exit(0);
//   writetofile_worker();

//   Preconditioning
  for(t=1; t<nsmooth; t++) {
    smooth(workers_mpi.start, workers_mpi.end);
    if(boundary_worker) {
      apply_boundary_conditions(taskid);
    }
    mpiexchange_left_right(taskid);
    mpiexchange_top_bottom(taskid);
    if (DIMENSION == 3) {
      mpiexchange_front_back(taskid);
    }
  }
//   printf("Finished smoothing\n");
//   exit(0);
  
  if (!WRITEHDF5) {
    if ((ASCII == 0)) {
      writetofile_mpi_binary(gridinfo_w, argv, 0 + STARTTIME);
    } else {
      writetofile_mpi(gridinfo_w, argv, 0 + STARTTIME);
    }
  } else {
    writetofile_mpi_hdf5(gridinfo_w, argv, 0 + STARTTIME);
  }
  printf("Finished smoothing\n");
//   exit(0);

  
  for(t=1;t<=ntimesteps;t++) {
    mpiexchange_left_right(taskid);
    mpiexchange_top_bottom(taskid);
    if (DIMENSION == 3) {
      mpiexchange_front_back(taskid);
    }
    
    solverloop_phasefield(workers_mpi.start, workers_mpi.end);
    
    if(boundary_worker) {
      apply_boundary_conditions(taskid);
    }
    mpiexchange_left_right(taskid);
    mpiexchange_top_bottom(taskid);
    if (DIMENSION == 3) {
      mpiexchange_front_back(taskid);
    }
    
    if ((FUNCTION_F != 5) && (!GRAIN_GROWTH)) {
      solverloop_concentration(workers_mpi.start,workers_mpi.end);
    }
    
    if (TEMPGRADY) {
      BASE_POS    = (temperature_gradientY.gradient_OFFSET/deltay) - shift_OFFSET + ((temperature_gradientY.velocity/deltay)*(t*deltat));
      GRADIENT    = (temperature_gradientY.DeltaT)*deltay/(temperature_gradientY.Distance);
      temp_bottom = temperature_gradientY.base_temp - BASE_POS*GRADIENT + (workers_mpi.offset[Y]-workers_mpi.offset_y)*GRADIENT;
      apply_temperature_gradientY(gridinfo_w, shift_OFFSET, t);
    }
    
    if (ELASTICITY) {
      for(iter=1; iter < MAX_ITERATIONS; iter++) {		//elasticity solver
		     mpiexchange_top_bottom_stress(taskid);
		     mpiexchange_left_right_stress(taskid);
         if (DIMENSION ==3) {
           mpiexchange_front_back_stress(taskid);
         }
		     iterative_stress_solver(workers_mpi.start, workers_mpi.end);
         if (boundary_worker) {
		        apply_boundary_conditions_stress(taskid);
         }
		     
// 		     if ((iter%100)==0) {
//            error = 0.0;
// 		       for(x=workers_mpi.start[X]; x<=workers_mpi.end[X]; x++) {
// 		         compute_error(x, &error);
// 		       }
// 		       printf("error=%le\n", error);
// 		       MPI_Reduce(&error,  &global_error,   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
// 		       if (fabs(global_error) < tolerance) {
// 		        break;
// 		       }
//          }
      }
    }
  
    if (t%100 == 0) {
      if(SHIFT) {
        MPI_Iallreduce(&workers_max_min.INTERFACE_POS_MAX,  &INTERFACE_POS_GLOBAL,  1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD, &request);
        shift_ON = 0;
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        if(INTERFACE_POS_GLOBAL > shiftj) {
          shift_ON = 1;
        }
        if (shift_ON) {
          apply_shiftY(gridinfo_w, INTERFACE_POS_GLOBAL); 
//           if (taskid == MASTER) {
          shift_OFFSET += (INTERFACE_POS_GLOBAL - shiftj);
//           }
          mpiexchange_top_bottom(taskid);
        }
      }
    }
    
    if(boundary_worker) {
      apply_boundary_conditions(taskid);
    }
    
    if(t%saveT == 0) {
      if (!WRITEHDF5) {
        if ((ASCII == 0)) {
          writetofile_mpi_binary(gridinfo_w, argv, t + STARTTIME);
        } else {
          writetofile_mpi(gridinfo_w, argv, t + STARTTIME);
        }
      } else {
        writetofile_mpi_hdf5(gridinfo_w, argv, t + STARTTIME);
      }
      if(SHIFT) {
        if (taskid == MASTER) {
          fp=fopen("DATA/shift.dat","a");
          fprintf(fp,"%ld %ld\n",t + STARTTIME, shift_OFFSET + shift_position);
          fclose(fp);
        }
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
        fprintf(stdout, "Time=%le\n", t*deltat + STARTTIME);
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
  MPI_Type_free(&MPI_gridinfo_vector_c);
  
  MPI_Type_free(&MPI_gridinfo);
  
  if (ELASTICITY) {
    MPI_Type_free(&MPI_gridinfo_vector_b_stress);
    MPI_Type_free(&MPI_gridinfo_vector_c_stress);
    
    MPI_Type_free(&MPI_iter_gridinfo);
  }
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



