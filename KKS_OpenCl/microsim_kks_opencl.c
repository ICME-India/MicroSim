#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h> 
#include <stdbool.h>
#include <sys/stat.h>
#include <mpi.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_roots.h>
#include "tdbs/Thermo.h"
#include "tdbs/Thermo.c"
#include "solverloop/defines.h"
#include "solverloop/defines1.h"
#include "functions/global_vars.h"
#include "functions/CL_global_vars.h" 
#include "functions/CL_initialize_variables.h"
#include "functions/CL_device_kernel_build.h"
#include "functions/CL_buffer_allocation.h"
#include "functions/CL_create_kernel_args.h"
#include "functions/CL_Initialize_domain.h"
#include "functions/CL_kernel_init_temperature.h"
#include "functions/CL_DeviceToHost.h"
#include "solverloop/FunctionDefines_CL.h"
#include "solverloop/Function_Thermo_calls.h"
#include "solverloop/CL_Update_Temperature.h"
#include "solverloop/CL_Solve_phi_com.h"
#include "solverloop/CL_Solve_phi_com_Function_F_2.h"
#include "solverloop/CL_Solve_phi_com_Function_F_3.h"
#include "solverloop/CL_Solve_phi_com_Function_F_4.h"
#include "functions/functions.h"
#include "functions/matrix.h"
#include "functions/utility_functions.h"
#include "functions/functionH.h"
#include "functions/functionF_01.h"
#include "functions/functionF_02.h"
#include "functions/functionF_03.h"
#include "functions/CL_functionF_03_HostUpdate.h"
#include "functions/functionF_04.h"
#include "functions/CL_functionF_04_HostUpdate.h"
#include "solverloop/FunctionF_4_SplineCPU.h"
#include "functions/filling.h"
#include "functions/reading_input_parameters.h"
#include "functions/read_boundary_conditions.h"
#include "functions/initialize_variables.h"
#include "functions/free_variables.h"
#include "functions/fill_domain.h"
#include "functions/Temperature_gradient.h"
#include "functions/shift.h"
#include "functions/CL_Shift.h"
#include "solverloop/serialinfo_xy.h"
#include "solverloop/mpi_xy.h"
#include "solverloop/initialize_functions_solverloop.h"
#include "solverloop/file_writer.h" 

int main(int argc, char * argv[]) { 
  
  // MPI_Comm comm=MPI_COMM_WORLD;


  MPI_Init(&argc,&argv);
  MPI_Comm_size(comm,&numtasks);
  MPI_Comm_rank(comm,&rank);

  printf("No. of processors in execution %d\n", numtasks);
  
  reading_input_parameters(argv);

  if ( RESTART ) {
    if ( (numworkers != numtasks)  ) {
      printf("ERROR:\n");
      printf("No. of numworkers = %d (processors) in input file are not same as no. of processors executed = %d (np)\n", numworkers, numtasks);
      exit(1);
    }
  }
  
  initialize_variables();
  
  //serialinfo_xy();
  
  initialize_functions_solverloop();
  
  read_boundary_conditions(argv);

  if (!(FUNCTION_F == 2)) {
    init_propertymatrices(T);
    if ( FUNCTION_F == 3 ) { 
      propf3Hostupdate(&propf3);
    }
    else if ( FUNCTION_F == 4 ) { 
      propf4Hostupdate(&propf4);
    }
  }


  serialinfo_xy();
  if ((STARTTIME == 0) && (RESTART ==0)) {
    fill_domain(argv);
  }

  if ( rank == MASTER ) { 
    mkdir("DATA",0777);
    if (!WRITEHDF5){
     for (n=0; n < numtasks; n++) {
       sprintf(dirname,"DATA/Processor_%d",n);
       mkdir(dirname, 0777);
     }
    }
    else {
      printf("Writing VTK format\n");
      for (n=0; n < numtasks; n++) {
        sprintf(dirname,"DATA/Processor_%d",n);
        mkdir(dirname, 0777);
      }
    }

    //print_input_parameters(argv);
    //print_boundary_conditions(argv);
    
    //serialinfo_xy();

    //if ((STARTTIME == 0) && (RESTART ==0)) {
    //  fill_domain(argv);
    //}
  }

  populate_table_names();

  if(ELASTICITY) {
    for (b=0; b<NUMPHASES; b++) {
      stiffness_phase_n[b].C11 = stiffness_phase[b].C11/stiffness_phase[NUMPHASES-1].C44;
      stiffness_phase_n[b].C12 = stiffness_phase[b].C12/stiffness_phase[NUMPHASES-1].C44;
      stiffness_phase_n[b].C44 = stiffness_phase[b].C44/stiffness_phase[NUMPHASES-1].C44;
    }
  }

  mpi_xy(rank, argv);
  
  if ((STARTTIME == 0) && (RESTART == 0)) {
    printf("Starting a new simulation\n");
    //fill_domain(argv);
  }
  else {
    printf("Restarting a simulation\n");
    if (ASCII) {
      readfromfile_serialmpi(gridinfomN, argv, STARTTIME);
    } else {
      readfromfile_serialmpi_binary(gridinfomN, argv, STARTTIME);
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
        temperature_gradientY.gradient_OFFSET = (temperature_gradientY.gradient_OFFSET) + floor((temperature_gradientY.velocity)*(STARTTIME*deltat));
      }
    }
  }
    if (TEMPGRADY) {
      t = STARTTIME;
      shift_OFFSET = 0;
      BASE_POS    = (temperature_gradientY.gradient_OFFSET/deltay) - shift_OFFSET + ((temperature_gradientY.velocity/deltay)*(t*deltat));
      GRADIENT    = (temperature_gradientY.DeltaT)*deltay/(temperature_gradientY.Distance);
      temp_bottom = temperature_gradientY.base_temp - BASE_POS*GRADIENT;
      apply_temperature_gradientY(gridinfomN, shift_OFFSET, t);
    }

  CL_initialize_variables();

  CL_device_kernel_build();
  
  CL_Initialize_domain();
  
  CL_buffer_allocation();

  CL_create_kernel_args();
  
  //CL_kernel_init_temperature();

  //writetofile_mpi(gridinfomN, argv, 0+STARTTIME);

  printf("All Initializations are done\n");

  if((RESTART == 0) || (STARTTIME ==0)) {
    if (ASCII == 0) {
      writetofile_mpi_binary(gridinfomN, argv, 0);
    } else {
      writetofile_mpi(gridinfomN, argv, 0);
    }
  }
  
  //Time-loop
  for(t=1;t<=ntimesteps;t++) {
    
    if (rank==MASTER)
    printf("Timestep=%ld\n",t);

    CL_Solve_phi_com_Function();
    
    CL_Update_Temperature(t + STARTTIME);
    
    if (SHIFT) {
      if (t%100 == 0) {
        CL_Shift();
      }
    }
	
    if(t%saveT == 0) {
      CL_DeviceToHost();
      if (ASCII == 0) {
        writetofile_mpi_binary(gridinfomN, argv, t + STARTTIME);
      } else {
        writetofile_mpi(gridinfomN, argv, t + STARTTIME);
      }
      fp=fopen("DATA/shift.dat","a");
      fprintf(fp,"%ld %ld\n", t + STARTTIME, shift_OFFSET + shift_position);
      fclose(fp);
      if (rank==MASTER)
      printf("Written time step = %ld\n", t + STARTTIME);
    }
    if (t%time_output == 0 ) {
     CL_Global_Max_Min();
      for (b=0; b<NUMPHASES; b++) {
        MPI_Reduce(&global_max_min.phi_max[b],        &global_max_min1.phi_max[b],         1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&global_max_min.phi_min[b],        &global_max_min1.phi_min[b],         1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&global_max_min.rel_change_phi[b], &global_max_min1.rel_change_phi[b],  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      }
      for (k=0; k<(NUMCOMPONENTS-1); k++) {
        MPI_Reduce(&global_max_min.com_max[k],         &global_max_min1.com_max[k],          1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&global_max_min.com_min[k],         &global_max_min1.com_min[k],          1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&global_max_min.rel_change_com[k],  &global_max_min1.rel_change_com[k],   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      }
      if (rank == MASTER) {
     fprintf(stdout, "Time=%le\n", t*deltat);
     for (b=0; b<NUMPHASES; b++) {
       fprintf(stdout, "%*s, Max = %le, Min = %le, Relative_Change=%le\n", max_length, Phases[b], global_max_min1.phi_max[b], global_max_min1.phi_min[b], sqrt(global_max_min1.rel_change_phi[b]));
     }
     for (k=0; k<NUMCOMPONENTS-1; k++) {
       fprintf(stdout, "%*s, Max = %le, Min = %le, Relative_Change=%le\n", max_length, Components[k], global_max_min1.com_max[k], global_max_min1.com_min[k], sqrt(global_max_min1.rel_change_com[k]));
     }
     fprintf(stdout, "\n");
    }
    }
  }
  free_variables();
  if(rank==MASTER) {

  printf("\nCompleted simulations!!!\n\n");
  }
  MPI_Finalize();
}



