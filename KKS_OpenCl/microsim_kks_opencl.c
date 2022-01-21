#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h> 
#include <stdbool.h>
#include <sys/stat.h>
#include "functions/global_vars.h"
#include "functions/CL_global_vars.h" 
#include "functions/CL_initialize_variables.h"
#include "functions/CL_device_kernel_build.h"
#include "functions/CL_buffer_allocation.h"
#include "functions/CL_create_kernel_args.h"
#include "functions/CL_Initialize_domain.h"
#include "functions/CL_kernel_init_temperature.h"
#include "functions/CL_DeviceToHost.h"
#include "solverloop/CL_Update_Temperature.h"
#include "solverloop/CL_Solve_phi_com.h"
#include "functions/functions.h"
#include "functions/matrix.h"
#include "functions/utility_functions.h"
#include "functions/filling.h"
#include "functions/reading_input_parameters.h"
#include "functions/read_boundary_conditions.h"
#include "functions/free_variables.h"
#include "functions/fill_domain.h"
#include "functions/Temperature_gradient.h"
#include "functions/shift.h"
#include "functions/CL_Shift.h"
#include "solverloop/serialinfo_xy.h"
#include "solverloop/file_writer.h" 

int main(int argc, char * argv[]) { 
  
  mkdir("DATA",0777);

  reading_input_parameters(argv);
  
  serialinfo_xy();
  
  read_boundary_conditions(argv);
  
  if ((STARTTIME == 0) && (RESTART ==0)) {
    fill_domain(argv);
    if (TEMPGRADY) {
      t = STARTTIME;
      shift_OFFSET = 0;
      BASE_POS    = (temperature_gradientY.gradient_OFFSET/deltay) - shift_OFFSET + ((temperature_gradientY.velocity/deltay)*(t*deltat));
      GRADIENT    = (temperature_gradientY.DeltaT)*deltay/(temperature_gradientY.Distance);
      temp_bottom = temperature_gradientY.base_temp - BASE_POS*GRADIENT;
      apply_temperature_gradientY(gridinfo, shift_OFFSET, t);
    }
  } else {
    if (ASCII) {
      readfromfile_serial2D(gridinfo, argv, STARTTIME);
    } else {
      readfromfile_serial2D_binary(gridinfo, argv, STARTTIME);
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

  CL_initialize_variables();

  CL_device_kernel_build();
  
  CL_Initialize_domain();
  
  CL_buffer_allocation();

  CL_create_kernel_args();
  
  CL_kernel_init_temperature();

  if((RESTART == 0) || (STARTTIME ==0)) {
    if (ASCII == 0) {
      writetofile_serial2D_binary(gridinfo, argv, 0);
    } else {
      writetofile_serial2D(gridinfo, argv, 0);
    }
  } 
  
  //Time-loop
  for(t=1;t<=ntimesteps;t++) {
    
    printf("Timestep=%ld\n",t);

    CL_Solve_phi_com();
    
    CL_Update_Temperature(t + STARTTIME);
  
    if (t%100 == 0) {
      CL_Shift();
    }
	
    if(t%saveT == 0) {
      CL_DeviceToHost();
      if (ASCII == 0) {
        writetofile_serial2D_binary(gridinfo, argv, t + STARTTIME);
      } else {
        writetofile_serial2D(gridinfo, argv, t + STARTTIME);
      }
      fp=fopen("DATA/shift.dat","a");
      fprintf(fp,"%ld %ld\n", t + STARTTIME, shift_OFFSET + shift_position);
      fclose(fp);
      printf("Written time step = %ld\n", t + STARTTIME);
    }
    if (t%time_output == 0) {
      fprintf(stdout, "Time=%le\n", t*deltat);
      CL_Global_Max_Min();
      for (b=0; b<NUMPHASES-1; b++) {
        fprintf(stdout, "%*s, Max = %le, Min = %le, Relative_Change=%le\n", max_length, Phases[b], global_max_min.phi_max[b], global_max_min.phi_min[b], sqrt(global_max_min.rel_change_phi[b]));
      }
      for (k=0; k<NUMCOMPONENTS-1; k++) {
        fprintf(stdout, "%*s, Max = %le, Min = %le, Relative_Change=%le\n", max_length, Components[k], global_max_min.com_max[k], global_max_min.com_min[k], sqrt(global_max_min.rel_change_com[k]));
      }
      fprintf(stdout, "\n");
    }

  }
  free_variables();
}



