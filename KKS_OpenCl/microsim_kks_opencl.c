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
#include "solverloop/serialinfo_xy.h"
#include "solverloop/file_writer.h" 

int main(int argc, char * argv[]) { 

  reading_input_parameters(argv);
  
  serialinfo_xy();
  
  read_boundary_conditions(argv);
  
  fill_domain(argv);

  CL_initialize_variables();

  CL_device_kernel_build();
  
  CL_Initialize_domain();
  
  CL_buffer_allocation();

  CL_create_kernel_args();
  
  CL_kernel_init_temperature();
  
  mkdir("DATA",0777);
  writetofile_serial2D(gridinfo, argv, 0);
  
  //Time-loop
  for(t=1;t<=ntimesteps;t++) {
    
    printf("Timestep=%d\n",t);
    
    CL_Update_Temperature();

    CL_Solve_phi_com();
	
    if(t%saveT == 0) {
      CL_DeviceToHost();
      if (ASCII == 0) {
        writetofile_serial2D_binary(gridinfo, argv, t);
      } else {
        writetofile_serial2D(gridinfo, argv, t);
      }
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



