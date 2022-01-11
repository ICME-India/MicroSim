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
#include "functions/functionW_01.h"
#include "functions/functionW_02.h"
#include "functions/function_A_00.h"
#include "functions/function_A_01.h"
#include "functions/anisotropy_01.h"
#include "functions/functionTau.h"
#include "functions/functionD.h"
#include "functions/filling.h"
#include "functions/reading_input_parameters.h"
#include "functions/read_boundary_conditions.h"
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
#include "solverloop/initialize_functions_solverloop.h"
#include "solverloop/solverloop.h"
#include "solverloop/boundary_serial.h"
#include "solverloop/file_writer.h"

int main(int argc, char * argv[]) {
  mkdir("DATA",0777);
  reading_input_parameters(argv);
  initialize_variables();
  serialinfo_xy();
  initialize_functions_solverloop();
  read_boundary_conditions(argv);
  if ((STARTTIME == 0) && (RESTART ==0)) {
    fill_domain(argv);
  } else {
    if (ASCII) {
      readfromfile_serial2D(gridinfo, argv, STARTTIME);
      init_propertymatrices(T);
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
    } else {
      readfromfile_serial2D_binary(gridinfo, argv, STARTTIME);
      init_propertymatrices(T);
    }
  }
  apply_boundary_conditions(0);

  if((RESTART == 0) || (STARTTIME ==0)) {
    if (ASCII == 0) {
      writetofile_serial2D_binary(gridinfo, argv, 0);
    } else {
      writetofile_serial2D(gridinfo, argv, 0);
    }
  } 
//   else {
//     writetofile_serial2D_binary(gridinfo, argv, STARTTIME +1);
//   }
  //Preconditioning
  for(t=1; t<nsmooth; t++) {
    smooth(start, end);
    apply_boundary_conditions(0);
  }
  //Time-loop
  for(t=1;t<=ntimesteps;t++) {
    solverloop_phasefield(start, end);
    apply_boundary_conditions(0);
    solverloop_concentration(start,end);
    apply_boundary_conditions(0);
    
    if (TEMPGRADY) {
      BASE_POS    = (temperature_gradientY.gradient_OFFSET/deltay) - shift_OFFSET + ((temperature_gradientY.velocity/deltay)*(t*deltat));
      GRADIENT    = (temperature_gradientY.DeltaT)*deltay/(temperature_gradientY.Distance);
      temp_bottom = temperature_gradientY.base_temp - BASE_POS*GRADIENT;
      apply_temperature_gradientY(gridinfo, shift_OFFSET, t);
    }
  
    if (t%100 == 0) {
      if(SHIFT) {
        if(MAX_INTERFACE_POS > shiftj) {
          shift_ON = 1;
        }
        if (shift_ON) {
          apply_shiftY(gridinfo, MAX_INTERFACE_POS); 
          shift_OFFSET += (MAX_INTERFACE_POS - shiftj);
          shift_ON=0;
        }
      }
    }
    
    
    if(t%saveT == 0) {
      if (ASCII == 0) {
        writetofile_serial2D_binary(gridinfo, argv, t + STARTTIME);
      } else {
        writetofile_serial2D(gridinfo, argv, t + STARTTIME);
      }
      fp=fopen("DATA/shift.dat","a");
      fprintf(fp,"%ld %ld\n",t + STARTTIME, shift_OFFSET + shift_position);
      fclose(fp);
    }
    if (t%time_output == 0) {
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
  free_variables();
}



