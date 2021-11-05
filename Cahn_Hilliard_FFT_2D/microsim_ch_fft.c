#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h> 
#include <stdbool.h>
#include <sys/stat.h>
#include <complex.h>
#include <gsl/gsl_rng.h>
#include "fftw3.h"
#include "functions/global_vars.h"
#include "functions/functions.h"
#include "functions/matrix.h"
#include "functions/utility_functions.h"
#include "functions/filling.h"
#include "functions/reading_input_parameters.h"
#include "functions/free_variables.h"
#include "functions/fill_domain.h"
#include "solverloop/serialinfo_xy.h"
#include "solverloop/functions_fftw.h"
#include "solverloop/file_writer.h"

int main(int argc, char * argv[]) {

  reading_input_parameters(argv);

  serialinfo_xy();

  fill_domain(argv);

  prepfftw();
  
  mkdir("DATA",0777);
  writetofile_serial2D(gridinfo, argv, 0);

  //Time-loop
  for(t=1;t<=ntimesteps;t++) {

    printf("Time step = %d\n", t);

    evolve_fftw();

    if(t%saveT == 0) {
      printf("Writing time step = %d\n", t);
      if (ASCII == 0) {
        writetofile_serial2D_binary(gridinfo, argv, t);
      } else {
        writetofile_serial2D(gridinfo, argv, t);
      }
    }
    if (t%time_output == 0) {
      fprintf(stdout, "Time=%le\n", t*deltat);
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



