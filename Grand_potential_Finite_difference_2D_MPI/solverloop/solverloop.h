#ifndef SOLVERLOOP_H_
#define SOLVERLOOP_H_

void solverloop_phasefield_tdb(long *start, long *end);
void solverloop_concentration_tdb(long *start, long *end);

// void solverloop(long *start, long *end) {
//   long x;
//   
//   calculate_gradients(0, gradient);
//   swaplayers();
//   calculate_gradients(1, gradient);
//   swaplayers();
//   
//   for(x=start[X]-2; x<=end[X]+2; x++) {
//     calculate_gradients(x+1, gradient);    
//     calculate_divergence_phasefield(x, gradient);
//       
//   //Updating concentrations for the layer x-1
//     if (x > 1) {
//       calculate_fluxes_concentration(x-1, gradient);
//       calculate_divergence_concentration(x-1, gradient);
//     }
//     swaplayers();
//   }
// }
// void solverloop_phasefield_F_01(long *start, long *end) {
//   long x;
//   
//   calculate_gradients_phasefield(0, gradient, 1);
//   swaplayers();
//   calculate_gradients_phasefield(1, gradient, 1);
//   swaplayers();
//   
//   for(x=start[X]-2; x<=end[X]+2; x++) {
//     calculate_gradients_phasefield(x+1, gradient, 1);
//     calculate_divergence_phasefield(x, gradient);
//     swaplayers();
//   }
// }
// void solverloop_concentration_F_01(long *start, long *end) {
//   long b,k;
//   long x;
//   long INTERFACE_POS;
//   workers_max_min.INTERFACE_POS_MAX = 0;
//   
//   for (b=0; b < NUMPHASES; b++) {
//     workers_max_min.rel_change_phi[b] = 0.0;
//   }
//   for (k=0; k < NUMCOMPONENTS-1; k++) {
//     workers_max_min.rel_change_mu[k] = 0.0;
//   }
//   
//   calculate_gradients_concentration(0, gradient);
//   swaplayers();
//   calculate_gradients_concentration(1, gradient);
//   swaplayers();
//   
//   for(x=start[X]-2; x<=end[X]+2; x++) {
//     calculate_gradients_concentration(   x+1, gradient);
//     calculate_gradients_phasefield(      x+1, gradient);      
//   //Updating concentrations for the layer x-1
//     if (x > 1) {
//       calculate_fluxes_concentration(    x-1, gradient);
//       calculate_divergence_concentration(x-1, gradient);
// /***********************************************************************************/
//       if(SHIFT) {
//         //Check condition for the shift only for the lowest level
//         INTERFACE_POS = check_SHIFT(x-1);
//         if (INTERFACE_POS > workers_max_min.INTERFACE_POS_MAX) {
//           workers_max_min.INTERFACE_POS_MAX = INTERFACE_POS;
//         }
//       }
// /*********************************************************************************/
//     }
//     swaplayers();
//   }
// }
void iterative_stress_solver(long *start, long *end) {
  calculate_gradients_stress(0, gradient);
  swaplayers();
  calculate_gradients_stress(1, gradient);
  swaplayers();
  long x;
  for(x=start[X]-2; x<=end[X]+2; x++) {
    calculate_gradients_stress (x+1, gradient);
    calculate_divergence_stress(x,   gradient);
    swaplayers();
  }
}
void solverloop_phasefield_tdb(long *start, long *end) {
  long x, b;
  long INTERFACE_POS;
  workers_max_min.INTERFACE_POS_MAX = 0;
  
  for (b=0; b < NUMPHASES; b++) {
    workers_max_min.rel_change_phi[b] = 0.0;
  }
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    workers_max_min.rel_change_mu[k] = 0.0;
  }
  
//   printf("rank=%d, rank_x=%d, rank_y=%d, rank_z=%d\n", taskid, workers_mpi.rank_x, workers_mpi.rank_y, workers_mpi.rank_z);
  if ((FUNCTION_F !=5) && (!GRAIN_GROWTH)) {
    calculate_diffusion_potential( 0, gradient);
  }
  calculate_gradients_phasefield(0, gradient, 0);
  swaplayers();
  
  if ((FUNCTION_F !=5) && (!GRAIN_GROWTH)) {
    calculate_diffusion_potential( 1, gradient);
  }
  calculate_gradients_phasefield(1, gradient, 0);
  swaplayers();
  
  for(x=start[X]-2; x<=end[X]+2; x++) {
    if ((FUNCTION_F != 5) && (!GRAIN_GROWTH)) {
      calculate_diffusion_potential( x+1, gradient);
    }
    calculate_gradients_phasefield(x+1, gradient, 0);
    calculate_divergence_phasefield(x, gradient);
    if ((FUNCTION_F == 5)) {
      if(SHIFT) {
        //Check condition for the shift only for the lowest level
        INTERFACE_POS = check_SHIFT(x-1);
        if (INTERFACE_POS > workers_max_min.INTERFACE_POS_MAX) {
          workers_max_min.INTERFACE_POS_MAX = INTERFACE_POS;
        }
      }
    }
    swaplayers();
  }
}
void solverloop_concentration_tdb(long *start, long *end) {
  long b,k;
  long x;
  long INTERFACE_POS;
  workers_max_min.INTERFACE_POS_MAX = 0;
  
  for (b=0; b < NUMPHASES; b++) {
    workers_max_min.rel_change_phi[b] = 0.0;
  }
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    workers_max_min.rel_change_mu[k] = 0.0;
  }
//   for (b=0; b < NUMPHASES; b++) {
//     workers_max_min.rel_change_phi[b] = 0.0;
//   }
//   for (k=0; k < NUMCOMPONENTS-1; k++) {
//     workers_max_min.rel_change_mu[k] = 0.0;
//   }
  
  calculate_gradients_phasefield(   0, gradient, 1); 
  swaplayers();
  calculate_gradients_phasefield(   1, gradient, 1); 
  calculate_gradients_concentration(0, gradient);
//   calculate_gradients_concentration(1, gradient);
  swaplayers();
  
  for(x=start[X]-2; x<=end[X]+2; x++) {
    calculate_gradients_phasefield(      x+1, gradient, 1); 
    calculate_gradients_concentration(   x,   gradient);
  //Updating concentrations for the layer x-1
    if (x > 1) {
      calculate_fluxes_concentration(    x-1, gradient);
      calculate_divergence_concentration(x-1, gradient);
/***********************************************************************************/
      if(SHIFT) {
        //Check condition for the shift only for the lowest level
        INTERFACE_POS = check_SHIFT(x-1);
        if (INTERFACE_POS > workers_max_min.INTERFACE_POS_MAX) {
          workers_max_min.INTERFACE_POS_MAX = INTERFACE_POS;
        }
      }
/*********************************************************************************/
    }
    swaplayers();
  }
}
void  smooth(long *start, long *end) {
  long x;
//   calculate_diffusion_potential( 0, gradient);
  calculate_gradients_phasefield(   0, gradient, 1); 
  swaplayers();
//   calculate_diffusion_potential( 1, gradient);
  calculate_gradients_phasefield(   1, gradient, 1);
  if ((FUNCTION_F!=5) && (!GRAIN_GROWTH)) {
    calculate_gradients_concentration(0, gradient);
  }
//   calculate_gradients_concentration(1, gradient);
  swaplayers();
  
  for(x=start[X]-2; x<=end[X]+2; x++) {
//     calculate_gradients(x+1, gradient);
//     calculate_divergence_phasefield_smooth(x, gradient);
//     calculate_diffusion_potential( x+1, gradient);
    calculate_gradients_phasefield(        x+1, gradient, 1);
    if ((FUNCTION_F !=5) && (!GRAIN_GROWTH)) {
      calculate_gradients_concentration(     x,   gradient);
    }
    calculate_divergence_phasefield_smooth(x,   gradient);  
  //Updating concentrations for the layer x-1
    if (x > 1) {
      if ((FUNCTION_F !=5) && (!GRAIN_GROWTH)) {
        calculate_divergence_concentration_smooth(x-1, gradient);
      }
    }
    swaplayers();
  }
}
// void  smooth_concentration(long *start, long *end) {
//   long x;
//   
//   calculate_gradients(0, gradient);
//   swaplayers();
//   calculate_gradients(1, gradient);
//   swaplayers();
//   
//   for(x=start[X]-2; x<=end[X]+2; x++) {
//     calculate_gradients(x+1, gradient);      
//   //Updating concentrations for the layer x-1
//     if (x > 1) {
//       calculate_divergence_concentration_smooth_concentration(x-1, gradient);
//     }
//     swaplayers();
//   }
// }
#endif
