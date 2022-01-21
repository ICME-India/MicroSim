void solverloop(long *start, long *end) {
  long x;
  
  calculate_gradients(0, gradient);
  swaplayers();
  calculate_gradients(1, gradient);
  swaplayers();
  
  for(x=start[X]-2; x<=end[X]+2; x++) {
    calculate_gradients(x+1, gradient);    
    calculate_divergence_phasefield(x, gradient);
      
  //Updating concentrations for the layer x-1
    if (x > 1) {
      calculate_fluxes_concentration(x-1, gradient);
      calculate_divergence_concentration(x-1, gradient);
    }
    swaplayers();
  }
}
void solverloop_phasefield(long *start, long *end) {
  long x;
  
  calculate_gradients_phasefield(0, gradient);
  swaplayers();
  calculate_gradients_phasefield(1, gradient);
  swaplayers();
  
  for(x=start[X]-2; x<=end[X]+2; x++) {
    calculate_gradients_phasefield(x+1, gradient);
    calculate_divergence_phasefield(x, gradient);
    swaplayers();
  }
}
void solverloop_concentration(long *start, long *end) {
  long b,k;
  long x;
  long INTERFACE_POS;
  MAX_INTERFACE_POS = 0;
  
  for (b=0; b < NUMPHASES; b++) {
    global_max_min.rel_change_phi[b] = 0.0;
  }
  for (k=0; k < NUMCOMPONENTS-1; k++) {
    global_max_min.rel_change_mu[k] = 0.0;
  }
  
  calculate_gradients_concentration(0, gradient);
  swaplayers();
  calculate_gradients_concentration(1, gradient);
  swaplayers();
  
  for(x=start[X]-2; x<=end[X]+2; x++) {
    calculate_gradients_concentration(   x+1, gradient);
    calculate_gradients_phasefield(      x+1, gradient);      
  //Updating concentrations for the layer x-1
    if (x > 1) {
      calculate_fluxes_concentration(    x-1, gradient);
      calculate_divergence_concentration(x-1, gradient);
/***********************************************************************************/
      if (SHIFT) {
        INTERFACE_POS = check_SHIFT(x-1);
        if (INTERFACE_POS > MAX_INTERFACE_POS) {
          MAX_INTERFACE_POS = INTERFACE_POS;
        }
      }
/*********************************************************************************/
    }
    swaplayers();
  }
}
void  smooth(long *start, long *end) {
  long x;
  
  calculate_gradients(0, gradient);
  swaplayers();
  calculate_gradients(1, gradient);
  swaplayers();
  
  for(x=start[X]-2; x<=end[X]+2; x++) {
    calculate_gradients(x+1, gradient);
    calculate_divergence_phasefield_smooth(x, gradient);
      
  //Updating concentrations for the layer x-1
    if (x > 1) {
      calculate_divergence_concentration_smooth(x-1, gradient);
    }
    swaplayers();
  }
}
void  smooth_concentration(long *start, long *end) {
  long x;
  
  calculate_gradients(0, gradient);
  swaplayers();
  calculate_gradients(1, gradient);
  swaplayers();
  
  for(x=start[X]-2; x<=end[X]+2; x++) {
    calculate_gradients(x+1, gradient);      
  //Updating concentrations for the layer x-1
    if (x > 1) {
      calculate_divergence_concentration_smooth_concentration(x-1, gradient);
    }
    swaplayers();
  }
}
