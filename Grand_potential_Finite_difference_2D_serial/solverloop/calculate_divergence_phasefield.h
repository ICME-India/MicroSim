#ifndef CALCULATE_DIVERGENCE_PHASEFIELD_H_
#define CALCULATE_DIVERGENCE_PHASEFIELD_H_

void calculate_divergence_phasefield_2D(long x, struct gradlayer **gradient) {
  for (gidy=1; gidy <= (end[Y]+2); gidy++) {
    center        =  gidy + x*layer_size;
    
    grad          =  &gradient[1][gidy];
    grad_left     =  &gradient[1][gidy-1];
    grad_back     =  &gradient[0][gidy];
    
    
    for (a=0; a < NUMPHASES; a++) {
      divphi[a]   =  (grad->gradphi[X][a] - grad_back->gradphi[X][a])/deltax;
      divphi[a]  +=  (grad->gradphi[Y][a] - grad_left->gradphi[Y][a])/deltay;
    }
    
  //Time step iteration, phase-field
  //Check if divergence of the phase-field exists.
    //For the free energy description
    if (!ISOTHERMAL) {
      T = gridinfo[center].temperature;
//       init_propertymatrices(T);
    }
    sum_lambdaphi=0.0;
    active_phases=0.0;
    for (a=0; a < NUMPHASES; a++) {
      if (fabs(divphi[a]) > 0.0) {      
        init_propertymatrices(T);
        lambda_phi[a] =  epsilon*(-dAdphi(gridinfo[center].phia, gradient, gidy, a) + 
        divdAdgradphi(gradient, center, gidy, a)/*- 0.0*d2gradphi(gridinfo[center].phia, gradient, center, gidy, a)*/) - (1.0/epsilon)*dwdphi(gridinfo[center].phia, divphi,gradient, gidy, a) - dpsi(gridinfo[center].compi,T, gridinfo[center].phia,a);
        
        if(NOISE_PHASEFIELD) {
          lambda_phi[a] += AMP_NOISE_PHASE*((double)(rand())/(double)RAND_MAX)*gridinfo[center].phia[a]*(1.0-gridinfo[center].phia[a]); 
        }
        sum_lambdaphi += lambda_phi[a];
        active_phases++;
      }
    }
    
    if (active_phases) {
      sum_lambdaphi /= active_phases;
    }
    
    for (a=0; a < NUMPHASES; a++) {
      if (fabs(divphi[a]) > 0.0) {
        grad->deltaphi[a] =  deltat*(lambda_phi[a] - sum_lambdaphi)/(FunctionTau(gridinfo[center].phia)*epsilon);
      } else {
        grad->deltaphi[a] = 0.0;
      }
    }
    
    //Gibbs-simplex back projection: Check if the change in the phase-field 
    //larger than the phase-field value or the change makes the value 
    //larger than one
    projection_on_simplex(divphi);
//     projection_on_simplex_without_weights();
    for (a=0; a < NUMPHASES; a++) {
      gridinfo[center].deltaphi[a] = grad->deltaphi[a];
    }
  }
}
void calculate_divergence_phasefield_smooth_2D(long x, struct gradlayer **gradient) {
  long a, b;
 for (gidy=1; gidy <=(end[Y]+2); gidy++) {
  center        =  gidy + x*layer_size;
  
  grad          =  &gradient[1][gidy];
  grad_left     =  &gradient[1][gidy-1];
  grad_back     =  &gradient[0][gidy];
  
  
  for (a=0; a < NUMPHASES; a++) {
    divphi[a]   =  (grad->gradphi[X][a] - grad_back->gradphi[X][a])/deltax;
    divphi[a]  +=  (grad->gradphi[Y][a] - grad_left->gradphi[Y][a])/deltay;
  }
  
//Time step iteration, phase-field
//Check if divergence of the phase-field exists.
  //For the free energy description
  if (!ISOTHERMAL) {
    T = gridinfo[center].temperature;
    init_propertymatrices(T);
  }
  sum_lambdaphi=0.0;
  active_phases=0.0;
  for (a=0; a < NUMPHASES; a++) {
    if (fabs(divphi[a]) > 0.0) {      
      if (ANISOTROPY) {
        lambda_phi[a] =  epsilon*(-dAdphi_smooth(gridinfo[center].phia, gradient, gidy, a) + 
        divdAdgradphi_smooth(gradient, center, gidy, a)) 
        - (1.0/epsilon)*dwdphi_smooth(gridinfo[center].phia, divphi,gradient, gidy, a);
      } else {
      lambda_phi[a] =  epsilon*(-dAdphi(gridinfo[center].phia, gradient, gidy, a) + 
      divdAdgradphi(gradient, center, gidy, a)) 
      - (1.0/epsilon)*dwdphi(gridinfo[center].phia, divphi,gradient, gidy, a);
      }
      sum_lambdaphi += lambda_phi[a];
      active_phases++;
    }
  }
//   
  if (active_phases) {
    sum_lambdaphi /= active_phases;
  }
  
  for (a=0; a < NUMPHASES; a++) {
    if (fabs(divphi[a]) > 0.0) {
      grad->deltaphi[a] =  deltat*(lambda_phi[a] - sum_lambdaphi)/(FunctionTau(gridinfo[center].phia)*epsilon);
    } else {
      grad->deltaphi[a] = 0.0;
    }
  }
  projection_on_simplex_without_weights(divphi);
  //Gibbs-simplex back projection: Check if the change in the phase-field 
  //larger than the phase-field value or the change makes the value 
  //larger than one
//   Deltaphi = 0.0;
//   count_phases = 0;
//   
//   //Find the number of phases which are okay.
//   for (a=0; a < NUMPHASES; a++) {
//     if ((fabs(divphi[a]) > 0.0) && (gridinfo[center].phia[a] + grad->deltaphi[a]) > 0.0  && (gridinfo[center].phia[a] + grad->deltaphi[a]) < 1.0) {
//       count_phases++;
//     }
//   }
//   
//   for (a=0; a < NUMPHASES; a++) {
//     if (fabs(divphi[a]) > 0.0) {
//       if ((gridinfo[center].phia[a] + grad->deltaphi[a]) > 1.0) {
// // 				      Deltaphi =  (grad->deltaphi[a] - (1.0-gridinfo[center].phia[a]));
// // 				      grad->deltaphi[a] -= Deltaphi;
// 	grad->deltaphi[a]  = (1.0-gridinfo[center].phia[a]);
// 	
// 	//Correct all the other phases due to this correction,
// 	//If you are bulk, all other phases must go to zero.
// 	for (b=0; b < NUMPHASES; b++) {
// 	  if ((fabs(divphi[b]) > 0.0) && b!=a) {
// 	    grad->deltaphi[b] = -gridinfo[center].phia[b];
// 	  }
// 	}
// 	break;
//       }
//       if ((gridinfo[center].phia[a] + grad->deltaphi[a]) < 0.0) {
// 	Deltaphi = fabs(grad->deltaphi[a] + gridinfo[center].phia[a]);
// 	grad->deltaphi[a] += Deltaphi;
// 	
// 	//Correct all the phases due to the correction in the given phase
// 	  for (b=0; b < NUMPHASES; b++) {
// 	    if (((fabs(divphi[b]) > 0.0) && b!=a) && (gridinfo[center].phia[b] + grad->deltaphi[b]) > 0.0  && (gridinfo[center].phia[b] + grad->deltaphi[b]) < 1.0) {
// 	      grad->deltaphi[b] -= Deltaphi/(count_phases);
// 	    }
// 	  }
// 	}
//       }
//     }
  }
}
#endif