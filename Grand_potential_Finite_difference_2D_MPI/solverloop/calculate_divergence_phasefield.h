#ifndef CALCULATE_DIVERGENCE_PHASEFIELD_H_
#define CALCULATE_DIVERGENCE_PHASEFIELD_H_

// void calculate_divergence_phasefield_2D(long x, struct gradlayer **gradient) {
//   for (gidy=1; gidy <= (workers_mpi.end[Y]+2); gidy++) {
//     center        =  gidy + x*workers_mpi.layer_size;
//     
//     grad          =  &gradient[1][gidy];
//     grad_left     =  &gradient[1][gidy-1];
//     grad_back     =  &gradient[0][gidy];
//     
//     
//     for (a=0; a < NUMPHASES; a++) {
//       divphi[a]   =  (grad->gradphi[X][a] - grad_back->gradphi[X][a])/deltax;
//       divphi[a]  +=  (grad->gradphi[Y][a] - grad_left->gradphi[Y][a])/deltay;
//     }
//     
//   //Time step iteration, phase-field
//   //Check if divergence of the phase-field exists.
//     //For the free energy description
//     if (!ISOTHERMAL) {
// //       gridinfo_w[center].temperature  = temp_bottom + gidy*GRADIENT;
//       T = gridinfo_w[center].temperature;
//     }
//     sum_lambdaphi=0.0;
//     active_phases=0.0;
//     for (a=0; a < NUMPHASES; a++) {
//       if (fabs(divphi[a]) > 0.0) {      
//         if (!ISOTHERMAL) {
//           init_propertymatrices(T);
//         }
//         lambda_phi[a] =  epsilon*(-dAdphi(gridinfo_w[center].phia, gradient, gidy, a) + 
//         divdAdgradphi(gradient, center, gidy, a)/*- 0.0*d2gradphi(gridinfo_w[center].phia, gradient, center, gidy, a)*/) - (1.0/epsilon)*dwdphi(gridinfo_w[center].phia, divphi,gradient, gidy, a) - dpsi(gridinfo_w[center].compi,T, gridinfo_w[center].phia,a);
//         
//         if(NOISE_PHASEFIELD) {
//           lambda_phi[a] += AMP_NOISE_PHASE*((double)(rand())/(double)RAND_MAX)*gridinfo_w[center].phia[a]*(1.0-gridinfo_w[center].phia[a]); 
//         }
//         sum_lambdaphi += lambda_phi[a];
//         active_phases++;
//       }
//     }
//     
//     if (active_phases) {
//       sum_lambdaphi /= active_phases;
//     }
//     
//     for (a=0; a < NUMPHASES; a++) {
//       if (fabs(divphi[a]) > 0.0) {
//         grad->deltaphi[a] =  deltat*(lambda_phi[a] - sum_lambdaphi)/(FunctionTau(gridinfo_w[center].phia)*epsilon);
//       } else {
//         grad->deltaphi[a] = 0.0;
//       }
//     }
//     
//     //Gibbs-simplex back projection: Check if the change in the phase-field 
//     //larger than the phase-field value or the change makes the value 
//     //larger than one
//     projection_on_simplex(divphi);
// //     projection_on_simplex_without_weights();
//     for (a=0; a < NUMPHASES; a++) {
//       gridinfo_w[center].deltaphi[a] = grad->deltaphi[a];
//     }
//   }
// }
void calculate_divergence_phasefield_smooth_2D(long x, struct gradlayer **gradient) {
  long a, b;
 for (gidy=1; gidy <=(workers_mpi.end[Y]+2); gidy++) {
  center        =  gidy + x*workers_mpi.layer_size;
  
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
//   if (!ISOTHERMAL) {
//     T = gridinfo_w[center].temperature;
// //     init_propertymatrices(T);
//   }
  sum_lambdaphi=0.0;
  active_phases=0.0;
  for (a=0; a < NUMPHASES; a++) {
    if (fabs(divphi[a]) > 0.0) {      
      if (ANISOTROPY) {
        lambda_phi[a] =  epsilon*(-dAdphi_smooth(gridinfo_w[center].phia, gradient, gidy, a) + 
        divdAdgradphi_smooth(gradient, center, gidy, a)) 
        - (1.0/epsilon)*dwdphi_smooth(gridinfo_w[center].phia, divphi,gradient, gidy, a);
      } else {
      lambda_phi[a] =  epsilon*(-dAdphi(gridinfo_w[center].phia, gradient, gidy, a) + 
      divdAdgradphi(gradient, center, gidy, a)) 
      - (1.0/epsilon)*dwdphi(gridinfo_w[center].phia, divphi,gradient, gidy, a);
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
      grad->deltaphi[a] =  deltat*(lambda_phi[a] - sum_lambdaphi)/(FunctionTau(gridinfo_w[center].phia)*epsilon);
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
//     if ((fabs(divphi[a]) > 0.0) && (gridinfo_w[center].phia[a] + grad->deltaphi[a]) > 0.0  && (gridinfo_w[center].phia[a] + grad->deltaphi[a]) < 1.0) {
//       count_phases++;
//     }
//   }
//   
//   for (a=0; a < NUMPHASES; a++) {
//     if (fabs(divphi[a]) > 0.0) {
//       if ((gridinfo_w[center].phia[a] + grad->deltaphi[a]) > 1.0) {
// // 				      Deltaphi =  (grad->deltaphi[a] - (1.0-gridinfo_w[center].phia[a]));
// // 				      grad->deltaphi[a] -= Deltaphi;
// 	grad->deltaphi[a]  = (1.0-gridinfo_w[center].phia[a]);
// 	
// 	//Correct all the other phases due to this correction,
// 	//If you are bulk, all other phases must go to zero.
// 	for (b=0; b < NUMPHASES; b++) {
// 	  if ((fabs(divphi[b]) > 0.0) && b!=a) {
// 	    grad->deltaphi[b] = -gridinfo_w[center].phia[b];
// 	  }
// 	}
// 	break;
//       }
//       if ((gridinfo_w[center].phia[a] + grad->deltaphi[a]) < 0.0) {
// 	Deltaphi = fabs(grad->deltaphi[a] + gridinfo_w[center].phia[a]);
// 	grad->deltaphi[a] += Deltaphi;
// 	
// 	//Correct all the phases due to the correction in the given phase
// 	  for (b=0; b < NUMPHASES; b++) {
// 	    if (((fabs(divphi[b]) > 0.0) && b!=a) && (gridinfo_w[center].phia[b] + grad->deltaphi[b]) > 0.0  && (gridinfo_w[center].phia[b] + grad->deltaphi[b]) < 1.0) {
// 	      grad->deltaphi[b] -= Deltaphi/(count_phases);
// 	    }
// 	  }
// 	}
//       }
//     }
  }
}
void calculate_divergence_phasefield_smooth_3D(long x, struct gradlayer **gradient) {
  long a, b;
 for (gidy=workers_mpi.rows_y+1; gidy <= (workers_mpi.layer_size - workers_mpi.rows_y-2); gidy++) {
  center        =  gidy + x*workers_mpi.layer_size;
  
  grad          =  &gradient[1][gidy];
  grad_left     =  &gradient[1][gidy-1];
  grad_back     =  &gradient[0][gidy];

  if ((gidy-workers_mpi.rows_y) >= 0) {
    grad_bottom  = &gradient[1][gidy-workers_mpi.rows_y];
  }
  
  for (a=0; a < NUMPHASES; a++) {
    divphi[a]   =  (grad->gradphi[X][a] - grad_back->gradphi[X][a])/deltax;
    divphi[a]  +=  (grad->gradphi[Y][a] - grad_left->gradphi[Y][a])/deltay;
    divphi[a]  +=  (grad->gradphi[Z][a] - grad_bottom->gradphi[Z][a])/deltaz;
  }
  
//Time step iteration, phase-field
//Check if divergence of the phase-field exists.
  //For the free energy description
//   if (!ISOTHERMAL) {
//     T = gridinfo_w[center].temperature;
// //     init_propertymatrices(T);
//   }
  sum_lambdaphi=0.0;
  active_phases=0.0;
  for (a=0; a < NUMPHASES; a++) {
    if (fabs(divphi[a]) > 0.0) {      
      if (ANISOTROPY) {
        lambda_phi[a] =  epsilon*(-dAdphi_smooth(gridinfo_w[center].phia, gradient, gidy, a) + 
        divdAdgradphi_smooth(gradient, center, gidy, a)) 
        - (1.0/epsilon)*dwdphi_smooth(gridinfo_w[center].phia, divphi,gradient, gidy, a);
      } else {
      lambda_phi[a] =  epsilon*(-dAdphi(gridinfo_w[center].phia, gradient, gidy, a) + 
      divdAdgradphi(gradient, center, gidy, a)) 
      - (1.0/epsilon)*dwdphi(gridinfo_w[center].phia, divphi,gradient, gidy, a);
      }
      sum_lambdaphi += lambda_phi[a];
      active_phases++;
    }
  }
//   
  if (active_phases) {
    sum_lambdaphi /= active_phases;
  }
//   
  for (a=0; a < NUMPHASES; a++) {
    if (fabs(divphi[a]) > 0.0) {
      grad->deltaphi[a] =  deltat*(lambda_phi[a] - sum_lambdaphi)/(FunctionTau(gridinfo_w[center].phia)*epsilon);
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
//     if ((fabs(divphi[a]) > 0.0) && (gridinfo_w[center].phia[a] + grad->deltaphi[a]) > 0.0  && (gridinfo_w[center].phia[a] + grad->deltaphi[a]) < 1.0) {
//       count_phases++;
//     }
//   }
//   
//   for (a=0; a < NUMPHASES; a++) {
//     if (fabs(divphi[a]) > 0.0) {
//       if ((gridinfo_w[center].phia[a] + grad->deltaphi[a]) > 1.0) {
// // 				      Deltaphi =  (grad->deltaphi[a] - (1.0-gridinfo_w[center].phia[a]));
// // 				      grad->deltaphi[a] -= Deltaphi;
// 	grad->deltaphi[a]  = (1.0-gridinfo_w[center].phia[a]);
// 	
// 	//Correct all the other phases due to this correction,
// 	//If you are bulk, all other phases must go to zero.
// 	for (b=0; b < NUMPHASES; b++) {
// 	  if ((fabs(divphi[b]) > 0.0) && b!=a) {
// 	    grad->deltaphi[b] = -gridinfo_w[center].phia[b];
// 	  }
// 	}
// 	break;
//       }
//       if ((gridinfo_w[center].phia[a] + grad->deltaphi[a]) < 0.0) {
// 	Deltaphi = fabs(grad->deltaphi[a] + gridinfo_w[center].phia[a]);
// 	grad->deltaphi[a] += Deltaphi;
// 	
// 	//Correct all the phases due to the correction in the given phase
// 	  for (b=0; b < NUMPHASES; b++) {
// 	    if (((fabs(divphi[b]) > 0.0) && b!=a) && (gridinfo_w[center].phia[b] + grad->deltaphi[b]) > 0.0  && (gridinfo_w[center].phia[b] + grad->deltaphi[b]) < 1.0) {
// 	      grad->deltaphi[b] -= Deltaphi/(count_phases);
// 	    }
// 	  }
// 	}
//       }
//     }
  }
}

void calculate_divergence_phasefield_2D(long x, struct gradlayer **gradient) {
  double noise;
  struct symmetric_tensor strain, sigma;
  
  for (gidy=1; gidy <= (workers_mpi.end[Y]+2); gidy++) {
    center        =  gidy + x*workers_mpi.layer_size;
    
    front         =  gidy   + (x+1)*workers_mpi.rows_y;				//front         =  gidy   + (x+1)*rows_y
    back          =  gidy   + (x-1)*workers_mpi.rows_y;				//back          =  gidy   + (x-1)*rows_y
    right         =  center + 1;
    left          =  center - 1;
    
    
    grad          =  &gradient[1][gidy];
    grad_left     =  &gradient[1][gidy-1];
    grad_back     =  &gradient[0][gidy];
    
    
    for (a=0; a < NUMPHASES; a++) {
      divphi[a]   =  (grad->gradphi[X][a] - grad_back->gradphi[X][a])/deltax;
      divphi[a]  +=  (grad->gradphi[Y][a] - grad_left->gradphi[Y][a])/deltay;
    }
    
    if (ELASTICITY) {
      grad->eigen_strain[X] = calculate_eigen_strain(center);
      
      grad->strain[X].xx = ((0.5)*(iter_gridinfo_w[front].disp[X][2] - iter_gridinfo_w[back].disp[X][2]) - grad->eigen_strain[X].xx);
      grad->strain[X].yy = ((0.5)*(iter_gridinfo_w[right].disp[Y][2] - iter_gridinfo_w[left].disp[Y][2]) - grad->eigen_strain[X].yy);
      grad->strain[X].xy = (0.25*((iter_gridinfo_w[right].disp[X][2] - iter_gridinfo_w[left].disp[X][2]) 
                        + (iter_gridinfo_w[front].disp[Y][2] - iter_gridinfo_w[back].disp[Y][2])));
      
      grad->stiffness_c[X] = calculate_stiffness(center);
  //     grad->stiffness_c[Y] = calculate_stiffness(center);
      
      
      sigma.xx  = grad->stiffness_c[X].C11*(grad->strain[X].xx)      + grad->stiffness_c[X].C12*(grad->strain[X].yy);
      sigma.yy  = grad->stiffness_c[X].C12*(grad->strain[X].xx)      + grad->stiffness_c[X].C11*(grad->strain[X].yy);
      sigma.xy  = 2.0*grad->stiffness_c[X].C44*(grad->strain[X].xy);
    }
    
    
  //Time step iteration, phase-field
  //Check if divergence of the phase-field exists.
    //For the free energy description
    if (!ISOTHERMAL) {
//       gridinfo_w[center].temperature  = temp_bottom + gidy*GRADIENT;
      T = gridinfo_w[center].temperature;
    }
    sum_lambdaphi=0.0;
    active_phases=0.0;
    
//     interface = 1;
//     
//     for (a=0; a < NUMPHASES; a++) {
//       if (gridinfo_w[center].phia[a] == 1.0) {
// 	bulk_phase=a;
// 	interface = 0;
// 	break;
//       }
//     }
//     noise = 4.0*AMP_NOISE_PHASE*((double)(rand())/(double)RAND_MAX);
    noise = 4.0*AMP_NOISE_PHASE*drand48();
    for (a=0; a < NUMPHASES; a++) {
      if (fabs(divphi[a]) > 0.0) {      
        lambda_phi[a] =  epsilon*(-dAdphi(gridinfo_w[center].phia, gradient, gidy, a) + 
        divdAdgradphi(gradient, center, gidy, a)/*- 0.0*d2gradphi(gridinfo_w[center].phia, gradient, center, gidy, a)*/) - (1.0/epsilon)*dwdphi(gridinfo_w[center].phia, divphi,gradient, gidy, a) - dpsi(gridinfo_w[center].compi, grad->phase_comp, T, gridinfo_w[center].phia,a)/V;
        
        if(NOISE_PHASEFIELD) {
//           lambda_phi[a] += AMP_NOISE_PHASE*((double)(rand())/(double)RAND_MAX)*gridinfo_w[center].phia[a]*(1.0-gridinfo_w[center].phia[a]); 
          lambda_phi[a] += noise*gridinfo_w[center].phia[a]*(1.0-gridinfo_w[center].phia[a])*dpsi(gridinfo_w[center].compi, grad->phase_comp, T, gridinfo_w[center].phia,a)/V;
        }
        if (ELASTICITY) {
          lambda_phi[a] -= df_elast(grad, sigma, gridinfo_w[center].phia, a);
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
        grad->deltaphi[a] =  deltat*(lambda_phi[a] - sum_lambdaphi)/(FunctionTau(gridinfo_w[center].phia)*epsilon);
      } else {
        grad->deltaphi[a] = 0.0;
      }
    }
    
    //Gibbs-simplex back projection: Check if the change in the phase-field 
    //larger than the phase-field value or the change makes the value 
    //larger than one
    projection_on_simplex(gridinfo_w[center].phia, grad->deltaphi, divphi);
//     projection_on_simplex_without_weights();
    for (a=0; a < NUMPHASES; a++) {
      gridinfo_w[center].deltaphi[a] = grad->deltaphi[a];
    }
  }
}
void calculate_divergence_phasefield_3D(long x, struct gradlayer **gradient) {
  double noise;
  struct symmetric_tensor strain, sigma;
  for (gidy=workers_mpi.rows_y+1; gidy <=(workers_mpi.layer_size - workers_mpi.rows_y-2); gidy++) {
    center        =  gidy + x*workers_mpi.layer_size;
    
    right         = center + 1;
    left          = center - 1;
    
    top           = center + workers_mpi.rows_y;
    bottom        = center - workers_mpi.rows_y;
    
    front         = center + workers_mpi.layer_size;
    back          = center - workers_mpi.layer_size;
    
    grad          =  &gradient[1][gidy];
    grad_left     =  &gradient[1][gidy-1];
    grad_back     =  &gradient[0][gidy];
    
    if (((gidy-workers_mpi.rows_y)/workers_mpi.rows_y) >= 0) {
      grad_bottom  = &gradient[1][gidy-workers_mpi.rows_y];
    }
    
    for (a=0; a < NUMPHASES; a++) {
      divphi[a]   =  (grad->gradphi[X][a] - grad_back->gradphi[X][a])/deltax;
      divphi[a]  +=  (grad->gradphi[Y][a] - grad_left->gradphi[Y][a])/deltay;
      divphi[a]  +=  (grad->gradphi[Z][a] - grad_bottom->gradphi[Z][a])/deltaz;
    }
    
    if (ELASTICITY) {
      grad->eigen_strain[X] = calculate_eigen_strain(center);
      
      grad->strain[X].xx = ((0.5)*(iter_gridinfo_w[front].disp[X][2] - iter_gridinfo_w[back].disp[X][2])   - grad->eigen_strain[X].xx);
      grad->strain[X].yy = ((0.5)*(iter_gridinfo_w[right].disp[Y][2] - iter_gridinfo_w[left].disp[Y][2])   - grad->eigen_strain[X].yy);
      grad->strain[X].zz = ((0.5)*(iter_gridinfo_w[top].disp[Z][2]   - iter_gridinfo_w[bottom].disp[Z][2]) - grad->eigen_strain[X].zz);
      grad->strain[X].xy = (0.25*((iter_gridinfo_w[right].disp[X][2] - iter_gridinfo_w[left].disp[X][2]) 
                         + (iter_gridinfo_w[front].disp[Y][2] - iter_gridinfo_w[back].disp[Y][2])));
      
      grad->strain[X].xz = (0.25*((iter_gridinfo_w[top].disp[X][2] - iter_gridinfo_w[bottom].disp[X][2]) 
                         + (iter_gridinfo_w[front].disp[Z][2] - iter_gridinfo_w[back].disp[Z][2])));
      
      grad->strain[X].yz = (0.25*((iter_gridinfo_w[right].disp[Z][2] - iter_gridinfo_w[left].disp[Z][2]) 
                         + (iter_gridinfo_w[top].disp[Y][2] - iter_gridinfo_w[bottom].disp[Y][2])));
      
      
      grad->stiffness_c[X] = calculate_stiffness(center);
  //     grad->stiffness_c[Y] = calculate_stiffness(center);
      
      
      sigma.xx  = grad->stiffness_c[X].C11*(grad->strain[X].xx) + grad->stiffness_c[X].C12*(grad->strain[X].yy + grad->strain[X].zz);
    
      sigma.yy  = grad->stiffness_c[X].C12*(grad->strain[X].xx  + grad->strain[X].zz) + grad->stiffness_c[X].C11*(grad->strain[X].yy);
     
      sigma.zz  = grad->stiffness_c[X].C12*(grad->strain[X].xx  + grad->strain[X].yy) + grad->stiffness_c[X].C11*(grad->strain[X].zz);
      
      sigma.xy  = 2.0*grad->stiffness_c[X].C44*(grad->strain[X].xy);
      
      sigma.xz  = 2.0*grad->stiffness_c[X].C44*(grad->strain[X].xz);
      
      sigma.yz  = 2.0*grad->stiffness_c[X].C44*(grad->strain[X].yz);
    }
    
    
  //Time step iteration, phase-field
  //Check if divergence of the phase-field exists.
    //For the free energy description
    if (!ISOTHERMAL) {
//       gridinfo_w[center].temperature  = temp_bottom + gidy*GRADIENT;
      T = gridinfo_w[center].temperature;
    }
    sum_lambdaphi=0.0;
    active_phases=0.0;
    
//     interface = 1;
//     
//     for (a=0; a < NUMPHASES; a++) {
//       if (gridinfo_w[center].phia[a] == 1.0) {
// 	bulk_phase=a;
// 	interface = 0;
// 	break;
//       }
//     }
//     noise = 4.0*AMP_NOISE_PHASE*((double)(rand())/(double)RAND_MAX);
    noise = 4.0*AMP_NOISE_PHASE*drand48();
    for (a=0; a < NUMPHASES; a++) {
      if (fabs(divphi[a]) > 0.0) {      
        lambda_phi[a] =  epsilon*(-dAdphi(gridinfo_w[center].phia, gradient, gidy, a) + 
        divdAdgradphi(gradient, center, gidy, a)/*- 0.0*d2gradphi(gridinfo_w[center].phia, gradient, center, gidy, a)*/) - (1.0/epsilon)*dwdphi(gridinfo_w[center].phia, divphi,gradient, gidy, a) - dpsi(gridinfo_w[center].compi, grad->phase_comp, T, gridinfo_w[center].phia,a)/V;
        
        if(NOISE_PHASEFIELD) {
//           lambda_phi[a] += AMP_NOISE_PHASE*((double)(rand())/(double)RAND_MAX)*gridinfo_w[center].phia[a]*(1.0-gridinfo_w[center].phia[a]); 
          lambda_phi[a] += noise*gridinfo_w[center].phia[a]*(1.0-gridinfo_w[center].phia[a])*dpsi(gridinfo_w[center].compi, grad->phase_comp, T, gridinfo_w[center].phia,a)/V;
        }
        if (ELASTICITY) {
          lambda_phi[a] -= df_elast(grad, sigma, gridinfo_w[center].phia, a);
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
        grad->deltaphi[a] =  deltat*(lambda_phi[a] - sum_lambdaphi)/(FunctionTau(gridinfo_w[center].phia)*epsilon);
      } else {
        grad->deltaphi[a] = 0.0;
      }
    }
    
    //Gibbs-simplex back projection: Check if the change in the phase-field 
    //larger than the phase-field value or the change makes the value 
    //larger than one
    projection_on_simplex(gridinfo_w[center].phia, grad->deltaphi, divphi);
//     projection_on_simplex_without_weights();
    for (a=0; a < NUMPHASES; a++) {
      gridinfo_w[center].deltaphi[a] = grad->deltaphi[a];
    }
  }
}
#endif
