#ifndef CALCULATE_DIVERGENCE_CONCENTRATION_H_
#define CALCULATE_DIVERGENCE_CONCENTRATION_H_

void calculate_divergence_concentration_2D(long x, struct gradlayer **gradient) {
  long k, l, a;
  for (gidy=1; gidy <=(workers_mpi.end[Y]+2); gidy++) {
    
    center        =  gidy   + (x)*workers_mpi.layer_size;
    front         =  gidy   + (x+1)*workers_mpi.layer_size;
    right         =  center + 1;
    left          =  center - 1;

    grad          =  &gradient[0][gidy];
    grad_left     =  &gradient[0][gidy-1];
    grad_back     =  &gradient[-1][gidy];
    
    if (!ISOTHERMAL) {
      T = gridinfo_w[center].temperature;
      init_propertymatrices(T);
    }
    
    double sum_dcbdT[NUMCOMPONENTS-1];
    double mu;
    
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      sum[k]        = 0.0;
      sum_dcbdT[k]  = 0.0;
      for (l=0; l < NUMCOMPONENTS-1; l++) {
	dcdmu[k][l] = 0.0;
      }
    }
    for (a=0; a < NUMPHASES; a++) {
      sum_dhphi = 0.0;
      for (b=0; b < NUMPHASES; b++) {
// 	sum_dhphi += dhphi(gridinfo_w[center].phia, a, b)*grad->deltaphi[b];
        sum_dhphi += dhphi(gridinfo_w[center].phia, a, b)*gridinfo_w[center].deltaphi[b];
      }
      for (k=0; k < NUMCOMPONENTS-1; k++) {
	sum[k]       += c_mu(gridinfo_w[center].compi, T, a, k)*sum_dhphi;
        sum_dcbdT[k] += dcbdT_phase[a][k]*hphi(gridinfo_w[center].phia, a);
//         sum_dcbdT[k] += func_dcbdT_phase(gridinfo_w[center].compi, T, a, k)*hphi(gridinfo_w[center].phia, a);
      }
    }
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      divflux[k] = 0.0;
      for (l=0; l < NUMCOMPONENTS-1; l++) {
	divflux[k] += (grad->Dmid[X][k][l]*grad->gradchempot[X][l] - grad_back->Dmid[X][k][l]*grad_back->gradchempot[X][l])/deltax;

		
	divflux[k] += (grad->Dmid[Y][k][l]*grad->gradchempot[Y][l] - grad_left->Dmid[Y][k][l]*grad_left->gradchempot[Y][l])/deltay;
      }
      if (OBSTACLE) {
        divjat[k]  = (grad->jat[X][k]-grad_back->jat[X][k])*(-0.25*M_PI*epsilon)/deltax;
        divjat[k] += (grad->jat[Y][k]-grad_left->jat[Y][k])*(-0.25*M_PI*epsilon)/deltay;
      }
      if (WELL) {
        divjat[k]  = (grad->jat[X][k]-grad_back->jat[X][k])*(-1.0*epsilon)/deltax;
        divjat[k] += (grad->jat[Y][k]-grad_left->jat[Y][k])*(-1.0*epsilon)/deltay;
      }
    }
    //Time step iteration, concentration
    //Check whether you are inside the interface
    interface = 1;
    for (a=0; a < NUMPHASES; a++) {
      if (gridinfo_w[center].phia[a] == 1.0) {
	bulk_phase=a;
	interface = 0;
	break;
      }
    }
    if (interface==0) {
      for (k=0; k < NUMCOMPONENTS-1; k++) {
	c_old[k] = c_mu(gridinfo_w[center].compi, T, bulk_phase,k);
	c_new[k] = c_old[k] + deltat*(divflux[k])-divjat[k];
      }
      for (k=0; k < NUMCOMPONENTS-1; k++) {
        deltamu[k] = Mu(c_new, T, bulk_phase, k) - gridinfo_w[center].compi[k];
	gridinfo_w[center].compi[k] = Mu(c_new, T, bulk_phase, k);
      }
    } else {
      for (k=0; k < NUMCOMPONENTS-1; k++) {
        if (ISOTHERMAL) {
	  deltac[k] = deltat*(divflux[k])-divjat[k]-sum[k];
        } else {
          deltac[k] = deltat*(divflux[k])-divjat[k]-sum[k] -sum_dcbdT[k]*deltat*(-temperature_gradientY.GRADIENT*temperature_gradientY.velocity);
        }
	for (l=0; l < NUMCOMPONENTS-1; l++) {
	  for (a=0; a < NUMPHASES; a++) {
	    dcdmu[k][l] += dc_dmu(gridinfo_w[center].compi, T, a, k, l)*hphi(gridinfo_w[center].phia, a);
	  }
	}
      }
      if (BINARY) {
        inv_dcdmu[0][0]             = 1.0/dcdmu[0][0];
        deltamu[0]                  = deltac[0]*inv_dcdmu[0][0];
        gridinfo_w[center].compi[0] += deltamu[0];
//         if ((fabs(gridinfo_w[center].compi[0]) > 1.0)) {
//           printf("taskid=%d, gridinfo_w[center].compi=%le, gidy=%ld\n",taskid, gridinfo_w[center].compi[0], gidy);
//         }
      }
      
      if (TERNARY) {
        DET = dcdmu[0][0]*dcdmu[1][1] - dcdmu[0][1]*dcdmu[1][0];
        inv_dcdmu[0][0]  = dcdmu[1][1]/DET;
        inv_dcdmu[1][1]  = dcdmu[0][0]/DET;
        inv_dcdmu[0][1]  = -dcdmu[0][1]/DET;
        inv_dcdmu[1][0]  = -dcdmu[1][0]/DET;
        
        for (k=0; k < NUMCOMPONENTS-1; k++ ) {
          deltamu[k] = 0.0;
          for (l=0; l < NUMCOMPONENTS-1; l++) {
            deltamu[k] += inv_dcdmu[k][l]*deltac[l];
          }
        }
        vectorsum(deltamu, gridinfo_w[center].compi, gridinfo_w[center].compi,NUMCOMPONENTS-1);
      }
      if (DILUTE) {
        for (k=0; k < NUMCOMPONENTS-1; k++) {
          inv_dcdmu[k][k]             =  1.0/dcdmu[k][k];
          deltamu[k]                  =  deltac[k]*inv_dcdmu[k][k];
          gridinfo_w[center].compi[k] += deltamu[k];
        }
      }
      if (!(DILUTE || BINARY || TERNARY)) {
       matinvnew(dcdmu,inv_dcdmu,NUMCOMPONENTS-1);
       multiply(inv_dcdmu,deltac,deltamu,NUMCOMPONENTS-1);
       vectorsum(deltamu,gridinfo_w[center].compi,gridinfo_w[center].compi,NUMCOMPONENTS-1);
      }
    }
    //Update of phase-field
    for (b=0; b < NUMPHASES; b++) {
//       gridinfo_w[center].phia[b] = gridinfo_w[center].phia[b] + grad->deltaphi[b];
       gridinfo_w[center].phia[b] = gridinfo_w[center].phia[b] + gridinfo_w[center].deltaphi[b];
       //FINDING MAXIMUM AND MINIMUM VALUES 
       if (gridinfo_w[center].phia[b] > workers_max_min.phi_max[b]) {
         workers_max_min.phi_max[b] = gridinfo_w[center].phia[b];
       }
       if (gridinfo_w[center].phia[b] < workers_max_min.phi_min[b]) {
         workers_max_min.phi_min[b] = gridinfo_w[center].phia[b];
       }
       workers_max_min.rel_change_phi[b] += gridinfo_w[center].deltaphi[b]*gridinfo_w[center].deltaphi[b];
       //FINISHED
    }
    for (k=0; k<NUMCOMPONENTS-1; k++) {
        //FINDING MAXIMUM AND MINIMUM VALUES 
       if (gridinfo_w[center].compi[k] > workers_max_min.mu_max[k]) {
         workers_max_min.mu_max[k] = gridinfo_w[center].compi[k];
       }
       if (gridinfo_w[center].compi[k] < workers_max_min.mu_min[k]) {
         workers_max_min.mu_min[k] = gridinfo_w[center].compi[k];
       }
       workers_max_min.rel_change_mu[k] += deltamu[k]*deltamu[k];
       //FINISHED
    }
  }
}
void calculate_divergence_concentration_smooth_2D(long x, struct gradlayer **gradient) {
  long k, l, a;
  for (gidy=1; gidy <=(workers_mpi.end[Y]+2); gidy++) {
    
    center        =  gidy   + (x)*workers_mpi.layer_size;
    front         =  gidy   + (x+1)*workers_mpi.rows_y;
    right         =  center + 1;
    left          =  center - 1;

    grad          =  &gradient[0][gidy];
    grad_left     =  &gradient[0][gidy-1];
    grad_back     =  &gradient[-1][gidy];
    
    if (!ISOTHERMAL) {
      T = gridinfo_w[center].temperature;
      init_propertymatrices(T);
    }
    
    for (k=0; k < NUMCOMPONENTS-1; k++) {
//       sum[k] = 0.0;
      for (l=0; l < NUMCOMPONENTS-1; l++) {
	dcdmu[k][l] = 0.0;
      }
    }
//     for (a=0; a < NUMPHASES; a++) {
//       sum_dhphi = 0.0;
//       for (b=0; b < NUMPHASES; b++) {
// 	sum_dhphi += dhphi(gridinfo_w[center].phia, a, b)*grad->deltaphi[b];
//       }
//       for (k=0; k < NUMCOMPONENTS-1; k++) {
// 	sum[k] += c_mu(gridinfo_w[center].compi, T, a, k)*sum_dhphi;
//       }
//     }
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      divflux[k] = 0.0;
      for (l=0; l < NUMCOMPONENTS-1; l++) {
	divflux[k] += (grad->Dmid[X][k][l]*grad->gradchempot[X][l] - grad_back->Dmid[X][k][l]*grad_back->gradchempot[X][l])/deltax;

		
	divflux[k] += (grad->Dmid[Y][k][l]*grad->gradchempot[Y][l] - grad_left->Dmid[Y][k][l]*grad_left->gradchempot[Y][l])/deltay;
	
      }
//   #ifdef OBSTACLE
//       divjat[k]  = (grad->jat[X][k]-grad_back->jat[X][k])*(-0.25*M_PI*epsilon)/deltax;
//       divjat[k] += (grad->jat[Y][k]-grad_left->jat[Y][k])*(-0.25*M_PI*epsilon)/deltay;
//   #endif				      
//   #ifdef WELL
// 
//       divjat[k]  = (grad->jat[X][k]-grad_back->jat[X][k])*(-1.0*epsilon)/deltax;
//       divjat[k] += (grad->jat[Y][k]-grad_left->jat[Y][k])*(-1.0*epsilon)/deltay;
//   #endif
    }
    //Time step iteration, concentration
    //Check whether you are inside the interface
    interface = 1;
    for (a=0; a < NUMPHASES; a++) {
      if (gridinfo_w[center].phia[a] == 1.0) {
	bulk_phase=a;
	interface = 0;
	break;
      }
    }
    if (interface==0) {
      for (k=0; k < NUMCOMPONENTS-1; k++) {
	c_old[k] = c_mu(gridinfo_w[center].compi, T, bulk_phase,k);
	c_new[k] = c_old[k] + deltat*(divflux[k]);
      }
      for (k=0; k < NUMCOMPONENTS-1; k++) {
	gridinfo_w[center].compi[k] = Mu(c_new, T, bulk_phase, k);
      }
    } else {
      for (k=0; k < NUMCOMPONENTS-1; k++) {
	deltac[k] = deltat*(divflux[k]);
	for (l=0; l < NUMCOMPONENTS-1; l++) {
	  for (a=0; a < NUMPHASES; a++) {
	    dcdmu[k][l] += dc_dmu(gridinfo_w[center].compi, T, a, k, l)*hphi(gridinfo_w[center].phia, a);
	  }
	}
      }
      matinvnew(dcdmu,inv_dcdmu,NUMCOMPONENTS-1);
      multiply(inv_dcdmu,deltac,deltamu,NUMCOMPONENTS-1);
      vectorsum(deltamu,gridinfo_w[center].compi,gridinfo_w[center].compi,NUMCOMPONENTS-1);
    }
    //Update of phase-field
    for (b=0; b < NUMPHASES; b++) {
      gridinfo_w[center].phia[b] = gridinfo_w[center].phia[b] + grad->deltaphi[b];
    }
    //The shift condition for growth of solid in
//   #ifdef SHIFT
//     if ((gridinfo_w[center-1].phia[NUMPHASES-1]-(1.0-gridinfo_w[center-1].phia[NUMPHASES-1]) < 0.0) 
//       && (gridinfo_w[center].phia[NUMPHASES-1]-(1.0-gridinfo_w[center].phia[NUMPHASES-1]) > 0.0) 
//       && (gidy > shiftj)) {
//       SHIFT_Y = TRUE;
//       if (gidy > MAX_INTERFACE_POS) {
// 	MAX_INTERFACE_POS = gidy;
//       }
//     }
//   #endif
  }
}
void calculate_divergence_concentration_smooth_concentration_2D(long x, struct gradlayer **gradient) {
  long k, l, a;
  for (gidy=1; gidy <=(workers_mpi.end[Y]+2); gidy++) {
    
    center        =  gidy   + (x)*workers_mpi.layer_size;
    front         =  gidy   + (x+1)*workers_mpi.layer_size;
    right         =  center + 1;
    left          =  center - 1;

    grad          =  &gradient[0][gidy];
    grad_left     =  &gradient[0][gidy-1];
    grad_back     =  &gradient[-1][gidy];
    
    if (!ISOTHERMAL) {
      T = gridinfo_w[center].temperature;
      init_propertymatrices(T);
    }
    
    for (k=0; k < NUMCOMPONENTS-1; k++) {
//       sum[k] = 0.0;
      for (l=0; l < NUMCOMPONENTS-1; l++) {
	dcdmu[k][l] = 0.0;
      }
    }
//     for (a=0; a < NUMPHASES; a++) {
//       sum_dhphi = 0.0;
//       for (b=0; b < NUMPHASES; b++) {
// 	sum_dhphi += dhphi(gridinfo_w[center].phia, a, b)*grad->deltaphi[b];
//       }
//       for (k=0; k < NUMCOMPONENTS-1; k++) {
// 	sum[k] += c_mu(gridinfo_w[center].compi, T, a, k)*sum_dhphi;
//       }
//     }
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      divflux[k] = 0.0;
      for (l=0; l < NUMCOMPONENTS-1; l++) {
	divflux[k] += (grad->Dmid[X][k][l]*grad->gradchempot[X][l] - grad_back->Dmid[X][k][l]*grad_back->gradchempot[X][l])/deltax;

		
	divflux[k] += (grad->Dmid[Y][k][l]*grad->gradchempot[Y][l] - grad_left->Dmid[Y][k][l]*grad_left->gradchempot[Y][l])/deltay;
	
      }
//   #ifdef OBSTACLE
//       divjat[k]  = (grad->jat[X][k]-grad_back->jat[X][k])*(-0.25*M_PI*epsilon)/deltax;
//       divjat[k] += (grad->jat[Y][k]-grad_left->jat[Y][k])*(-0.25*M_PI*epsilon)/deltay;
//   #endif				      
//   #ifdef WELL
// 
//       divjat[k]  = (grad->jat[X][k]-grad_back->jat[X][k])*(-1.0*epsilon)/deltax;
//       divjat[k] += (grad->jat[Y][k]-grad_left->jat[Y][k])*(-1.0*epsilon)/deltay;
//   #endif
    }
    //Time step iteration, concentration
    //Check whether you are inside the interface
    interface = 1;
    for (a=0; a < NUMPHASES; a++) {
      if (gridinfo_w[center].phia[a] == 1.0) {
	bulk_phase=a;
	interface = 0;
	break;
      }
    }
    if (interface==0) {
      for (k=0; k < NUMCOMPONENTS-1; k++) {
	c_old[k] = c_mu(gridinfo_w[center].compi, T, bulk_phase,k);
	c_new[k] = c_old[k] + deltat*(divflux[k]);
      }
      for (k=0; k < NUMCOMPONENTS-1; k++) {
	gridinfo_w[center].compi[k] = Mu(c_new, T, bulk_phase, k);
      }
    } else {
      for (k=0; k < NUMCOMPONENTS-1; k++) {
	deltac[k] = deltat*(divflux[k]);
	for (l=0; l < NUMCOMPONENTS-1; l++) {
	  for (a=0; a < NUMPHASES; a++) {
	    dcdmu[k][l] += dc_dmu(gridinfo_w[center].compi, T, a, k, l)*hphi(gridinfo_w[center].phia, a);
	  }
	}
      }
      matinvnew(dcdmu,inv_dcdmu,NUMCOMPONENTS-1);
      multiply(inv_dcdmu,deltac,deltamu,NUMCOMPONENTS-1);
      vectorsum(deltamu,gridinfo_w[center].compi,gridinfo_w[center].compi,NUMCOMPONENTS-1);
    }
    //Update of phase-field
//     for (b=0; b < NUMPHASES; b++) {
//       gridinfo_w[center].phia[b] = gridinfo_w[center].phia[b] + grad->deltaphi[b];
//     }
    //The shift condition for growth of solid in
//   #ifdef SHIFT
//     if ((gridinfo_w[center-1].phia[NUMPHASES-1]-(1.0-gridinfo_w[center-1].phia[NUMPHASES-1]) < 0.0) 
//       && (gridinfo_w[center].phia[NUMPHASES-1]-(1.0-gridinfo_w[center].phia[NUMPHASES-1]) > 0.0) 
//       && (gidy > shiftj)) {
//       SHIFT_Y = TRUE;
//       if (gidy > MAX_INTERFACE_POS) {
// 	MAX_INTERFACE_POS = gidy;
//       }
//     }
//   #endif
  }
}
#endif
// double func_dcbdT_phase(double *mu, double T, long a, long k) {
//   return dc_dmu(mu, T, a, k, k)*(-func_dBbdT(T, a, k));
// }
// double func_dBbdT(double T, long a, long k) {
//   double dercomp_l, dercomp_s;
//   double comp_l, comp_s;
//   dercomp_l = calculate_der_equilibrium_composition(T, k, NUMPHASES-1);
//   dercomp_s = calculate_der_equilibrium_composition(T, k, a);
//   
//   comp_l    = calculate_equilibrium_composition(T, k, NUMPHASES-1);
//   comp_s    = calculate_equilibrium_composition(T, k, a);
//   
//   return T*(dercomp_l/comp_l - dercomp_s/comp_s);
// }