#ifndef CALCULATE_GRADIENTS_H_
#define CALCULATE_GRADIENTS_H_

void calculate_diffusion_potential(long x, struct gradlayer **gradient);
struct Stiffness_cubic calculate_stiffness(long index);
struct symmetric_tensor calculate_eigen_strain(long index);
struct Stiffness_cubic calculate_stiffness_n(long index);
// void calculate_gradients_phasefield_tdb_2D(long x, struct gradlayer **gradient, int CALCULATE_COMPOSITION);
// void calculate_gradients_concentration_tdb_2D(long x, struct gradlayer **gradient);
// void calculate_gradients_tdb_2D(long x, struct gradlayer **gradient);

// void calculate_gradients_2D(long x, struct gradlayer **gradient) {
// //Calculation relevant to  phase-field
//  long k, l, a;
//   for (gidy=0; gidy <= workers_mpi.end[Y]+2; gidy++) {
// 	  
//     grad          =  &gradient[2][gidy];
//     grad_back     =  &gradient[1][gidy];
//     
//     center        =  gidy   + (x)*workers_mpi.layer_size;
//     front         =  gidy   + (x+1)*workers_mpi.layer_size;
//     right         =  center + 1;
//     left          =  center - 1;
//     if (DIMENSION != 2) {
//       top         =  center + workers_mpi.rows_y;
//       bottom      =  center - workers_mpi.rows_y;
//     }
//     for (a=0; a < NUMPHASES; a++) {
//       if (x <=(workers_mpi.end[X]+2)) {
// 	grad->gradphi[X][a]  = (gridinfo_w[front].phia[a]     - gridinfo_w[center].phia[a])/deltax;
// 	grad->phistagg[X][a] = 0.5*(gridinfo_w[front].phia[a] + gridinfo_w[center].phia[a]);
//       } else {
// 	grad->gradphi[X][a]  = 0.0;
// 	grad->phistagg[X][a] = (gridinfo_w[center].phia[a]);
//       }
//       grad->gradphi[Y][a]    = (gridinfo_w[right].phia[a]     - gridinfo_w[center].phia[a])/deltay;
//       grad->phistagg[Y][a]   = 0.5*(gridinfo_w[right].phia[a] + gridinfo_w[center].phia[a]);
//     }
// //    
//     if(ISOTHERMAL) {
//       for (k=0; k < NUMCOMPONENTS-1; k++) {
//         if (x <=(workers_mpi.end[X]+2)) {
//           grad->gradchempot[X][k] = (gridinfo_w[front].compi[k] - gridinfo_w[center].compi[k])/deltax;
//         } else {
//           grad->gradchempot[X][k] = 0.0;
//         }
//         grad->gradchempot[Y][k] = (gridinfo_w[right].compi[k] - gridinfo_w[center].compi[k])/deltay;
//         for (l=0; l < NUMCOMPONENTS-1; l++) {
//         if (x <=(workers_mpi.end[X]+2)) {
//           grad->Dmid[X][k][l] = (D(gridinfo_w, T, front, k, l) + D(gridinfo_w, T, center, k, l))*0.5;
//         } else {
//             grad->Dmid[X][k][l] = D(gridinfo_w, T, center, k, l);
//         }
//           grad->Dmid[Y][k][l] = (D(gridinfo_w, T, right, k, l) + D(gridinfo_w, T, center, k, l))*0.5;
//         }
//       }
//     } else {
//       for (k=0; k < NUMCOMPONENTS-1; k++) {
//         if(x<=(workers_mpi.end[X]+2)) {
//           grad->gradchempot[X][k] = (gridinfo_w[front].compi[k] - gridinfo_w[center].compi[k])/deltax;
//         } else {
//           grad->gradchempot[X][k] = 0.0;
//         }
//         grad->gradchempot[Y][k] = (gridinfo_w[right].compi[k] - gridinfo_w[center].compi[k])/deltay;
//         for (l=0; l < NUMCOMPONENTS-1; l++) {  
//           if(x<=(workers_mpi.end[X]+2)) {
//             grad->Dmid[X][k][l] = (D(gridinfo_w, gridinfo_w[front].temperature, front, k, l) + D(gridinfo_w, gridinfo_w[center].temperature, center, k, l))*0.5;
//           } else {
//             grad->Dmid[X][k][l] = D(gridinfo_w, gridinfo_w[center].temperature, center, k, l);
//           }
//           grad->Dmid[Y][k][l] = (D(gridinfo_w, gridinfo_w[right].temperature, right, k, l) + D(gridinfo_w, gridinfo_w[center].temperature, center, k, l))*0.5;
//         }
//       }
//    }
//     for (a=0; a < NUMPHASES; a++) {
//       if (gidy > 0 && gidy <=(workers_mpi.end[Y]+2)) {
//         grad->gradphi_c[Y][a] = (gridinfo_w[right].phia[a] - gridinfo_w[left].phia[a])/(2.0*deltay);
//       } else {
//         grad->gradphi_c[Y][a] = grad->gradphi[Y][a];
//       }
//       if (x > 0) {
//         grad->gradphi_c[X][a]  = 0.5*(grad->gradphi[X][a] + grad_back->gradphi[X][a]);
// //         grad->gradphi_c[X][a]  = grad->gradphi[X][a];
//       }
//     }
//   }
// //   if(taskid==1) {
// //     printf("taskid=%d,Coming here; x=%ld, y=%ld, end[X]=%ld, end[Y]=%ld\n", taskid,x,gidy,workers_mpi.end[X],workers_mpi.end[Y]);
// //   }
//   
// //   for (gidy=0; gidy <= (end[Y]+2); gidy++) {
// //     grad          = &gradient[2][gidy];
// //     grad_back     = &gradient[1][gidy];
// // //       grad_right    = &gradient[2][gidy+1];
// //     grad_left     = &gradient[2][gidy-1];
// //   
// //     for(a=0; a<NUMPHASES; a++) {
// //       if((x > 0) && (x<=(end[X]+2)) && (gidy > 0) && (gidy <= (end[Y]+2))) {
// // 	grad->d2gradphi[a]  = (grad->gradphi[X][a] -  grad_back->gradphi[X][a])/deltax;
// // 	grad->d2gradphi[a] += (grad->gradphi[Y][a] -  grad_left->gradphi[Y][a])/deltay;
// //       }
// //     }
// //   }
// }
// void calculate_gradients_tdb_2D(long x, struct gradlayer **gradient) {
// //Calculation relevant to  phase-field
//  long k, l, a;
//   for (gidy=0; gidy <= workers_mpi.end[Y]+2; gidy++) {
// 	  
//     grad          =  &gradient[2][gidy];
//     grad_back     =  &gradient[1][gidy];
//     
//     center        =  gidy   + (x)*workers_mpi.layer_size;
//     front         =  gidy   + (x+1)*workers_mpi.layer_size;
//     right         =  center + 1;
//     left          =  center - 1;
//     if (DIMENSION != 2) {
//       top         =  center + workers_mpi.rows_y;
//       bottom      =  center - workers_mpi.rows_y;
//     }
//     for (a=0; a < NUMPHASES; a++) {
//       if (x <=(workers_mpi.end[X]+2)) {
// 	grad->gradphi[X][a]  = (gridinfo_w[front].phia[a]     - gridinfo_w[center].phia[a])/deltax;
// 	grad->phistagg[X][a] = 0.5*(gridinfo_w[front].phia[a] + gridinfo_w[center].phia[a]);
//       } else {
// 	grad->gradphi[X][a]  = 0.0;
// 	grad->phistagg[X][a] = (gridinfo_w[center].phia[a]);
//       }
//       grad->gradphi[Y][a]    = (gridinfo_w[right].phia[a]     - gridinfo_w[center].phia[a])/deltay;
//       grad->phistagg[Y][a]   = 0.5*(gridinfo_w[right].phia[a] + gridinfo_w[center].phia[a]);
//     }
//     
// //     interface = 1;
// //     
// //     for (a=0; a < NUMPHASES; a++) {
// //       if (gridinfo_w[center].phia[a] == 1.0) {
// // 	bulk_phase=a;
// // 	interface = 0;
// // 	break;
// //       }
// //     }
// //     
// //     if (interface) {
// //       if (!ISOTHERMAL) {
// //         T = gridinfo_w[center].temperature;
// //       }
// //       for (a=0; a < NUMPHASES; a++) {
// //         function_F_02_getPhaseComposition(gridinfo_w[center].compi, grad->phase_comp[a], T, a, ceq[a][a]);
// //       }
// //       for (a=0; a < NUMPHASES; a++) {
// //         tdb_dc_dmu(gridinfo_w[center].compi, grad->phase_comp[a], T, a, grad->dcdmu_phase[a]);
// //       }
// //     } else {
// //       for (k=0; k < NUMCOMPONENTS-1; k++) {
// //         grad->phase_comp[bulk_phase][k] = gridinfo_w[center].composition[k];
// //       }
// //       tdb_dc_dmu(gridinfo_w[center].compi, grad->phase_comp[bulk_phase], T, bulk_phase, grad->dcdmu_phase[bulk_phase]);
// //     }
// // //    
// //     if(ISOTHERMAL) {
// //       for (k=0; k < NUMCOMPONENTS-1; k++) {
// //         if (x <=(workers_mpi.end[X]+2)) {
// //           grad->gradchempot[X][k] = (gridinfo_w[front].compi[k] - gridinfo_w[center].compi[k])/deltax;
// //         } else {
// //           grad->gradchempot[X][k] = 0.0;
// //         }
// //         grad->gradchempot[Y][k] = (gridinfo_w[right].compi[k] - gridinfo_w[center].compi[k])/deltay;
// //         for (l=0; l < NUMCOMPONENTS-1; l++) {
// //         if (x <=(workers_mpi.end[X]+2)) {
// //           grad->Dmid[X][k][l] = (D(gridinfo_w, T, front, k, l) + D(gridinfo_w, T, center, k, l))*0.5;
// //         } else {
// //             grad->Dmid[X][k][l] = D(gridinfo_w, T, center, k, l);
// //         }
// //           grad->Dmid[Y][k][l] = (D(gridinfo_w, T, right, k, l) + D(gridinfo_w, T, center, k, l))*0.5;
// //         }
// //       }
// //     } else {
// //       for (k=0; k < NUMCOMPONENTS-1; k++) {
// //         if(x<=(workers_mpi.end[X]+2)) {
// //           grad->gradchempot[X][k] = (gridinfo_w[front].compi[k] - gridinfo_w[center].compi[k])/deltax;
// //         } else {
// //           grad->gradchempot[X][k] = 0.0;
// //         }
// //         grad->gradchempot[Y][k] = (gridinfo_w[right].compi[k] - gridinfo_w[center].compi[k])/deltay;
// //         for (l=0; l < NUMCOMPONENTS-1; l++) {  
// //           if(x<=(workers_mpi.end[X]+2)) {
// //             grad->Dmid[X][k][l] = (D(gridinfo_w, gridinfo_w[front].temperature, front, k, l) + D(gridinfo_w, gridinfo_w[center].temperature, center, k, l))*0.5;
// //           } else {
// //             grad->Dmid[X][k][l] = D(gridinfo_w, gridinfo_w[center].temperature, center, k, l);
// //           }
// //           grad->Dmid[Y][k][l] = (D(gridinfo_w, gridinfo_w[right].temperature, right, k, l) + D(gridinfo_w, gridinfo_w[center].temperature, center, k, l))*0.5;
// //         }
// //       }
// //    }
//     for (a=0; a < NUMPHASES; a++) {
//       if (gidy > 0 && gidy <=(workers_mpi.end[Y]+2)) {
//         grad->gradphi_c[Y][a] = (gridinfo_w[right].phia[a] - gridinfo_w[left].phia[a])/(2.0*deltay);
//       } else {
//         grad->gradphi_c[Y][a] = grad->gradphi[Y][a];
//       }
//       if (x > 0) {
//         grad->gradphi_c[X][a]  = 0.5*(grad->gradphi[X][a] + grad_back->gradphi[X][a]);
// //         grad->gradphi_c[X][a]  = grad->gradphi[X][a];
//       }
//     }
//   }
// //   if(taskid==1) {
// //     printf("taskid=%d,Coming here; x=%ld, y=%ld, end[X]=%ld, end[Y]=%ld\n", taskid,x,gidy,workers_mpi.end[X],workers_mpi.end[Y]);
// //   }
//   
// //   for (gidy=0; gidy <= (end[Y]+2); gidy++) {
// //     grad          = &gradient[2][gidy];
// //     grad_back     = &gradient[1][gidy];
// // //       grad_right    = &gradient[2][gidy+1];
// //     grad_left     = &gradient[2][gidy-1];
// //   
// //     for(a=0; a<NUMPHASES; a++) {
// //       if((x > 0) && (x<=(end[X]+2)) && (gidy > 0) && (gidy <= (end[Y]+2))) {
// // 	grad->d2gradphi[a]  = (grad->gradphi[X][a] -  grad_back->gradphi[X][a])/deltax;
// // 	grad->d2gradphi[a] += (grad->gradphi[Y][a] -  grad_left->gradphi[Y][a])/deltay;
// //       }
// //     }
// //   }
// }
// void calculate_gradients_phasefield_2D(long x, struct gradlayer **gradient) {
//   long k, l, a, gidy;
//   for (gidy=0; gidy <= (workers_mpi.end[Y]+2); gidy++) {
//           
//     grad          =  &gradient[2][gidy];
//     grad_back     =  &gradient[1][gidy];
//     
//     center        =  gidy   + (x)*workers_mpi.layer_size;
//     front         =  gidy   + (x+1)*workers_mpi.layer_size;
//     right         =  center + 1;
//     left          =  center - 1;
//     if (DIMENSION != 2) {
//       top         =  center + workers_mpi.rows_y;
//       bottom      =  center - workers_mpi.rows_y;
//     }
//     
//     for (a=0; a < NUMPHASES; a++) {
//       if (x <=(workers_mpi.end[X]+2)) {
//         grad->gradphi[X][a]  = (gridinfo_w[front].phia[a]     - gridinfo_w[center].phia[a])/deltax;
//         grad->phistagg[X][a] = 0.5*(gridinfo_w[front].phia[a] + gridinfo_w[center].phia[a]);
//       } else {
//         grad->gradphi[X][a]  = 0.0;
//         grad->phistagg[X][a] = (gridinfo_w[center].phia[a]);
//       }
//       grad->gradphi[Y][a]  = (gridinfo_w[right].phia[a] - gridinfo_w[center].phia[a])/deltay;
//       grad->phistagg[Y][a] = 0.5*(gridinfo_w[right].phia[a] + gridinfo_w[center].phia[a]);
//     }
//     
//     for (a=0; a < NUMPHASES; a++) {
//       if (gidy > 0 && gidy <=(workers_mpi.end[Y]+2)) {
//         grad->gradphi_c[Y][a] = (gridinfo_w[right].phia[a] - gridinfo_w[left].phia[a])/(2.0*deltay);
//       } else {
//         grad->gradphi_c[Y][a] = grad->gradphi[Y][a];
//       }
//       if (x > 0) {
//         grad->gradphi_c[X][a] = 0.5*(grad->gradphi[X][a] + grad_back->gradphi[X][a]);
//       }
//     }
//   }
// }
// void calculate_gradients_concentration_2D(long x, struct gradlayer **gradient) {
//   long k, l, a, gidy;
//   for (gidy=0; gidy <= (workers_mpi.end[Y]+2); gidy++) {
//           
//     grad          =  &gradient[2][gidy];
//     grad_back     =  &gradient[1][gidy];
//     
//     center        =  gidy   + (x)*workers_mpi.layer_size;
//     front         =  gidy   + (x+1)*workers_mpi.layer_size;
//     right         =  center + 1;
//     left          =  center - 1;
//     if (DIMENSION != 2) {
//       top         =  center + workers_mpi.rows_y;
//       bottom      =  center - workers_mpi.rows_y;
//     }
//     
//     if(ISOTHERMAL) {
//       for (k=0; k < NUMCOMPONENTS-1; k++) {
//         if (x <=(workers_mpi.end[X]+2)) {
//           grad->gradchempot[X][k]  = (gridinfo_w[front].compi[k] - gridinfo_w[center].compi[k])/deltax;
//         } else {
//           grad->gradchempot[X][k] = 0.0;
//         }
//         grad->gradchempot[Y][k]   = (gridinfo_w[right].compi[k]  - gridinfo_w[center].compi[k])/deltay;
//         
//         for (l=0; l < NUMCOMPONENTS-1; l++) {
//         if (x <=(workers_mpi.end[X]+2)) {
//           grad->Dmid[X][k][l] = (D(gridinfo_w, T, front, k, l) + D(gridinfo_w, T, center, k, l))*0.5;
//           } else {
//             grad->Dmid[X][k][l] = D(gridinfo_w, T, center, k, l);
//           }
//           grad->Dmid[Y][k][l] = (D(gridinfo_w, T, right, k, l) + D(gridinfo_w, T, center, k, l))*0.5;
//         }
//       }
//     } else {
//       for (k=0; k < NUMCOMPONENTS-1; k++) {
//         if(x<=(workers_mpi.end[X]+2)) {
//           grad->gradchempot[X][k] = (gridinfo_w[front].compi[k] - gridinfo_w[center].compi[k])/deltax;
//         } else {
//           grad->gradchempot[X][k] = 0.0;
//         }
//         grad->gradchempot[Y][k]   = (gridinfo_w[right].compi[k] - gridinfo_w[center].compi[k])/deltay;
//         for (l=0; l < NUMCOMPONENTS-1; l++) {  
//           if(x<=(workers_mpi.end[X]+2)) {
//             grad->Dmid[X][k][l]   = (D(gridinfo_w, gridinfo_w[front].temperature, front, k, l) + D(gridinfo_w, gridinfo_w[center].temperature, center, k, l))*0.5;
//           } else {
//             grad->Dmid[X][k][l]   = D(gridinfo_w, gridinfo_w[center].temperature, center, k, l);
//           }
//           grad->Dmid[Y][k][l]     = (D(gridinfo_w, gridinfo_w[right].temperature, right, k, l) + D(gridinfo_w, gridinfo_w[center].temperature, center, k, l))*0.5;
//         }
//       }
//     } 
//   }
// }
void calculate_gradients_phasefield_2D(long x, struct gradlayer **gradient, int CALCULATE_COMPOSITION) {
  long k, l, a, gidy;
  for (gidy=0; gidy <= (workers_mpi.end[Y]+2); gidy++) {
          
    grad          =  &gradient[2][gidy];
    grad_back     =  &gradient[1][gidy];
    
    center        =  gidy   + (x)*workers_mpi.layer_size;
    front         =  gidy   + (x+1)*workers_mpi.layer_size;
    right         =  center + 1;
    left          =  center - 1;
    if (DIMENSION != 2) {
      top         =  center + workers_mpi.rows_y;
      bottom      =  center - workers_mpi.rows_y;
    }
    
    grad->interface = 1;
    grad->bulk_phase = 0;
    
    for (a=0; a < NUMPHASES; a++) {
      if (gridinfo_w[center].phia[a] == 1.0) {
        grad->bulk_phase=a;
        grad->interface = 0;
        break;
      }
    }
    if (!ISOTHERMAL) {
      T = gridinfo_w[center].temperature;
      DELTAT = deltat*(-temperature_gradientY.GRADIENT*temperature_gradientY.velocity);
    }
    
    for (a=0; a < NUMPHASES; a++) {
      if (x <=(workers_mpi.end[X]+2)) {
        grad->gradphi[X][a]  = (gridinfo_w[front].phia[a]     - gridinfo_w[center].phia[a])/deltax;
        grad->phistagg[X][a] = 0.5*(gridinfo_w[front].phia[a] + gridinfo_w[center].phia[a]);
      } else {
        grad->gradphi[X][a]  = 0.0;
        grad->phistagg[X][a] = (gridinfo_w[center].phia[a]);
      }
      grad->gradphi[Y][a]  = (gridinfo_w[right].phia[a] - gridinfo_w[center].phia[a])/deltay;
      grad->phistagg[Y][a] = 0.5*(gridinfo_w[right].phia[a] + gridinfo_w[center].phia[a]);
    }
    
    for (a=0; a < NUMPHASES; a++) {
      if (gidy > 0 && gidy <=(workers_mpi.end[Y]+2)) {
        grad->gradphi_c[Y][a] = (gridinfo_w[right].phia[a] - gridinfo_w[left].phia[a])/(2.0*deltay);
      } else {
        grad->gradphi_c[Y][a] = grad->gradphi[Y][a];
      }
      if (x > 0) {
        grad->gradphi_c[X][a] = 0.5*(grad->gradphi[X][a] + grad_back->gradphi[X][a]);
      }
    }
    if (CALCULATE_COMPOSITION) {
      if (grad->interface) {
        for (a=0; a < NUMPHASES; a++) {
          c_mu(gridinfo_w[center].compi, grad->phase_comp[a], T, a, c_guess[a][a]);
//           if ((a==0) || (a==NUMPHASES-1)) {
//             c_mu(gridinfo_w[center].compi, grad->phase_comp[a], T, a, c_guess[a][a]);
//           }
          if (!ISOTHERMAL) {
            c_mu(gridinfo_w[center].compi, c_tdt, T+DELTAT, a, ceq[a][a]);
            for (k=0; k < NUMCOMPONENTS-1; k++) {
              grad->dcbdT_phase[a][k] = c_tdt[k] - grad->phase_comp[a][k];
            }
          }
        }
//         for (a=1; a<NUMPHASES-1; a++) {
//           for (k=0; k< NUMCOMPONENTS-1; k++) {
//             grad->phase_comp[a][k] = grad->phase_comp[0][k];
//           }
//         }
        for (a=0; a < NUMPHASES; a++) {
          dc_dmu(gridinfo_w[center].compi, grad->phase_comp[a], T, a, grad->dcdmu_phase[a]);
        }
      } else {
        for (k=0; k < NUMCOMPONENTS-1; k++) {
          grad->phase_comp[grad->bulk_phase][k] = gridinfo_w[center].composition[k];
        }
        dc_dmu(gridinfo_w[center].compi, grad->phase_comp[grad->bulk_phase], T, grad->bulk_phase, grad->dcdmu_phase[grad->bulk_phase]);
      }
    }
  }
}
void calculate_gradients_phasefield_3D(long x, struct gradlayer **gradient, int CALCULATE_COMPOSITION) {
  long k, l, a, gidy;
  double temp;
  for (gidy=1; gidy < workers_mpi.layer_size-1; gidy++) {
          
    grad          =  &gradient[2][gidy];
    grad_back     =  &gradient[1][gidy];
    
    center        =  gidy   + (x)*workers_mpi.layer_size;
    front         =  gidy   + (x+1)*workers_mpi.layer_size;
    right         =  center + 1;
    left          =  center - 1;
  
    if (DIMENSION != 2) {
      top         =  center + workers_mpi.rows_y;
      bottom      =  center - workers_mpi.rows_y;
    }
    
    grad->interface = 1;
    grad->bulk_phase = 0;
    
    for (a=0; a < NUMPHASES; a++) {
      if (gridinfo_w[center].phia[a] == 1.0) {
        grad->bulk_phase=a;
        grad->interface = 0;
        break;
      }
    }
    if (!ISOTHERMAL) {
      T      = gridinfo_w[center].temperature;
      DELTAT = deltat*(-temperature_gradientY.GRADIENT*temperature_gradientY.velocity);
    }
    
//     printf("workers_mpi.start[X]=%ld, workers_mpi.end[X]=%ld, front=%ld, center=%ld, layer_size=%ld\n", workers_mpi.start[X], workers_mpi.end[X], front, center, workers_mpi.layer_size);

    
//     printf("grad->gradphi[X][0]=%p,grad->gradphi[X][1]=%p, NUMPHASES=%d\n", &grad->gradphi[X][0], &grad->gradphi[X][1], NUMPHASES);
    for (a=0; a < NUMPHASES; a++) {
      if (x <=(workers_mpi.end[X]+2)) {
        grad->gradphi[X][a]  = (gridinfo_w[front].phia[a]          - gridinfo_w[center].phia[a])/deltax;
        grad->phistagg[X][a] = 0.5*(gridinfo_w[front].phia[a]     + gridinfo_w[center].phia[a]);
      } else {
        grad->gradphi[X][a]  = 0.0;
        grad->phistagg[X][a] = (gridinfo_w[center].phia[a]);
      }
      
      if (((gidy + workers_mpi.rows_y)/workers_mpi.rows_y) < workers_mpi.rows_z) {
        grad->gradphi[Z][a]  =     (gridinfo_w[top].phia[a] - gridinfo_w[center].phia[a])/deltaz;
        grad->phistagg[Z][a] = 0.5*(gridinfo_w[top].phia[a] + gridinfo_w[center].phia[a]);
      } 
      
      grad->gradphi[Y][a]  =     (gridinfo_w[right].phia[a] - gridinfo_w[center].phia[a])/deltay;
      grad->phistagg[Y][a] = 0.5*(gridinfo_w[right].phia[a] + gridinfo_w[center].phia[a]);
    }
    
    for (a=0; a < NUMPHASES; a++) {
      if (gidy > 0 && gidy < (workers_mpi.layer_size-1)) {
        grad->gradphi_c[Y][a] = (gridinfo_w[right].phia[a] - gridinfo_w[left].phia[a])/(2.0*deltay);
      } else {
        grad->gradphi_c[Y][a] = grad->gradphi[Y][a];
      }
      
      if ((((gidy-workers_mpi.rows_y)/workers_mpi.rows_y) > 0) && (((gidy+workers_mpi.rows_y)/workers_mpi.rows_y) < (workers_mpi.rows_z))) {
        grad->gradphi_c[Z][a] = (gridinfo_w[top].phia[a]   - gridinfo_w[bottom].phia[a])/(2.0*deltaz);
      } else {
        grad->gradphi_c[Z][a] = grad->gradphi[Z][a]; 
      }
      if(x>0) {
        grad_back             =  &gradient[1][gidy];
        grad->gradphi_c[X][a] = 0.5*(grad->gradphi[X][a] + grad_back->gradphi[X][a]);
      } else {
        grad->gradphi_c[X][a] = grad->gradphi[X][a];
      }
    }
    if (CALCULATE_COMPOSITION) {
      if (grad->interface) {
        for (a=0; a < NUMPHASES; a++) {
          c_mu(gridinfo_w[center].compi, grad->phase_comp[a], T, a, c_guess[a][a]);
          if (!ISOTHERMAL) {
            c_mu(gridinfo_w[center].compi, c_tdt, T+DELTAT, a, ceq[a][a]);
            for (k=0; k < NUMCOMPONENTS-1; k++) {
              grad->dcbdT_phase[a][k] = c_tdt[k] - grad->phase_comp[a][k];
            }
          }
        }
        for (a=0; a < NUMPHASES; a++) {
          dc_dmu(gridinfo_w[center].compi, grad->phase_comp[a], T, a, grad->dcdmu_phase[a]);
        }
      } else {
        for (k=0; k < NUMCOMPONENTS-1; k++) {
          grad->phase_comp[grad->bulk_phase][k] = gridinfo_w[center].composition[k];
        }
        dc_dmu(gridinfo_w[center].compi, grad->phase_comp[grad->bulk_phase], T, grad->bulk_phase, grad->dcdmu_phase[grad->bulk_phase]);
      }
    }
  }
}
void calculate_gradients_stress_2D(long x, struct gradlayer **gradient) {
    long a,k;
    long index;
    
    for (gidy=0; gidy <= (workers_mpi.end[Y]+2); gidy++) {				//gidy=0 end[Y]+2
    
      grad          =  &gradient[2][gidy];
    
      center        =  gidy   + (x)*(workers_mpi.rows_y);				//center        =  gidy   + (x)*rows_y
      front         =  gidy   + (x+1)*workers_mpi.rows_y;				//front         =  gidy   + (x+1)*rows_y
      back          =  gidy   + (x-1)*workers_mpi.rows_y;				//back          =  gidy   + (x-1)*rows_y
      right         =  center + 1;
      left          =  center - 1;
        
      grad->eigen_strain[X]    = calculate_eigen_strain(center);
      grad->eigen_strain[Y]    = calculate_eigen_strain(center);
      
//       printf("grad->eigen_strain[X].xx=%le,grad->eigen_strain[Y].xx=%le,grad->eigen_strain[X].yy=%le,grad->eigen_strain[Y].yy=%le\n",grad->eigen_strain[X].xx,grad->eigen_strain[Y].xx,grad->eigen_strain[X].yy, grad->eigen_strain[Y].yy);
      

      if (x>0) {
//         grad->strain[X].xx = ((0.5/deltax)*(iter_gridinfo_w[front].disp[X][2] - iter_gridinfo_w[back].disp[X][2]))
//                               - grad->eigen_strain[X].xx;
         grad->strain[X].xx = ((0.5)*(iter_gridinfo_w[front].disp[X][2] - iter_gridinfo_w[back].disp[X][2]))
                              - grad->eigen_strain[X].xx;
      }
      if (gidy > 0) {  
//         grad->strain[X].yy = (0.5/deltay)*((iter_gridinfo_w[right].disp[Y][2] - iter_gridinfo_w[left].disp[Y][2]))  
//                                  - grad->eigen_strain[X].yy; 
//         grad->strain[Y].yy = ((0.5/deltay)*(iter_gridinfo_w[right].disp[Y][2] - iter_gridinfo_w[left].disp[Y][2])) 
//                               - grad->eigen_strain[Y].yy;
        grad->strain[X].yy = (0.5)*((iter_gridinfo_w[right].disp[Y][2] - iter_gridinfo_w[left].disp[Y][2]))  
                                 - grad->eigen_strain[X].yy; 
        grad->strain[Y].yy = ((0.5)*(iter_gridinfo_w[right].disp[Y][2] - iter_gridinfo_w[left].disp[Y][2])) 
                              - grad->eigen_strain[Y].yy;                    
                              
      }
      
      if ((x > 0) && (gidy >0)) {
//         grad->strain[Y].xx = ((0.5/deltax)*(iter_gridinfo_w[front].disp[X][2] - iter_gridinfo_w[back].disp[X][2]))
//                                           - grad->eigen_strain[Y].xx;
//         grad->strain[X].xy =  (0.25)*((iter_gridinfo_w[right].disp[X][2]      - iter_gridinfo_w[left].disp[X][2])/deltay  
//                                     + (iter_gridinfo_w[front].disp[Y][2]      - iter_gridinfo_w[back].disp[Y][2])/deltax); 
        
        grad->strain[Y].xx =  ((0.5)*(iter_gridinfo_w[front].disp[X][2] - iter_gridinfo_w[back].disp[X][2]))
                                     - grad->eigen_strain[Y].xx;
        grad->strain[X].xy =  (0.25)*((iter_gridinfo_w[right].disp[X][2]      - iter_gridinfo_w[left].disp[X][2]) 
                                    + (iter_gridinfo_w[front].disp[Y][2]      - iter_gridinfo_w[back].disp[Y][2])); 
     
        grad->strain[Y].xy = grad->strain[X].xy;
      }
    }
// 
    for (gidy=0; gidy <= (workers_mpi.end[Y]+2); gidy++) {				//gidy=0  end[Y]+2
     grad          =  &gradient[2][gidy];
     grad_back     =  &gradient[1][gidy];
    
     center        =  gidy    + (x)*(workers_mpi.rows_y);
  
     grad->stiffness_c[X]= calculate_stiffness_n(center);
     
     grad->stiffness_c[Y] = calculate_stiffness_n(center);
   }
}
void calculate_gradients_stress_3D(long x, struct gradlayer **gradient) {
    long a,k;
    long index;
    
    for (gidy=0; gidy <= (workers_mpi.layer_size-1); gidy++) {				//gidy=0 end[Y]+2
    
      grad          =  &gradient[2][gidy];
    
      center        =  gidy   + (x)*(workers_mpi.layer_size);				//center        =  gidy   + (x)*rows_y
      front         =  gidy   + (x+1)*workers_mpi.layer_size;				//front         =  gidy   + (x+1)*rows_y
      back          =  gidy   + (x-1)*workers_mpi.layer_size;				//back          =  gidy   + (x-1)*rows_y
      top           =  center + workers_mpi.rows_y;
      bottom        =  center - workers_mpi.rows_y;
      right         =  center + 1;
      left          =  center - 1;
        
      grad->eigen_strain[X]    = calculate_eigen_strain(center);
//       grad->eigen_strain[Y]    = calculate_eigen_strain(center);
      grad->eigen_strain[Y]    = grad->eigen_strain[X];
      grad->eigen_strain[Z]    = grad->eigen_strain[X];
      
      
  //For normal strains
      if (x>0) {
        grad->strain[X].xx = ((0.5)*(iter_gridinfo_w[front].disp[X][2] - iter_gridinfo_w[back].disp[X][2]))
                              - grad->eigen_strain[X].xx;
        grad->strain[Y].xx =  grad->strain[X].xx;
        grad->strain[Z].xx =  grad->strain[X].xx;
      }
    
      if (gidy > 0) {
        grad->strain[X].yy = (0.5)*((iter_gridinfo_w[right].disp[Y][2] - iter_gridinfo_w[left].disp[Y][2]))  
//                                  + (iter_gridinfo_w[rightfront].disp[Y][2] - iter_gridinfo_w[leftfront].disp[Y][2])))
                                 - grad->eigen_strain[X].yy; 
        grad->strain[Y].yy = grad->strain[X].yy;
        grad->strain[Z].yy = grad->strain[X].yy;
      }
      
      if ((((gidy-workers_mpi.rows_y)/workers_mpi.rows_y) > 0) && (((gidy+workers_mpi.rows_y)/workers_mpi.rows_y) < (workers_mpi.rows_z))) {
        grad->strain[X].zz = (0.5)*((iter_gridinfo_w[top].disp[Z][2] - iter_gridinfo_w[bottom].disp[Z][2]))  
                                 - grad->eigen_strain[X].zz; 
        grad->strain[Y].zz = grad->strain[X].zz;
        grad->strain[Z].zz = grad->strain[X].zz;                         
      }
      
  //For shear strains
      if ((x > 0) && (gidy > 0) ) {
        grad->strain[X].xy =  (0.25)*((iter_gridinfo_w[right].disp[X][2]      - iter_gridinfo_w[left].disp[X][2])  
                                    + (iter_gridinfo_w[front].disp[Y][2]      - iter_gridinfo_w[back].disp[Y][2]));
        grad->strain[Y].xy =  grad->strain[X].xy;
        grad->strain[Z].xy =  grad->strain[X].xy;
      }
      if ((x > 0) && (((gidy-workers_mpi.rows_y)/workers_mpi.rows_y) > 0) && (((gidy+workers_mpi.rows_y)/workers_mpi.rows_y) < (workers_mpi.rows_z)) ) {
        grad->strain[X].xz =  (0.25)*((iter_gridinfo_w[top].disp[X][2]        - iter_gridinfo_w[bottom].disp[X][2])  
                                    + (iter_gridinfo_w[front].disp[Z][2]      - iter_gridinfo_w[back].disp[Z][2]));
        grad->strain[Y].xz =  grad->strain[X].xz;
        grad->strain[Z].xz =  grad->strain[X].xz;
      }
      
      if ((gidy > 0) && (((gidy-workers_mpi.rows_y)/workers_mpi.rows_y) > 0) && (((gidy+workers_mpi.rows_y)/workers_mpi.rows_y) < (workers_mpi.rows_z)) ) {
        grad->strain[X].yz =  (0.25)*((iter_gridinfo_w[right].disp[Z][2]      - iter_gridinfo_w[left].disp[Z][2])  
                                    + (iter_gridinfo_w[top].disp[Y][2]        - iter_gridinfo_w[bottom].disp[Y][2]));
        grad->strain[Y].yz =  grad->strain[X].yz;
        grad->strain[Z].yz =  grad->strain[X].yz;
      }
    }

  
// For stiffness
    for (gidy=1; gidy <= (workers_mpi.layer_size-1); gidy++) {				//gidy=0  end[Y]+2
     grad          =  &gradient[2][gidy];
     grad_back     =  &gradient[1][gidy];
    
     center        =  gidy    + (x)*(workers_mpi.layer_size);
  
     grad->stiffness_c[X] = calculate_stiffness_n(center);
     grad->stiffness_c[Y] = grad->stiffness_c[X];
     grad->stiffness_c[Z] = grad->stiffness_c[X];
   }
}

void calculate_diffusion_potential(long x, struct gradlayer **gradient) {
  long k, l, a, gidy;
  for (gidy=0; gidy < (workers_mpi.layer_size); gidy++) {
          
    grad          =  &gradient[2][gidy];
    center        =  gidy   + (x)*workers_mpi.layer_size;
 
    interface = 1;
    
    for (a=0; a < NUMPHASES; a++) {
      if (gridinfo_w[center].phia[a] == 1.0) {
        bulk_phase=a;
        interface = 0;
        break;
      }
    }
    if(!ISOTHERMAL){
      T = gridinfo_w[center].temperature;
    }
    
    if (interface) {
      if (FUNCTION_F==2) {
    //    function_F_02_getMu(gridinfo_w[center].compi, gridinfo_w[center].composition, gridinfo_w[center].phia, T, grad);
//          c_mu(gridinfo_w[center].compi, grad->phase_comp[a], T, a, ceq[a][a]); 
        for (a=0; a < NUMPHASES; a++) {
//           if ((a==0) || (a==NUMPHASES-1)) {
            c_mu(gridinfo_w[center].compi, grad->phase_comp[a], T, a, c_guess[a][a]);
//           }
        }
//         for (a=1; a<NUMPHASES-1; a++) {
//           for (k=0; k< NUMCOMPONENTS-1; k++) {
//             grad->phase_comp[a][k] = grad->phase_comp[0][k];
//           }
//         }
      } else {
        for (a=0; a < NUMPHASES; a++) {
//           if ((a==0) || (a==NUMPHASES-1)) {
            c_mu(gridinfo_w[center].compi, grad->phase_comp[a], T, a, c_guess[a][a]);
//           }
        }
//         for (a=1; a<NUMPHASES-1; a++) {
//           for (k=0; k< NUMCOMPONENTS-1; k++) {
//             grad->phase_comp[a][k] = grad->phase_comp[0][k];
//           }
//         }
      }
    }
  }
}
void calculate_gradients_concentration_2D(long x, struct gradlayer **gradient) {
  long k, l, a, gidy;
  for (gidy=0; gidy <= (workers_mpi.end[Y]+2); gidy++) {
          
    grad          =  &gradient[1][gidy];
    grad_back     =  &gradient[0][gidy];
    grad_front    =  &gradient[2][gidy];
    grad_right    =  &gradient[1][gidy+1];
    
    center        =  gidy   + (x)*workers_mpi.layer_size;
    front         =  gidy   + (x+1)*workers_mpi.layer_size;
    right         =  center + 1;
    left          =  center - 1;
    if (DIMENSION != 2) {
      top         =  center + workers_mpi.rows_y;
      bottom      =  center - workers_mpi.rows_y;
    }
    
//     interface = 1;
    
//     for (a=0; a < NUMPHASES; a++) {
//       if (gridinfo_w[center].phia[a] == 1.0) {
// 	bulk_phase=a;
// 	interface = 0;
// 	break;
//       }
//     }
    
    if(ISOTHERMAL) {
      for (k=0; k < NUMCOMPONENTS-1; k++) {
        if (x <=(workers_mpi.end[X]+2)) {
          grad->gradchempot[X][k]  = (gridinfo_w[front].compi[k] - gridinfo_w[center].compi[k])/deltax;
        } else {
          grad->gradchempot[X][k] = 0.0;
        }
        grad->gradchempot[Y][k]   = (gridinfo_w[right].compi[k]  - gridinfo_w[center].compi[k])/deltay;
        
        for (l=0; l < NUMCOMPONENTS-1; l++) {
        if (x <=(workers_mpi.end[X]+2)) {
          grad->Dmid[X][k][l] = (D(gridinfo_w, grad_front, T, front, k, l) + D(gridinfo_w, grad, T, center, k, l))*0.5;
          } else {
            grad->Dmid[X][k][l] = D(gridinfo_w, grad, T, center, k, l);
          }
          grad->Dmid[Y][k][l] = (D(gridinfo_w, grad_right, T, right, k, l) + D(gridinfo_w, grad, T, center, k, l))*0.5;
        }
      }
    } else {
      for (k=0; k < NUMCOMPONENTS-1; k++) {
        if(x<=(workers_mpi.end[X]+2)) {
          grad->gradchempot[X][k] = (gridinfo_w[front].compi[k] - gridinfo_w[center].compi[k])/deltax;
        } else {
          grad->gradchempot[X][k] = 0.0;
        }
        grad->gradchempot[Y][k]   = (gridinfo_w[right].compi[k] - gridinfo_w[center].compi[k])/deltay;
        for (l=0; l < NUMCOMPONENTS-1; l++) {  
          if(x<=(workers_mpi.end[X]+2)) {
            grad->Dmid[X][k][l]   = (D(gridinfo_w, grad_front, gridinfo_w[front].temperature, front, k, l) + D(gridinfo_w, grad, gridinfo_w[center].temperature, center, k, l))*0.5;
          } else {
            grad->Dmid[X][k][l]   = D(gridinfo_w, grad, gridinfo_w[center].temperature, center, k, l);
          }
          grad->Dmid[Y][k][l]     = (D(gridinfo_w, grad_right, gridinfo_w[right].temperature, right, k, l) + D(gridinfo_w, grad, gridinfo_w[center].temperature, center, k, l))*0.5;
        }
      }
    } 
//     if (grad->interface) {
// //      if (interface) {
//       if (!ISOTHERMAL) {
//         T = gridinfo_w[center].temperature;
//         DELTAT = deltat*(-temperature_gradientY.GRADIENT*temperature_gradientY.velocity);
//         for (a=0; a<NUMPHASES; a++) {
//           c_mu(gridinfo_w[center].compi, c_tdt, T+DELTAT, a, ceq[a][a]);
//           for (k=0; k < NUMCOMPONENTS-1; k++) {
//             grad->dcbdT_phase[a][k] = c_tdt[k] - grad->phase_comp[a][k];
//           }
//         }
//       }
//     }
  }
}
void calculate_gradients_concentration_3D(long x, struct gradlayer **gradient) {
  long k, l, a, gidy;
  
  for (gidy=1; gidy < (workers_mpi.layer_size-1); gidy++) {
    
    grad          =  &gradient[1][gidy];
    grad_back     =  &gradient[0][gidy];
    grad_front    =  &gradient[2][gidy];
    grad_right    =  &gradient[1][gidy+1];
    if ((gidy+workers_mpi.rows_y) < workers_mpi.layer_size) {
      grad_top      =  &gradient[1][gidy+workers_mpi.rows_y];
    }
    
    center        =  gidy   + (x)*workers_mpi.layer_size;
    front         =  gidy   + (x+1)*workers_mpi.layer_size;
    right         =  center + 1;
    left          =  center - 1;
    top           =  center + workers_mpi.rows_y;
    bottom        =  center - workers_mpi.rows_y;
    
        
    if(ISOTHERMAL) {
      for (k=0; k < NUMCOMPONENTS-1; k++) {
        if (x <=(workers_mpi.end[X]+2)) {
          grad->gradchempot[X][k]  = (gridinfo_w[front].compi[k] - gridinfo_w[center].compi[k])/deltax;
        } else {
          grad->gradchempot[X][k] = 0.0;
        }
        if (((gidy + workers_mpi.rows_y)/workers_mpi.rows_y) < workers_mpi.rows_z) {
          grad->gradchempot[Z][k]  = (gridinfo_w[top].compi[k] - gridinfo_w[center].compi[k])/deltaz;
        }
        grad->gradchempot[Y][k]   = (gridinfo_w[right].compi[k]  - gridinfo_w[center].compi[k])/deltay;
        
        for (l=0; l < NUMCOMPONENTS-1; l++) {
          if (x <=(workers_mpi.end[X]+2)) {
            grad->Dmid[X][k][l] = (D(gridinfo_w, grad_front, T, front, k, l) + D(gridinfo_w, grad, T, center, k, l))*0.5;
          } else {
            grad->Dmid[X][k][l] = D(gridinfo_w, grad, T, center, k, l);
          }
          if (((gidy + workers_mpi.rows_y)/workers_mpi.rows_y) < workers_mpi.rows_z) {
            grad->Dmid[Z][k][l] = (D(gridinfo_w, grad_top, T, top, k, l) + D(gridinfo_w, grad, T, center, k, l))*0.5;
          }
          
          grad->Dmid[Y][k][l] = (D(gridinfo_w, grad_right, T, right, k, l) + D(gridinfo_w, grad, T, center, k, l))*0.5;
        }
      }
    } else {
      for (k=0; k < NUMCOMPONENTS-1; k++) {
        if(x<=(workers_mpi.end[X]+2)) {
          grad->gradchempot[X][k] = (gridinfo_w[front].compi[k] - gridinfo_w[center].compi[k])/deltax;
        } else {
          grad->gradchempot[X][k] = 0.0;
        }
        if (((gidy + workers_mpi.rows_y)/workers_mpi.rows_y) < workers_mpi.rows_z) {
          grad->gradchempot[Z][k]  = (gridinfo_w[top].compi[k] - gridinfo_w[center].compi[k])/deltaz;
        }
        grad->gradchempot[Y][k]   = (gridinfo_w[right].compi[k] - gridinfo_w[center].compi[k])/deltay;
        
        for (l=0; l < NUMCOMPONENTS-1; l++) {  
          if(x<=(workers_mpi.end[X]+2)) {
            grad->Dmid[X][k][l]   = (D(gridinfo_w, grad_front, gridinfo_w[front].temperature, front, k, l) + D(gridinfo_w, grad, gridinfo_w[center].temperature, center, k, l))*0.5;
          } else {
            grad->Dmid[X][k][l]   = D(gridinfo_w, grad, gridinfo_w[center].temperature, center, k, l);
          }
          
          if (((gidy + workers_mpi.rows_y)/workers_mpi.rows_y) < workers_mpi.rows_z) {
            grad->Dmid[Z][k][l] = (D(gridinfo_w, grad_top, gridinfo_w[top].temperature, top, k, l) + D(gridinfo_w, grad, gridinfo_w[center].temperature, center, k, l))*0.5;
          }
          grad->Dmid[Y][k][l]     = (D(gridinfo_w, grad_right, gridinfo_w[right].temperature, right, k, l) + D(gridinfo_w, grad, gridinfo_w[center].temperature, center, k, l))*0.5;
        }
      }
    } 
  }
}
struct Stiffness_cubic calculate_stiffness_n(long index) {
  long a;
  struct Stiffness_cubic stiffness;
  stiffness.C11 = 0.0;
  stiffness.C12 = 0.0;
  stiffness.C44 = 0.0; 
  for (a=0; a<NUMPHASES; a++) {
//     stiffness.C11 += (stiffness_phase[a].C11)*hphi(gridinfo_w[index].phia,a);
//     stiffness.C12 += (stiffness_phase[a].C12)*hphi(gridinfo_w[index].phia,a);
//     stiffness.C44 += (stiffness_phase[a].C44)*hphi(gridinfo_w[index].phia,a);
    stiffness.C11 += (stiffness_phase_n[a].C11)*gridinfo_w[index].phia[a];
    stiffness.C12 += (stiffness_phase_n[a].C12)*gridinfo_w[index].phia[a];
    stiffness.C44 += (stiffness_phase_n[a].C44)*gridinfo_w[index].phia[a];
  }
  return stiffness;
}
struct Stiffness_cubic calculate_stiffness(long index) {
  long a;
  struct Stiffness_cubic stiffness;
  stiffness.C11 = 0.0;
  stiffness.C12 = 0.0;
  stiffness.C44 = 0.0; 
  for (a=0; a<NUMPHASES; a++) {
//     stiffness.C11 += (stiffness_phase[a].C11)*hphi(gridinfo_w[index].phia,a);
//     stiffness.C12 += (stiffness_phase[a].C12)*hphi(gridinfo_w[index].phia,a);
//     stiffness.C44 += (stiffness_phase[a].C44)*hphi(gridinfo_w[index].phia,a);
    stiffness.C11 += (stiffness_phase[a].C11)*gridinfo_w[index].phia[a];
    stiffness.C12 += (stiffness_phase[a].C12)*gridinfo_w[index].phia[a];
    stiffness.C44 += (stiffness_phase[a].C44)*gridinfo_w[index].phia[a];
  }
  return stiffness;
}

struct symmetric_tensor calculate_eigen_strain(long index) {
  long a,b;
//   double sum=0;
  struct symmetric_tensor eigen_strain;
  eigen_strain.xx = 0.0;
  eigen_strain.yy = 0.0;
  eigen_strain.zz = 0.0;
  eigen_strain.yz = 0.0;
  eigen_strain.xz = 0.0;
  eigen_strain.xy = 0.0;
  for (a=0; a<NUMPHASES; a++) {
//     eigen_strain.xx += eigen_strain_phase[a].xx*hphi(gridinfo_w[index].phia,a);
//     eigen_strain.yy += eigen_strain_phase[a].yy*hphi(gridinfo_w[index].phia,a);
//     eigen_strain.zz += eigen_strain_phase[a].zz*hphi(gridinfo_w[index].phia,a);
//     eigen_strain.yz += eigen_strain_phase[a].yz*hphi(gridinfo_w[index].phia,a);
//     eigen_strain.xz += eigen_strain_phase[a].xz*hphi(gridinfo_w[index].phia,a);
//     eigen_strain.xy += eigen_strain_phase[a].xy*hphi(gridinfo_w[index].phia,a);
    eigen_strain.xx += eigen_strain_phase[a].xx*gridinfo_w[index].phia[a];
    eigen_strain.yy += eigen_strain_phase[a].yy*gridinfo_w[index].phia[a];
    eigen_strain.zz += eigen_strain_phase[a].zz*gridinfo_w[index].phia[a];
    eigen_strain.yz += eigen_strain_phase[a].yz*gridinfo_w[index].phia[a];
    eigen_strain.xz += eigen_strain_phase[a].xz*gridinfo_w[index].phia[a];
    eigen_strain.xy += eigen_strain_phase[a].xy*gridinfo_w[index].phia[a];
  }
  return eigen_strain;
}

#endif
