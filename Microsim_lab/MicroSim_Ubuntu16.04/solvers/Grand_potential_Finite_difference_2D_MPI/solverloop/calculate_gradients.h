#ifndef CALCULATE_GRADIENTS_H_
#define CALCULATE_GRADIENTS_H_

void calculate_gradients_2D(long x, struct gradlayer **gradient) {
//Calculation relevant to  phase-field
 long k, l, a;
  for (gidy=0; gidy <= workers_mpi.end[Y]+2; gidy++) {
	  
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
    for (a=0; a < NUMPHASES; a++) {
      if (x <=(workers_mpi.end[X]+2)) {
	grad->gradphi[X][a]  = (gridinfo_w[front].phia[a]     - gridinfo_w[center].phia[a])/deltax;
	grad->phistagg[X][a] = 0.5*(gridinfo_w[front].phia[a] + gridinfo_w[center].phia[a]);
      } else {
	grad->gradphi[X][a]  = 0.0;
	grad->phistagg[X][a] = (gridinfo_w[center].phia[a]);
      }
      grad->gradphi[Y][a]    = (gridinfo_w[right].phia[a]     - gridinfo_w[center].phia[a])/deltay;
      grad->phistagg[Y][a]   = 0.5*(gridinfo_w[right].phia[a] + gridinfo_w[center].phia[a]);
    }
//    
    if(ISOTHERMAL) {
      for (k=0; k < NUMCOMPONENTS-1; k++) {
        if (x <=(workers_mpi.end[X]+2)) {
          grad->gradchempot[X][k] = (gridinfo_w[front].compi[k] - gridinfo_w[center].compi[k])/deltax;
        } else {
          grad->gradchempot[X][k] = 0.0;
        }
        grad->gradchempot[Y][k] = (gridinfo_w[right].compi[k] - gridinfo_w[center].compi[k])/deltay;
        for (l=0; l < NUMCOMPONENTS-1; l++) {
        if (x <=(workers_mpi.end[X]+2)) {
          grad->Dmid[X][k][l] = (D(gridinfo_w, T, front, k, l) + D(gridinfo_w, T, center, k, l))*0.5;
        } else {
            grad->Dmid[X][k][l] = D(gridinfo_w, T, center, k, l);
        }
          grad->Dmid[Y][k][l] = (D(gridinfo_w, T, right, k, l) + D(gridinfo_w, T, center, k, l))*0.5;
        }
      }
    } else {
      for (k=0; k < NUMCOMPONENTS-1; k++) {
        if(x<=(workers_mpi.end[X]+2)) {
          grad->gradchempot[X][k] = (gridinfo_w[front].compi[k] - gridinfo_w[center].compi[k])/deltax;
        } else {
          grad->gradchempot[X][k] = 0.0;
        }
        grad->gradchempot[Y][k] = (gridinfo_w[right].compi[k] - gridinfo_w[center].compi[k])/deltay;
        for (l=0; l < NUMCOMPONENTS-1; l++) {  
          if(x<=(workers_mpi.end[X]+2)) {
            grad->Dmid[X][k][l] = (D(gridinfo_w, gridinfo_w[front].temperature, front, k, l) + D(gridinfo_w, gridinfo_w[center].temperature, center, k, l))*0.5;
          } else {
            grad->Dmid[X][k][l] = D(gridinfo_w, gridinfo_w[center].temperature, center, k, l);
          }
          grad->Dmid[Y][k][l] = (D(gridinfo_w, gridinfo_w[right].temperature, right, k, l) + D(gridinfo_w, gridinfo_w[center].temperature, center, k, l))*0.5;
        }
      }
   }
    for (a=0; a < NUMPHASES; a++) {
      if (gidy > 0 && gidy <=(workers_mpi.end[Y]+2)) {
        grad->gradphi_c[Y][a] = (gridinfo_w[right].phia[a] - gridinfo_w[left].phia[a])/(2.0*deltay);
      } else {
        grad->gradphi_c[Y][a] = grad->gradphi[Y][a];
      }
      if (x > 0) {
        grad->gradphi_c[X][a]  = 0.5*(grad->gradphi[X][a] + grad_back->gradphi[X][a]);
//         grad->gradphi_c[X][a]  = grad->gradphi[X][a];
      }
    }
  }
//   if(taskid==1) {
//     printf("taskid=%d,Coming here; x=%ld, y=%ld, end[X]=%ld, end[Y]=%ld\n", taskid,x,gidy,workers_mpi.end[X],workers_mpi.end[Y]);
//   }
  
//   for (gidy=0; gidy <= (end[Y]+2); gidy++) {
//     grad          = &gradient[2][gidy];
//     grad_back     = &gradient[1][gidy];
// //       grad_right    = &gradient[2][gidy+1];
//     grad_left     = &gradient[2][gidy-1];
//   
//     for(a=0; a<NUMPHASES; a++) {
//       if((x > 0) && (x<=(end[X]+2)) && (gidy > 0) && (gidy <= (end[Y]+2))) {
// 	grad->d2gradphi[a]  = (grad->gradphi[X][a] -  grad_back->gradphi[X][a])/deltax;
// 	grad->d2gradphi[a] += (grad->gradphi[Y][a] -  grad_left->gradphi[Y][a])/deltay;
//       }
//     }
//   }
}
void calculate_gradients_phasefield_2D(long x, struct gradlayer **gradient) {
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
  }
}
void calculate_gradients_concentration_2D(long x, struct gradlayer **gradient) {
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
          grad->Dmid[X][k][l] = (D(gridinfo_w, T, front, k, l) + D(gridinfo_w, T, center, k, l))*0.5;
          } else {
            grad->Dmid[X][k][l] = D(gridinfo_w, T, center, k, l);
          }
          grad->Dmid[Y][k][l] = (D(gridinfo_w, T, right, k, l) + D(gridinfo_w, T, center, k, l))*0.5;
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
            grad->Dmid[X][k][l]   = (D(gridinfo_w, gridinfo_w[front].temperature, front, k, l) + D(gridinfo_w, gridinfo_w[center].temperature, center, k, l))*0.5;
          } else {
            grad->Dmid[X][k][l]   = D(gridinfo_w, gridinfo_w[center].temperature, center, k, l);
          }
          grad->Dmid[Y][k][l]     = (D(gridinfo_w, gridinfo_w[right].temperature, right, k, l) + D(gridinfo_w, gridinfo_w[center].temperature, center, k, l))*0.5;
        }
      }
    } 
  }
}
#endif