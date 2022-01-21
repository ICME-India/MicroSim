#ifndef CALCULATE_GRADIENTS_H_
#define CALCULATE_GRADIENTS_H_

void calculate_gradients_2D(long x, struct gradlayer **gradient) {
//Calculation relevant to  phase-field
 long k, l, a;
  for (gidy=0; gidy <= end[Y]+2; gidy++) {
	  
    grad          =  &gradient[2][gidy];
    grad_back     =  &gradient[1][gidy];
    
    center        =  gidy   + (x)*layer_size;
    front         =  gidy   + (x+1)*layer_size;
    right         =  center + 1;
    left          =  center - 1;
    if (DIMENSION != 2) {
      top         =  center + rows_y;
      bottom      =  center - rows_y;
    }
    
    
    for (a=0; a < NUMPHASES; a++) {
      if (x <=(end[X]+2)) {
	grad->gradphi[X][a]  = (gridinfo[front].phia[a]     - gridinfo[center].phia[a])/deltax;
	grad->phistagg[X][a] = 0.5*(gridinfo[front].phia[a] + gridinfo[center].phia[a]);
      } else {
	grad->gradphi[X][a]  = 0.0;
	grad->phistagg[X][a] = (gridinfo[center].phia[a]);
      }
      grad->gradphi[Y][a]    = (gridinfo[right].phia[a]     - gridinfo[center].phia[a])/deltay;
      grad->phistagg[Y][a]   = 0.5*(gridinfo[right].phia[a] + gridinfo[center].phia[a]);
    }
      
    if(ISOTHERMAL) {
      for (k=0; k < NUMCOMPONENTS-1; k++) {
        if (x <=(end[X]+2)) {
          grad->gradchempot[X][k] = (gridinfo[front].compi[k] - gridinfo[center].compi[k])/deltax;
        } else {
          grad->gradchempot[X][k] = 0.0;
        }
        grad->gradchempot[Y][k] = (gridinfo[right].compi[k] - gridinfo[center].compi[k])/deltay;
        for (l=0; l < NUMCOMPONENTS-1; l++) {
        if (x <=(end[X]+2)) {
          grad->Dmid[X][k][l] = (D(gridinfo, T, front, k, l) + D(gridinfo, T, center, k, l))*0.5;
        } else {
            grad->Dmid[X][k][l] = D(gridinfo, T, center, k, l);
        }
          grad->Dmid[Y][k][l] = (D(gridinfo, T, right, k, l) + D(gridinfo, T, center, k, l))*0.5;
        }
      }
    } else {
      for (k=0; k < NUMCOMPONENTS-1; k++) {
        if(x<=(end[X]+2)) {
          grad->gradchempot[X][k] = (gridinfo[front].compi[k] - gridinfo[center].compi[k])/deltax;
        } else {
          grad->gradchempot[X][k] = 0.0;
        }
        grad->gradchempot[Y][k] = (gridinfo[right].compi[k] - gridinfo[center].compi[k])/deltay;
        for (l=0; l < NUMCOMPONENTS-1; l++) {  
          if(x<=(end[X]+2)) {
            grad->Dmid[X][k][l] = (D(gridinfo, gridinfo[front].temperature, front, k, l) + D(gridinfo, gridinfo[center].temperature, center, k, l))*0.5;
          } else {
            grad->Dmid[X][k][l] = D(gridinfo, gridinfo[center].temperature, center, k, l);
          }
          grad->Dmid[Y][k][l] = (D(gridinfo, gridinfo[right].temperature, right, k, l) + D(gridinfo, gridinfo[center].temperature, center, k, l))*0.5;
        }
      }
   }
    for (a=0; a < NUMPHASES; a++) {
      if (gidy > 0 && gidy <=(end[Y]+2)) {
        grad->gradphi_c[Y][a] = (gridinfo[right].phia[a] - gridinfo[left].phia[a])/(2.0*deltay);
      } else {
        grad->gradphi_c[Y][a] = grad->gradphi[Y][a];
      }
      grad->gradphi_c[X][a]  = 0.5*(grad->gradphi[X][a] + grad_back->gradphi[X][a]);
    }
  }
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
  for (gidy=0; gidy <= (end[Y]+2); gidy++) {
          
    grad          =  &gradient[2][gidy];
    grad_back     =  &gradient[1][gidy];
    
    center        =  gidy   + (x)*layer_size;
    front         =  gidy   + (x+1)*layer_size;
    right         =  center + 1;
    left          =  center - 1;
    if (DIMENSION != 2) {
      top         =  center + rows_y;
      bottom      =  center - rows_y;
    }
    
    for (a=0; a < NUMPHASES; a++) {
      if (x <=(end[X]+2)) {
        grad->gradphi[X][a]  = (gridinfo[front].phia[a]     - gridinfo[center].phia[a])/deltax;
        grad->phistagg[X][a] = 0.5*(gridinfo[front].phia[a] + gridinfo[center].phia[a]);
      } else {
        grad->gradphi[X][a]  = 0.0;
        grad->phistagg[X][a] = (gridinfo[center].phia[a]);
      }
      grad->gradphi[Y][a]  = (gridinfo[right].phia[a] - gridinfo[center].phia[a])/deltay;
      grad->phistagg[Y][a] = 0.5*(gridinfo[right].phia[a] + gridinfo[center].phia[a]);
    }
    
    for (a=0; a < NUMPHASES; a++) {
      if (gidy > 0 && gidy <=(end[Y]+2)) {
        grad->gradphi_c[Y][a] = (gridinfo[right].phia[a] - gridinfo[left].phia[a])/(2.0*deltay);
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
  for (gidy=0; gidy <= (end[Y]+2); gidy++) {
          
    grad          =  &gradient[2][gidy];
    grad_back     =  &gradient[1][gidy];
    
    center        =  gidy   + (x)*layer_size;
    front         =  gidy   + (x+1)*layer_size;
    right         =  center + 1;
    left          =  center - 1;
    if (DIMENSION != 2) {
      top         =  center + rows_y;
      bottom      =  center - rows_y;
    }
    
    if(ISOTHERMAL) {
      for (k=0; k < NUMCOMPONENTS-1; k++) {
        if (x <=(end[X]+2)) {
          grad->gradchempot[X][k]  = (gridinfo[front].compi[k] - gridinfo[center].compi[k])/deltax;
        } else {
          grad->gradchempot[X][k] = 0.0;
        }
        grad->gradchempot[Y][k]   = (gridinfo[right].compi[k]     - gridinfo[center].compi[k])/deltay;
        
        for (l=0; l < NUMCOMPONENTS-1; l++) {
        if (x <=(end[X]+2)) {
          grad->Dmid[X][k][l] = (D(gridinfo, T, front, k, l) + D(gridinfo, T, center, k, l))*0.5;
        } else {
          grad->Dmid[X][k][l] = D(gridinfo, T, center, k, l);
        }
          grad->Dmid[Y][k][l] = (D(gridinfo, T, right, k, l) + D(gridinfo, T, center, k, l))*0.5;
        }
      }
    } else {
      for (k=0; k < NUMCOMPONENTS-1; k++) {
        if(x<=(end[X]+2)) {
          grad->gradchempot[X][k] = (gridinfo[front].compi[k] - gridinfo[center].compi[k])/deltax;
        } else {
          grad->gradchempot[X][k] = 0.0;
        }
        grad->gradchempot[Y][k]   = (gridinfo[right].compi[k] - gridinfo[center].compi[k])/deltay;
        for (l=0; l < NUMCOMPONENTS-1; l++) {  
          if(x<=(end[X]+2)) {
            grad->Dmid[X][k][l]   = (D(gridinfo, gridinfo[front].temperature, front, k, l) + D(gridinfo, gridinfo[center].temperature, center, k, l))*0.5;
          } else {
            grad->Dmid[X][k][l]   = D(gridinfo, gridinfo[center].temperature, center, k, l);
          }
          grad->Dmid[Y][k][l]     = (D(gridinfo, gridinfo[right].temperature, right, k, l) + D(gridinfo, gridinfo[center].temperature, center, k, l))*0.5;
        }
      }
    } 
  }
}
#endif