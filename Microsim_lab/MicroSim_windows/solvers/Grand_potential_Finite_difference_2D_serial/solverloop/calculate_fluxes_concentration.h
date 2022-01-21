#ifndef CALCULATE_FLUXES_CONCENTRATION_H_
#define CALCULATE_FLUXES_CONCENTRATION_H_

void calculate_fluxes_concentration_2D(long x, struct gradlayer **gradient) {
 double cl_front, cl_center, cs_front, cs_center, cl_right, cs_right;
 double DET;
 for (gidy=0; gidy <=(end[Y]+2); gidy++) {
    
    grad          =  &gradient[0][gidy];
    grad_right    =  &gradient[0][gidy+1];
    grad_front    =  &gradient[1][gidy];
    
    center        =  gidy   + (x)*layer_size;
    front         =  gidy   + (x+1)*layer_size;
    right         =  center + 1;
    left          =  center - 1;
    
    if(!ISOTHERMAL) {
      T = gridinfo[center].temperature;
//       init_propertymatrices(T);
    }
    
    //For multi-phase antitrapping. The last phase is by default the liquid phase.
    //Generating the normal vector corresponding to the liquid contour
    gradphix_l      = 0.5*(grad->gradphi_c[X][NUMPHASES-1] + grad_right->gradphi_c[X][NUMPHASES-1]);
    normgradphiy_l  = sqrt(gradphix_l*gradphix_l + grad->gradphi[Y][NUMPHASES-1]*grad->gradphi[Y][NUMPHASES-1]);
    
    gradphiy_l      = 0.5*(grad->gradphi_c[Y][NUMPHASES-1] + grad_front->gradphi_c[Y][NUMPHASES-1]);
    normgradphix_l  = sqrt(gradphiy_l*gradphiy_l + grad->gradphi[X][NUMPHASES-1]*grad->gradphi[X][NUMPHASES-1]);
    
    //Initializing the jats to zero
    
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      grad->jat[X][k] = 0.0;
      grad->jat[Y][k] = 0.0;
    }
    for (a=0; a < NUMPHASES-1; a++) {
    
      gradphix      = 0.5*(grad->gradphi_c[X][a] + grad_right->gradphi_c[X][a]);
      normgradphi   = sqrt(gradphix*gradphix + grad->gradphi[Y][a]*grad->gradphi[Y][a]);
      
      if (OBSTACLE) {		  
        if (gridinfo[right].phia[a]*(1.0-gridinfo[right].phia[a]) > 0.0) {
          s_phi_right = gridinfo[right].phia[a]*(1.0-hphi(gridinfo[right].phia,a))/(sqrt(gridinfo[right].phia[a]*(1.0-gridinfo[right].phia[a])));
        } else {
          s_phi_right = 0.0;
        }
                                          
        if (gridinfo[center].phia[a]*(1.0-gridinfo[center].phia[a]) > 0.0) {
          s_phi_center = gridinfo[center].phia[a]*(1.0-hphi(gridinfo[center].phia,a))/(sqrt(gridinfo[center].phia[a]*(1.0-gridinfo[center].phia[a])));
        } else {
          s_phi_center = 0.0;
        }
        
        if (gridinfo[front].phia[a]*(1.0-gridinfo[front].phia[a]) > 0.0) {
          s_phi_front = gridinfo[front].phia[a]*(1.0-hphi(gridinfo[front].phia,a))/(sqrt(gridinfo[front].phia[a]*(1.0-gridinfo[front].phia[a])));
        } else {
          s_phi_front = 0.0;
        }
      }
      
      for (k=0; k < NUMCOMPONENTS-1; k++) {
// 					      grad->jat[Y][k] =  0.0;
	if (normgradphi != 0.0) {
	  if(gidy!=0 && gidy!=(MESH_Y-2)) {
            if (ISOTHERMAL) {
              cl_right  = c_mu(gridinfo[right].compi,  T, NUMPHASES-1, k);
              cl_center = c_mu(gridinfo[center].compi, T, NUMPHASES-1, k);
      
              cs_right  = c_mu(gridinfo[right].compi,  T, a, k);
              cs_center = c_mu(gridinfo[center].compi, T, a, k);
            } else {
              cl_right  = c_mu(gridinfo[right].compi,  gridinfo[right].temperature,  NUMPHASES-1, k);
              cl_center = c_mu(gridinfo[center].compi, gridinfo[center].temperature, NUMPHASES-1, k);
      
              cs_right  = c_mu(gridinfo[right].compi,  gridinfo[right].temperature,  a, k);
              cs_center = c_mu(gridinfo[center].compi, gridinfo[center].temperature, a, k);
            }
	    
	    //Computing the projection of the phase contour along the liquid contour
	    scalprod  = grad->gradphi[Y][a]*grad->gradphi[Y][NUMPHASES-1];
	    scalprod += gradphix*gradphix_l;
	    
	    if (normgradphiy_l > 0.0) {
	      scalprod /= (normgradphiy_l*normgradphi);
	    }
	    if (WELL) { 
// 	    grad->jat[Y][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*(0.5/sqrt(2.0))*((cl_right-cs_right)*grad_right->deltaphi[a] + (cl_center-cs_center)*grad->deltaphi[a])*grad->gradphi[Y][a]/normgradphi;
              grad->jat[Y][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*(0.5/sqrt(2.0))*((cl_right-cs_right)*gridinfo[right].deltaphi[a] + (cl_center-cs_center)*gridinfo[center].deltaphi[a])*grad->gradphi[Y][a]/normgradphi;
            }
            if (OBSTACLE) {
// 	    grad->jat[Y][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*0.5*(s_phi_right*(cl_right-cs_right)*grad_right->deltaphi[a] + s_phi_center*(cl_center-cs_center)*grad->deltaphi[a])*grad->gradphi[Y][a]/normgradphi;
              grad->jat[Y][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*0.5*(s_phi_right*(cl_right-cs_right)*gridinfo[right].deltaphi[a] + s_phi_center*(cl_center-cs_center)*gridinfo[center].deltaphi[a])*grad->gradphi[Y][a]/normgradphi;
            }
	  }
	}
      }
      gradphiy    = 0.5*(grad->gradphi_c[Y][a] + grad_front->gradphi_c[Y][a]);
      normgradphi = sqrt(gradphiy*gradphiy + grad->gradphi[X][a]*grad->gradphi[X][a]);
      
      for (k=0; k < NUMCOMPONENTS-1; k++) {
// 					      grad->jat[X][k] = 0.0;
	if (normgradphi != 0.0) {
          if (ISOTHERMAL) {
            cl_front  = c_mu(gridinfo[front].compi,  T, NUMPHASES-1, k);
            cl_center = c_mu(gridinfo[center].compi, T, NUMPHASES-1, k);
            
            cs_front  = c_mu(gridinfo[front].compi,  T, a, k);
            cs_center = c_mu(gridinfo[center].compi, T, a, k);
          } else {
            cl_front  = c_mu(gridinfo[front].compi,  gridinfo[front].temperature,  NUMPHASES-1, k);
            cl_center = c_mu(gridinfo[center].compi, gridinfo[center].temperature, NUMPHASES-1, k);
            
            cs_front  = c_mu(gridinfo[front].compi,  gridinfo[front].temperature,  a, k);
            cs_center = c_mu(gridinfo[center].compi, gridinfo[center].temperature, a, k);
          }
	  
	  //Computing the projection of the phase contour along the liquid contour
	  scalprod  = grad->gradphi[X][a]*grad->gradphi[X][NUMPHASES-1];
	  scalprod += gradphiy*gradphiy_l;
	  
	  if (normgradphix_l > 0.0) {
	    scalprod /= (normgradphix_l*normgradphi);
          }
          
          if (WELL) {
// 	  grad->jat[X][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*(0.5/sqrt(2.0))*((cl_front-cs_front)*grad_front->deltaphi[a] + (cl_center-cs_center)*grad->deltaphi[a])*(grad->gradphi[X][a])/normgradphi;
            grad->jat[X][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*(0.5/sqrt(2.0))*((cl_front-cs_front)*gridinfo[front].deltaphi[a] + (cl_center-cs_center)*gridinfo[center].deltaphi[a])*(grad->gradphi[X][a])/normgradphi;
          }
          if (OBSTACLE) {
// 	  grad->jat[X][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*0.5*(s_phi_front*(cl_front-cs_front)*grad_front->deltaphi[a] + s_phi_center*(cl_center-cs_center)*grad->deltaphi[a])*(grad->gradphi[X][a])/normgradphi;
//           if (grad_front->deltaphi[a] != 0.0) {
//             printf("%le %ld %ld\n", grad_front->deltaphi[a], front, a);
//           }
            grad->jat[X][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*0.5*(s_phi_front*(cl_front-cs_front)*gridinfo[front].deltaphi[a] + s_phi_center*(cl_center-cs_center)*gridinfo[center].deltaphi[a])*(grad->gradphi[X][a])/normgradphi;
//            if (gridinfo[front].deltaphi[a] != 0.0) {
//              printf("%le %ld %ld\n", gridinfo[front].deltaphi[a], front, a);
//            }
          }
	}
      }
    }
    //End of multi-phase antitrapping current
  }
}
#endif