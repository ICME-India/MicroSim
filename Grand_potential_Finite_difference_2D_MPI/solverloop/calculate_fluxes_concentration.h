#ifndef CALCULATE_FLUXES_CONCENTRATION_H_
#define CALCULATE_FLUXES_CONCENTRATION_H_

// void calculate_fluxes_concentration_2D(long x, struct gradlayer **gradient) {
//  double cl_front, cl_center, cs_front, cs_center, cl_right, cs_right;
//  double DET;
//  for (gidy=0; gidy <=(workers_mpi.end[Y]+2); gidy++) {
//     
//     grad          =  &gradient[0][gidy];
//     grad_right    =  &gradient[0][gidy+1];
//     grad_front    =  &gradient[1][gidy];
//     
//     center        =  gidy   + (x)*workers_mpi.layer_size;
//     front         =  gidy   + (x+1)*workers_mpi.layer_size;
//     right         =  center + 1;
//     left          =  center - 1;
//     
//     if(!ISOTHERMAL) {
//       //gridinfo_w[center].temperature  = temp_bottom + gidy*GRADIENT;
//       T = gridinfo_w[center].temperature;
// //       init_propertymatrices(T);
//     }
//     
//     //For multi-phase antitrapping. The last phase is by default the liquid phase.
//     //Generating the normal vector corresponding to the liquid contour
//     gradphix_l      = 0.5*(grad->gradphi_c[X][NUMPHASES-1] + grad_right->gradphi_c[X][NUMPHASES-1]);
//     normgradphiy_l  = sqrt(gradphix_l*gradphix_l + grad->gradphi[Y][NUMPHASES-1]*grad->gradphi[Y][NUMPHASES-1]);
//     
//     gradphiy_l      = 0.5*(grad->gradphi_c[Y][NUMPHASES-1] + grad_front->gradphi_c[Y][NUMPHASES-1]);
//     normgradphix_l  = sqrt(gradphiy_l*gradphiy_l + grad->gradphi[X][NUMPHASES-1]*grad->gradphi[X][NUMPHASES-1]);
//     
//     //Initializing the jats to zero
//     
//     for (k=0; k < NUMCOMPONENTS-1; k++) {
//       grad->jat[X][k] = 0.0;
//       grad->jat[Y][k] = 0.0;
//     }
//     for (a=0; a < NUMPHASES-1; a++) {
//     
//       gradphix      = 0.5*(grad->gradphi_c[X][a] + grad_right->gradphi_c[X][a]);
//       normgradphi   = sqrt(gradphix*gradphix + grad->gradphi[Y][a]*grad->gradphi[Y][a]);
//       
//       if (OBSTACLE) {		  
//         if (gridinfo_w[right].phia[a]*(1.0-gridinfo_w[right].phia[a]) > 0.0) {
//           s_phi_right = gridinfo_w[right].phia[a]*(1.0-hphi(gridinfo_w[right].phia,a))/(sqrt(gridinfo_w[right].phia[a]*(1.0-gridinfo_w[right].phia[a])));
//         } else {
//           s_phi_right = 0.0;
//         }
//                                           
//         if (gridinfo_w[center].phia[a]*(1.0-gridinfo_w[center].phia[a]) > 0.0) {
//           s_phi_center = gridinfo_w[center].phia[a]*(1.0-hphi(gridinfo_w[center].phia,a))/(sqrt(gridinfo_w[center].phia[a]*(1.0-gridinfo_w[center].phia[a])));
//         } else {
//           s_phi_center = 0.0;
//         }
//         
//         if (gridinfo_w[front].phia[a]*(1.0-gridinfo_w[front].phia[a]) > 0.0) {
//           s_phi_front = gridinfo_w[front].phia[a]*(1.0-hphi(gridinfo_w[front].phia,a))/(sqrt(gridinfo_w[front].phia[a]*(1.0-gridinfo_w[front].phia[a])));
//         } else {
//           s_phi_front = 0.0;
//         }
//       }
//       
//       for (k=0; k < NUMCOMPONENTS-1; k++) {
// // 					      grad->jat[Y][k] =  0.0;
// 	if (normgradphi != 0.0) {
// 	  if(gidy!=0 && gidy!=(workers_mpi.end[Y]+2)) {
//             if (ISOTHERMAL) {
//               cl_right  = c_mu(gridinfo_w[right].compi,  T, NUMPHASES-1, k);
//               cl_center = c_mu(gridinfo_w[center].compi, T, NUMPHASES-1, k);
//       
//               cs_right  = c_mu(gridinfo_w[right].compi,  T, a, k);
//               cs_center = c_mu(gridinfo_w[center].compi, T, a, k);
//             } else {
//               cl_right  = c_mu(gridinfo_w[right].compi,  gridinfo_w[right].temperature,  NUMPHASES-1, k);
//               cl_center = c_mu(gridinfo_w[center].compi, gridinfo_w[center].temperature, NUMPHASES-1, k);
//       
//               cs_right  = c_mu(gridinfo_w[right].compi,  gridinfo_w[right].temperature,  a, k);
//               cs_center = c_mu(gridinfo_w[center].compi, gridinfo_w[center].temperature, a, k);
//             }
// 	    
// 	    //Computing the projection of the phase contour along the liquid contour
// 	    scalprod  = grad->gradphi[Y][a]*grad->gradphi[Y][NUMPHASES-1];
// 	    scalprod += gradphix*gradphix_l;
// 	    
// 	    if (normgradphiy_l > 0.0) {
// 	      scalprod /= (normgradphiy_l*normgradphi);
// 	    }
// 	    if (WELL) { 
// // 	    grad->jat[Y][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*(0.5/sqrt(2.0))*((cl_right-cs_right)*grad_right->deltaphi[a] + (cl_center-cs_center)*grad->deltaphi[a])*grad->gradphi[Y][a]/normgradphi;
//               grad->jat[Y][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*(0.5/sqrt(2.0))*((cl_right-cs_right)*gridinfo_w[right].deltaphi[a] + (cl_center-cs_center)*gridinfo_w[center].deltaphi[a])*grad->gradphi[Y][a]/normgradphi;
//             }
//             if (OBSTACLE) {
// // 	    grad->jat[Y][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*0.5*(s_phi_right*(cl_right-cs_right)*grad_right->deltaphi[a] + s_phi_center*(cl_center-cs_center)*grad->deltaphi[a])*grad->gradphi[Y][a]/normgradphi;
//               grad->jat[Y][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*0.5*(s_phi_right*(cl_right-cs_right)*gridinfo_w[right].deltaphi[a] + s_phi_center*(cl_center-cs_center)*gridinfo_w[center].deltaphi[a])*grad->gradphi[Y][a]/normgradphi;
//             }
// 	  }
// 	}
//       }
//       gradphiy    = 0.5*(grad->gradphi_c[Y][a] + grad_front->gradphi_c[Y][a]);
//       normgradphi = sqrt(gradphiy*gradphiy + grad->gradphi[X][a]*grad->gradphi[X][a]);
//       
//       for (k=0; k < NUMCOMPONENTS-1; k++) {
// // 					      grad->jat[X][k] = 0.0;
// 	if (normgradphi != 0.0) {
//           if (ISOTHERMAL) {
//             cl_front  = c_mu(gridinfo_w[front].compi,  T, NUMPHASES-1, k);
//             cl_center = c_mu(gridinfo_w[center].compi, T, NUMPHASES-1, k);
//             
//             cs_front  = c_mu(gridinfo_w[front].compi,  T, a, k);
//             cs_center = c_mu(gridinfo_w[center].compi, T, a, k);
//           } else {
//             cl_front  = c_mu(gridinfo_w[front].compi,  gridinfo_w[front].temperature,  NUMPHASES-1, k);
//             cl_center = c_mu(gridinfo_w[center].compi, gridinfo_w[center].temperature, NUMPHASES-1, k);
//             
//             cs_front  = c_mu(gridinfo_w[front].compi,  gridinfo_w[front].temperature,  a, k);
//             cs_center = c_mu(gridinfo_w[center].compi, gridinfo_w[center].temperature, a, k);
//           }
// 	  
// 	  //Computing the projection of the phase contour along the liquid contour
// 	  scalprod  = grad->gradphi[X][a]*grad->gradphi[X][NUMPHASES-1];
// 	  scalprod += gradphiy*gradphiy_l;
// 	  
// 	  if (normgradphix_l > 0.0) {
// 	    scalprod /= (normgradphix_l*normgradphi);
//           }
//           
//           if (WELL) {
// // 	  grad->jat[X][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*(0.5/sqrt(2.0))*((cl_front-cs_front)*grad_front->deltaphi[a] + (cl_center-cs_center)*grad->deltaphi[a])*(grad->gradphi[X][a])/normgradphi;
//             grad->jat[X][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*(0.5/sqrt(2.0))*((cl_front-cs_front)*gridinfo_w[front].deltaphi[a] + (cl_center-cs_center)*gridinfo_w[center].deltaphi[a])*(grad->gradphi[X][a])/normgradphi;
//           }
//           if (OBSTACLE) {
// // 	  grad->jat[X][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*0.5*(s_phi_front*(cl_front-cs_front)*grad_front->deltaphi[a] + s_phi_center*(cl_center-cs_center)*grad->deltaphi[a])*(grad->gradphi[X][a])/normgradphi;
// //           if (grad_front->deltaphi[a] != 0.0) {
// //             printf("%le %ld %ld\n", grad_front->deltaphi[a], front, a);
// //           }
//             grad->jat[X][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*0.5*(s_phi_front*(cl_front-cs_front)*gridinfo_w[front].deltaphi[a] + s_phi_center*(cl_center-cs_center)*gridinfo_w[center].deltaphi[a])*(grad->gradphi[X][a])/normgradphi;
// //            if (gridinfo_w[front].deltaphi[a] != 0.0) {
// //              printf("%le %ld %ld\n", gridinfo_w[front].deltaphi[a], front, a);
// //            }
//           }
// 	}
//       }
//     }
//     //End of multi-phase antitrapping current
//   }
// }
void calculate_fluxes_concentration_2D(long x, struct gradlayer **gradient) {
 double cl_front, cl_center, cs_front, cs_center, cl_right, cs_right;
 double DET;
 for (gidy=0; gidy <=(workers_mpi.end[Y]+2); gidy++) {
    
    grad          =  &gradient[0][gidy];
    grad_right    =  &gradient[0][gidy+1];
    grad_front    =  &gradient[1][gidy];
    
    center        =  gidy   + (x)*workers_mpi.layer_size;
    front         =  gidy   + (x+1)*workers_mpi.layer_size;
    right         =  center + 1;
    left          =  center - 1;
    
    if(!ISOTHERMAL) {
      //gridinfo_w[center].temperature  = temp_bottom + gidy*GRADIENT;
      T = gridinfo_w[center].temperature;
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
        if (gridinfo_w[right].phia[a]*(1.0-gridinfo_w[right].phia[a]) > 0.0) {
          s_phi_right = gridinfo_w[right].phia[a]*(1.0-hphi(gridinfo_w[right].phia,a))/(sqrt(gridinfo_w[right].phia[a]*(1.0-gridinfo_w[right].phia[a])));
        } else {
          s_phi_right = 0.0;
        }
                                          
        if (gridinfo_w[center].phia[a]*(1.0-gridinfo_w[center].phia[a]) > 0.0) {
          s_phi_center = gridinfo_w[center].phia[a]*(1.0-hphi(gridinfo_w[center].phia,a))/(sqrt(gridinfo_w[center].phia[a]*(1.0-gridinfo_w[center].phia[a])));
        } else {
          s_phi_center = 0.0;
        }
        
        if (gridinfo_w[front].phia[a]*(1.0-gridinfo_w[front].phia[a]) > 0.0) {
          s_phi_front = gridinfo_w[front].phia[a]*(1.0-hphi(gridinfo_w[front].phia,a))/(sqrt(gridinfo_w[front].phia[a]*(1.0-gridinfo_w[front].phia[a])));
        } else {
          s_phi_front = 0.0;
        }
      }
      
      for (k=0; k < NUMCOMPONENTS-1; k++) {
// 					      grad->jat[Y][k] =  0.0;
        if (normgradphi != 0.0) {
          if(gidy!=0 && gidy!=(workers_mpi.end[Y]+2)) {
            cl_right  = grad_right->phase_comp[NUMPHASES-1][k];
            cl_center = grad->phase_comp[NUMPHASES-1][k];

            cs_right  = grad_right->phase_comp[a][k];
            cs_center = grad->phase_comp[a][k];
      //             if (ISOTHERMAL) {
      //               cl_right  = c_mu(gridinfo_w[right].compi,  T, NUMPHASES-1, k);
      //               cl_center = c_mu(gridinfo_w[center].compi, T, NUMPHASES-1, k);
      //       
      //               cs_right  = c_mu(gridinfo_w[right].compi,  T, a, k);
      //               cs_center = c_mu(gridinfo_w[center].compi, T, a, k);
      //             } else {
      //               cl_right  = c_mu(gridinfo_w[right].compi,  gridinfo_w[right].temperature,  NUMPHASES-1, k);
      //               cl_center = c_mu(gridinfo_w[center].compi, gridinfo_w[center].temperature, NUMPHASES-1, k);
      //       
      //               cs_right  = c_mu(gridinfo_w[right].compi,  gridinfo_w[right].temperature,  a, k);
      //               cs_center = c_mu(gridinfo_w[center].compi, gridinfo_w[center].temperature, a, k);
      //             }
            
            //Computing the projection of the phase contour along the liquid contour
            scalprod  = grad->gradphi[Y][a]*grad->gradphi[Y][NUMPHASES-1];
            scalprod += gradphix*gradphix_l;
            
            if (normgradphiy_l > 0.0) {
              scalprod /= (normgradphiy_l*normgradphi);
            }
            if (WELL) { 
      // 	    grad->jat[Y][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*(0.5/sqrt(2.0))*((cl_right-cs_right)*grad_right->deltaphi[a] + (cl_center-cs_center)*grad->deltaphi[a])*grad->gradphi[Y][a]/normgradphi;
              grad->jat[Y][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*(0.5/sqrt(2.0))*((cl_right-cs_right)*gridinfo_w[right].deltaphi[a] + (cl_center-cs_center)*gridinfo_w[center].deltaphi[a])*grad->gradphi[Y][a]/normgradphi;
            }
            if (OBSTACLE) {
      // 	    grad->jat[Y][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*0.5*(s_phi_right*(cl_right-cs_right)*grad_right->deltaphi[a] + s_phi_center*(cl_center-cs_center)*grad->deltaphi[a])*grad->gradphi[Y][a]/normgradphi;
              grad->jat[Y][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*0.5*(s_phi_right*(cl_right-cs_right)*gridinfo_w[right].deltaphi[a] + s_phi_center*(cl_center-cs_center)*gridinfo_w[center].deltaphi[a])*grad->gradphi[Y][a]/normgradphi;
            }
          }
        }
      }
      gradphiy    = 0.5*(grad->gradphi_c[Y][a] + grad_front->gradphi_c[Y][a]);
      normgradphi = sqrt(gradphiy*gradphiy + grad->gradphi[X][a]*grad->gradphi[X][a]);
      
      for (k=0; k < NUMCOMPONENTS-1; k++) {
// 					      grad->jat[X][k] = 0.0;
        if (normgradphi != 0.0) {
          cl_front  = grad_front->phase_comp[NUMPHASES-1][k];
          cl_center = grad->phase_comp[NUMPHASES-1][k];
          
          cs_front  = grad_front->phase_comp[a][k];
          cs_center = grad->phase_comp[a][k];
      //           if (ISOTHERMAL) {
      //             cl_front  = c_mu(gridinfo_w[front].compi,  T, NUMPHASES-1, k);
      //             cl_center = c_mu(gridinfo_w[center].compi, T, NUMPHASES-1, k);
      //             
      //             cs_front  = c_mu(gridinfo_w[front].compi,  T, a, k);
      //             cs_center = c_mu(gridinfo_w[center].compi, T, a, k);
      //           } else {
      //             cl_front  = c_mu(gridinfo_w[front].compi,  gridinfo_w[front].temperature,  NUMPHASES-1, k);
      //             cl_center = c_mu(gridinfo_w[center].compi, gridinfo_w[center].temperature, NUMPHASES-1, k);
      //             
      //             cs_front  = c_mu(gridinfo_w[front].compi,  gridinfo_w[front].temperature,  a, k);
      //             cs_center = c_mu(gridinfo_w[center].compi, gridinfo_w[center].temperature, a, k);
      //           }
          
          //Computing the projection of the phase contour along the liquid contour
          scalprod  = grad->gradphi[X][a]*grad->gradphi[X][NUMPHASES-1];
          scalprod += gradphiy*gradphiy_l;
          
          if (normgradphix_l > 0.0) {
            scalprod /= (normgradphix_l*normgradphi);
          }
                
          if (WELL) {
// 	  grad->jat[X][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*(0.5/sqrt(2.0))*((cl_front-cs_front)*grad_front->deltaphi[a] + (cl_center-cs_center)*grad->deltaphi[a])*(grad->gradphi[X][a])/normgradphi;
            grad->jat[X][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*(0.5/sqrt(2.0))*((cl_front-cs_front)*gridinfo_w[front].deltaphi[a] + (cl_center-cs_center)*gridinfo_w[center].deltaphi[a])*(grad->gradphi[X][a])/normgradphi;
          }
          if (OBSTACLE) {
// 	  grad->jat[X][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*0.5*(s_phi_front*(cl_front-cs_front)*grad_front->deltaphi[a] + s_phi_center*(cl_center-cs_center)*grad->deltaphi[a])*(grad->gradphi[X][a])/normgradphi;
//           if (grad_front->deltaphi[a] != 0.0) {
//             printf("%le %ld %ld\n", grad_front->deltaphi[a], front, a);
//           }
            grad->jat[X][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*0.5*(s_phi_front*(cl_front-cs_front)*gridinfo_w[front].deltaphi[a] + s_phi_center*(cl_center-cs_center)*gridinfo_w[center].deltaphi[a])*(grad->gradphi[X][a])/normgradphi;
//            if (gridinfo_w[front].deltaphi[a] != 0.0) {
//              printf("%le %ld %ld\n", gridinfo_w[front].deltaphi[a], front, a);
//            }
          }
        }
      }
    }
    //End of multi-phase antitrapping current
  }
}
void calculate_fluxes_concentration_3D(long x, struct gradlayer **gradient) {
 double cl_front, cl_center, cs_front, cs_center, cl_right, cs_right, cl_top, cl_bottom, cs_top, cs_bottom;
 double DET;
 for (gidy=0; gidy < (workers_mpi.layer_size); gidy++) {
    
    grad          =  &gradient[0][gidy];
    grad_front    =  &gradient[1][gidy];
    
    if ((gidy%workers_mpi.rows_y) < (workers_mpi.rows_y-1)) {
      grad_right   =  &gradient[0][gidy+1];
    }
    
    if ((gidy+workers_mpi.rows_y) < workers_mpi.layer_size) {
      grad_top     =  &gradient[0][gidy+workers_mpi.rows_y];
    }
    
    
    center        =  gidy   + (x)*workers_mpi.layer_size;
    front         =  gidy   + (x+1)*workers_mpi.layer_size;
    right         =  center + 1;
    left          =  center - 1;
    top           =  center + workers_mpi.rows_y;
    bottom        =  center - workers_mpi.rows_y;
    
    if(!ISOTHERMAL) {
      T = gridinfo_w[center].temperature;
    }
    
    //For multi-phase antitrapping. The last phase is by default the liquid phase.
    //Generating the normal vector corresponding to the liquid contour    
    gradphiy_l_x      = 0.5*(grad->gradphi_c[Y][NUMPHASES-1] + grad_front->gradphi_c[Y][NUMPHASES-1]);
    gradphiz_l_x      = 0.5*(grad->gradphi_c[Z][NUMPHASES-1] + grad_front->gradphi_c[Z][NUMPHASES-1]);
    normgradphix_l    = sqrt(gradphiy_l_x*gradphiy_l_x + gradphiz_l_x*gradphiz_l_x + grad->gradphi[X][NUMPHASES-1]*grad->gradphi[X][NUMPHASES-1]);
    
    if ((((gidy)%workers_mpi.rows_y) < (workers_mpi.rows_y-1))) {
      gradphix_l_y      = 0.5*(grad->gradphi_c[X][NUMPHASES-1] + grad_right->gradphi_c[X][NUMPHASES-1]);
      gradphiz_l_y      = 0.5*(grad->gradphi_c[Z][NUMPHASES-1] + grad_right->gradphi_c[Z][NUMPHASES-1]);
      normgradphiy_l    = sqrt(gradphiz_l_y*gradphiz_l_y + gradphix_l_y*gradphix_l_y + grad->gradphi[Y][NUMPHASES-1]*grad->gradphi[Y][NUMPHASES-1]);
    } else {
      gradphix_l_y      = 0.0;
      gradphiz_l_y      = 0.0;
      normgradphiy_l    = sqrt(gradphiz_l_y*gradphiz_l_y + gradphix_l_y*gradphix_l_y + grad->gradphi[Y][NUMPHASES-1]*grad->gradphi[Y][NUMPHASES-1]);
    }
        
    if ((((gidy+workers_mpi.rows_y)/workers_mpi.rows_y) < (workers_mpi.rows_z))) {
      gradphix_l_z    = 0.5*(grad->gradphi_c[X][NUMPHASES-1] + grad_top->gradphi_c[X][NUMPHASES-1]);
      gradphiy_l_z    = 0.5*(grad->gradphi_c[Y][NUMPHASES-1] + grad_top->gradphi_c[Y][NUMPHASES-1]);
      normgradphiz_l  = sqrt(gradphix_l_z*gradphix_l_z + gradphiy_l_z*gradphiy_l_z + grad->gradphi[Z][NUMPHASES-1]*grad->gradphi[Z][NUMPHASES-1]);
    } else {
      gradphix_l_z    = 0.0;
      gradphiy_l_z    = 0.0;
      normgradphiz_l  = sqrt(gradphix_l_z*gradphix_l_z + gradphiy_l_z*gradphiy_l_z + grad->gradphi[Z][NUMPHASES-1]*grad->gradphi[Z][NUMPHASES-1]);
    }
    
    //Initializing the jats to zero
    
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      grad->jat[X][k] = 0.0;
      grad->jat[Y][k] = 0.0;
      grad->jat[Z][k] = 0.0;
    }
    for (a=0; a < NUMPHASES-1; a++) {
      if (((gidy)%workers_mpi.rows_y) < (workers_mpi.rows_y-1)) {
        gradphix      = 0.5*(grad->gradphi_c[X][a] + grad_right->gradphi_c[X][a]);
        gradphiz      = 0.5*(grad->gradphi_c[Z][a] + grad_right->gradphi_c[Z][a]);
      } else {
        gradphix      = 0.5*(grad->gradphi_c[X][a] + 0.0);
        gradphiz      = 0.5*(grad->gradphi_c[Z][a] + 0.0);
      }
      
      normgradphi   = sqrt(gradphix*gradphix + gradphiz*gradphiz + grad->gradphi[Y][a]*grad->gradphi[Y][a]);
      
      if (OBSTACLE) {  
        if ((gidy%workers_mpi.rows_y) < (workers_mpi.rows_y-1)) {
          if (gridinfo_w[right].phia[a]*(1.0-gridinfo_w[right].phia[a]) > 0.0) {
            s_phi_right = gridinfo_w[right].phia[a]*(1.0-hphi(gridinfo_w[right].phia,a))/(sqrt(gridinfo_w[right].phia[a]*(1.0-gridinfo_w[right].phia[a])));
          } else {
            s_phi_right = 0.0;
          }
        } else {
          s_phi_right = 0.0;
        }
                                          
        if (gridinfo_w[center].phia[a]*(1.0-gridinfo_w[center].phia[a]) > 0.0) {
          s_phi_center = gridinfo_w[center].phia[a]*(1.0-hphi(gridinfo_w[center].phia,a))/(sqrt(gridinfo_w[center].phia[a]*(1.0-gridinfo_w[center].phia[a])));
        } else {
          s_phi_center = 0.0;
        }
        
        if (gridinfo_w[front].phia[a]*(1.0-gridinfo_w[front].phia[a]) > 0.0) {
          s_phi_front = gridinfo_w[front].phia[a]*(1.0-hphi(gridinfo_w[front].phia,a))/(sqrt(gridinfo_w[front].phia[a]*(1.0-gridinfo_w[front].phia[a])));
        } else {
          s_phi_front = 0.0;
        }
        
        if ((((gidy+workers_mpi.rows_y)/workers_mpi.rows_y) < (workers_mpi.rows_z))) {
          if (gridinfo_w[top].phia[a]*(1.0-gridinfo_w[top].phia[a]) > 0.0) {
            s_phi_top = gridinfo_w[top].phia[a]*(1.0-hphi(gridinfo_w[top].phia,a))/(sqrt(gridinfo_w[top].phia[a]*(1.0-gridinfo_w[top].phia[a])));
          } else {
            s_phi_top = 0.0;
          }
        } else {
          s_phi_top = 0.0;
        }
      }
      
      for (k=0; k < NUMCOMPONENTS-1; k++) {
        if (normgradphi != 0.0) {
          if((gidy > 0) && (gidy < workers_mpi.layer_size)) {
            cl_right  = grad_right->phase_comp[NUMPHASES-1][k];
            cl_center = grad->phase_comp[NUMPHASES-1][k];
      
            cs_right  = grad_right->phase_comp[a][k];
            cs_center = grad->phase_comp[a][k];
          
	    //Computing the projection of the phase contour along the liquid contour
            scalprod  = grad->gradphi[Y][a]*grad->gradphi[Y][NUMPHASES-1];
            scalprod += gradphix*gradphix_l_y;
            scalprod += gradphiz*gradphiz_l_y;
            
            if (normgradphiy_l > 0.0) {
              scalprod /= (normgradphiy_l*normgradphi);
            }
            if (WELL) { 
              grad->jat[Y][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*(0.5/sqrt(2.0))*((cl_right-cs_right)*gridinfo_w[right].deltaphi[a] + (cl_center-cs_center)*gridinfo_w[center].deltaphi[a])*grad->gradphi[Y][a]/normgradphi;
            }
            if (OBSTACLE) {
              grad->jat[Y][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*0.5*(s_phi_right*(cl_right-cs_right)*gridinfo_w[right].deltaphi[a] + s_phi_center*(cl_center-cs_center)*gridinfo_w[center].deltaphi[a])*grad->gradphi[Y][a]/normgradphi;
            }
          }
        }
      }     
      gradphiy    = 0.5*(grad->gradphi_c[Y][a] + grad_front->gradphi_c[Y][a]);
      gradphiz    = 0.5*(grad->gradphi_c[Z][a] + grad_front->gradphi_c[Z][a]);
      normgradphi = sqrt(gradphiy*gradphiy + gradphiz*gradphiz + grad->gradphi[X][a]*grad->gradphi[X][a]);
            
      for (k=0; k < NUMCOMPONENTS-1; k++) {
        if (normgradphi != 0.0) {
          cl_front  = grad_front->phase_comp[NUMPHASES-1][k];
          cl_center = grad->phase_comp[NUMPHASES-1][k];
          
          cs_front  = grad_front->phase_comp[a][k];
          cs_center = grad->phase_comp[a][k];
          
          //Computing the projection of the phase contour along the liquid contour
          scalprod  = grad->gradphi[X][a]*grad->gradphi[X][NUMPHASES-1];
          scalprod += gradphiy*gradphiy_l_x;
          scalprod += gradphiz*gradphiz_l_x;
          
          if (normgradphix_l > 0.0) {
            scalprod /= (normgradphix_l*normgradphi);
          }
          
          if (WELL) {
            grad->jat[X][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*(0.5/sqrt(2.0))*((cl_front-cs_front)*gridinfo_w[front].deltaphi[a] + (cl_center-cs_center)*gridinfo_w[center].deltaphi[a])*(grad->gradphi[X][a])/normgradphi;
          }
          if (OBSTACLE) {
            grad->jat[X][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*0.5*(s_phi_front*(cl_front-cs_front)*gridinfo_w[front].deltaphi[a] + s_phi_center*(cl_center-cs_center)*gridinfo_w[center].deltaphi[a])*(grad->gradphi[X][a])/normgradphi;

          }
	      }  
      }
      
      if (((gidy+workers_mpi.rows_y)/workers_mpi.rows_y) < workers_mpi.rows_z) {
        gradphix    = 0.5*(grad->gradphi_c[X][a] + grad_top->gradphi_c[X][a]);
        gradphiy    = 0.5*(grad->gradphi_c[Y][a] + grad_top->gradphi_c[Y][a]);
      } else {
        gradphix    = 0.5*(grad->gradphi_c[X][a] + 0.0);
        gradphiy    = 0.5*(grad->gradphi_c[Y][a] + 0.0);
      }
      
      normgradphi = sqrt(gradphiy*gradphiy + gradphix*gradphix + grad->gradphi[Z][a]*grad->gradphi[Z][a]);
      
      for (k=0; k < NUMCOMPONENTS-1; k++) {
        if (normgradphi != 0.0) {
          cl_top    = grad_top->phase_comp[NUMPHASES-1][k];
          cl_center = grad->phase_comp[NUMPHASES-1][k];
          
          cs_top    = grad_top->phase_comp[a][k];
          cs_center = grad->phase_comp[a][k];
          
          //Computing the projection of the phase contour along the liquid contour
          scalprod  = grad->gradphi[Z][a]*grad->gradphi[Z][NUMPHASES-1];
          scalprod += gradphix*gradphix_l_z;
          scalprod += gradphiy*gradphiy_l_z;
          
          if (normgradphiz_l > 0.0) {
            scalprod /= (normgradphiz_l*normgradphi);
          }
          
          if (WELL) {
            grad->jat[Z][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*(0.5/sqrt(2.0))*((cl_top-cs_top)*gridinfo_w[top].deltaphi[a] + (cl_center-cs_center)*gridinfo_w[center].deltaphi[a])*(grad->gradphi[Z][a])/normgradphi;
          }
          if (OBSTACLE) {
            grad->jat[Z][k] += (1.0-Diffusivity[a][k][k]/Diffusivity[NUMPHASES-1][k][k])*fabs(scalprod)*0.5*(s_phi_top*(cl_top-cs_top)*gridinfo_w[top].deltaphi[a] + s_phi_center*(cl_center-cs_center)*gridinfo_w[center].deltaphi[a])*(grad->gradphi[Z][a])/normgradphi;
          }
	      }  
      }
    }
    //End of multi-phase antitrapping current
  }
}
#endif
