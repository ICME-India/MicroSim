#ifndef FUNCTION_W_01_H_
#define FUNCTION_W_01_H_


// #ifdef ISOTROPY
double function_W_01_dwdphi(double *phi, double *divphi, struct gradlayer **gradient, long gidy, long a) {
  long b,c;
  double sum=0.0;
  double phibphic;
  //The normal obstacle potential
  for (b=0; b < NUMPHASES; b++) {
    if (b!=a) {
      if (fabs(divphi[b]) > 0.0) {
        if ((phi[a]*phi[b])>=0.00000) {
          sum += Gamma[a][b]*(phi[b]);
        } else {
          sum += -Gamma[a][b]*(phi[b]);
        }
      }
    }
  }
  sum *= 16.0/(M_PI*M_PI);
  //The third order potential
  for (b=0; b < NUMPHASES; b++) {
    for (c=0; c < NUMPHASES; c++) {
      if (b!=a && c!=a && b < c) {
        if (fabs(divphi[b]) > 0.0 && fabs(divphi[c]) > 0.0) {
          phibphic = phi[b]*phi[c];
          if (phi[a]*phibphic >= 0.00000) {
            sum += Gamma_abc[a][b][c]*phibphic;
          } else {
            sum += -Gamma_abc[a][b][c]*phibphic;
          }
        }
      }
    }
  }
  return sum;
}

double function_W_01_dwdphi_smooth(double *phi, double *divphi, struct gradlayer **gradient, long gidy, long a) {
  long b,c;
  double sum=0.0;
  //The normal obstacle potential
  for (b=0; b < NUMPHASES; b++) {
    if (b!=a) {
      if (fabs(divphi[b]) > 0.0) {
        sum += Gamma[a][b]*phi[b];
      }
    }
  }
  sum *= 16.0/(M_PI*M_PI);
  //The third order potential
  for (b=0; b < NUMPHASES; b++) {
    for (c=0; c < NUMPHASES; c++) {
      if (b!=a && c!=a && b < c) {
        if (fabs(divphi[b]) > 0.0 && fabs(divphi[c]) > 0.0) {
          sum += Gamma_abc[a][b][c]*phi[b]*phi[c];
        }
      }
    }
  }
  return sum;
}


// #endif



// #ifdef ANISOTROPY
// double function_W_01_dwdphi(double *phi, double *divphi, struct gradlayer **gradient, long gidy, long a) {
// //double dwdphi(double *phi, double *divphi, long a) {//double dwdphi(double *phi, double *divphi, long a) {
//   long b,c;
//   double sum=0.0, sum1=0.0;
//   double divergence =0.0;
//   
//   long dir;
//   double dadq[2];
//   double qab[2];
//   double ac1, ac2, qab2;
//   
//   double dadq_front[2], dadq_back[2], dadq_left[2], dadq_right[2];
//   double qab_front[2],  qab_back[2],  qab_left[2],  qab_right[2];
//   double rotated_dadq_front[2], rotated_dadq_back[2], rotated_dadq_left[2], rotated_dadq_right[2];
//   
//   double ac_front,  ac_back,   ac_right,  ac_left;
//   double ac_front2, ac_back2,  ac_right2, ac_left2;
//   double qab_front2,qab_back2, qab_right2,qab_left2;
//   
//   double dqdphi[2], Rotated_dqdphi[2];
//   
//   struct gradlayer *grad1;
//   struct gradlayer *grad1_left;
//   struct gradlayer *grad1_back;
//   struct gradlayer *grad1_right;
//   struct gradlayer *grad1_front;
//   
//   grad1       =  &gradient[1][gidy];
//   grad1_left  =  &gradient[1][gidy-1];
//   grad1_back  =  &gradient[0][gidy];
//   grad1_right =  &gradient[1][gidy+1];
//   grad1_front =  &gradient[2][gidy];
//   
//   double phibphic;
// 
// #ifdef OBSTACLE
//   //The normal obstacle potential
//   for (b=0; b < NUMPHASES; b++) {
//     if (b!=a) {
//       
// #ifdef ANISOTROPY_GRADIENT
// //     if (fabs(divphi[b]) > 0.0) {
// //       sum += Gamma[a][b]*phi[b];
// //     }
//   
//   //The normal obstacle potential
//   if (fabs(divphi[b]) > 0.0) {
//     if ((phi[a]*phi[b]) >= 0.00000) {
//       sum +=  Gamma[a][b]*(phi[b]);
//     } else {
//       sum += -Gamma[a][b]*(phi[b]);
//     }
//   }
// #endif
// #ifdef ANISOTROPY_POTENTIAL
//     if (fabs(divphi[b]) > 0.0) {
//       q_dadphi(phi, grad1, a,  b , qab);
//       dAdq(qab, dadq, a, b);
//       ac1  = function_ac(qab, a, b);
//       ac2  = ac1*ac1;
//       
//       sum1  = phi[b]*ac2;
//       
//       dqdphi[X] = grad1->gradphi_c[X][b];
//       dqdphi[Y] = grad1->gradphi_c[Y][b];
//       
//       multiply(Rotation_matrix[a][b], dqdphi, Rotated_dqdphi, DIMENSION);
//       
//       for(dir=0; dir<DIMENSION; dir++) {
//         dqdphi[dir] = Rotated_dqdphi[dir];
//       }
//       
//       for (dir=0; dir < 2; dir++)  {
//         sum1 += 2.0*ac1*dadq[dir]*dqdphi[dir]*phi[a]*phi[b];
//       }
//       sum1 *= Gamma[a][b];
//       
//       q_divx (grad1, grad1_front, a,  b, qab_front);
//       ac_front = function_ac(qab_front, a, b);
//       dAdq(qab_front, dadq_front, a, b);
//       rotate_vector(dadq_front, rotated_dadq_front, a, b);
//       
//       q_divx (grad1_back, grad1, a,  b, qab_back);
//       ac_back   = function_ac(qab_back, a, b);
//       dAdq(qab_back, dadq_back, a, b);
//       rotate_vector(dadq_back, rotated_dadq_back, a, b);
//       
//       q_divy (grad1, grad1_right, a,  b, qab_right);
//       ac_right   = function_ac(qab_right, a, b);
//       dAdq(qab_right, dadq_right, a, b);
//       rotate_vector(dadq_right, rotated_dadq_right, a, b);
//      
//       q_divy (grad1_left, grad1, a,  b, qab_left);
//       ac_left = function_ac(qab_left, a, b);
//       dAdq(qab_left, dadq_left, a, b);
//       rotate_vector(dadq_left, rotated_dadq_left, a, b);
//       
//       divergence  =  2.0*(ac_front*rotated_dadq_front[X]*grad1->phistagg[X][b]*grad1->phistagg[X][a]*grad1->phistagg[X][b] - ac_back*rotated_dadq_back[X]*(grad1_back->phistagg[X][b])*grad1_back->phistagg[X][a]*grad1_back->phistagg[X][b])/deltax;
//       
//       divergence +=  2.0*(ac_right*rotated_dadq_right[Y]*(grad1->phistagg[Y][b])*grad1->phistagg[Y][a]*grad1->phistagg[Y][b] - ac_left*rotated_dadq_left[Y]*grad1_left->phistagg[Y][b]*grad1_left->phistagg[Y][a]*grad1_left->phistagg[Y][b])/deltay;
//       
//       divergence *= Gamma[a][b];
//       
//       sum += (divergence + sum1);
//     }
// #endif
// #ifdef ANISOTROPY_GRADIENT_POTENTIAL
//      if (fabs(divphi[b]) > 0.0) {
//       q_dadphi(phi, grad1, a,  b , qab);
//       dAdq(qab, dadq, a, b);
//       ac1  = function_ac(qab, a, b);
//       
//       sum1  = phi[b]*ac1;
//       
//       dqdphi[X] = grad1->gradphi_c[X][b];
//       dqdphi[Y] = grad1->gradphi_c[Y][b];
//       
//       multiply(Rotation_matrix[a][b], dqdphi, Rotated_dqdphi, DIMENSION);
//       
//       for(dir=0; dir<DIMENSION; dir++) {
//         dqdphi[dir] = Rotated_dqdphi[dir];
//       }
// 
//       for (dir=0; dir < 2; dir++)  {
//         sum1 += dadq[dir]*dqdphi[dir]*phi[a]*phi[b];
//       }
//       
//       sum1 *= Gamma[a][b];
//       
//       q_divx (grad1, grad1_front, a,  b, qab_front);
//       dAdq(qab_front, dadq_front, a, b);
//       rotate_vector(dadq_front, rotated_dadq_front, a, b);
//       
//       q_divx (grad1_back, grad1, a,  b, qab_back);
//       dAdq(qab_back, dadq_back, a, b);
//       rotate_vector(dadq_back, rotated_dadq_back, a, b);
//       
//       q_divy (grad1, grad1_right, a,  b, qab_right);
//       dAdq(qab_right, dadq_right, a, b);
//       rotate_vector(dadq_right, rotated_dadq_right,a , b);
//      
//       q_divy (grad1_left, grad1, a,  b, qab_left);
//       dAdq(qab_left, dadq_left, a, b);
//       rotate_vector(dadq_left, rotated_dadq_left, a, b);
//       
//       divergence  =  (rotated_dadq_front[X]*(grad1->phistagg[X][b])*grad1->phistagg[X][a]*grad1->phistagg[X][b] - rotated_dadq_back[X]*(grad1_back->phistagg[X][b])*grad1_back->phistagg[X][a]*grad1_back->phistagg[X][b])/deltax;
//       
//       divergence +=  (rotated_dadq_right[Y]*(grad1->phistagg[Y][b])*grad1->phistagg[Y][a]*grad1->phistagg[Y][b] - rotated_dadq_left[Y]*(grad1_left->phistagg[Y][b])*grad1_left->phistagg[Y][a]*grad1_left->phistagg[Y][b])/deltay;
//       
//       divergence *=   Gamma[a][b];
//       
//       sum += (divergence + sum1);
//      }
// #endif
//     }
//   }
//   sum *= 16.0/(M_PI*M_PI);
// #endif
// #ifdef WELL
//   //The well potential
//   for (b=0; b < NUMPHASES; b++) {
//     if (b!=a) {
//       if (fabs(divphi[b]) > 0.0) {
//         sum += 2.0*Gamma[a][b]*phi[b]*phi[b]*phi[a];
//       }
//     }
//   }
//   sum *= 9.0;
// #endif
//   //The third order potential
//   for (b=0; b < NUMPHASES; b++) {
//     for (c=0; c < NUMPHASES; c++) {
//       if (b!=a && c!=a && b < c) {
//         if (fabs(divphi[b]) > 0.0 && fabs(divphi[c]) > 0.0) {
//           phibphic = phi[b]*phi[c];
//           if (phi[a]*phibphic >= 0.00000) {
//             sum += Gamma_abc[a][b][c]*phibphic;
//           } else {
//             sum += -Gamma_abc[a][b][c]*phibphic;
//           }
//         }
//       }
//     }
//   }
//   return sum;
// }
// #endif


#endif