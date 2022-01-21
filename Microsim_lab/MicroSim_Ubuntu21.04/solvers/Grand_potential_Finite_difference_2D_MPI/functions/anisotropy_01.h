#ifndef ANISOTROPY_01_H_
#define ANISOTROPY_01_H_

// #include<stdio.h>
// #include<math.h>
// #include<stdlib.h>


void anisotropy_01_dAdq (double *qab, double* dadq, long a, long b) {
  double qx3 = qab[X]*qab[X]*qab[X];
  double qy3 = qab[Y]*qab[Y]*qab[Y];
  double qz3 = qab[Z]*qab[Z]*qab[Z];
  
  double q2  = qab[X]*qab[X] + qab[Y]*qab[Y] + qab[Z]*qab[Z];
  double q22 = q2*q2;
  double q23 = q2*q2*q2;
  double q4  = qx3*qab[X] + qy3*qab[Y] + qz3*qab[Z];
  
  if (fabs(q2) > 1.0e-15) {
   dadq[X] = 16.0*dab[a][b]*(qx3/q22 - qab[X]*q4/q23);
   dadq[Y] = 16.0*dab[a][b]*(qy3/q22 - qab[Y]*q4/q23);
   dadq[Z] = 16.0*dab[a][b]*(qz3/q22 - qab[Z]*q4/q23);
  } else {
    dadq[X] = 0.0;
    dadq[Y] = 0.0;
    dadq[Z] = 0.0;
  }
}

double anisotropy_01_function_ac(double *qab, long a, long b) {
  double qx2 = qab[X]*qab[X];
  double qx4 = qx2*qx2;
  
  double qy2 = qab[Y]*qab[Y];
  double qy4 = qy2*qy2;
  
  double qz2 = qab[Z]*qab[Z];
  double qz4 = qz2*qz2;
  
  double q2  = qx2 + qy2 + qz2;
  double ac;
  if (fabs(q2) > 1.0e-15) {
   ac = 1.0 - dab[a][b]*(3.0 - 4.0*(qx4 + qy4 + qz4)/(q2*q2));
  } else {
   ac = 1.0;
  }
  return ac;
}

// #ifdef FOLD2
// void dAdq (double *qab, double* dadq, long a, long b) {
//   //double qx2 = qab[X]*qab[X];
//   //double qy2 = qab[Y]*qab[Y];
//   double q2  = qab[X]*qab[X] + qab[Y]*qab[Y];
//   double q22 = q2*q2;
//   double q4  = qab[X]*qab[X] - qab[Y]*qab[Y];
//   //double qx3 = qab[X]*qab[X]*qab[X];
//   //double qy3 = qab[Y]*qab[Y]*qab[Y];
//   //double q4  = qx3*qab[X]    + qy3*qab[Y];
//   if (fabs(q2) > 1.0e-15) {
//    dadq[X] = 2.0*dab[a][b]*(qab[X]/q2 - (qab[X]*q4/q22));
//    dadq[Y] = -2.0*dab[a][b]*(qab[Y]/q2 + (qab[Y]*q4/q22));
//   } else {
//     dadq[X] = 0.0;
//     dadq[Y] = 0.0;
//   }
// }
// 
// double function_ac(double *qab, long a, long b) {
//   double qx2 = qab[X]*qab[X];
//   double qy2 = qab[Y]*qab[Y];
//   double q2  = qx2 + qy2;
//   double ac;
//   
//   if (fabs(q2) > 1.0e-15) {
//    ac = 1.0 + dab[a][b]*((qx2-qy2)/q2);
//   } else {
//    ac = 1.0;
//   }
//   return ac;
// }
// #endif
// 
// #ifdef COMPLEX
// void dAdq (double *qab, double* dadq, long a, long b) {
//   double qx3 = qab[X]*qab[X]*qab[X];
//   double qy3 = qab[Y]*qab[Y]*qab[Y];
//   double q2  = qab[X]*qab[X] + qab[Y]*qab[Y];
//   double q22 = q2*q2;
//   double q23 = q2*q2*q2;
//   double q4  = qx3*qab[X]    + qy3*qab[Y];
//   double qm  = qab[X]*qab[X] - qab[Y]*qab[Y];
//   double qd  = atan(qab[Y]/qab[X]);
//   double qw  = qd/w;
//   double term2 = exp(-qw*qw);
//   double term1 = term2*(-(2.0*qd)/(w*w));
//   
//   if (fabs(q2) > 1.0e-15) {
//    dadq[X] = - ec[a][b]*term1*(-qab[Y]/q2) - 2.0*e2[a][b]*qab[X]/q2 + 2*e2[a][b]*qab[X]*qm/q22 - 16.0*e4[a][b]*qx3/q22 + 16.0*e4[a][b]*qab[X]*q4/q23;
//    dadq[Y] = - ec[a][b]*term1*(qab[X]/q2)  + 2.0*e2[a][b]*qab[Y]/q2 + 2*e2[a][b]*qab[Y]*qm/q22 - 16.0*e4[a][b]*qy3/q22 + 16.0*e4[a][b]*qab[Y]*q4/q23;
//   } else {
//     dadq[X] = 0.0;
//     dadq[Y] = 0.0;
//   }
//   
//   //double stiffness;
//   //if(t==timebreak) {
//     //stiffness = 1.0 - ec[a][b]*term2 - (2*ec[a][b]/(w*w))*term2*(1-2*qw*qw)+3*e2[a][b]*qm/q2+15*e4[a][b]*(4*q4/q22-3);
//     //printf("%lf\t",stiffness);
//   //}
// }
// 
// double function_ac(double *qab, long a, long b) {
//   double qx2 = qab[X]*qab[X];
//   double qx4 = qx2*qx2;
//   double qy2 = qab[Y]*qab[Y];
//   double qy4 = qy2*qy2;
//   double q2  = qx2+qy2;
//   double q22 = q2*q2;
//   double x   = atan(qab[Y]/qab[X]);
//   double pw  = x/w;
//   double ac;
//   if (fabs(q2) > 1.0e-15) {
//    ac = 1.0 - ec[a][b] * exp(-pw*pw) - e2[a][b] * ((qx2-qy2)/q2) - e4[a][b] * (4.0 *((qx4+qy4)/q22) - 3.0);
//   } else {
//    ac = 1.0;
//   }
//   return ac;
// }
// #endif
// 
// #ifdef MIXED
// void dAdq (double *qab, double* dadq, long a, long b) {
//   double qx3 = qab[X]*qab[X]*qab[X];
//   double qy3 = qab[Y]*qab[Y]*qab[Y];
//   double q2  = qab[X]*qab[X] + qab[Y]*qab[Y];
//   double q22 = q2*q2;
//   double q23 = q2*q2*q2;
//   double q4  = qx3*qab[X]    + qy3*qab[Y];
//   double q4_2  = qab[X]*qab[X] - qab[Y]*qab[Y];
//   if (fabs(q2) > 1.0e-15) {
//    dadq[X] = 16.0*dab[a][b]*FAC_FOLD_4*(qx3/q22 - qab[X]*q4/q23);
//    
//    dadq[X] += 2.0*dab[a][b]*FAC_FOLD_2*(qab[X]/q2 - (qab[X]*q4_2/q22));
//    
//    dadq[Y] = 16.0*dab[a][b]*FAC_FOLD_4*(qy3/q22 - qab[Y]*q4/q23);
//    
//    dadq[Y] += -2.0*dab[a][b]*FAC_FOLD_2*(qab[Y]/q2 + (qab[Y]*q4_2/q22));
//   } else {
//     dadq[X] = 0.0;
//     dadq[Y] = 0.0;
//   }
// }
// 
// double function_ac(double *qab, long a, long b) {
//   double qx2 = qab[X]*qab[X];
//   double qx4 = qx2*qx2;
//   double qy2 = qab[Y]*qab[Y];
//   double qy4 = qy2*qy2;
//   double q2  = qx2+qy2;
//   double ac;
//   if (fabs(q2) > 1.0e-15) {
//    ac = 1.0 - dab[a][b]*FAC_FOLD_4*(3.0 - 4.0*(qx4+qy4)/(q2*q2)) + FAC_FOLD_2*dab[a][b]*((qx2-qy2)/q2);
//   } else {
//    ac = 1.0;
//   }
//   return ac;
// }
// #endif

#endif
