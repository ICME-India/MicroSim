#include "functionW_01.cuh"

/*
 * __device__ double calcDoubleObstacleDerivative
 *
 * Calculate g'(\phi)
 *
 * Arguments:
 *              1. cufftDoubleComplex **phi -> all the phase volume-fraction values
 *              2. long phase -> differentiate wrt to this phase
 *              3. double *theta_i   -> coefficients for theta_i, one per phase
 *              4. double *theta_ij  -> coefficients for theta_ij, one per pair of phases
 *              5. double *theta_ijk -> coefficients for theta_ijk, one per triplet of phases
 *              6. long idx -> position of cell in 1D
 *              7. long NUMPHASES -> number of phases
 * Return:
 *              numerical evaluation of derivative of interpolation polynomial, as a double datatype
 */
extern __device__
double calcDoubleObstacleDerivative(double **phi, long phase,
                                    double *theta_i, double *theta_ij, double *theta_ijk,
                                    long idx, long NUMPHASES)
{
    if (NUMPHASES < 2)
        return 0.0;

  long b,c;
  double sum=0.0;
//   double phibphic;
//   //The normal obstacle potential
//   for (b=0; b < NUMPHASES; b++) {
//     if (b!=phase) {
//       if (fabs(divphi[b]) > 0.0) {
//         if ((phi[phase]*phi[b])>=0.00000) {
//           sum += Gamma[phase][b]*(phi[b]);
//         } else {
//           sum += -Gamma[phase][b]*(phi[b]);
//         }
//       }
//     }
//   }
//   sum *= 16.0/(M_PI*M_PI);
//   //The third order potential
//   for (b=0; b < NUMPHASES; b++) {
//     for (c=0; c < NUMPHASES; c++) {
//       if (b!=a && c!=a && b < c) {
//         if (fabs(divphi[b]) > 0.0 && fabs(divphi[c]) > 0.0) {
//           phibphic = phi[b]*phi[c];
// //           if (phi[a]*phibphic >= 0.00000) {
//             sum += Gamma_abc[a][b][c]*phibphic;
// //           } else {
// //             sum += -Gamma_abc[a][b][c]*phibphic;
// //           }
//         }
//       }
//     }
//   }
  return sum;
}
