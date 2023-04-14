#ifndef FUNCTION_F_ELAST_H_
#define FUNCTION_F_ELAST_H_
#include "global_vars.h"
double df_elast_2D(struct gradlayer *grad, struct symmetric_tensor sigma, double *phi, long a) {
  double delast;
  long b;
  double sum=0.0;
  struct symmetric_tensor sigma_phase;
  
//   sigma.xx  = grad->stiffness_c[X].C11*(grad->strain[X].xx)      + grad->stiffness_c[X].C12*(grad->strain[X].yy);
//   sigma.yy  = grad->stiffness_c[Y].C12*(grad->strain[Y].xx)      + grad->stiffness_c[Y].C11*(grad->strain[Y].yy);
//   sigma.xy  = 2.0*grad->stiffness_c[X].C44*(grad->strain[X].xy);
  
//   for(b=0; b < NUMPHASES; b++) {
//     delast         = -(sigma.xx*eigen_strain_phase[b].xx + sigma.yy*eigen_strain_phase[b].yy + 2.0*sigma.xy*eigen_strain_phase[b].xy);
//     sigma_phase.xx = stiffness_phase[b].C11*grad->strain[X].xx + stiffness_phase[b].C12*grad->strain[X].yy;
//     sigma_phase.yy = stiffness_phase[b].C12*grad->strain[X].xx + stiffness_phase[b].C11*grad->strain[X].yy;
//     sigma_phase.xy = 2.0*stiffness_phase[b].C44*grad->strain[X].xy;
//     delast        += 0.5*(sigma_phase.xx*grad->strain[X].xx + sigma_phase.yy*grad->strain[X].yy + 2.0*sigma_phase.xy*grad->strain[X].xy);
//     delast        *= dhphi(phi, b, a);
//     sum           += delast;
    
    delast         = -(sigma.xx*eigen_strain_phase[a].xx + sigma.yy*eigen_strain_phase[a].yy + 2.0*sigma.xy*eigen_strain_phase[a].xy);
    sigma_phase.xx =   stiffness_phase[a].C11*grad->strain[X].xx + stiffness_phase[a].C12*grad->strain[X].yy;
    sigma_phase.yy =   stiffness_phase[a].C12*grad->strain[X].xx + stiffness_phase[a].C11*grad->strain[X].yy;
    sigma_phase.xy =   2.0*stiffness_phase[a].C44*grad->strain[X].xy;
    delast        +=   0.5*(sigma_phase.xx*grad->strain[X].xx + sigma_phase.yy*grad->strain[X].yy + 2.0*sigma_phase.xy*grad->strain[X].xy);
//     delast        *= dhphi(phi, b, a);
    sum            =   delast;
    
//   }
  return sum;
}
double df_elast_3D(struct gradlayer *gradient, struct symmetric_tensor sigma, double *phi, long a) {
  double delast;
  long b;
  double sum=0.0;
  struct symmetric_tensor sigma_phase;
//   for(b=0; b < NUMPHASES; b++) {
    delast         = -(sigma.xx*eigen_strain_phase[a].xx + sigma.yy*eigen_strain_phase[a].yy + sigma.zz*eigen_strain_phase[a].zz + 2.0*sigma.xy*eigen_strain_phase[a].xy + 2.0*sigma.xz*eigen_strain_phase[a].xz + 2.0*sigma.yz*eigen_strain_phase[a].yz);
    
    sigma_phase.xx = stiffness_phase[a].C11*grad->strain[X].xx  + stiffness_phase[a].C12*(grad->strain[X].yy + grad->strain[X].zz);
    sigma_phase.yy = stiffness_phase[a].C12*(grad->strain[X].xx + grad->strain[X].zz) + stiffness_phase[a].C11*grad->strain[X].yy;
    sigma_phase.zz = stiffness_phase[a].C12*(grad->strain[X].xx + grad->strain[X].yy) + stiffness_phase[a].C11*grad->strain[X].zz;
    
    
    sigma_phase.xy = 2.0*stiffness_phase[a].C44*grad->strain[X].xy;
    sigma_phase.xz = 2.0*stiffness_phase[a].C44*grad->strain[X].xz;
    sigma_phase.yz = 2.0*stiffness_phase[a].C44*grad->strain[X].yz;
    
    delast        += 0.5*(sigma_phase.xx*grad->strain[X].xx + sigma_phase.yy*grad->strain[X].yy + sigma_phase.zz*grad->strain[X].zz + 2.0*sigma_phase.xy*grad->strain[X].xy + 2.0*sigma_phase.xz*grad->strain[X].xz + 2.0*sigma_phase.yz*grad->strain[X].yz);
//     delast        *= dhphi(phi, b, a);
    sum           = delast;
//   }
  return sum;
}
#endif
