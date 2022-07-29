#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#include "Thermo3.c"
#include "Thermo3.h"

void ges(double T, double *y, double *Ge);
void mus(double T, double *y, double *mus);
void dmudcs(double T, double *y, double *dmu);

void ges(double T, double *y, double *Ge) {
  int ip, is, is1, is2;
  double tmp0;
  double cies[npha][nsol+1], retG[npha][1]; 
  
  for ( ip = 0; ip < npha; ip++ ) { 
    for ( is = 0; is < nsol+1; is++ ) { 
      cies[ip][is] = y[ip*(nsol+1)+is]; 
    }
  }

  GE(T, cies[0], cies[1], retG[0], retG[1]);
  
  for ( ip = 0; ip < npha; ip++ ) { 
    Ge[ip] = retG[ip][0];
  }

}

void mus(double T, double *y, double *Mu) {
  int ip, is, is1, is2;
  double tmp0;
  double cies[npha][nsol+1], retmu[npha][nsol]; 
  
  for ( ip = 0; ip < npha; ip++ ) { 
    for ( is = 0; is < nsol+1; is++ ) { 
      cies[ip][is] = y[ip*(nsol+1)+is]; 
    }
  }
  
  mu(T, cies[0], cies[1], retmu[0], retmu[1]);
  
  for ( ip = 0; ip < npha; ip++ ) { 
    for ( is = 0; is < nsol; is++ ) { 
      Mu[ip*nsol+is] = retmu[ip][is];
    }
  }

}

void dmudcs(double T, double *y, double *dmu) {
  int ip, is, is1, is2;
  double tmp0;
  double cies[npha][nsol+1], retdmu[npha][nsol*nsol]; 
  
  for ( ip = 0; ip < npha; ip++ ) { 
    for ( is = 0; is < nsol+1; is++ ) { 
      cies[ip][is] = y[ip*(nsol+1)+is]; 
    }
  }
  
  dmudc(T, cies[0], cies[1], retdmu[0], retdmu[1]);
  
  for ( ip = 0; ip < npha; ip++ ) { 
    for ( is1 = 0; is1 < nsol; is1++ ) { 
      for ( is2 = 0; is2 < nsol; is2++ ) { 
        dmu[(ip*nsol+is1)*nsol+is2] = retdmu[ip][is1*nsol+is2];
      }
    }
  }

}

