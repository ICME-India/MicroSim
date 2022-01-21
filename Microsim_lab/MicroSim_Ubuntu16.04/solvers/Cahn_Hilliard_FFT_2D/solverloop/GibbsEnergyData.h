#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#include "GibbsEnergyFunctions/GSol.h"
#include "GibbsEnergyFunctions/GSol.c"
#include "GibbsEnergyFunctions/GSol_m.h"
#include "GibbsEnergyFunctions/GSol_m.c"

double GSOL(double T, double X1S);
double dGSOLdX1S(double T, double X1S);
double GSOL_m(double T, double X1L); 
double dGSOL_mdX1S(double T, double X1L); 

double GSOL(double T, double X1S) {

  double xS[2];
  double GS[1];
  double gs;

  xS[0] = X1S;
  xS[1] = 1.0 - X1S; 

  GES(T, xS, GS);

  gs = GS[0];

  return gs;
}

double dGSOLdX1S(double T, double X1S) {

  double xS[2];
  double dGS[2];
  double dgsdxs;

  xS[0] = X1S;
  xS[1] = 1.0 - X1S; 

  dGES(T, xS, dGS);

  //dgsdxs = 2.0*A_fm[0]*X1S*(1.0-X1S)*(1.0-2.0*X1S);//dGS[0] - dGS[1];
  dgsdxs = dGS[0] - dGS[1];

  return dgsdxs;
}


double GSOL_m(double T, double X1L) {
  double xL[2];
  double GL[1];
  double gl;

  xL[0] = X1L;
  xL[1] = 1.0 - X1L; 

  GEL(T, xL, GL);

  gl = GL[0];

  return gl;
}

double dGSOL_mdX1S(double T, double X1L) {

  double xL[2];
  double dGL[2];
  double dgldxl;

  xL[0] = X1L;
  xL[1] = 1.0 - X1L; 

  dGEL(T, xL, dGL);

  dgldxl = dGL[0] - dGL[1];

  return dgldxl;
}