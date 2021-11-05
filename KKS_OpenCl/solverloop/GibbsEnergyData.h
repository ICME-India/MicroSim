#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#include "GibbsEnergyFunctions/GLiq.h"
#include "GibbsEnergyFunctions/GSol.h"
#include "GibbsEnergyFunctions/GLiq.c"
#include "GibbsEnergyFunctions/GSol.c"

double GLIQ(double T, double X1L);
double GSOL(double T, double X1S);

double dGLIQdX1L(double T, double X1L);
double dGSOLdX1S(double T, double X1S);

double ddGLIQdX1LdX1L(double T, double X1L);
double ddGSOLdX1SdX1S(double T, double X1S);
 
double GLIQ(double T, double X1L) {
  double xL[2];
  double GL[1];
  double gl;

  xL[0] = X1L;
  xL[1] = 1.0 - X1L; 

  GEL(T, xL, GL);

  gl = GL[0];

  return gl;
}

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

double dGLIQdX1L(double T, double X1L) {

  double xL[2];
  double dGL[2];
  double dgldxl;

  xL[0] = X1L;
  xL[1] = 1.0 - X1L; 

  dGEL(T, xL, dGL);

  dgldxl = dGL[0] - dGL[1];

  return dgldxl;
}

double dGSOLdX1S(double T, double X1S) {

  double xS[2];
  double dGS[2];
  double dgsdxs;

  xS[0] = X1S;
  xS[1] = 1.0 - X1S; 

  dGES(T, xS, dGS);

  dgsdxs = dGS[0] - dGS[1];

  return dgsdxs;
}

double ddGLIQdX1LdX1L(double T, double X1L) {

  double xL[2];
  double ddGL[2];
  double ddgldxldxl;

  xL[0] = X1L;
  xL[1] = 1.0 - X1L; 

  ddGEL(T, xL, ddGL);

  ddgldxldxl = ddGL[0] - ddGL[1];

  return ddgldxldxl;
}

double ddGSOLdX1SdX1S(double T, double X1S) {

  double xS[2];
  double ddGS[2];
  double ddgsdxsdxs;

  xS[0] = X1S;
  xS[1] = 1.0 - X1S; 

  ddGES(T, xS, ddGS);

  ddgsdxsdxs = ddGS[0] - ddGS[1];

  return ddgsdxsdxs;
}
