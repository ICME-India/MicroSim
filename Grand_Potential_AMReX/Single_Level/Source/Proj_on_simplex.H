#ifndef PROJ_ON_SIMPLEX
#define PROJ_ON_SIMPLEX

#include<AMReX_Utility.H>
using namespace amrex;

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void projection_on_simplex(int i, int j, int k, amrex::Array4<amrex::Real> const& phiold, amrex::Array4<amrex::Real> const& phinew, Array1D<Real,0,phasecount-1> &deltaphi, Array1D<Real,0,phasecount-1> &div, int numph) {
  long a, b;
  Real Deltaphi = 0.0;
  Real count_phases = 0;
  double Deltaphi_alpha = 0.0;
  //Find the number of phases which are okay.
  double sum_phib = 0.0;
  for (a=0; a < numph; a++) {
    
    if ((fabs(div(a)) > 0.0) && (phiold(i,j,k,a)+deltaphi(a)) > 0.0  && (phiold(i,j,k,a)+deltaphi(a)) < 1.0) {
      count_phases++;
      sum_phib += phiold(i,j,k,a);
    }
  }
  
  for (a=0; a < numph; a++) {
    if ((phiold(i,j,k,a)+deltaphi(a)) < 0.0) {
      Deltaphi_alpha = fabs(phiold(i,j,k,a)+deltaphi(a));
      deltaphi(a)   += Deltaphi_alpha;
      Deltaphi      += Deltaphi_alpha;
    }
  }
  if (Deltaphi > 0.0) {
    for (b=0; b < numph; b++) {
      
      if (((fabs(div(b)) > 0.0)) && ((phiold(i,j,k,b)+deltaphi(b)) > 0.0)  && ((phiold(i,j,k,b)+deltaphi(b)) < 1.0)) {
        if (fabs(sum_phib) > 0.0) {
          deltaphi(b) -= Deltaphi*(phiold(i,j,k,b))/(sum_phib);
        } else {
          deltaphi(b) -= Deltaphi/(count_phases);
        }
      }
    }
  }
  for (a=0; a < numph; a++) {

    if (fabs(div(a)) > 0.0) {
      if ((phiold(i,j,k,a)+deltaphi(a)) > 1.0) {
        deltaphi(a)  = (1.0-phiold(i,j,k,a));
        //Correct all the other phases due to this correction,
        //If you are bulk, all other phases must go to zero.
        for (b=0; b < numph; b++) {
            
          if ((fabs(div(b)) > 0.0) && b!=a) {
            deltaphi(b) = -phiold(i,j,k,b);
          }
        }
        Deltaphi = 0.0;
        break;
      }
    }
  }
}



#endif