#ifndef _FUNCTION_Q_H_
#define _FUNCTION_Q_H_

// 0 = center 
// 1 = i+1/2
// 2 = i-1/2
// 3 = j+1/2
// 4 = j-1/2 

#include <AMReX_Utility.H>
using namespace amrex;

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Function_Q(int i, int j, int k, amrex::Array4<Real const> const& phi, GpuArray<Real,AMREX_SPACEDIM> delta, int a, int b, Array2D<Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1, Order::C> &qab){

    qab(cent,X) = phi(i,j,k,a)*(phi(i+1,j,k,b)-phi(i-1,j,k,b))/(2*delta[0]) -phi(i,j,k,b)*(phi(i+1,j,k,a)-phi(i-1,j,k,a))/(2*delta[0]);
	qab(iph,X) = 0.5*(phi(i,j,k,a)+phi(i+1,j,k,a))*(phi(i+1,j,k,b)-phi(i,j,k,b))/(delta[0]) -0.5*(phi(i,j,k,b)+phi(i+1,j,k,b))*(phi(i+1,j,k,a)-phi(i,j,k,a))/(delta[0]);
	qab(imh,X) = 0.5*(phi(i,j,k,a)+phi(i-1,j,k,a))*(phi(i,j,k,b)-phi(i-1,j,k,b))/(delta[0]) -0.5*(phi(i,j,k,b)+phi(i-1,j,k,b))*(phi(i,j,k,a)-phi(i-1,j,k,a))/(delta[0]);
	qab(jph,X) = 0.5*(phi(i,j+1,k,a)+phi(i,j,k,a))*(-phi(i-1,j+1,k,b)+phi(i+1,j+1,k,b)-phi(i-1,j,k,b)+phi(i+1,j,k,b))/(4*delta[0]) -0.5*(phi(i,j+1,k,b)+phi(i,j,k,b))*(-phi(i-1,j+1,k,a)+phi(i+1,j+1,k,a)-phi(i-1,j,k,a)+phi(i+1,j,k,a))/(4*delta[0]);
	qab(jmh,X) = 0.5*(phi(i,j,k,a)+phi(i,j-1,k,a))*(-phi(i-1,j,k,b)+phi(i+1,j,k,b)-phi(i-1,j-1,k,b)+phi(i+1,j-1,k,b))/(4*delta[0]) -0.5*(phi(i,j,k,b)+phi(i,j-1,k,b))*(-phi(i-1,j,k,a)+phi(i+1,j,k,a)-phi(i-1,j-1,k,a)+phi(i+1,j-1,k,a))/(4*delta[0]);
#if(AMREX_SPACEDIM>2)
	qab(kph,X) = 0.5*(phi(i,j,k+1,a)+phi(i,j,k,a))*(-phi(i-1,j,k+1,b)+phi(i+1,j,k+1,b)-phi(i-1,j,k,b)+phi(i+1,j,k,b))/(4*delta[0]) -0.5*(phi(i,j,k+1,b)+phi(i,j,k,b))*(-phi(i-1,j,k+1,a)+phi(i+1,j,k+1,a)-phi(i-1,j,k,a)+phi(i+1,j,k,a))/(4*delta[0]);
	qab(kmh,X) = 0.5*(phi(i,j,k,a)+phi(i,j,k-1,a))*(-phi(i-1,j,k,b)+phi(i+1,j,k,b)-phi(i-1,j,k-1,b)+phi(i+1,j,k-1,b))/(4*delta[0]) -0.5*(phi(i,j,k,b)+phi(i,j,k-1,b))*(-phi(i-1,j,k,a)+phi(i+1,j,k,a)-phi(i-1,j,k-1,a)+phi(i+1,j,k-1,a))/(4*delta[0]);
#endif

	qab(cent,Y) = phi(i,j,k,a)*(phi(i,j+1,k,b)-phi(i,j-1,k,b))/(2*delta[1]) -phi(i,j,k,b)*(phi(i,j+1,k,a)-phi(i,j-1,k,a))/(2*delta[1]);
	qab(iph,Y) = 0.5*(phi(i,j,k,a)+phi(i+1,j,k,a))*(phi(i+1,j+1,k,b)-phi(i+1,j-1,k,b)+phi(i,j+1,k,b)-phi(i,j-1,k,b))/(4*delta[1]) -0.5*(phi(i,j,k,b)+phi(i+1,j,k,b))*(phi(i+1,j+1,k,a)-phi(i+1,j-1,k,a)+phi(i,j+1,k,a)-phi(i,j-1,k,a))/(4*delta[1]);
	qab(imh,Y) = 0.5*(phi(i,j,k,a)+phi(i-1,j,k,a))*(phi(i,j+1,k,b)-phi(i,j-1,k,b)+phi(i-1,j+1,k,b)-phi(i-1,j-1,k,b))/(4*delta[1]) -0.5*(phi(i,j,k,b)+phi(i-1,j,k,b))*(phi(i,j+1,k,a)-phi(i,j-1,k,a)+phi(i-1,j+1,k,a)-phi(i-1,j-1,k,a))/(4*delta[1]);
	qab(jph,Y) = 0.5*(phi(i,j+1,k,a)+phi(i,j,k,a))*(phi(i,j+1,k,b)-phi(i,j,k,b))/(delta[1]) -0.5*(phi(i,j+1,k,b)+phi(i,j,k,b))*(phi(i,j+1,k,a)-phi(i,j,k,a))/(delta[1]);
	qab(jmh,Y) = 0.5*(phi(i,j,k,a)+phi(i,j-1,k,a))*(phi(i,j,k,b)-phi(i,j-1,k,b))/(delta[1]) -0.5*(phi(i,j,k,b)+phi(i,j-1,k,b))*(phi(i,j,k,a)-phi(i,j-1,k,a))/(delta[1]);
#if(AMREX_SPACEDIM>2)
	qab(kph,Y) = 0.5*(phi(i,j,k+1,a)+phi(i,j,k,a))*(-phi(i,j-1,k+1,b)+phi(i,j+1,k+1,b)-phi(i,j-1,k,b)+phi(i,j+1,k,b))/(4*delta[1]) -0.5*(phi(i,j,k+1,b)+phi(i,j,k,b))*(-phi(i,j-1,k+1,a)+phi(i,j+1,k+1,a)-phi(i,j-1,k,a)+phi(i,j+1,k,a))/(4*delta[1]);
	qab(kmh,Y) = 0.5*(phi(i,j,k,a)+phi(i,j,k-1,a))*(-phi(i,j-1,k,b)+phi(i,j+1,k,b)-phi(i,j-1,k-1,b)+phi(i,j+1,k-1,b))/(4*delta[1]) -0.5*(phi(i,j,k,b)+phi(i,j,k-1,b))*(-phi(i,j-1,k,a)+phi(i,j+1,k,a)-phi(i,j-1,k-1,a)+phi(i,j+1,k-1,a))/(4*delta[1]);
#endif

#if (AMREX_SPACEDIM > 2)	

	qab(cent,Z) = phi(i,j,k,a)*(phi(i,j,k+1,b)-phi(i,j,k-1,b))/(2*delta[2]) -phi(i,j,k,b)*(phi(i,j,k+1,a)-phi(i,j,k-1,a))/(2*delta[2]);
	qab(iph,Z) = 0.5*(phi(i,j,k,a)+phi(i+1,j,k,a))*(phi(i+1,j,k+1,b)-phi(i+1,j,k-1,b)+phi(i,j,k+1,b)-phi(i,j,k-1,b))/(4*delta[2]) -0.5*(phi(i,j,k,b)+phi(i+1,j,k,b))*(phi(i+1,j,k+1,a)-phi(i+1,j,k-1,a)+phi(i,j,k+1,a)-phi(i,j,k-1,a))/(4*delta[2]);
	qab(imh,Z) = 0.5*(phi(i,j,k,a)+phi(i-1,j,k,a))*(phi(i,j,k+1,b)-phi(i,j,k-1,b)+phi(i-1,j,k+1,b)-phi(i-1,j,k-1,b))/(4*delta[2]) -0.5*(phi(i,j,k,b)+phi(i-1,j,k,b))*(phi(i,j,k+1,a)-phi(i,j,k-1,a)+phi(i-1,j,k+1,a)-phi(i-1,j,k-1,a))/(4*delta[2]);
	qab(jph,Z) = 0.5*(phi(i,j,k,a)+phi(i,j+1,k,a))*(phi(i,j+1,k+1,b)-phi(i,j+1,k-1,b)+phi(i,j,k+1,b)-phi(i,j,k-1,b))/(4*delta[2]) -0.5*(phi(i,j,k,b)+phi(i,j+1,k,b))*(phi(i,j+1,k+1,a)-phi(i,j+1,k-1,a)+phi(i,j,k+1,a)-phi(i,j,k-1,a))/(4*delta[2]);
	qab(jmh,Z) = 0.5*(phi(i,j,k,a)+phi(i,j-1,k,a))*(phi(i,j,k+1,b)-phi(i,j,k-1,b)+phi(i,j-1,k+1,b)-phi(i,j-1,k-1,b))/(4*delta[2]) -0.5*(phi(i,j,k,b)+phi(i,j-1,k,b))*(phi(i,j,k+1,a)-phi(i,j,k-1,a)+phi(i,j-1,k+1,a)-phi(i,j-1,k-1,a))/(4*delta[2]);
	qab(kph,Z) = 0.5*(phi(i,j,k,a)+phi(i,j,k+1,a))*(phi(i,j,k+1,b)-phi(i,j,k,b))/(delta[2]) -0.5*(phi(i,j,k,b)+phi(i,j,k+1,b))*(phi(i,j,k+1,a)-phi(i,j,k,a))/(delta[2]);
	qab(kmh,Z) = 0.5*(phi(i,j,k,a)+phi(i,j,k-1,a))*(phi(i,j,k,b)-phi(i,j,k-1,b))/(delta[2]) -0.5*(phi(i,j,k,b)+phi(i,j,k-1,b))*(phi(i,j,k,a)-phi(i,j,k-1,a))/(delta[2]);

#endif

}


#endif