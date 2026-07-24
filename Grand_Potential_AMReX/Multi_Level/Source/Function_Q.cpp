#include <AMReX_Utility.H>
#include "AmrCoreGP.H"
using namespace amrex;
//####################################################################################################################
void 
AmrCoreGP::Function_Q(int i, int j, int k, amrex::Array4<Real const> const& phi, 
			GpuArray<Real,AMREX_SPACEDIM> delta, int a, int b, Array2D<Real,0,
			AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1, Order::C> &qab)
{

    qab(cent,X) = phi(i,j,k,a)*(phi(i+1,j,k,b)-phi(i-1,j,k,b))/(2*delta[X]) 
					-phi(i,j,k,b)*(phi(i+1,j,k,a)-phi(i-1,j,k,a))/(2*delta[X]);
	
	qab(iph,X) = 0.5*(phi(i,j,k,a)+phi(i+1,j,k,a))*(phi(i+1,j,k,b)-phi(i,j,k,b))/(delta[X]) 
					-0.5*(phi(i,j,k,b)+phi(i+1,j,k,b))*(phi(i+1,j,k,a)-phi(i,j,k,a))/(delta[X]);
	
	qab(imh,X) = 0.5*(phi(i,j,k,a)+phi(i-1,j,k,a))*(phi(i,j,k,b)-phi(i-1,j,k,b))/(delta[X]) 
					-0.5*(phi(i,j,k,b)+phi(i-1,j,k,b))*(phi(i,j,k,a)-phi(i-1,j,k,a))/(delta[X]);
	
	qab(jph,X) = 0.5*(phi(i,j+1,k,a)+phi(i,j,k,a))*(-phi(i-1,j+1,k,b)+phi(i+1,j+1,k,b)-phi(i-1,j,k,b)+phi(i+1,j,k,b))/(4*delta[X]) 
					-0.5*(phi(i,j+1,k,b)+phi(i,j,k,b))*(-phi(i-1,j+1,k,a)+phi(i+1,j+1,k,a)-phi(i-1,j,k,a)+phi(i+1,j,k,a))/(4*delta[X]);
	
	qab(jmh,X) = 0.5*(phi(i,j,k,a)+phi(i,j-1,k,a))*(-phi(i-1,j,k,b)+phi(i+1,j,k,b)-phi(i-1,j-1,k,b)+phi(i+1,j-1,k,b))/(4*delta[X]) 
					-0.5*(phi(i,j,k,b)+phi(i,j-1,k,b))*(-phi(i-1,j,k,a)+phi(i+1,j,k,a)-phi(i-1,j-1,k,a)+phi(i+1,j-1,k,a))/(4*delta[X]);
#if(AMREX_SPACEDIM>2)
	qab(kph,X) = 0.5*(phi(i,j,k+1,a)+phi(i,j,k,a))*(-phi(i-1,j,k+1,b)+phi(i+1,j,k+1,b)-phi(i-1,j,k,b)+phi(i+1,j,k,b))/(4*delta[X]) 
					-0.5*(phi(i,j,k+1,b)+phi(i,j,k,b))*(-phi(i-1,j,k+1,a)+phi(i+1,j,k+1,a)-phi(i-1,j,k,a)+phi(i+1,j,k,a))/(4*delta[X]);
	
	qab(kmh,X) = 0.5*(phi(i,j,k,a)+phi(i,j,k-1,a))*(-phi(i-1,j,k,b)+phi(i+1,j,k,b)-phi(i-1,j,k-1,b)+phi(i+1,j,k-1,b))/(4*delta[X]) 
					-0.5*(phi(i,j,k,b)+phi(i,j,k-1,b))*(-phi(i-1,j,k,a)+phi(i+1,j,k,a)-phi(i-1,j,k-1,a)+phi(i+1,j,k-1,a))/(4*delta[X]);
#endif

	qab(cent,Y) = phi(i,j,k,a)*(phi(i,j+1,k,b)-phi(i,j-1,k,b))/(2*delta[Y]) 
					-phi(i,j,k,b)*(phi(i,j+1,k,a)-phi(i,j-1,k,a))/(2*delta[Y]);
	
	qab(iph,Y) = 0.5*(phi(i,j,k,a)+phi(i+1,j,k,a))*(phi(i+1,j+1,k,b)-phi(i+1,j-1,k,b)+phi(i,j+1,k,b)-phi(i,j-1,k,b))/(4*delta[Y]) 
					-0.5*(phi(i,j,k,b)+phi(i+1,j,k,b))*(phi(i+1,j+1,k,a)-phi(i+1,j-1,k,a)+phi(i,j+1,k,a)-phi(i,j-1,k,a))/(4*delta[Y]);
	
	qab(imh,Y) = 0.5*(phi(i,j,k,a)+phi(i-1,j,k,a))*(phi(i,j+1,k,b)-phi(i,j-1,k,b)+phi(i-1,j+1,k,b)-phi(i-1,j-1,k,b))/(4*delta[Y]) 
					-0.5*(phi(i,j,k,b)+phi(i-1,j,k,b))*(phi(i,j+1,k,a)-phi(i,j-1,k,a)+phi(i-1,j+1,k,a)-phi(i-1,j-1,k,a))/(4*delta[Y]);
	
	qab(jph,Y) = 0.5*(phi(i,j+1,k,a)+phi(i,j,k,a))*(phi(i,j+1,k,b)-phi(i,j,k,b))/(delta[Y]) 
					-0.5*(phi(i,j+1,k,b)+phi(i,j,k,b))*(phi(i,j+1,k,a)-phi(i,j,k,a))/(delta[Y]);
	
	qab(jmh,Y) = 0.5*(phi(i,j,k,a)+phi(i,j-1,k,a))*(phi(i,j,k,b)-phi(i,j-1,k,b))/(delta[Y]) 
					-0.5*(phi(i,j,k,b)+phi(i,j-1,k,b))*(phi(i,j,k,a)-phi(i,j-1,k,a))/(delta[Y]);
#if(AMREX_SPACEDIM>2)
	qab(kph,Y) = 0.5*(phi(i,j,k+1,a)+phi(i,j,k,a))*(-phi(i,j-1,k+1,b)+phi(i,j+1,k+1,b)-phi(i,j-1,k,b)+phi(i,j+1,k,b))/(4*delta[Y]) 
					-0.5*(phi(i,j,k+1,b)+phi(i,j,k,b))*(-phi(i,j-1,k+1,a)+phi(i,j+1,k+1,a)-phi(i,j-1,k,a)+phi(i,j+1,k,a))/(4*delta[Y]);
	
	qab(kmh,Y) = 0.5*(phi(i,j,k,a)+phi(i,j,k-1,a))*(-phi(i,j-1,k,b)+phi(i,j+1,k,b)-phi(i,j-1,k-1,b)+phi(i,j+1,k-1,b))/(4*delta[Y]) 
					-0.5*(phi(i,j,k,b)+phi(i,j,k-1,b))*(-phi(i,j-1,k,a)+phi(i,j+1,k,a)-phi(i,j-1,k-1,a)+phi(i,j+1,k-1,a))/(4*delta[Y]);
#endif

#if (AMREX_SPACEDIM > 2)	

	qab(cent,Z) = phi(i,j,k,a)*(phi(i,j,k+1,b)-phi(i,j,k-1,b))/(2*delta[Z]) 
					-phi(i,j,k,b)*(phi(i,j,k+1,a)-phi(i,j,k-1,a))/(2*delta[Z]);
	
	qab(iph,Z) = 0.5*(phi(i,j,k,a)+phi(i+1,j,k,a))*(phi(i+1,j,k+1,b)-phi(i+1,j,k-1,b)+phi(i,j,k+1,b)-phi(i,j,k-1,b))/(4*delta[Z]) 
					-0.5*(phi(i,j,k,b)+phi(i+1,j,k,b))*(phi(i+1,j,k+1,a)-phi(i+1,j,k-1,a)+phi(i,j,k+1,a)-phi(i,j,k-1,a))/(4*delta[Z]);
	
	qab(imh,Z) = 0.5*(phi(i,j,k,a)+phi(i-1,j,k,a))*(phi(i,j,k+1,b)-phi(i,j,k-1,b)+phi(i-1,j,k+1,b)-phi(i-1,j,k-1,b))/(4*delta[Z]) 
					-0.5*(phi(i,j,k,b)+phi(i-1,j,k,b))*(phi(i,j,k+1,a)-phi(i,j,k-1,a)+phi(i-1,j,k+1,a)-phi(i-1,j,k-1,a))/(4*delta[Z]);
	
	qab(jph,Z) = 0.5*(phi(i,j,k,a)+phi(i,j+1,k,a))*(phi(i,j+1,k+1,b)-phi(i,j+1,k-1,b)+phi(i,j,k+1,b)-phi(i,j,k-1,b))/(4*delta[Z]) 
					-0.5*(phi(i,j,k,b)+phi(i,j+1,k,b))*(phi(i,j+1,k+1,a)-phi(i,j+1,k-1,a)+phi(i,j,k+1,a)-phi(i,j,k-1,a))/(4*delta[Z]);
	
	qab(jmh,Z) = 0.5*(phi(i,j,k,a)+phi(i,j-1,k,a))*(phi(i,j,k+1,b)-phi(i,j,k-1,b)+phi(i,j-1,k+1,b)-phi(i,j-1,k-1,b))/(4*delta[Z]) 
					-0.5*(phi(i,j,k,b)+phi(i,j-1,k,b))*(phi(i,j,k+1,a)-phi(i,j,k-1,a)+phi(i,j-1,k+1,a)-phi(i,j-1,k-1,a))/(4*delta[Z]);
	
	qab(kph,Z) = 0.5*(phi(i,j,k,a)+phi(i,j,k+1,a))*(phi(i,j,k+1,b)-phi(i,j,k,b))/(delta[Z]) 
					-0.5*(phi(i,j,k,b)+phi(i,j,k+1,b))*(phi(i,j,k+1,a)-phi(i,j,k,a))/(delta[Z]);
	
	qab(kmh,Z) = 0.5*(phi(i,j,k,a)+phi(i,j,k-1,a))*(phi(i,j,k,b)-phi(i,j,k-1,b))/(delta[Z]) 
					-0.5*(phi(i,j,k,b)+phi(i,j,k-1,b))*(phi(i,j,k,a)-phi(i,j,k-1,a))/(delta[Z]);

#endif
}
//####################################################################################################################