#ifndef MY_KERNEL_H_
#define MY_KERNEL_H_

#include <AMReX_FArrayBox.H>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <AMReX_ParallelDescriptor.H>
#include <time.h>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void init_phi (int i, int j, int k,
               amrex::Array4<amrex::Real> const& phi,
               amrex::Real c0)
{
    using amrex::Real;;
	Real noise = 0.02;

	
    using namespace std;

    phi(i,j,k) = c0 + noise*(0.5 - amrex::Random());
    //phi(i,j,k) = c0 + noise*(0.5 - ((double) rand() / (RAND_MAX)));
    //std::cout<<phi(i,j,k)<<"\n";
}



AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void update_phi (int i, int j, int k,
                 amrex::Array4<amrex::Real const> const& phiold,
                 amrex::Array4<amrex::Real      > const& phinew,
		 amrex::Array4<amrex::Real      > const& df,
		 amrex::Real M,
                 amrex::Real dt)
{
    phinew(i,j,k) = phiold(i,j,k) + (dt * M * df(i,j,k));
    if (phinew(i,j,k) >= 0.9999)
    {
      phinew(i,j,k) = 0.9999;
    }
    if (phinew(i,j,k) < 0.00001)
    {
      phinew(i,j,k) = 0.00001;
    }
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void derivative(int i, int j, int k,
		amrex::Array4<amrex::Real const> const& phiOld,
		amrex::Array4<amrex::Real      > const& deriv)
{
		
            	deriv(i,j,k) = (2.0*phiOld(i,j,k)*(1.0-phiOld(i,j,k))*(1.0-phiOld(i,j,k)) - 2.0*phiOld(i,j,k)*phiOld(i,j,k)*(1.0 - phiOld(i,j,k)));
        	

}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void laplacian(int i, int j, int k,
		amrex::Array4<amrex::Real const> const& Phiold,
		amrex::Array4<amrex::Real      > const& lap,
		Geometry const& geom)
{
	GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();
	//std::cout<<dx[0]<<"   "<<dx[1]<<"\n";
	
            	lap(i,j,k) = (Phiold(i-1,j,k) + Phiold(i+1,j,k) + Phiold(i,j-1,k) + Phiold(i,j+1,k) - 4.0 * Phiold(i,j,k))/(dx[0]*dx[1]);
        		
	
}










/*
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void energy1cal(int i, int j, int k,
		amrex::Array4<amrex::Real const> const& ma,
		amrex::Real const& energ)
{
	energ = energ + ma(i,j,k)*ma(i,j,k)*(1.0 - ma(i,j,k))*(1.0 - ma(i,j,k));

}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void energy2cal(int i, int j, int k,
		amrex::Array4<amrex::Real const> const& ma,
		amrex::Real const& energ,
		amrex::Real const& K,
		Geometry const& geom)
{
	const Box& domain_bx = geom.Domain();
	const auto dom_hi = ubound(domain_bx);
	
	if (i<dom_hi.x && j<dom_hi.y)
			{
				energ = energ + (K/2)*((ma(i+1,j,k)-ma(i,j,k))*(ma(i+1,j,k)-ma(i,j,k)) + (ma(i,j+1,k)-ma(i,j,k))*(ma(i,j+1,k)-ma(i,j,k)));
			}

}


*/


#endif