#ifndef INTERP_H_
#define INTERP_H_

#include <AMReX_Utility.H>
#include "Function_F4.H"

using namespace amrex;

void Func_H(MultiFab& h_phi,MultiFab& dh_dphi, MultiFab& phi_old){

	BL_PROFILE("Func_H()");

    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();
		Array4<Real> const& phiold = phi_old.array(mfi);
        Array4<Real> const& arrh_phi = h_phi.array(mfi);
        Array4<Real> const& arrdh_dphi = dh_dphi.array(mfi);
	
		amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // if(phiold(i,j,k) == 0)
            // {
            //     arrh_phi(i,j,k) = 0;
            // }
            // else if(phiold(i,j,k) == 1){
                
            //     arrh_phi(i,j,k) = 1;
            // }
            //else{
            arrh_phi(i,j,k) = pow(phiold(i,j,k),3)*(10.0 - 15.0*phiold(i,j,k) + 6.0*pow(phiold(i,j,k),2));
            //}

            //if(phiold(i,j,k) == 0 || phiold(i,j,k) == 1)
            //{
            //    arrdh_dphi(i,j,k) = 0;
            //}
            //else{
            arrdh_dphi(i,j,k) = 30.0*pow(phiold(i,j,k),2)*pow((1.0-phiold(i,j,k)),2);
            //}
		});
	
	}


}



#endif
