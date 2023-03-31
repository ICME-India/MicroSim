#ifndef _FUNCTIONW_H_
#define _FUNCTIONW_H_

using namespace amrex;

void computeterm2(amrex::MultiFab& term2, amrex::MultiFab& phi_old)
{
	BL_PROFILE("computeterm2()");
	
	for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();
		Array4<Real> const& phiold = phi_old.array(mfi);
		Array4<Real> const& term = term2.array(mfi);
		Real GAM = gammaa;
	
		amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
			term(i,j,k) = 9.0*GAM*2.0*phiold(i,j,k)*(1.0-phiold(i,j,k))*(1.0 - 2.0*phiold(i,j,k));
		});
	
	}
}


#endif
