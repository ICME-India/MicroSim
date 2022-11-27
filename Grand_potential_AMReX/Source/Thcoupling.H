#ifndef THCOUPLING_H_
#define THCOUPLING_H_

using namespace amrex;

void computeterm3(amrex::MultiFab& mu_old, amrex::MultiFab& term3, amrex::MultiFab& phi_old, amrex::MultiFab& dh_dphi, amrex::MultiFab& psi)
{
	
	int nc = psi.nComp();
	for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();
		Array4<Real> const& phiold = phi_old.array(mfi);
		Array4<Real> const& term = term3.array(mfi);
		Array4<Real> const& mu = mu_old.array(mfi);
		Array4<Real> const& psii = psi.array(mfi);
		Array4<Real> const& arr_dh_dphi = dh_dphi.array(mfi);
		Real numphase = nump;
		Array1D <Real,0,1> BB;
		BB(0)=B[0];
		BB(1)=B[1];
		
		Array1D <Real,0,1> DD;
		DD(0)=D[0];
		DD(1)=D[1];
		
		Array1D <Real,0,1> der_cmu;
		der_cmu(0) = dcdmu[0];
		der_cmu(1) = dcdmu[1];
	
		amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {	
			
			for(int a=0; a<numphase; a++){
			
				psii(i,j,k,a) = -pow((mu(i,j,k) - BB(a)),2)*der_cmu(a)*0.5 + DD(a);
			
			}
			
			term(i,j,k) = (arr_dh_dphi(i,j,k))*(psii(i,j,k,0) - psii(i,j,k,1));
			
		});
	
	}
}

#endif
