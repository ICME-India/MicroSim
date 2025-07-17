#ifndef MERGE_INIT_H_
#define MERGE_INIT_H_

#include <AMReX_Utility.H>
#include "Variables.h"

using namespace std;

void merge_initial(MultiFab& merge_phi, MultiFab& phi_new, int numphase){
    
	#ifdef AMREX_USE_OMP
	#pragma omp parallel if (Gpu::notInLaunchRegion())
	#endif

	for ( MFIter mfi(phi_new); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();									//Defining the box for iteration space
		Array4<Real> const& phiNew = phi_new.array(mfi);					//Taking the Multifabs as arrays
		Array4<Real> const& mp = merge_phi.array(mfi);						//Taking the Multifabs as arrays
		
        int numph = numphase;
		
		amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {	
			for(int a = 0; a<numph;a++){
				mp(i,j,k) += (a+1)*phiNew(i,j,k,a);
			}

		});

	}
}

#endif