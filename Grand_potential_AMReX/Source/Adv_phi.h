#ifndef ADV_PHI_H
#define ADV_PHI_H

#include <AMReX_Utility.H>
#include "Proj_on_simplex.H"

using namespace amrex;

void update_phi(MultiFab& phi_new, MultiFab& phi_old, MultiFab& term1, MultiFab& term2, MultiFab& term3,MultiFab& lambad, Geometry const& geom)
{	
	BL_PROFILE("update_phi()");	

	#ifdef AMREX_USE_OMP
		#pragma omp parallel if (Gpu::notInLaunchRegion())
	#endif

	for (MFIter mfi(phi_old); mfi.isValid(); ++mfi)
	{
		const Box& dbx = mfi.validbox();
		Array4<Real> const& fin_term1 = term1.array(mfi);
		Array4<Real> const& fin_term2 = term2.array(mfi);
		Array4<Real> const& fin_term3 = term3.array(mfi);
		Array4<Real> const& phiNew = phi_new.array(mfi);
		Array4<Real> const& phiOld = phi_old.array(mfi);
		Array4<Real> const& lamb = lambad.array(mfi);
		
		// Array2D <Real, 0, phasecount, 0, phasecount> tauab{};
		// for(int a =0; a<nump; a++){
		// 	for(int b=0; b<nump; b++){
		// 		tauab(a,b) = tau_ab[a][b];
		// 	}
		// }

		Real Tauu = tau;
		Real time_step = dt;
		Real epsilon = eps;
		Real molar_vol = Vm;
		int numphase = nump; 
		int dimsn = dim;
        GpuArray<Real,AMREX_SPACEDIM> delta = geom.CellSizeArray();
				

		amrex::ParallelFor( dbx, 
		[=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
		{	

			Real sum_lambda{0.0};
			Real active_phase{0.0};
			Array1D<Real,0,phasecount-1> deltaphi{};
			Array1D<Real,0,phasecount-1> div{};

			//Calculate divergence of phi for 2D -----------------------------------------------
			if(dimsn == 2){
				for(int a=0; a<numphase; a++){
				div(a) = (phiOld(i+1,j,k,a)-2.0*phiOld(i,j,k,a)+phiOld(i-1,j,k,a))/(delta[X]*delta[X])+(phiOld(i,j+1,k,a)-2.0*phiOld(i,j,k,a)+phiOld(i,j-1,k,a))/(delta[Y]*delta[Y]);
			}
			}

			//Calculate divergence of phi for 3D -----------------------------------------------
			if(dimsn == 3){
				for(int a=0; a<numphase; a++){
				div(a) = (phiOld(i+1,j,k,a)-2.0*phiOld(i,j,k,a)+phiOld(i-1,j,k,a))/(delta[X]*delta[X])+(phiOld(i,j+1,k,a)-2.0*phiOld(i,j,k,a)+phiOld(i,j-1,k,a))/(delta[Y]*delta[Y])+(phiOld(i,j,k+1,a)-2.0*phiOld(i,j,k,a)+phiOld(i,j,k-1,a))/(delta[Z]*delta[Z]);
			}
			}

			for(int a=0; a<numphase; a++){

				if(fabs(div(a))>0.0){
					
					lamb(i,j,k,a) = epsilon*fin_term1(i,j,k,a)-fin_term2(i,j,k,a)/epsilon-fin_term3(i,j,k,a)/molar_vol;

					sum_lambda += lamb(i,j,k,a); 
					active_phase++;
				}
			}

			 if (active_phase) {
      			sum_lambda /= active_phase;
    			}
				
			//Calculating delta phi ----------------------------------------------------------------------------
			for(int a=0; a<numphase; a++){

				if(fabs(div(a))>0.0){
					deltaphi(a) = (time_step/(epsilon*FunctionTau(i,j,k,numphase,Tauu,phiOld)))*(lamb(i,j,k,a)-sum_lambda);
				}
				else{
					deltaphi(a) = 0.0;
				}
			}

			//Keeping values between bounds --------------------------------------------------
			projection_on_simplex(i,j,k,phiOld,phiNew,deltaphi,div,numphase);

			//Updating phi -------------------------------------------------
			for(int a = 0; a<numphase;a++){
				phiNew(i,j,k,a) = phiOld(i,j,k,a) + deltaphi(a);
			}

		});
	}
}

#endif