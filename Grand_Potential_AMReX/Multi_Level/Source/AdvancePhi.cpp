#include <AMReX_MultiFabUtil.H>
#include "AmrCoreGP.H"
#include "Proj_on_simplex.H"

using namespace amrex;
//####################################################################################################################
void 
AmrCoreGP::update_phi(amrex::MultiFab& phi_new, amrex::MultiFab& phi_old,amrex::MultiFab& term1, amrex::MultiFab& term2, amrex::MultiFab& term3,amrex::MultiFab& term4,amrex::MultiFab& lambad, amrex::GpuArray<Real,AMREX_SPACEDIM> dx, amrex::Real dt)
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
		Array4<Real> const& fin_term4 = term4.array(mfi);
		Array4<Real> const& phiNew = phi_new.array(mfi);
		Array4<Real> const& phiOld = phi_old.array(mfi);
		Array4<Real> const& lamb = lambad.array(mfi);
		Real Tauu = tau;
		Real time_step = dt;
		Real epsilon = eps;
		Real molar_vol = Vm;
		int numphase = nump; 
		int dimsn = dim;
		
        GpuArray<Real,AMREX_SPACEDIM> delta = dx;
		amrex::ParallelFor( dbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
		{	
            Real sum_lambda{0.0};
			Real active_phase{0.0};
			Array1D<Real,0,phasecount-1> deltaphi{};
			Array1D<Real,0,phasecount-1> div{};
			//Calculate divergence of phi for 2D -----------------------------------------------
			if(dimsn == 2)
            {
				for(int a=0; a<numphase; a++)
                {
				    div(a) = (phiOld(i+1,j,k,a)-2.0*phiOld(i,j,k,a)+phiOld(i-1,j,k,a))/(delta[X]*delta[X])+(phiOld(i,j+1,k,a)-2.0*phiOld(i,j,k,a)+phiOld(i,j-1,k,a))/(delta[Y]*delta[Y]);
			    }
			}
			//Calculate divergence of phi for 3D -----------------------------------------------
			if(dimsn == 3)
            {
				for(int a=0; a<numphase; a++)
                {
				    div(a) = (phiOld(i+1,j,k,a)-2.0*phiOld(i,j,k,a)+phiOld(i-1,j,k,a))/(delta[X]*delta[X])+(phiOld(i,j+1,k,a)-2.0*phiOld(i,j,k,a)+phiOld(i,j-1,k,a))/(delta[Y]*delta[Y])+(phiOld(i,j,k+1,a)-2.0*phiOld(i,j,k,a)+phiOld(i,j,k-1,a))/(delta[Z]*delta[Z]);
			    }
			}

			for(int a=0; a<numphase; a++)
            {
                if(fabs(div(a))>0.0)
                {
					lamb(i,j,k,a) = epsilon*fin_term1(i,j,k,a)-fin_term2(i,j,k,a)/epsilon-fin_term3(i,j,k,a)/molar_vol;
					if(ELASTICITY==1)
                    {
						lamb(i,j,k,a) = lamb(i,j,k,a) - fin_term4(i,j,k,a);
					}
					sum_lambda += lamb(i,j,k,a); 
					active_phase++;
				}
			}

			if (active_phase) 
            {
      			sum_lambda /= active_phase;
    		}
			//Calculating delta phi ----------------------------------------------------------------------------
			for(int a=0; a<numphase; a++)
            {
				if(fabs(div(a))>0.0)
                {
					deltaphi(a) = (time_step/(epsilon*FunctionTau(i,j,k,numphase,Tauu,phiOld)))*(lamb(i,j,k,a)-sum_lambda);
				}
				else
                {
					deltaphi(a) = 0.0;
				}
			}
			//Keeping values between bounds --------------------------------------------------
			projection_on_simplex(i,j,k,phiOld,phiNew,deltaphi,div,numphase);
			//Updating phi -------------------------------------------------
			for(int a = 0; a<numphase;a++)
            {
				phiNew(i,j,k,a) = phiOld(i,j,k,a) + deltaphi(a);
			}
		});
	}
}

//####################################################################################################################
void 
AmrCoreGP::AdvancePhiAtLevel(amrex::Real time,amrex::Real dt, int lev)
{		
	GpuArray<Real,AMREX_SPACEDIM> dx = geom[lev].CellSizeArray();
    MultiFab temp_phi(grids[lev], dmap[lev], phi_new[lev].nComp(),1);
	FillPatch("P",lev, time, temp_phi, 0, temp_phi.nComp());
	std::swap(phi_old[lev], phi_new[lev]);
    MultiFab& term1_state = term1[lev];
    if(funcANI == 0)
    {
        if(dim == 2)
        {                                                                                                                                                   
            function_A_00_iso_2D(term1_state, temp_phi, dx);
        }
        if(dim == 3)
        {
            function_A_00_iso_3D(term1_state, temp_phi, dx);
        }
    }
    if(funcANI == 1)
    {
        if(dim == 2)
        {
            function_A_01_ani_2D(term1_state, temp_phi, dx);
        }
        if(dim == 3)
        {
            function_A_01_ani_3D(term1_state, temp_phi, dx);
        }
    }
    MultiFab& term2_state = term2[lev];
    if(funcW == 1)
    {    
		function_W_01_dwdphi(term2_state, temp_phi, dx);
    }
    if(funcW == 2)
    {
        function_W_02_dwdphi(term2_state, temp_phi, dx);
    }
    MultiFab temp_mu(grids[lev], dmap[lev], mu_new[lev].nComp(),1);
	FillPatch("M",lev, time, temp_mu, 0, temp_mu.nComp());
	std::swap(mu_old[lev], mu_new[lev]);
    MultiFab& term3_state = term3[lev];
    MultiFab& psi_state = psi[lev];
    dpsi(temp_mu, term3_state, temp_phi, psi_state, dx);

	MultiFab temp_term1(grids[lev], dmap[lev], term1[lev].nComp(),1);
    MultiFab temp_term2(grids[lev], dmap[lev], term2[lev].nComp(),1);
    MultiFab temp_term3(grids[lev], dmap[lev], term3[lev].nComp(),1);
    MultiFab temp_term4(grids[lev], dmap[lev], term4[lev].nComp(),1);
    FillPatch("T1",lev, time, temp_term1, 0, temp_term1.nComp());
    FillPatch("T2",lev, time, temp_term2, 0, temp_term2.nComp());
    FillPatch("T3",lev, time, temp_term3, 0, temp_term3.nComp());
    FillPatch("T4",lev, time, temp_term4, 0, temp_term4.nComp());
    MultiFab& phi_state = phi_new[lev];
    MultiFab& lambad_state = lambad[lev];
    update_phi(phi_state, temp_phi, temp_term1, temp_term2, temp_term3,temp_term4,lambad_state,dx, dt);

	MultiFab temp_phi1(grids[lev], dmap[lev], phi_new[lev].nComp(),1);
	MultiFab temp_comp(grids[lev], dmap[lev], comp_new[lev].nComp(),1);
	FillPatch("P",lev, time, temp_phi1, 0, temp_phi1.nComp());
	FillPatch("C",lev, time, temp_comp, 0, temp_comp.nComp());
	std::swap(comp_old[lev], comp_new[lev]);
		
	MultiFab& mu_state = mu_new[lev];
	MultiFab& comp_state = comp_new[lev];
    if(dim == 2)
	{
		dmudt_2D(mu_state, temp_mu, temp_phi1, temp_phi, comp_state, temp_comp,dx, dt);
    }
    if(dim == 3)
	{
		dmudt_3D(mu_state, temp_mu, temp_phi1, temp_phi, comp_state, temp_comp,dx, dt);
    }
}

//####################################################################################################################
void 
AmrCoreGP::AdvancePhiAllLevels(amrex::Real time,amrex::Real dt)
{		
    for (int lev = 0; lev <= finest_level; lev++)
    {
        AdvancePhiAtLevel(time,dt,lev);
    }
}

//####################################################################################################################