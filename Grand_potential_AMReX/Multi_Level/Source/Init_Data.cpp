#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Utility.H>
#include <AMReX_Print.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "AmrCoreGP.H"

using namespace amrex;
using namespace std;

//####################################################################################################################

void 
AmrCoreGP::Init_liq(MultiFab& phi_new)
{
	for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
	{
		const Box& wbx = mfi.validbox();
		Array4<Real> const& phiNew = phi_new.array(mfi);
		Real numphase = nump;
		amrex::ParallelFor( wbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
		{		
			Real sum{0.0};
			for(int a=0; a<numphase-1; a++)
			{
				sum = sum + phiNew(i,j,k,a);
			}

			phiNew(i,j,k,numphase-1) = 1.0 - sum;
		});
	}
}
//####################################################################################################################

void
AmrCoreGP::Init_mu(amrex::MultiFab& mu_new, amrex::MultiFab& phi_new)
{ 	
	for (MFIter mfi(mu_new); mfi.isValid(); ++mfi)
	{
		const Box& wbx = mfi.validbox();
		Array4<Real> const& muNew = mu_new.array(mfi);
		Array4<Real> const& phiNew = phi_new.array(mfi);
		Array2D<Real,0,compcount-2,0,compcount-2> A_liq{};
		for(int l=0; l<numcom-1; l++)
		{
			for(int m=0; m<numcom-1;m++)
			{
				A_liq(l,m) = A[nump-1][l][m];
			}
		}

		Array2D<Real,0,compcount-2,0,compcount-2> A_sol{};
		for(int l=0; l<numcom-1; l++)
		{
			for(int m=0; m<numcom-1;m++)
			{
				A_sol(l,m) = A[0][l][m];
			}
		}

		Array1D<Real,0,compcount-2> ceq_liq{};
		for(int l=0; l<numcom-1; l++)
		{
			ceq_liq(l) = conceq[nump-1][l];
		}

		Array1D<Real,0,compcount-2> ceq_sol{};
		for(int l=0; l<numcom-1; l++)
		{
			ceq_sol(l) = conceq[0][l];
		}

		Array2D <Real,0,phasecount-1,0,compcount-2, Order::C> BB{};
		for(int a=0; a<nump; a++)
		{
			for(int l=0; l<numcom-1; l++)
			{
				BB(a,l) = B[a][l];
			}
		}

		Real numcomp = numcom;
		amrex::ParallelFor( wbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
		{	
			double sum =0.0;
			if(phiNew(i,j,k,0)==1)
			{
				for(int l=0; l<numcomp-1; l++)
				{
					for(int m=0; m<numcomp-1; m++)
					{
						if(l==m)
						{
							sum += 2.0*A_sol(l,m)*ceq_sol(m);
						}
						else
						{
							sum += A_sol(l,m)*ceq_sol(m);
						}
					}
					muNew(i,j,k,l) = sum + BB(0,l);
					sum=0.0;
				}
			}
			else
			{
				for(int l=0; l<numcomp-1; l++)
				{
					for(int m=0; m<numcomp-1; m++)
					{
						if(l==m)
						{
							sum += 2.0*A_liq(l,m)*ceq_liq(m);
						}
						else
						{
							sum += A_liq(l,m)*ceq_liq(m);
						}
					}
					muNew(i,j,k,l) = sum;
					sum=0.0;
				}
			}
		});
	}
}	

//####################################################################################################################
void 
AmrCoreGP::Init_comp(MultiFab& phi_new, MultiFab& comp_new)
{
	for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
	{	
		const Box& pbx = mfi.validbox();
		Array4<Real> const& phiNew = phi_new.array(mfi);
		Array4<Real> const& compNew = comp_new.array(mfi);
		int numphase = nump;
		int numcomp = numcom;

		Array2D<Real,0,phasecount-1,0,compcount-2,Order::C> co{};
		for(int a=0; a<nump; a++){
			for(int l=0; l<numcom-1; l++){
				co(a,l) = conc[a][l];
			}
		}

		Array2D<Real,0,phasecount-1,0,compcount-2,Order::C> coeq{};
		for(int a=0; a<nump; a++){
			for(int l=0; l<numcom-1; l++){
				coeq(a,l) = conceq[a][l];
			}
		} 
		
		amrex::ParallelFor( pbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
		{		
				Real sum{0.0};
				int val{0};

				for(int a=0; a<numphase-1; a++){
					
					if(phiNew(i,j,k,a)==1.0){
						for(int l=0; l<numcomp-1; l++){
						compNew(i,j,k,l) = coeq(a,l);
					}
					val=1;
					break;
					}
				}

				if(val==0){
					for(int l=0; l<numcomp-1; l++){
						compNew(i,j,k,l) = co(nump-1,l);;
					}
				}

		});
		}
}

//####################################################################################################################
void
AmrCoreGP::Init_phi_cyl(amrex::MultiFab& phi_new)
{
	for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
	{
		const Box& wbx = mfi.validbox();
		Array4<Real> const& phiNew = phi_new.array(mfi);
		for(int p=0; p<cylinder.size(); p++)
        {
			Real cyl_comp = cylinder[p][0];
			Real cyl_X_cent = cylinder[p][1];
			Real cyl_Y_cent = cylinder[p][2];
			Real cyl_Z_strt = cylinder[p][3];
			Real cyl_Z_end = cylinder[p][4];
			Real cyl_rad = cylinder[p][5];
		    amrex::ParallelFor(wbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
		    {	
			    #if (AMREX_SPACEDIM <= 2)
			    if(((i-cyl_X_cent)*(i-cyl_X_cent) + (j-cyl_Y_cent)*(j-cyl_Y_cent)) < cyl_rad*cyl_rad)
			    {
			    	phiNew(i,j,k,cyl_comp) = 1.0;
			    }
			    #endif
			    #if (AMREX_SPACEDIM > 2)
			    if(((i-cyl_X_cent)*(i-cyl_X_cent) + (j-cyl_Y_cent)*(j-cyl_Y_cent)) < cyl_rad*cyl_rad && cyl_Z_strt<=k<=cyl_Z_end)
			    {
				    phiNew(i,j,k,cyl_comp) = 1.0;
			    }
			    #endif 
		    });
		}
	}
}

//####################################################################################################################

void
AmrCoreGP:: Init_phi_sph(MultiFab& phi_new)
{
	for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
	{
		const Box& wbx = mfi.validbox();
		auto const& phiNew = phi_new.array(mfi);
		
		for(int p=0; p<sphere.size(); p++){
			Real sph_phase = sphere[p][0];
			Real sph_X_cent = sphere[p][1];
			Real sph_Y_cent = sphere[p][2];
			Real sph_Z_cent = sphere[p][3];
			Real sph_rad = sphere[p][4];

		amrex::ParallelFor( wbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
		{	
			#if (AMREX_SPACEDIM>2)
			if(((i-sph_X_cent)*(i-sph_X_cent) + (j-sph_Y_cent)*(j-sph_Y_cent)+(k-sph_Z_cent)*(k-sph_Z_cent)) < sph_rad*sph_rad)
			{
				phiNew(i,j,k,sph_phase) = 1.0;
			}
			#endif
			
		});
	}
	}
}

//####################################################################################################################
