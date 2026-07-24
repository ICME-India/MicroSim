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
AmrCoreGP::function_W_01_dwdphi(amrex::MultiFab& term2, amrex::MultiFab& phi_old, amrex::GpuArray<Real,AMREX_SPACEDIM> dx)
{
	BL_PROFILE("computeterm2()");
	
	#ifdef AMREX_USE_OMP
	#pragma omp parallel if (Gpu::notInLaunchRegion())
	#endif

	for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();									//Defining the box for iteration space
		Array4<Real> const& phiOld = phi_old.array(mfi);					//Taking the Multifabs as arrays
		Array4<Real> const& term = term2.array(mfi);						//Taking the Multifabs as arrays
		//Redefining variables in GPU space	--------------------------------------------------------------------------
		int numphase = nump;
		int dimsn = dim;
		//Turning the vector gamma to a GPU compatible array----------------------------------------------------------
		Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> gamm{};
		for(int a=0; a<numphase; a++)
        {
			for(int b=a; b<numphase; b++)
			{
				gamm(a,b) = gam[a][b];
				gamm(b,a) = gam[a][b];
			}
		}
		//Turning the vector gamma_abc to a GPU compatible array----------------------------------------------------------
		Array3D <Real,0,phasecount-1,0,phasecount-1,0,phasecount-1,Order::C> gamm_abc{};
		if(nump>2)
        {
		    for(int a=0; a<numphase; a++)
            {
			    for(int b=0; b<numphase; b++)
                {
				    for(int c=0; c<numphase;c++)
                    {
					    gamm_abc(a,b,c) = gam_abc[a][b][c];
				    }
			    }
		    }
		}
		//delta stores dx, dy and dz -----------------------------------------------------------
		GpuArray<Real,AMREX_SPACEDIM> delta = dx;
		amrex::ParallelFor(vbx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {	
			//Declaring variables ---------------------------------------------------------------
			Real sum{0.0};
			Real phi_phi{0.0};
			Array1D<Real,0,phasecount-1> div{};
			//Calculating divergence for 2-D ----------------------------------------------------------
			if(dimsn == 2)
            {
				for(int a=0; a<numphase; a++)
                {
				    div(a) = (phiOld(i+1,j,k,a)-2.0*phiOld(i,j,k,a)+phiOld(i-1,j,k,a))/(delta[X]*delta[X])
						    +(phiOld(i,j+1,k,a)-2.0*phiOld(i,j,k,a)+phiOld(i,j-1,k,a))/(delta[Y]*delta[Y]);
			    }
			}
			
			//Calculating divergence for 3-D ----------------------------------------------------------
			if(dimsn == 3)
            {
				for(int a=0; a<numphase; a++)
                {
				    div(a) = (phiOld(i+1,j,k,a)-2.0*phiOld(i,j,k,a)+phiOld(i-1,j,k,a))/(delta[X]*delta[X])
						    +(phiOld(i,j+1,k,a)-2.0*phiOld(i,j,k,a)+phiOld(i,j-1,k,a))/(delta[Y]*delta[Y])
						    +(phiOld(i,j,k+1,a)-2.0*phiOld(i,j,k,a)+phiOld(i,j,k-1,a))/(delta[Z]*delta[Z]);
			    }
			}
			
			//Calculating \frac{16}{\pi^{2} }\cdot \phi_{a} \sum_{a<b}^{N}\gamma_{ab} \cdot \phi_{b}+ \sum_{a<b<c}^{N}\gamma_{abc}\cdot \phi_{b} \cdot \phi_{c}  ------------------------------------------
			for(int a =0; a<numphase; a++)
            {

				for(int b=0; b<numphase;b++)
                {
					if(b!=a)
                    {
						if(fabs(div(b)) > 0.0)
                        {
							if(phiOld(i,j,k,a)*phiOld(i,j,k,b)>=0.00)
                            {
								sum+= gamm(a,b)*phiOld(i,j,k,b);	
							}
							else
                            {
								sum+= -gamm(a,b)*phiOld(i,j,k,b);
							}
							//sum = sum + gamm(a,b)*phiOld(i,j,k,b);
						}
					}
				}
				sum*= 16/pow(M_PI ,2);
				for(int b =0; b<numphase; b++)
                {
					for(int c = 0; c<numphase; c++)
                    {

						if(b!=a && c!=a && b<c)
                        {

							if(fabs(div(b))>0.0 && fabs(div(c))>0.0)
                            { 
								phi_phi = phiOld(i,j,k,b)*phiOld(i,j,k,c);
								sum += gamm_abc(a,b,c)*phi_phi;
							}
						}
					}
				}
			    term(i,j,k,a) = sum;
			    phi_phi=0.0;
			    sum=0.0;
			}

		});

	}
}

//####################################################################################################################
void
AmrCoreGP::function_W_02_dwdphi(amrex::MultiFab& term2, amrex::MultiFab& phi_old, amrex::GpuArray<Real,AMREX_SPACEDIM> dx)
{
	BL_PROFILE("computeterm2()");
	#ifdef AMREX_USE_OMP
	#pragma omp parallel if (Gpu::notInLaunchRegion())
	#endif
	for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();									//Defining the box for iteration space
		Array4<Real> const& phiOld = phi_old.array(mfi);					//Taking the Multifabs as arrays
		Array4<Real> const& term = term2.array(mfi);						//Taking the Multifabs as arrays
		//Redefining variables in GPU space	--------------------------------------------------------------------------
		int numphase = nump;
		int dimsn = dim;
		//Turning the vector gamma to a GPU compatible array----------------------------------------------------------
		Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> gamm{};
		for(int a=0; a<numphase; a++)
        {
			for(int b=a; b<numphase; b++)
            {
				gamm(a,b) = gam[a][b];
				gamm(b,a) = gam[a][b];
			}
		}
		//Turning the vector gamma_abc to a GPU compatible array----------------------------------------------------------
		Array3D <Real,0,phasecount-1,0,phasecount-1,0,phasecount-1,Order::C> gamm_abc{};
		if(nump>2)
        {
		    for(int a=0; a<numphase; a++)
            {
			    for(int b=0; b<numphase; b++)
                {
				    for(int c=0; c<numphase;c++)
                    {
					    gamm_abc(a,b,c) = gam_abc[a][b][c];
				    }
			    }
		    }
		}

		//delta stores dx, dy and dz -----------------------------------------------------------
		GpuArray<Real,AMREX_SPACEDIM> delta = dx;
		amrex::ParallelFor(vbx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {	
			//Declaring variables ---------------------------------------------------------------
			Real sum{0.0};
			Real phi_phi{0.0};
			Array1D<Real,0,phasecount-1> div{};
			//Calculating divergence for 2-D ----------------------------------------------------------
			if(dimsn == 2)
            {
				for(int a=0; a<numphase; a++)
                {
				    div(a) = (phiOld(i+1,j,k,a)-2.0*phiOld(i,j,k,a)+phiOld(i-1,j,k,a))/(delta[X]*delta[X])
						    +(phiOld(i,j+1,k,a)-2.0*phiOld(i,j,k,a)+phiOld(i,j-1,k,a))/(delta[Y]*delta[Y]);
			    }
			}
			//Calculating divergence for 3-D ----------------------------------------------------------
			if(dimsn == 3)
            {
				for(int a=0; a<numphase; a++)
                {
				    div(a) = (phiOld(i+1,j,k,a)-2.0*phiOld(i,j,k,a)+phiOld(i-1,j,k,a))/(delta[X]*delta[X])
						    +(phiOld(i,j+1,k,a)-2.0*phiOld(i,j,k,a)+phiOld(i,j-1,k,a))/(delta[Y]*delta[Y])
						    +(phiOld(i,j,k+1,a)-2.0*phiOld(i,j,k,a)+phiOld(i,j,k-1,a))/(delta[Z]*delta[Z]);
			    }
			}
			//Calculating 18.0 \cdot \phi_{a} \sum_{a<b}^{N}\gamma_{ab} \cdot \phi_{b}^{2} + \sum_{a<b<c}^{N}\gamma_{abc}\cdot \phi_{b} \cdot \phi_{c} ------------------------------------------
			for(int a =0; a<numphase; a++)
            {

				for(int b=0; b<numphase;b++)
                {
					if(b!=a)
                    {
						if(fabs(div(b)) > 0.0)
                        {
							sum = sum + 2.0*gamm(a,b)*phiOld(i,j,k,a)*phiOld(i,j,k,b)*phiOld(i,j,k,b);
						}
					}
				}
				sum*=9.0;
				if(numphase>2)
                {
					for(int b =0; b<numphase; b++)
                    {
						for(int c = 0; c<numphase; c++)
                        {
							if(b!=a && c!=a && b<c && a<b)
                            {
								if(fabs(div(b))>0.0 && fabs(div(c))>0.0)
                                { 
									phi_phi = phiOld(i,j,k,b)*phiOld(i,j,k,c);
									sum += gamm_abc(a,b,c)*phi_phi;
								}
							}
						}
					}
				}
			    term(i,j,k,a) = sum;
			    phi_phi=0.0;
			    sum=0.0;
			}
		});
	}
}

//####################################################################################################################