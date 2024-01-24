#ifndef ADV_CHEM_POT_H_
#define ADV_CHEM_POT_H_

#include <Function_H.h>
#include <Cal_MOI.h>

using namespace amrex;

//Code for evolution of chemical potential in 2-D -----------------------------------------------------------------------------------------------------
void dmudt_2D(MultiFab& mu_new, MultiFab& mu_old, MultiFab& phi_new, MultiFab& phi_old, MultiFab& comp_new, MultiFab& comp_old, Geometry const& geom)
{
	BL_PROFILE("dmudt()");

	#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif

	for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();								//Defining the box for iteration space
		Array4<Real> const& phiOld = phi_old.array(mfi);				//Taking the Multifabs as array
		Array4<Real> const& phiNew = phi_new.array(mfi);				//Taking the Multifabs as array
		Array4<Real> const& mun = mu_new.array(mfi);					//Taking the Multifabs as array
		Array4<Real> const& muo = mu_old.array(mfi);					//Taking the Multifabs as array
		Array4<Real> const& compn = comp_new.array(mfi);
		Array4<Real> const& compo = comp_old.array(mfi);
		
		//Turning the vector B to a GPU compatible array----------------------------------------------------------
		Array2D <Real,0,phasecount-1,0,compcount-2, Order::C> BB{};
		for(int a=0; a<nump; a++){
			for(int l=0; l<numcom-1; l++){
				BB(a,l) = B[a][l];
			}
		}

		//Turning the vector A to a GPU compatible array-----------------------------------------------------------
		Array3D <Real,0,phasecount-1,0,compcount-2,0,compcount-2, Order::C> AA{};
		for(int a=0; a<nump; a++){
			for(int l=0; l<numcom-1; l++){
				for(int m=0; m<numcom-1; m++){
					AA(a,l,m) = A[a][l][m];
				}
			}
		}

		//Turning the vector dcdmu to a GPU compatible array----------------------------------------------------------
		Array3D <Real,0, phasecount-1, 0, compcount-1, 0, compcount-1, Order::C> der_cmu;
		for(int a=0; a<nump; a++){
			for(int l=0; l<numcom-1; l++){
				for(int m=0; m<numcom-1; m++){
					der_cmu(a,l,m) = dcdmu[a][l][m];
				}
			}
		}

		//Turning the vector diff to a GPU compatible array----------------------------------------------------------
		Array3D <Real,0, phasecount-1, 0, compcount-1, 0, compcount-1, Order::C> diffs;
		for(int a=0; a<nump; a++){
		for(int l=0; l<numcom-1; l++){
				for(int m=0; m<numcom-1; m++){
					diffs(a,l,m) = diff[a][l][m];
				}
			}
		}

		//Redefining variables in GPU space	--------------------------------------------------------------------------
		Real time_step = dt;
		Real epsilon = eps;
		int numphase = nump;
		int numcomp = numcom;
		
		//delta stores dx, dy and dz -----------------------------------------------------------
		GpuArray<Real,AMREX_SPACEDIM> delta = geom.CellSizeArray();

		//Iteration ----------------------------------------------------------------------------
		amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {	
			//Declaring variables --------------------------------------------------------------
			Real sum=0;
			Array2D <Real, 0, phasecount-1, 0, AMREX_SPACEDIM*2, Order::C> norm_x{};
			Array2D <Real, 0, phasecount-1, 0, AMREX_SPACEDIM*2, Order::C> norm_y{};
			Array1D <Real, 0, AMREX_SPACEDIM*2> mod{};
			//Array2D <Real, 0, phasecount-1, 0, AMREX_SPACEDIM*2, Order::C> conct{};
			Array2D <Real, 0, AMREX_SPACEDIM*2, 0, AMREX_SPACEDIM-1, Order::C> gradphi{};
			Array2D <Real, 0, AMREX_SPACEDIM*2, 0, compcount-2, Order::C> dmu{};
			Array3D <Real, 0, AMREX_SPACEDIM*2, 0, compcount-2, 0, compcount-2, Order::C> der_M{};
			Array2D <Real, 0, AMREX_SPACEDIM*2, 0, compcount-2, Order::C> jat{};
			Array1D <Real, 0, compcount-2> divflux{};
			Array1D <Real, 0, compcount-2> divjat{};
			Array1D <Real, 0, compcount-2> sum_fin{};
			Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> c{};
			Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> c_ip{};
			Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> c_im{};
			Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> c_jp{};
			Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> c_jm{};
			Array2D<Real,0,compcount-2,0,compcount-2, Order::C> denom{};
			Array2D<Real,0,compcount-2,0,compcount-2, Order::C> inv_denom{};

			//Computing \frac{\partial \mu}{\partial x} and \frac{\partial \mu}{\partial y} ------------------------
			for(int l=0; l< numcomp-1; l++){
			dmu(cent,l) = 0.0;
			dmu(iph,l) = (muo(i+1,j,k,l)-muo(i,j,k,l))/(delta[X]);
			dmu(imh,l) = (muo(i,j,k,l)-muo(i-1,j,k,l))/(delta[X]);
			dmu(jph,l) = (muo(i,j+1,k,l)-muo(i,j,k,l))/(delta[Y]);
			dmu(jmh,l) = (muo(i,j,k,l)-muo(i,j-1,k,l))/(delta[Y]);
			}

			//Computing M_{ij} (\phi) = \sum_{\alpha =1}^{N-1} M_{ij}^{\alpha} g_{a}(\phi) ------------------------------------
			Cal_derM(i,j,k,der_M, phiOld,der_cmu,diffs,numcomp,numphase);

			//Computing \nabla \cdot \left ( \sum_{j=1}^{K-1} M_{ij}^{\alpha}(\phi) \nabla \mu_{j}\right) -----------------------------------------
			for(int l=0; l<numcomp-1; l++){
				for(int m=0; m<numcomp-1; m++){
					divflux(l) += ((0.5*(der_M(1,l,m)+der_M(0,l,m))*dmu(iph,m)) - (0.5*(der_M(0,l,m)+der_M(2,l,m))*dmu(imh,m)))/delta[X];
					divflux(l) += ((0.5*(der_M(3,l,m)+der_M(0,l,m))*dmu(jph,m)) - (0.5*(der_M(0,l,m)+der_M(4,l,m))*dmu(jmh,m)))/delta[Y];
				}
			} 

			
			//Computing \sum_{a}^{N} c^{a} \cdot \frac{\partial h_{a}}{\partial t} -------------------------------------------------------------------------
			for(int a=0; a<numphase; a++){
				for(int b=0; b<numphase; b++){
					
						sum += dhphi(i,j,k,phiOld,a,b,numphase)*(phiNew(i,j,k,b)-phiOld(i,j,k,b))/(time_step);

				}

				//Computing c^{a} at i,j,k ---------------------------------- 
				c_mu(i,j,k,muo,c,BB,AA,numcomp,a);

					for (int l=0; l<numcomp-1; l++) {
						
						sum_fin(l)  += c(a,l)*sum;
	
					}
				sum=0.0;
			}

			//Calculating \frac{\nabla\phi_{a}}{\left |\nabla\phi_{a} \right |} ---------------------------------------------
			for(int a=0; a<numphase;a++){
				
				gradphi(cent,X) = 0.0;
				gradphi(iph,X)=(phiOld(i+1,j,k,a)-phiOld(i,j,k,a))/delta[X];
            	gradphi(imh,X)=(phiOld(i,j,k,a)-phiOld(i-1,j,k,a))/delta[X];
            	gradphi(jph,X)=(phiOld(i+1,j+1,k,a)-phiOld(i-1,j+1,k,a)+phiOld(i+1,j,k,a)-phiOld(i-1,j,k,a))/(4*delta[X]);
            	gradphi(jmh,X)=(phiOld(i+1,j,k,a)-phiOld(i-1,j,k,a)+phiOld(i+1,j-1,k,a)-phiOld(i-1,j-1,k,a))/(4*delta[X]);

				gradphi(cent,Y) = 0.0;
            	gradphi(iph,Y)=(phiOld(i+1,j+1,k,a)-phiOld(i+1,j-1,k,a)+phiOld(i,j+1,k,a)-phiOld(i,j-1,k,a))/(4*delta[Y]);
            	gradphi(imh,Y)=(phiOld(i,j+1,k,a)-phiOld(i,j-1,k,a)+phiOld(i-1,j+1,k,a)-phiOld(i-1,j-1,k,a))/(4*delta[Y]);
            	gradphi(jph,Y)=(phiOld(i,j+1,k,a)-phiOld(i,j,k,a))/delta[Y];
            	gradphi(jmh,Y)=(phiOld(i,j,k,a)-phiOld(i,j-1,k,a))/delta[Y];


				for(int p=0; p<mod.len();p++){
					mod(p) = sqrt(pow(gradphi(p,X),2)+pow(gradphi(p,Y),2)); 
				}			

				for(int p=0; p<AMREX_SPACEDIM*2+1;p++){
					if(mod(p)>1e-15 && p!=0){
						norm_x(a,p) = gradphi(p,X)/(mod(p));
						norm_y(a,p) = gradphi(p,Y)/(mod(p));
					}
					else{
						norm_x(a,p) = 0.0;
						norm_y(a,p) = 0.0;
					}	
				}
			}

			//Calculating \sum_{a=1}^{N} \left (  j_{at}^{a\rightarrow l}\right )_{i}\left ( -\frac{\nabla\phi_{a}}{\left | \nabla\phi_{a} \right |} \cdot \frac{\nabla\phi_{liq}}{\left | \nabla\phi_{liq} \right |} \right ) ------------------------------------------------------------------------------------
			for(int a=0; a<numphase-1;a++){
					
					//Computing c^{a} at i+1/2, i-1/2, j+1/2, j-1/2 ----------------------
					c_mu(i+1,j,k,muo,c_ip,BB,AA,numcomp,a);
					c_mu(i-1,j,k,muo,c_im,BB,AA,numcomp,a);
					c_mu(i,j+1,k,muo,c_jp,BB,AA,numcomp,a);
					c_mu(i,j-1,k,muo,c_jm,BB,AA,numcomp,a);

					//Computing c^{liq} at i+1/2, i-1/2, j+1/2, j-1/2 ----------------------
					c_mu(i+1,j,k,muo,c_ip,BB,AA,numcomp,numphase-1);
					c_mu(i-1,j,k,muo,c_im,BB,AA,numcomp,numphase-1);
					c_mu(i,j+1,k,muo,c_jp,BB,AA,numcomp,numphase-1);
					c_mu(i,j-1,k,muo,c_jm,BB,AA,numcomp,numphase-1);

				
				for(int l=0; l<numcomp-1; l++){
					
					jat(iph,l) += (1.0-diffs(a,l,l)/diffs(numphase-1,l,l))*(0.5/sqrt(2))*0.5*((c_ip(numphase-1,l)-c_ip(a,l))*(phiNew(i+1,j,k,a)-phiOld(i+1,j,k,a)) + (c(numphase-1,l)-c(a,l))*(phiNew(i,j,k,a)-phiOld(i,j,k,a)))*(norm_x(a,iph))*fabs(norm_x(a,iph)*norm_x(numphase-1,iph)+norm_y(a,iph)*norm_y(numphase-1,iph));
					jat(imh,l) += (1.0-diffs(a,l,l)/diffs(numphase-1,l,l))*(0.5/sqrt(2))*0.5*((c(numphase-1,l)-c(a,l))*(phiNew(i,j,k,a)-phiOld(i,j,k,a)) + (c_im(numphase-1,l)-c_im(a,l))*(phiNew(i-1,j,k,a)-phiOld(i-1,j,k,a)))*(norm_x(a,imh))*fabs(norm_x(a,imh)*norm_x(numphase-1,imh)+norm_y(a,imh)*norm_y(numphase-1,imh));
						
					jat(jph,l) += (1.0-diffs(a,l,l)/diffs(numphase-1,l,l))*(0.5/sqrt(2))*0.5*((c_jp(numphase-1,l)-c_jp(a,l))*(phiNew(i,j+1,k,a)-phiOld(i,j+1,k,a)) + (c(numphase-1,l)-c(a,l))*(phiNew(i,j,k,a)-phiOld(i,j,k,a)))*(norm_y(a,jph))*fabs(norm_x(a,jph)*norm_x(numphase-1,jph)+norm_y(a,jph)*norm_y(numphase-1,jph));
					jat(jmh,l) += (1.0-diffs(a,l,l)/diffs(numphase-1,l,l))*(0.5/sqrt(2))*0.5*((c(numphase-1,l)-c(a,l))*(phiNew(i,j,k,a)-phiOld(i,j,k,a)) + (c_jm(numphase-1,l)-c_jm(a,l))*(phiNew(i,j-1,k,a)-phiOld(i,j-1,k,a)))*(norm_y(a,jmh))*fabs(norm_x(a,jmh)*norm_x(numphase-1,jmh)+norm_y(a,jmh)*norm_y(numphase-1,jmh));
				}
			}

			//Calculating \frac{-\epsilon}{\Delta t} \left (\frac{\partial j_{at}^{x}}{\partial x}  + \frac{\partial j_{at}^{y}}{\partial y} \right ) ---------------------------
			for(int l=0; l<numcomp-1;l++){
				divjat(l) = (jat(iph,l)-jat(imh,l))*(-1.0*epsilon)/(time_step*delta[X]) + (jat(jph,l)-jat(jmh,l))*(-1.0*epsilon)/(time_step*delta[Y]);
			}
			

			//Calculating \sum_{a}^{N} h_{a}(\phi) \cdot \frac{\partial c_{a}}{\partial \mu} -------------------------------------------------------------
			for(int l=0; l<numcomp-1; l++){
				for(int m=0; m<numcomp-1; m++){
					for(int a=0; a<numphase; a++){
							denom(l,m) += der_cmu(a,l,m)*hphi(i,j,k,phiOld,a,numphase);
					}
				}
			}
			
			if(numcomp==2){
				inv_denom(0,0) = 1.0/denom(0,0); 

				mun(i,j,k,0) = muo(i,j,k,0) + time_step*(divflux(0) -sum_fin(0)- divjat(0))*inv_denom(0,0);

				//compn(i,j,k,0) = compo(i,j,k,0) + time_step*(divflux(0)-divjat(0));
			}

			if(numcomp==3){
				Real DET = denom(0,0)*denom(1,1) - denom(0,1)*denom(1,0);
				inv_denom(0,0)  = denom(1,1)/DET;
				inv_denom(1,1)  = denom(0,0)/DET;
				inv_denom(0,1)  = -denom(0,1)/DET;
				inv_denom(1,0)  = -denom(1,0)/DET;

				for (int l=0; l < numcomp-1; l++ ) {
					//compn(i,j,k,l) = compo(i,j,k,l) + time_step*(divflux(l)-divjat(l));
          			for (int m=0; m < numcomp-1; m++) {
            			mun(i,j,k,l) = muo(i,j,k,l) + time_step*(divflux(m) - sum_fin(m) - divjat(m))*inv_denom(l,m);
          			}
        		}
			}
		});
	}
}


//Code for evolution of chemical potential in 3-D -------------------------------------------------------------------------------------------------------
void dmudt_3D(MultiFab& mu_new, MultiFab& mu_old, MultiFab& phi_new, MultiFab& phi_old, MultiFab& comp_new, MultiFab& comp_old, Geometry const& geom)
{
	BL_PROFILE("dmudt()");

	#ifdef AMREX_USE_OMP
		#pragma omp parallel if (Gpu::notInLaunchRegion())
	#endif

	for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();										//Defining the box for iteration space
		Array4<Real> const& phiOld = phi_old.array(mfi);						//Taking the Multifabs as array
		Array4<Real> const& phiNew = phi_new.array(mfi);						//Taking the Multifabs as array
		Array4<Real> const& mun = mu_new.array(mfi);							//Taking the Multifabs as array
		Array4<Real> const& muo = mu_old.array(mfi);							//Taking the Multifabs as array
		//Array4<Real> const& compn = comp_new.array(mfi);
		//Array4<Real> const& compo = comp_old.array(mfi);
		
		//Turning the vector B to a GPU compatible array----------------------------------------------------------
		Array2D <Real,0,phasecount-1,0,compcount-2, Order::C> BB{};
		for(int a=0; a<nump; a++){
			for(int l=0; l<numcom-1; l++){
				BB(a,l) = B[a][l];
			}
		}
		
		//Turning the vector A to a GPU compatible array-----------------------------------------------------------
		Array3D <Real,0,phasecount-1,0,compcount-2,0,compcount-2, Order::C> AA{};
		for(int a=0; a<nump; a++){
			for(int l=0; l<numcom-1; l++){
				for(int m=0; m<numcom-1; m++){
					AA(a,l,m) = A[a][l][m];
				}
			}
		}

		//Turning the vector dcdmu to a GPU compatible array----------------------------------------------------------
		Array3D <Real,0, phasecount-1, 0, compcount-1, 0, compcount-1, Order::C> der_cmu;
		for(int a=0; a<nump; a++){
			for(int l=0; l<numcom-1; l++){
				for(int m=0; m<numcom-1; m++){
					der_cmu(a,l,m) = dcdmu[a][l][m];
				}
			}
		}

		//Turning the vector diff to a GPU compatible array----------------------------------------------------------
		Array3D <Real,0, phasecount-1, 0, compcount-1, 0, compcount-1, Order::C> diffs;
		for(int a=0; a<nump; a++){
		for(int l=0; l<numcom-1; l++){
				for(int m=0; m<numcom-1; m++){
					diffs(a,l,m) = diff[a][l][m];
				}
			}
		}

		//Redefining variables in GPU space	--------------------------------------------------------------------------
		Real time_step = dt;
		Real epsilon = eps;
		int numphase = nump;
		int numcomp = numcom;
		
		//delta stores dx, dy and dz -----------------------------------------------------------
		GpuArray<Real,AMREX_SPACEDIM> delta = geom.CellSizeArray();
	

		//Iteration --------------------------------------------------------------------------
		amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {	
			//Declaring variables ------------------------------------------------------------
			Real sum=0;
			Array2D <Real, 0, phasecount-1,0, AMREX_SPACEDIM*2,Order::C> norm_x{};
			Array2D <Real, 0, phasecount-1,0, AMREX_SPACEDIM*2,Order::C> norm_y{};
			Array2D <Real, 0, phasecount-1,0, AMREX_SPACEDIM*2,Order::C> norm_z{};
			Array1D <Real, 0, AMREX_SPACEDIM*2> mod{};
			Array2D <Real, 0, phasecount-1,0, AMREX_SPACEDIM*2,Order::C> conct{};
			Array2D <Real, 0, AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> gradphi{};
			Array2D <Real, 0, AMREX_SPACEDIM*2, 0, compcount-2, Order::C> dmu{};
			Array3D <Real, 0, AMREX_SPACEDIM*2, 0, compcount-2, 0, compcount-2, Order::C> der_M{};
			Array2D <Real, 0, AMREX_SPACEDIM*2, 0, compcount-2, Order::C> jat{};
			Array1D <Real, 0, compcount-2> divflux{};
			Array1D <Real, 0, compcount-2> divjat{};
			Array1D <Real, 0, compcount-2> sum_fin{};
			Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> c{};
			Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> c_ip{};
			Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> c_im{};
			Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> c_jp{};
			Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> c_jm{};
			Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> c_kp{};
			Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> c_km{};
			Array2D<Real,0,compcount-2,0,compcount-2, Order::C> denom{};
			Array2D<Real,0,compcount-2,0,compcount-2, Order::C> inv_denom{};

			//Computing \frac{\partial \mu}{\partial x}, \frac{\partial \mu}{\partial y} and \frac{\partial \mu}{\partial z} ------------------------
			for(int l=0; l< numcomp-1; l++){
			dmu(cent,l) = 0.0;
			dmu(iph,l) = (muo(i+1,j,k,l)-muo(i,j,k,l))/(delta[X]);
			dmu(imh,l) = (muo(i,j,k,l)-muo(i-1,j,k,l))/(delta[X]);
			dmu(jph,l) = (muo(i,j+1,k,l)-muo(i,j,k,l))/(delta[Y]);
			dmu(jmh,l) = (muo(i,j,k,l)-muo(i,j-1,k,l))/(delta[Y]);
			dmu(kph,l) = (muo(i,j,k+1,l)-muo(i,j,k,l))/(delta[Z]);
			dmu(kmh,l) = (muo(i,j,k,l)-muo(i,j,k-1,l))/(delta[Z]);
			}
			
			//Computing M_{ij} (\phi) = \sum_{\alpha =1}^{N-1} M_{ij}^{\alpha} g_{a}(\phi) ------------------------------------
			Cal_derM(i,j,k,der_M, phiOld,der_cmu,diffs,numcomp,numphase);

			//Computing \nabla \cdot \left ( \sum_{j=1}^{K-1} M_{ij}^{\alpha}(\phi) \nabla \mu_{j}\right) -----------------------------------------
			for(int l=0; l<numcomp-1; l++){
				for(int m=0; m<numcomp-1; m++){
					divflux(l) += ((0.5*(der_M(1,l,m)+der_M(0,l,m))*dmu(iph,m)) - (0.5*(der_M(0,l,m)+der_M(2,l,m))*dmu(imh,m)))/delta[X];
					divflux(l) += ((0.5*(der_M(3,l,m)+der_M(0,l,m))*dmu(jph,m)) - (0.5*(der_M(0,l,m)+der_M(4,l,m))*dmu(jmh,m)))/delta[Y];
					divflux(l) += ((0.5*(der_M(5,l,m)+der_M(0,l,m))*dmu(kph,m)) - (0.5*(der_M(0,l,m)+der_M(6,l,m))*dmu(kmh,m)))/delta[Z];
				}
			} 
			

			//Computing \sum_{a}^{N} c^{a} \cdot \frac{\partial h_{a}}{\partial t} -------------------------------------------------------------------------
			for(int a=0; a<numphase; a++){
				for(int b=0; b<numphase; b++){
					
						sum += dhphi(i,j,k,phiOld,a,b,numphase)*(phiNew(i,j,k,b)-phiOld(i,j,k,b))/(time_step);

				}

				//Computing c^{a} at i,j,k ---------------------------------- 
				c_mu(i,j,k,muo,c,BB,AA,numcomp,a);

					for (int l=0; l<numcomp-1; l++) {
						
						sum_fin(l)  += c(a,l)*sum;
	
					}
				sum=0.0;
			}


			//Calculating \frac{\nabla\phi_{a}}{\left |\nabla\phi_{a} \right |} ---------------------------------------------
			for(int a=0; a<numphase;a++){
				
				gradphi(cent,X) = 0.0;
				gradphi(iph,X)=(phiOld(i+1,j,k,a)-phiOld(i,j,k,a))/delta[X];
            	gradphi(imh,X)=(phiOld(i,j,k,a)-phiOld(i-1,j,k,a))/delta[X];
            	gradphi(jph,X)=(phiOld(i+1,j+1,k,a)-phiOld(i-1,j+1,k,a)+phiOld(i+1,j,k,a)-phiOld(i-1,j,k,a))/(4*delta[X]);
            	gradphi(jmh,X)=(phiOld(i+1,j,k,a)-phiOld(i-1,j,k,a)+phiOld(i+1,j-1,k,a)-phiOld(i-1,j-1,k,a))/(4*delta[X]);
				gradphi(kph,X)=(phiOld(i+1,j,k+1,a)-phiOld(i-1,j,k+1,a)+phiOld(i+1,j,k,a)-phiOld(i-1,j,k,a))/(4*delta[X]);
				gradphi(kmh,X)=(phiOld(i+1,j,k,a)-phiOld(i-1,j,k,a)+phiOld(i+1,j,k-1,a)-phiOld(i-1,j,k-1,a))/(4*delta[X]);

				gradphi(cent,Y) = 0.0;
            	gradphi(iph,Y)=(phiOld(i+1,j+1,k,a)-phiOld(i+1,j-1,k,a)+phiOld(i,j+1,k,a)-phiOld(i,j-1,k,a))/(4*delta[Y]);
            	gradphi(imh,Y)=(phiOld(i,j+1,k,a)-phiOld(i,j-1,k,a)+phiOld(i-1,j+1,k,a)-phiOld(i-1,j-1,k,a))/(4*delta[Y]);
            	gradphi(jph,Y)=(phiOld(i,j+1,k,a)-phiOld(i,j,k,a))/delta[Y];
            	gradphi(jmh,Y)=(phiOld(i,j,k,a)-phiOld(i,j-1,k,a))/delta[Y];
				gradphi(kph,Y)=(phiOld(i,j+1,k+1,a)-phiOld(i,j-1,k+1,a)+phiOld(i,j+1,k,a)-phiOld(i,j-1,k,a))/(4*delta[Y]);
				gradphi(kmh,Y)=(phiOld(i,j+1,k,a)-phiOld(i,j-1,k,a)+phiOld(i,j+1,k-1,a)-phiOld(i,j-1,k-1,a))/(4*delta[Y]);

				gradphi(cent,Z) = 0.0;
				gradphi(iph,Z)=(phiOld(i+1,j,k+1,a)-phiOld(i+1,j,k-1,a)+phiOld(i,j,k+1,a)-phiOld(i,j,k-1,a))/(4*delta[Z]);
				gradphi(imh,Z)=(phiOld(i,j,k+1,a)-phiOld(i,j,k-1,a)+phiOld(i-1,j,k+1,a)-phiOld(i-1,j,k-1,a))/(4*delta[Z]);
				gradphi(jph,Z)=(phiOld(i,j+1,k+1,a)-phiOld(i,j+1,k-1,a)+phiOld(i,j,k+1,a)-phiOld(i,j,k-1,a))/(4*delta[Z]);
				gradphi(jmh,Z)=(phiOld(i,j,k+1,a)-phiOld(i,j,k-1,a)+phiOld(i,j-1,k+1,a)-phiOld(i,j-1,k-1,a))/(4*delta[Z]);
				gradphi(kph,Z)=(phiOld(i,j,k+1,a)-phiOld(i,j,k,a))/delta[Z];
            	gradphi(kmh,Z)=(phiOld(i,j,k,a)-phiOld(i,j,k-1,a))/delta[Z];


				for(int p=0; p<mod.len();p++){
					#if(AMREX_SPACEDIM>2)
					mod(p) = sqrt(pow(gradphi(p,X),2)+pow(gradphi(p,Y),2)+pow(gradphi(p,Z),2));
					#else
					mod(p) = sqrt(pow(gradphi(p,X),2)+pow(gradphi(p,Y),2));
					#endif
				}			

				for(int p=0; p<AMREX_SPACEDIM*2+1;p++){
					if(mod(p)>1e-15 && p!=0){
						norm_x(a,p) = gradphi(p,X)/(mod(p));
						norm_y(a,p) = gradphi(p,Y)/(mod(p));
						norm_z(a,p) = gradphi(p,Z)/(mod(p));
					}
					else{
						norm_x(a,p) = 0.0;
						norm_y(a,p) = 0.0;
						norm_z(a,p) = 0.0;
					}	
				}
			}

			//Calculating \sum_{a=1}^{N} \left (  j_{at}^{a\rightarrow l}\right )_{i}\left ( -\frac{\nabla\phi_{a}}{\left | \nabla\phi_{a} \right |} \cdot \frac{\nabla\phi_{liq}}{\left | \nabla\phi_{liq} \right |} \right ) ------------------------------------------------------------------------------------
			for(int a=0; a<numphase-1;a++){
					
					//Computing c^{a} at i+1/2, i-1/2, j+1/2, j-1/2, k+1/2, k-1/2 ----------------------
					c_mu(i+1,j,k,muo,c_ip,BB,AA,numcomp,a);
					c_mu(i-1,j,k,muo,c_im,BB,AA,numcomp,a);
					c_mu(i,j+1,k,muo,c_jp,BB,AA,numcomp,a);
					c_mu(i,j-1,k,muo,c_jm,BB,AA,numcomp,a);
					c_mu(i,j,k+1,muo,c_kp,BB,AA,numcomp,a);
					c_mu(i,j,k-1,muo,c_km,BB,AA,numcomp,a);

					//Computing c^{liq} at i+1/2, i-1/2, j+1/2, j-1/2, k+1/2, k-1/2 ----------------------
					c_mu(i+1,j,k,muo,c_ip,BB,AA,numcomp,numphase-1);
					c_mu(i-1,j,k,muo,c_im,BB,AA,numcomp,numphase-1);
					c_mu(i,j+1,k,muo,c_jp,BB,AA,numcomp,numphase-1);
					c_mu(i,j-1,k,muo,c_jm,BB,AA,numcomp,numphase-1);
					c_mu(i,j,k+1,muo,c_kp,BB,AA,numcomp,numphase-1);
					c_mu(i,j,k-1,muo,c_km,BB,AA,numcomp,numphase-1);

				
				for(int l=0; l<numcomp-1; l++){
					
					jat(iph,l) += (1.0-diffs(a,l,l)/diffs(numphase-1,l,l))*(0.5/sqrt(2))*0.5*((c_ip(numphase-1,l)-c_ip(a,l))*(phiNew(i+1,j,k,a)-phiOld(i+1,j,k,a)) + (c(numphase-1,l)-c(a,l))*(phiNew(i,j,k,a)-phiOld(i,j,k,a)))*(norm_x(a,iph))*fabs(norm_x(a,iph)*norm_x(numphase-1,iph)+norm_y(a,iph)*norm_y(numphase-1,iph)+norm_z(a,iph)*norm_z(numphase-1,iph));
					jat(imh,l) += (1.0-diffs(a,l,l)/diffs(numphase-1,l,l))*(0.5/sqrt(2))*0.5*((c(numphase-1,l)-c(a,l))*(phiNew(i,j,k,a)-phiOld(i,j,k,a)) + (c_im(numphase-1,l)-c_im(a,l))*(phiNew(i-1,j,k,a)-phiOld(i-1,j,k,a)))*(norm_x(a,imh))*fabs(norm_x(a,imh)*norm_x(numphase-1,imh)+norm_y(a,imh)*norm_y(numphase-1,imh)+norm_z(a,imh)*norm_z(numphase-1,imh));
						
					jat(jph,l) += (1.0-diffs(a,l,l)/diffs(numphase-1,l,l))*(0.5/sqrt(2))*0.5*((c_jp(numphase-1,l)-c_jp(a,l))*(phiNew(i,j+1,k,a)-phiOld(i,j+1,k,a)) + (c(numphase-1,l)-c(a,l))*(phiNew(i,j,k,a)-phiOld(i,j,k,a)))*(norm_y(a,jph))*fabs(norm_x(a,jph)*norm_x(numphase-1,jph)+norm_y(a,jph)*norm_y(numphase-1,jph)+norm_z(a,jph)*norm_z(numphase-1,jph));
					jat(jmh,l) += (1.0-diffs(a,l,l)/diffs(numphase-1,l,l))*(0.5/sqrt(2))*0.5*((c(numphase-1,l)-c(a,l))*(phiNew(i,j,k,a)-phiOld(i,j,k,a)) + (c_jm(numphase-1,l)-c_jm(a,l))*(phiNew(i,j-1,k,a)-phiOld(i,j-1,k,a)))*(norm_y(a,jmh))*fabs(norm_x(a,jmh)*norm_x(numphase-1,jmh)+norm_y(a,jmh)*norm_y(numphase-1,jmh)+norm_z(a,jmh)*norm_z(numphase-1,jmh));

					jat(kph,l) += (1.0-diffs(a,l,l)/diffs(numphase-1,l,l))*(0.5/sqrt(2))*0.5*((c_kp(numphase-1,l)-c_kp(a,l))*(phiNew(i,j,k+1,a)-phiOld(i,j,k+1,a)) + (c(numphase-1,l)-c(a,l))*(phiNew(i,j,k,a)-phiOld(i,j,k,a)))*(norm_z(a,kph))*fabs(norm_x(a,kph)*norm_x(numphase-1,kph)+norm_y(a,kph)*norm_y(numphase-1,kph)+norm_z(a,kph)*norm_z(numphase-1,kph));
					jat(kmh,l) += (1.0-diffs(a,l,l)/diffs(numphase-1,l,l))*(0.5/sqrt(2))*0.5*((c(numphase-1,l)-c(a,l))*(phiNew(i,j,k,a)-phiOld(i,j,k,a)) + (c_km(numphase-1,l)-c_km(a,l))*(phiNew(i,j,k-1,a)-phiOld(i,j,k-1,a)))*(norm_z(a,kmh))*fabs(norm_x(a,kmh)*norm_x(numphase-1,kmh)+norm_y(a,kmh)*norm_y(numphase-1,kmh)+norm_z(a,kmh)*norm_z(numphase-1,kmh));
				}
			}

			//Calculating \frac{-\epsilon}{\Delta t} \left (\frac{\partial j_{at}^{x}}{\partial x}  + \frac{\partial j_{at}^{y}}{\partial y} \right ) ---------------------------
			for(int l=0; l<numcomp-1;l++){
				divjat(l) = (jat(iph,l)-jat(imh,l))*(-1.0*epsilon)/(time_step*delta[X]) + (jat(jph,l)-jat(jmh,l))*(-1.0*epsilon)/(time_step*delta[Y]) + (jat(kph,l)-jat(kmh,l))*(-1.0*epsilon)/(time_step*delta[Z]);
			}

			//Calculating \sum_{a}^{N} h_{a}(\phi) \cdot \frac{\partial c_{a}}{\partial \mu} -------------------------------------------------------------
			for(int l=0; l<numcomp-1; l++){
				for(int m=0; m<numcomp-1; m++){
					for(int a=0; a<numphase; a++){
							denom(l,m) += der_cmu(a,l,m)*hphi(i,j,k,phiOld,a,numphase);
					}
				}
			}

			if(numcomp==2){
				inv_denom(0,0) = 1.0/denom(0,0); 
			
				mun(i,j,k,0) = muo(i,j,k,0) + time_step*(divflux(0) -sum_fin(0)- divjat(0))*inv_denom(0,0);

				//compn(i,j,k,0) = compo(i,j,k,0) + time_step*(divflux(0)-divjat(0));
			}

			if(numcomp==3){
				Real DET = denom(0,0)*denom(1,1) - denom(0,1)*denom(1,0);
				inv_denom(0,0)  = denom(1,1)/DET;
				inv_denom(1,1)  = denom(0,0)/DET;
				inv_denom(0,1)  = -denom(0,1)/DET;
				inv_denom(1,0)  = -denom(1,0)/DET;

				for (int l=0; l < numcomp-1; l++ ) {
					//compn(i,j,k,l) = compo(i,j,k,l) + time_step*(divflux(l)-divjat(l));
          			for (int m=0; m < numcomp-1; m++) {
            			mun(i,j,k,l) = muo(i,j,k,l) + time_step*(divflux(m) - sum_fin(m) - divjat(m))*inv_denom(l,m);
          			}
        		}
			}
		});
	}
}


#endif