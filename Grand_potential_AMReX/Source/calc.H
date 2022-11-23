#ifndef _CALC_H_
#define _CALC_H_

#include<AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include "Tau.H"

#define iph 0
#define imh 1
#define jph 2
#define jmh 3
#define kph 4
#define kmh 5

using namespace amrex;

///////////////////////////////////////////////////////////////////////////////////////////////////////

//AMREX_GPU_DEVICE AMREX_FORCE_INLINE

void update_phi(MultiFab& phi_new, MultiFab& phi_old, MultiFab& term1, MultiFab& term2, MultiFab& term3)
{	

	for (MFIter mfi(phi_old); mfi.isValid(); ++mfi)
	{
		const Box& dbx = mfi.validbox();
		Array4<Real> const& fin_term1 = term1.array(mfi);
		Array4<Real> const& fin_term2 = term2.array(mfi);
		Array4<Real> const& fin_term3 = term3.array(mfi);
		Array4<Real> const& phiNew = phi_new.array(mfi);
		Array4<Real> const& phiOld = phi_old.array(mfi);
		Real tauu = tau_final;
		Real time_step = dt;
		Real epsilon = eps;
		Real molar_vol = Vm; 

		amrex::ParallelFor( dbx, 
		[=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
		{
		
			phiNew(i,j,k) = phiOld(i,j,k) + (time_step/tauu)*(fin_term1(i,j,k)-fin_term2(i,j,k)/(pow(epsilon,2)) - fin_term3(i,j,k)/(epsilon*molar_vol));
		 	
		});

		
	}
	
}
	
///////////////////////////////////////////////////////////////////////////////////////////////////////

void dmudt(MultiFab& mu_new, MultiFab& mu_old, MultiFab& phi_new, MultiFab& phi_old, MultiFab& derivx, MultiFab& derivy, MultiFab& h_phi, MultiFab& dh_dphi, Geometry const& geom)
{
	GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

	for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();
		Array4<Real> const& phiold = phi_old.array(mfi);
		Array4<Real> const& phinew = phi_new.array(mfi);
		Array4<Real> const& mun = mu_new.array(mfi);
		Array4<Real> const& muo = mu_old.array(mfi);
		Array4<Real> const& derivxx = derivx.array(mfi);
        Array4<Real> const& derivyy = derivy.array(mfi);
		Array4<Real> const& arr_h_phi = h_phi.array(mfi);
		Array4<Real> const& arr_dh_dphi = dh_dphi.array(mfi);
		Array1D <Real,0,1> BB;
		BB(0)=B[0];
		BB(1)=B[1];
		Array1D <Real,0,1> der_cmu;
		der_cmu(0) = dcdmu[0];
		der_cmu(1) = dcdmu[1];

		Array2D <Real,0,1,0,2,Order::C> flag;
		flag(0,0)=diff[0][0];
		flag(0,1)=diff[0][1];
		flag(0,2)=diff[0][2];
		flag(1,0)=diff[1][0];
		flag(1,1)=diff[1][1];
		flag(1,2)=diff[1][2];
		Real time_step = dt;
		Real epsilon = eps;
		
	
		amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {	
			if(phiold(i,j,k)==1 || phiold(i,j,k)==0){
					mun(i,j,k) = muo(i,j,k);
			}

			else{
			Real calpha = (muo(i,j,k) - (BB(0)))*(der_cmu(0));

			Real cbeta = muo(i,j,k)*(der_cmu(1));

			Real dmudx_iph = (muo(i+1,j,k)-muo(i,j,k))/(dx[0]);
			Real dmudx_imh = (muo(i,j,k)-muo(i-1,j,k))/(dx[0]);
			Real dmudy_jph = (muo(i,j+1,k)-muo(i,j,k))/(dx[1]);
			Real dmudy_jmh = (muo(i,j,k)-muo(i,j-1,k))/(dx[1]);

			Real modphisq_iph = pow(derivxx(i,j,k,0),2)+pow(derivyy(i,j,k,0),2);
			Real modphisq_imh = pow(derivxx(i,j,k,1),2)+pow(derivyy(i,j,k,1),2);
			Real modphisq_jph = pow(derivxx(i,j,k,2),2)+pow(derivyy(i,j,k,2),2);
			Real modphisq_jmh = pow(derivxx(i,j,k,3),2)+pow(derivyy(i,j,k,3),2);

			Real jatx_iph, jatx_imh, jaty_jph, jaty_jmh;

			if(modphisq_iph>1e-15){
			jatx_iph = (0.5/sqrt(2))*epsilon*(cbeta-calpha)*(((phinew(i+1,j,k)+phinew(i,j,k))*0.5 - (phiold(i+1,j,k)+phiold(i,j,k))*0.5)/time_step)*(derivxx(i,j,k,0)/sqrt(modphisq_iph));
			}
			else{
				jatx_iph = 0.0;
			}

			if(modphisq_imh>1e-15){
			jatx_imh = (0.5/sqrt(2))*epsilon*(cbeta-calpha)*(((phinew(i-1,j,k)+phinew(i,j,k))*0.5 - (phiold(i-1,j,k)+phiold(i,j,k))*0.5)/time_step)*(derivxx(i,j,k,1)/sqrt(modphisq_imh));
			}
			else{
				jatx_imh = 0.0;
			}

			if(modphisq_jph>1e-15){
			jaty_jph = (0.5/sqrt(2))*epsilon*(cbeta-calpha)*(((phinew(i,j+1,k)+phinew(i,j,k))*0.5 - (phiold(i,j+1,k)+phiold(i,j,k))*0.5)/time_step)*(derivyy(i,j,k,2)/sqrt(modphisq_jph));
			}
			else{
				jaty_jph = 0.0;
			}

			if(modphisq_jmh>1e-15){
			jaty_jmh = (0.5/sqrt(2))*epsilon*(cbeta-calpha)*(((phinew(i,j-1,k)+phinew(i,j,k))*0.5 - (phiold(i,j-1,k)+phiold(i,j,k))*0.5)/time_step)*(derivyy(i,j,k,3)/sqrt(modphisq_jmh));
			}
			else{
				jaty_jmh = 0.0;
			}

			Real dbdx = (((flag(0,2)*0.5*(phiold(i,j,k)+phiold(i+1,j,k))*der_cmu(0) + flag(1,2)*(1.0 - 0.5*(phiold(i,j,k)+phiold(i+1,j,k)))*der_cmu(1))*dmudx_iph - jatx_iph)
					-	 ((flag(0,2)*0.5*(phiold(i,j,k)+phiold(i-1,j,k))*der_cmu(0) + flag(1,2)*(1.0 - 0.5*(phiold(i,j,k)+phiold(i-1,j,k)))*der_cmu(1))*dmudx_imh - jatx_imh))/dx[0];

			Real dbdy = (((flag(0,2)*0.5*(phiold(i,j,k)+phiold(i,j+1,k))*der_cmu(0) + flag(1,2)*(1.0 - 0.5*(phiold(i,j,k)+phiold(i,j+1,k)))*der_cmu(1))*dmudy_jph - jaty_jph)
					-	 ((flag(0,2)*0.5*(phiold(i,j,k)+phiold(i,j-1,k))*der_cmu(0) + flag(1,2)*(1.0 - 0.5*(phiold(i,j,k)+phiold(i,j-1,k)))*der_cmu(1))*dmudy_jmh - jaty_jmh))/dx[1];


//			Real cdhdt = (calpha - cbeta)*(6.0*phiold(i,j,k)*(1.0-phiold(i,j,k)))*(phinew(i,j,k) - phiold(i,j,k))/dt;
			Real cdhdt = (calpha - cbeta)*(arr_dh_dphi(i,j,k))*(phinew(i,j,k) - phiold(i,j,k))/time_step;

//			Real coeffdmudt = pow(phiold(i,j,k),2)*(3.0 - 2.0*phiold(i,j,k))*dcdmu[0] + (1.0 - pow(phiold(i,j,k),2)*(3.0 - 2.0*phiold(i,j,k)))*dcdmu[1];
			Real coeffdmudt = arr_h_phi(i,j,k)*der_cmu(0) + (1.0 - arr_h_phi(i,j,k))*der_cmu(1);

			Real dmudt = (dbdx + dbdy - cdhdt)/coeffdmudt;

			mun(i,j,k) = muo(i,j,k) + dmudt*time_step;
			}
		});
	
	}
}





/*

//Function derivative
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void derivative_x(int i, int j, int k,
		amrex::Array4<amrex::Real const> const& phiOld,
		amrex::Array4<amrex::Real      > const& deriv,
		Geometry const& geom)
{
		GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();
		//Two sided difference
        deriv(i,j,k,iph) = (phiOld(i+1,j,k,0)-phiOld(i,j,k,0))/(dx[0]);
		deriv(i,j,k,imh) = (phiOld(i,j,k,0)-phiOld(i-1,j,k,0))/(dx[0]);
		deriv(i,j,k,jph) = 0.5*((phiOld(i+1,j+1,k,0)-phiOld(i-1,j+1,k,0))/(2.0*dx[0])+(phiOld(i+1,j,k,0)-phiOld(i-1,j,k,0))/(2.0*dx[0]));	//dPhi/dy - stored as component 1
		deriv(i,j,k,jmh) = 0.5*((phiOld(i+1,j,k,0)-phiOld(i-1,j,k,0))/(2.0*dx[0])+(phiOld(i+1,j-1,k,0)-phiOld(i-1,j-1,k,0))/(2.0*dx[0]));	//dPhi/dy - stored as component 1
		#if (AMREX_SPACEDIM > 2)
		deriv(i,j,k,kph) = 0.5*((phiOld(i+1,j,k+1,0)-phiOld(i-1,j,k+1,0))/(2.0*dx[0])+(phiOld(i+1,j,k,0)-phiOld(i-1,j,k,0))/(2.0*dx[0]));
		deriv(i,j,k,kmh) = 0.5*((phiOld(i+1,j,k,0)-phiOld(i-1,j,k,0))/(2.0*dx[0])+(phiOld(i+1,j,k-1,0)-phiOld(i-1,j,k-1,0))/(2.0*dx[0]));
		#endif

	//	Print(3)<<"("<<i<<","<<j<<","<<k<<") ";
	//	Print(3)<<deriv(i,j,k,iph)<<"\n";
		// AllPrint()<<deriv(i,j,k,imh)<<"\n";
		//Print(3)<<phiOld(i,j,k)<<"\n";
}

void derivative_y(int i, int j, int k,
		amrex::Array4<amrex::Real const> const& phiOld,
		amrex::Array4<amrex::Real      > const& deriv,
		Geometry const& geom)
{
		GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();
		//Two sided difference
        deriv(i,j,k,jph) = (phiOld(i,j+1,k,0)-phiOld(i,j,k,0))/(dx[1]);
		deriv(i,j,k,jmh) = (phiOld(i,j,k,0)-phiOld(i,j-1,k,0))/(dx[1]);
		deriv(i,j,k,iph) = 0.5*((phiOld(i+1,j+1,k,0)-phiOld(i+1,j-1,k,0))/(2.0*dx[1])+(phiOld(i,j+1,k,0)-phiOld(i,j-1,k,0))/(2.0*dx[1]));	//dPhi/dy - stored as component 1
		deriv(i,j,k,imh) = 0.5*((phiOld(i,j+1,k,0)-phiOld(i,j-1,k,0))/(2.0*dx[1])+(phiOld(i-1,j+1,k,0)-phiOld(i-1,j-1,k,0))/(2.0*dx[1]));	//dPhi/dy - stored as component 1
		#if (AMREX_SPACEDIM > 2)
		deriv(i,j,k,kph) = 0.5*((phiOld(i,j+1,k+1,0)-phiOld(i,j-1,k+1,0))/(2.0*dx[1])+(phiOld(i,j+1,k,0)-phiOld(i,j-1,k,0))/(2.0*dx[1]));
		deriv(i,j,k,kmh) = 0.5*((phiOld(i,j+1,k,0)-phiOld(i,j-1,k,0))/(2.0*dx[1])+(phiOld(i,j+1,k-1,0)-phiOld(i,j-1,k-1,0))/(2.0*dx[1]));
		#endif
}

#if (AMREX_SPACEDIM > 2)
void derivative_z(int i, int j, int k,
		amrex::Array4<amrex::Real const> const& phiOld,
		amrex::Array4<amrex::Real      > const& deriv,
		Geometry const& geom)
{
		GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();
		//Two sided difference
		deriv(i,j,k,jph) = 0.5*((phiOld(i,j+1,k+1,0)-phiOld(i,j+1,k-1,0))/(2.0*dx[2])+(phiOld(i,j,k+1,0)-phiOld(i,j,k-1,0))/(2.0*dx[2]));
		deriv(i,j,k,jmh) = 0.5*((phiOld(i,j,k+1,0)-phiOld(i,j,k-1,0))/(2.0*dx[2])+(phiOld(i,j-1,k+1,0)-phiOld(i,j-1,k-1,0))/(2.0*dx[2]));
		deriv(i,j,k,iph) = 0.5*((phiOld(i+1,j,k+1,0)-phiOld(i+1,j,k-1,0))/(2.0*dx[2])+(phiOld(i,j,k+1,0)-phiOld(i,j,k-1,0))/(2.0*dx[2]));	//dPhi/dy - stored as component 1
		deriv(i,j,k,imh) = 0.5*((phiOld(i,j,k+1,0)-phiOld(i,j,k-1,0))/(2.0*dx[2])+(phiOld(i-1,j,k+1,0)-phiOld(i-1,j,k-1,0))/(2.0*dx[2]));	//dPhi/dy - stored as component 1
		#if (AMREX_SPACEDIM > 2)
		deriv(i,j,k,kph) = (phiOld(i,j,k+1,0)-phiOld(i,j,k,0))/(dx[2]);
		deriv(i,j,k,kmh) = (phiOld(i,j,k,0)-phiOld(i,j,k-1,0))/(dx[2]);
		#endif
}
#endif
*/


#endif	
	

	
