#ifndef _CALC_H_
#define _CALC_H_

#include<AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include "Tau.H"
//#include "simplex.H"

#define iph 0
#define imh 1
#define jph 2
#define jmh 3
#define kph 4
#define kmh 5

using namespace amrex;

///////////////////////////////////////////////////////////////////////////////////////////////////////

//AMREX_GPU_DEVICE AMREX_FORCE_INLINE

void update_phi(MultiFab& phi_new, MultiFab& phi_old, MultiFab& term1, MultiFab& term2, MultiFab& term3,MultiFab& lambad, MultiFab& derivx, MultiFab& derivy, Geometry const& geom)
{	
	BL_PROFILE("update_phi()");	

	for (MFIter mfi(phi_old); mfi.isValid(); ++mfi)
	{
		const Box& dbx = mfi.validbox();
		Array4<Real> const& fin_term1 = term1.array(mfi);
		Array4<Real> const& fin_term2 = term2.array(mfi);
		Array4<Real> const& fin_term3 = term3.array(mfi);
		Array4<Real> const& phiNew = phi_new.array(mfi);
		Array4<Real> const& phiOld = phi_old.array(mfi);
		Array4<Real> const& lamb = lambad.array(mfi);
		Array4<Real> const& derivxx = derivx.array(mfi);
        Array4<Real> const& derivyy = derivy.array(mfi);
		Real tauu = tau_final;
		Real time_step = dt;
		Real epsilon = eps;
		Real molar_vol = Vm;
		Real numphase = nump; 
		GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();
		
		amrex::ParallelFor( dbx, 
		[=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
		{	
			Real sum_lambda{0.0};
			Real div{0.0};
			Real active_phase{0.0};

			for(int a=0; a<numphase; a++){

				div =  (phiOld(i+1,j,k,a)-2.0*phiOld(i,j,k,a)+phiOld(i-1,j,k,a))/(dx[0]*dx[0])+(phiOld(i,j+1,k,a)-2.0*phiOld(i,j,k,a)+phiOld(i,j-1,k,a))/(dx[1]*dx[1]); 

			if(fabs(div)>0.0){

				lamb(i,j,k,a) = epsilon*fin_term1(i,j,k,a)-fin_term2(i,j,k,a)/epsilon-fin_term3(i,j,k,a)/molar_vol;

				sum_lambda += lamb(i,j,k,a); 
				active_phase++;
			}
			}

			 if (active_phase) {
      			sum_lambda /= active_phase;
    			}
				
	
			for(int a=0; a<numphase; a++){
			 div =  (phiOld(i+1,j,k,a)-2.0*phiOld(i,j,k,a)+phiOld(i-1,j,k,a))/(dx[0]*dx[0])+(phiOld(i,j+1,k,a)-2.0*phiOld(i,j,k,a)+phiOld(i,j-1,k,a))/(dx[1]*dx[1]); 
			 //div =  (derivxx(i,j,k,4*a+0) - derivxx(i,j,k,4*a+1))/dx[0] + (derivyy(i,j,k,4*a+2) - derivyy(i,j,k,4*a+3))/dx[1]; 

			if(fabs(div)>0.0){
				phiNew(i,j,k,a) = phiOld(i,j,k,a) + (time_step/(epsilon*tauu))*(lamb(i,j,k,a)-sum_lambda);
				}
			else{
				phiNew(i,j,k,a) = phiOld(i,j,k,a);
				}
			}
		});
	}
}
	



void dmudt(MultiFab& mu_new, MultiFab& mu_old, MultiFab& phi_new, MultiFab& phi_old, MultiFab& derivx, MultiFab& derivy, MultiFab& h_phi, MultiFab& dh_dphi, Geometry const& geom)
{
	BL_PROFILE("dmudt()");

	GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

	for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();
		Array4<Real> const& phiOld = phi_old.array(mfi);
		Array4<Real> const& phiNew = phi_new.array(mfi);
		Array4<Real> const& mun = mu_new.array(mfi);
		Array4<Real> const& muo = mu_old.array(mfi);
		Array4<Real> const& derivxx = derivx.array(mfi);
        Array4<Real> const& derivyy = derivy.array(mfi);
		Array4<Real> const& arr_h_phi = h_phi.array(mfi);
		Array4<Real> const& arr_dh_dphi = dh_dphi.array(mfi);
		
		Array1D <Real,0,20> BB;
		for(int a=0; a<nump; a++){
		BB(a) = B[a];
		}
		
		Array1D <Real,0,20> der_cmu;
		for(int a=0; a<nump; a++){
		der_cmu(a) = dcdmu[a];
		}

		Array1D <Real,0,20> diffs;
		for(int a=0; a<nump; a++){
		diffs(a) = diff[a];
		}

		Real time_step = dt;
		Real epsilon = eps;
		Real numphase = nump;
		GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();
	
		amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {	
		
			Real sum=0;
			Array2D <Real,0,2,0,2> der_dhdphi{0.0};
			Array2D <Real,0,1,0,1> der_dhdphi_2p{0.0};
			Array2D <Real,0,1,0,3> jat_x{0.0};
			Array2D <Real,0,1,0,3> jat_y{0.0};
			Array2D <Real,0,2,0,3> norm_x{0.0};
			Array2D <Real,0,2,0,3> norm_y{0.0};
			Array1D <Real,0,4> mod{0.0};
			Array1D <Real,0,4> conct{0.0};
			Array1D <Real,0,4> phix{0.0};
			Array1D <Real,0,4> phiy{0.0};

			Real dmudx_iph = (muo(i+1,j,k)-muo(i,j,k))/(dx[0]);
			Real dmudx_imh = (muo(i,j,k)-muo(i-1,j,k))/(dx[0]);
			Real dmudy_jph = (muo(i,j+1,k)-muo(i,j,k))/(dx[1]);
			Real dmudy_jmh = (muo(i,j,k)-muo(i,j-1,k))/(dx[1]);

			Real M_iph{0.0}, M_imh{0.0}, M_jph{0.0}, M_jmh{0.0}, subterm1{0.0}, subterm2{0.0} ,subterm3{0.0}, denom{0.0};
			Real adn{0.0},adn2{0.0},adn3{0.0},adn4{0.0};

			for(int a=0; a<numphase;a++){
				conct(a) = (muo(i,j,k)-BB(a))*der_cmu(a);
			}

			for(int a=0; a<numphase; a++){
			
				M_iph = M_iph + diffs(a)*0.5*(phiOld(i+1,j,k,a)+phiOld(i,j,k,a))*der_cmu(a);
				M_imh = M_imh + diffs(a)*0.5*(phiOld(i,j,k,a)+phiOld(i-1,j,k,a))*der_cmu(a); 
				M_jph = M_jph + diffs(a)*0.5*(phiOld(i,j+1,k,a)+phiOld(i,j,k,a))*der_cmu(a);
				M_jmh = M_jmh + diffs(a)*0.5*(phiOld(i,j,k,a)+phiOld(i,j-1,k,a))*der_cmu(a);
			
			}

			subterm1 = (M_iph*dmudx_iph - M_imh*dmudx_imh)/dx[0] + (M_jph*dmudy_jph - M_jmh*dmudy_jmh)/dx[1];

			for(int a=0;a<numphase;a++){
				for(int b=0;b<numphase;b++){

						if(numphase>2){
							for(int c=0;c<numphase;c++){
								if(a!=c && b!=c && a!=b){
									der_dhdphi(a,b) = -15*0.5*(pow(phiOld(i,j,k,a),2)-pow(phiOld(i,j,k,a),3))*(phiOld(i,j,k,b)-phiOld(i,j,k,c));
						}
						}
						}
						else if(a!=b && numphase==2){
							der_dhdphi_2p(a,b) = -15*0.5*(pow(phiOld(i,j,k,a),2)-pow(phiOld(i,j,k,a),3))*(phiOld(i,j,k,b));
						}
				}
			}



			for(int a=0; a<numphase; a++){
				for(int b=0; b<numphase; b++){
					
					if(a==b){
						sum = sum + arr_dh_dphi(i,j,k,a)*(phiNew(i,j,k,a)-phiOld(i,j,k,a))/(time_step);
					}  
					
					else{
						if(numphase>2){
							sum = sum +der_dhdphi(a,b)*(phiNew(i,j,k,b)-phiOld(i,j,k,b))/(time_step);
						}

						else if(numphase==2){
							sum = sum +der_dhdphi_2p(a,b)*(phiNew(i,j,k,b)-phiOld(i,j,k,b))/(time_step);
						}
					}
				}

				subterm2 = subterm2 + (muo(i,j,k)-BB(a))*der_cmu(a)*sum;
				sum=0.0;
			}

			for(int a=0; a<numphase;a++){

				phix(0)=(phiOld(i+1,j,k,a)-phiOld(i,j,k,a))/dx[0];
            	phix(1)=(phiOld(i,j,k,a)-phiOld(i-1,j,k,a))/dx[0];
            	phix(2)=(phiOld(i+1,j+1,k,a)-phiOld(i-1,j+1,k,a)+phiOld(i+1,j,k,a)-phiOld(i-1,j,k,a))/(4*dx[0]);
            	phix(3)=(phiOld(i+1,j,k,a)-phiOld(i-1,j,k,a)+phiOld(i+1,j-1,k,a)-phiOld(i-1,j-1,k,a))/(4*dx[0]);

            	phiy(0)=(phiOld(i+1,j+1,k,a)-phiOld(i+1,j-1,k,a)+phiOld(i,j+1,k,a)-phiOld(i,j-1,k,a))/(4*dx[1]);
            	phiy(1)=(phiOld(i,j+1,k,a)-phiOld(i,j-1,k,a)+phiOld(i-1,j+1,k,a)-phiOld(i-1,j-1,k,a))/(4*dx[1]);
            	phiy(2)=(phiOld(i,j+1,k,a)-phiOld(i,j,k,a))/dx[1];
            	phiy(3)=(phiOld(i,j,k,a)-phiOld(i,j-1,k,a))/dx[1];
			
				mod(0) = sqrt(pow(phix(0),2)+pow(phiy(0),2));
				mod(1) = sqrt(pow(phix(1),2)+pow(phiy(1),2));
				mod(2) = sqrt(pow(phix(2),2)+pow(phiy(2),2));
				mod(3) = sqrt(pow(phix(3),2)+pow(phiy(3),2));

				if(mod(0)>1e-15){
				 	norm_x(a,0) = phix(0)/(mod(0));
					norm_y(a,0) = phiy(0)/(mod(0));
				}
				 else{
					norm_x(a,0) = 0.0;
					norm_y(a,0) = 0.0;
				}

				if(mod(1)>1e-15){
					norm_x(a,1) = phix(1)/(mod(1));
					norm_y(a,1) = phiy(1)/(mod(1));
				}
				else{
					norm_x(a,1) = 0.0;
					norm_y(a,1) = 0.0;
				}
					
				if(mod(2)>1e-15){
					norm_x(a,2) = phix(2)/(mod(2));
					norm_y(a,2) = phiy(2)/(mod(2));
				}
				else{
					norm_x(a,2) = 0.0;
					norm_y(a,2) = 0.0;
				}

				if(mod(3)>1e-15){
					norm_x(a,3) = phix(3)/(mod(3));
					norm_y(a,3) = phiy(3)/(mod(3));
				}
				else{
					norm_x(a,3) = 0.0;
					norm_y(a,3) = 0.0;
				}
			}

			for(int a=0; a<numphase-1;a++){

				jat_x(a,0) = -(0.5/sqrt(2))*epsilon*(conct(numphase-1)-conct(a))*(((phiNew(i+1,j,k)+phiNew(i,j,k))*0.5 - (phiOld(i+1,j,k)+phiOld(i,j,k))*0.5)/time_step)*(norm_x(a,0))*fabs(norm_x(a,0)*norm_x(numphase-1,0)+norm_y(a,0)*norm_y(numphase-1,0));
				jat_x(a,1) = -(0.5/sqrt(2))*epsilon*(conct(numphase-1)-conct(a))*(((phiNew(i-1,j,k)+phiNew(i,j,k))*0.5 - (phiOld(i-1,j,k)+phiOld(i,j,k))*0.5)/time_step)*(norm_x(a,1))*fabs(norm_x(a,1)*norm_x(numphase-1,1)+norm_y(a,1)*norm_y(numphase-1,1));
					
				jat_y(a,2) = -(0.5/sqrt(2))*epsilon*(conct(numphase-1)-conct(a))*(((phiNew(i,j+1,k)+phiNew(i,j,k))*0.5 - (phiOld(i,j+1,k)+phiOld(i,j,k))*0.5)/time_step)*(norm_y(a,2))*fabs(norm_x(a,2)*norm_x(numphase-1,2)+norm_y(a,2)*norm_y(numphase-1,2));
				jat_y(a,3) = -(0.5/sqrt(2))*epsilon*(conct(numphase-1)-conct(a))*(((phiNew(i,j-1,k)+phiNew(i,j,k))*0.5 - (phiOld(i,j-1,k)+phiOld(i,j,k))*0.5)/time_step)*(norm_y(a,3))*fabs(norm_x(a,3)*norm_x(numphase-1,3)+norm_y(a,3)*norm_y(numphase-1,3));
			}

			for(int a=0;a<numphase-1;a++){
				adn = adn + jat_x(a,0);
				adn2 = adn2 + jat_x(a,1);
				adn3 = adn3 + jat_y(a,2);
				adn4 = adn4 + jat_y(a,3);
			}

			subterm3 = (adn-adn2)/dx[0] + (adn3-adn4)/dx[1];


			for(int a=0; a<numphase; a++){

				denom = denom + arr_h_phi(i,j,k,a)*der_cmu(a);

			}

			//Print()<<"subterm3: "<<subterm3<<"\n";

			mun(i,j,k) = muo(i,j,k) + time_step*(subterm1 - subterm2 - subterm3)/denom;


		});
	
	}
}



///////////////////////////////////////////////////////////////////////////////////////////////////////

// void dmudt(MultiFab& mu_new, MultiFab& mu_old, MultiFab& phi_new, MultiFab& phi_old, MultiFab& derivx, MultiFab& derivy, MultiFab& h_phi, MultiFab& dh_dphi, Geometry const& geom)
// {
// 	BL_PROFILE("dmudt()");

// 	GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

// 	for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
//     {
//         const Box& vbx = mfi.validbox();
// 		Array4<Real> const& phiold = phi_old.array(mfi);
// 		Array4<Real> const& phinew = phi_new.array(mfi);
// 		Array4<Real> const& mun = mu_new.array(mfi);
// 		Array4<Real> const& muo = mu_old.array(mfi);
// 		Array4<Real> const& derivxx = derivx.array(mfi);
//         Array4<Real> const& derivyy = derivy.array(mfi);
// 		Array4<Real> const& arr_h_phi = h_phi.array(mfi);
// 		Array4<Real> const& arr_dh_dphi = dh_dphi.array(mfi);
		
// 		Array1D <Real,0,20> BB;
// 		for(int a=0; a<nump; a++){
// 		BB(a) = B[a];
// 		}
		
// 		Array1D <Real,0,20> der_cmu;
// 		for(int a=0; a<nump; a++){
// 		der_cmu(a) = dcdmu[a];
// 		}

// 		Array1D <Real,0,20> diffs;
// 		for(int a=0; a<nump; a++){
// 		diffs(a) = diff[a];
// 		}

		
// 		Array2D <Real,0,1,0,2,Order::C> flag;
// 		flag(0,0)=diffu[0][0];
// 		flag(0,1)=diffu[0][1];
// 		flag(0,2)=diffu[0][2];
// 		flag(1,0)=diffu[1][0];
// 		flag(1,1)=diffu[1][1];
// 		flag(1,2)=diffu[1][2];
// 		Real time_step = dt;
// 		Real epsilon = eps;
		
	
// 		amrex::ParallelFor(vbx,
//         [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//         {	
// 			//if(phiold(i,j,k)==1 || phiold(i,j,k)==0){
// 			//		mun(i,j,k) = muo(i,j,k);
// 			//}

// 			//else{
// 			Real calpha = (muo(i,j,k) - (BB(0)))*(der_cmu(0));

// 			Real cbeta = muo(i,j,k)*(der_cmu(1));

// 			Real dmudx_iph = (muo(i+1,j,k)-muo(i,j,k))/(dx[0]);
// 			Real dmudx_imh = (muo(i,j,k)-muo(i-1,j,k))/(dx[0]);
// 			Real dmudy_jph = (muo(i,j+1,k)-muo(i,j,k))/(dx[1]);
// 			Real dmudy_jmh = (muo(i,j,k)-muo(i,j-1,k))/(dx[1]);

// 			Real modphisq_iph = pow(derivxx(i,j,k,0),2)+pow(derivyy(i,j,k,0),2);
// 			Real modphisq_imh = pow(derivxx(i,j,k,1),2)+pow(derivyy(i,j,k,1),2);
// 			Real modphisq_jph = pow(derivxx(i,j,k,2),2)+pow(derivyy(i,j,k,2),2);
// 			Real modphisq_jmh = pow(derivxx(i,j,k,3),2)+pow(derivyy(i,j,k,3),2);

// 			Real jatx_iph, jatx_imh, jaty_jph, jaty_jmh;

// 			if(modphisq_iph>1e-15){
// 			jatx_iph = -(0.5/sqrt(2))*epsilon*(cbeta-calpha)*(((phinew(i+1,j,k)+phinew(i,j,k))*0.5 - (phiold(i+1,j,k)+phiold(i,j,k))*0.5)/time_step)*(derivxx(i,j,k,0)/sqrt(modphisq_iph));
// 			}
// 			else{
// 				jatx_iph = 0.0;
// 			}

// 			if(modphisq_imh>1e-15){
// 			jatx_imh = -(0.5/sqrt(2))*epsilon*(cbeta-calpha)*(((phinew(i-1,j,k)+phinew(i,j,k))*0.5 - (phiold(i-1,j,k)+phiold(i,j,k))*0.5)/time_step)*(derivxx(i,j,k,1)/sqrt(modphisq_imh));
// 			}
// 			else{
// 				jatx_imh = 0.0;
// 			}

// 			if(modphisq_jph>1e-15){
// 			jaty_jph = -(0.5/sqrt(2))*epsilon*(cbeta-calpha)*(((phinew(i,j+1,k)+phinew(i,j,k))*0.5 - (phiold(i,j+1,k)+phiold(i,j,k))*0.5)/time_step)*(derivyy(i,j,k,2)/sqrt(modphisq_jph));
// 			}
// 			else{
// 				jaty_jph = 0.0;
// 			}

// 			if(modphisq_jmh>1e-15){
// 			jaty_jmh = -(0.5/sqrt(2))*epsilon*(cbeta-calpha)*(((phinew(i,j-1,k)+phinew(i,j,k))*0.5 - (phiold(i,j-1,k)+phiold(i,j,k))*0.5)/time_step)*(derivyy(i,j,k,3)/sqrt(modphisq_jmh));
// 			}
// 			else{
// 				jaty_jmh = 0.0;
// 			}

// 			Real dbdx = (((flag(0,2)*0.5*(phiold(i,j,k)+phiold(i+1,j,k))*der_cmu(0) + flag(1,2)*(1.0 - 0.5*(phiold(i,j,k)+phiold(i+1,j,k)))*der_cmu(1))*dmudx_iph - jatx_iph)
// 					-	 ((flag(0,2)*0.5*(phiold(i,j,k)+phiold(i-1,j,k))*der_cmu(0) + flag(1,2)*(1.0 - 0.5*(phiold(i,j,k)+phiold(i-1,j,k)))*der_cmu(1))*dmudx_imh - jatx_imh))/dx[0];

// 			Real dbdy = (((flag(0,2)*0.5*(phiold(i,j,k)+phiold(i,j+1,k))*der_cmu(0) + flag(1,2)*(1.0 - 0.5*(phiold(i,j,k)+phiold(i,j+1,k)))*der_cmu(1))*dmudy_jph - jaty_jph)
// 					-	 ((flag(0,2)*0.5*(phiold(i,j,k)+phiold(i,j-1,k))*der_cmu(0) + flag(1,2)*(1.0 - 0.5*(phiold(i,j,k)+phiold(i,j-1,k)))*der_cmu(1))*dmudy_jmh - jaty_jmh))/dx[1];


// //			Real cdhdt = (calpha - cbeta)*(6.0*phiold(i,j,k)*(1.0-phiold(i,j,k)))*(phinew(i,j,k) - phiold(i,j,k))/dt;
// 			Real cdhdt = (calpha - cbeta)*(arr_dh_dphi(i,j,k))*(phinew(i,j,k) - phiold(i,j,k))/time_step;

// //			Real coeffdmudt = pow(phiold(i,j,k),2)*(3.0 - 2.0*phiold(i,j,k))*dcdmu[0] + (1.0 - pow(phiold(i,j,k),2)*(3.0 - 2.0*phiold(i,j,k)))*dcdmu[1];
// 			Real coeffdmudt = arr_h_phi(i,j,k)*der_cmu(0) + (1.0 - arr_h_phi(i,j,k))*der_cmu(1);

// 			Real dmudt = (dbdx + dbdy - cdhdt)/coeffdmudt;

// 			mun(i,j,k) = muo(i,j,k) + dmudt*time_step;
// 			//}
// 		});
	
// 	}
// }































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
	

	
