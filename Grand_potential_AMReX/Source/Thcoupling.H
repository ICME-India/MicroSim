#ifndef THCOUPLING_H_
#define THCOUPLING_H_

using namespace amrex;

void computeterm3(amrex::MultiFab& mu_old, amrex::MultiFab& term3, amrex::MultiFab& phi_old, amrex::MultiFab& dh_dphi, amrex::MultiFab& psi)
{

	BL_PROFILE("computeterm3()");
	
	for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();
		Array4<Real> const& phiOld = phi_old.array(mfi);
		Array4<Real> const& term = term3.array(mfi);
		Array4<Real> const& mu = mu_old.array(mfi);
		Array4<Real> const& psii = psi.array(mfi);
		Array4<Real> const& arr_dh_dphi = dh_dphi.array(mfi);
		int numphase = nump;
		
		Array1D <Real,0,3> BB;
		for(int a=0; a<nump; a++){
		BB(a) = B[a];
		}

		Array1D <Real,0,3> DD;
		for(int a=0; a<nump; a++){
		DD(a) = D[a];
		}
		
		Array1D <Real,0,3> der_cmu;
		for(int a=0; a<nump; a++){
		der_cmu(a) = dcdmu[a];
		}
	
		amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {	
			int sum = 0;
			Array2D <Real,0,2,0,2> der_dhdphi{0.0};
			Array2D <Real,0,1,0,1> der_dhdphi_2p{0.0};
			//int count = 0;
			
			for(int a=0; a<numphase; a++){
			
				psii(i,j,k,a) = -pow((mu(i,j,k) - BB(a)),2)*der_cmu(a)*0.5 + DD(a);
			
			}

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
					sum = sum + psii(i,j,k,a)*arr_dh_dphi(i,j,k,a);
					}  
					
					else{
						if(numphase>2){
							sum = sum + der_dhdphi(b,a)*psii(i,j,k,b);
						}
						if(numphase==2){
							sum = sum + der_dhdphi_2p(b,a)*psii(i,j,k,b);
						}
					}
				}
				term(i,j,k,a) = sum;
				sum=0.0;



				//term(i,j,k,0) = (psii(i,j,k,0)-psii(i,j,k,1))*arr_dh_dphi(i,j,k,0);
				//term(i,j,k,1) = (-psii(i,j,k,0)+psii(i,j,k,1))*arr_dh_dphi(i,j,k,1);


				//term(i,j,k,a) = psii(i,j,k,a)*arr_dh_dphi(i,j,k,a) - psii(i,j,k,nump-1)*0.5*15.0*(pow(phiOld(i,j,k,nump-1),2)-pow(phiOld(i,j,k,nump-1),3)*(val(a))); 

				// term(i,j,k,nump-1) = (-psii(i,j,k,a) + psii(i,j,k,nump-1))*arr_dh_dphi(i,j,k,nump-1);

				//term(i,j,k,nump-1) = psii(i,j,k,nump-1)*arr_dh_dphi(i,j,k,nump-1) - 0.5*15.0*(pow(phiOld(i,j,k,a),2)-pow(phiOld(i,j,k,a),3)*(val(nump-1)))*psii(i,j,k,a);

			}
			
			//term(i,j,k) = (arr_dh_dphi(i,j,k))*(psii(i,j,k,0) - psii(i,j,k,1));
			
		});
	
	}
}

#endif
