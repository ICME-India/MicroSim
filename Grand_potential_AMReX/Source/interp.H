#ifndef INTERP_H_
#define INTERP_H_

#include <AMReX_Utility.H>
#include "Function_F4.H"

using namespace amrex;

void Func_H(MultiFab& h_phi,MultiFab& dh_dphi, MultiFab& phi_old){

	BL_PROFILE("Func_H()");

    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();
		Array4<Real> const& phiOld = phi_old.array(mfi);
        Array4<Real> const& arrh_phi = h_phi.array(mfi);
        Array4<Real> const& arrdh_dphi = dh_dphi.array(mfi);
        int numphase=nump;
	
		amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {   
            Array1D <Real,0,3> sum{0.0};
            Real count=0;

            for(int a=0;a<numphase;a++){ 
                for(int b=0;b<numphase;b++){
                    if(numphase>2){
                        for(int c=0;c<numphase;c++){
                            if(a!=b && a!=c && b<c)
                            {
                                sum(count) = phiOld(i,j,k,b)-phiOld(i,j,k,c);
                                count++;
                            }
                        }
                    }
                    else if(a!=b && numphase==2){
                        
                        sum(count) = phiOld(i,j,k,b);
                        count++;
                    }
                }

                arrh_phi(i,j,k,a) = 0.25*pow(phiOld(i,j,k,a),2)*((15.0*(1-phiOld(i,j,k,a))*(1+phiOld(i,j,k,a)-pow(sum(a),2))) + phiOld(i,j,k,a)*(9.0*pow(phiOld(i,j,k,a),2)-5.0));

                arrdh_dphi(i,j,k,a) = (15.0/4.0)*(3.0*pow(phiOld(i,j,k,a),4)-4.0*pow(phiOld(i,j,k,a),3)-pow(phiOld(i,j,k,a),2)+2.0*phiOld(i,j,k,a)-pow(sum(a),2)*(2.0*phiOld(i,j,k,a)-3.0*pow(phiOld(i,j,k,a),2)));

            }
		});
	}
}





// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// double dhphi(int i, int j, int k, amrex::Array4<amrex::Real const> const& phi, long a, long b,int numphase) {
//   long e, f;
//   double val=0.0;
//   double ph_sum=0.0;
  
//   if(a==b){
//     if(numphase>2){ 
//     for(int m=0; m<numphase; m++){
//         if(a!=m && b!=m){
//             ph_sum = phi(i,j,k,m);
//             if(m==nump-1){
//                 val =  (15.0/4.0)*(3.0*pow(phi(i,j,k,a),4)-4.0*pow(phi(i,j,k,a),3)-pow(phi(i,j,k,a),2)+2.0*phi(i,j,k,a)-pow(ph_sum,2)*(2.0*phi(i,j,k,a)-3.0*pow(phi(i,j,k,a),2)));
//             }
//         }
//         ph_sum=0.0;
//     }
//     }
//     else{
//         val =  (15.0/4.0)*(3.0*pow(phi(i,j,k,a),4)-4.0*pow(phi(i,j,k,a),3)-pow(phi(i,j,k,a),2)+2.0*phi(i,j,k,a)-pow(phi(i,j,k,b),2)*(2.0*phi(i,j,k,a)-3.0*pow(phi(i,j,k,a),2)));
//     }
//   }
  
//   else{
//     if(numphase>2){
//     for(int m=0; m<numphase; m++){
//         if(a!=b && a!=m && b!=m){
//             val =  -15*0.5*(pow(phi(i,j,k,a),2)-pow(phi(i,j,k,a),3))*(phi(i,j,k,b)-phi(i,j,k,m));
//         }
//     }
//     }
//     else{
//             val = -15*0.5*(pow(phi(i,j,k,a),2)-pow(phi(i,j,k,a),3))*(phi(i,j,k,b));
//     }
//   }

//   return val;
// } 


// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// double hphi(int i, int j, int k, amrex::Array4<amrex::Real const> const& phi, long a) {
//   long b,c;
//   double sum1 = 0.0;
//   double val = 0.0;

//     for(int b=0;b<nump;b++){
//         if(nump>2){
//             for(int c=0;c<nump;c++){
//                 if(a!=b && a!=c && b<c)
//                 {
//                     sum1 = phi(i,j,k,b)-phi(i,j,k,c);
//                     val = 0.25*pow(phi(i,j,k,a),2)*((15.0*(1-phi(i,j,k,a))*(1+phi(i,j,k,a)-pow(sum1,2))) + phi(i,j,k,a)*(9.0*pow(phi(i,j,k,a),2)-5.0));
                    
//                 }
//                 sum1=0.0;
//             }
//         }
//         else{
//             if(a!=b){           
//                 val = 0.25*pow(phi(i,j,k,a),2)*((15.0*(1-phi(i,j,k,a))*(1+phi(i,j,k,a)-pow(phi(i,j,k,b),2))) + phi(i,j,k,a)*(9.0*pow(phi(i,j,k,a),2)-5.0));
//             }
//         }
//     }

//   return val;
// }





#endif
