#ifndef CALC_MOI_H
#define CALC_MOI_H

using namespace amrex;

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Cal_derM(int i, int j, int k, Array3D <Real, 0, AMREX_SPACEDIM*2, 0, compcount-2, 0, compcount-2, Order::C> &der_M, amrex::Array4<Real const> const& phi, Array3D <Real,0, phasecount-1, 0, compcount-1, 0, compcount-1, Order::C> der_cmu, Array3D <Real,0, phasecount-1, 0, compcount-1, 0, compcount-1, Order::C> diffs, int numcomp, int numphase){
    
    for(int l=0; l<numcomp-1; l++){
        for(int m=0; m< numcomp-1; m++){
            for (int b=0; b < numphase; b++) {
                for (int n=0; n < numcomp-1; n++) {
                    der_M(0,l,m) += diffs(b,l,n)*der_cmu(b,n,m)*phi(i,j,k,b);
                    der_M(1,l,m) += diffs(b,l,n)*der_cmu(b,n,m)*phi(i+1,j,k,b);
                    der_M(2,l,m) += diffs(b,l,n)*der_cmu(b,n,m)*phi(i-1,j,k,b);
                    der_M(3,l,m) += diffs(b,l,n)*der_cmu(b,n,m)*phi(i,j+1,k,b);
                    der_M(4,l,m) += diffs(b,l,n)*der_cmu(b,n,m)*phi(i,j-1,k,b);
                #if(AMREX_SPACEDIM>2)
                    der_M(5,l,m) += diffs(b,l,n)*der_cmu(b,n,m)*phi(i,j,k+1,b);
                    der_M(6,l,m) += diffs(b,l,n)*der_cmu(b,n,m)*phi(i,j,k-1,b);
                #endif
                }
            }
        }
    }

}


#endif