#ifndef _ANISOTROPY_01_H
#define _ANISOTROPY_01_H

#include <AMReX_Utility.H>
using namespace amrex;


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void dAdq(int i, int j , int k, Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> daby, Array2D<Real, 0, AMREX_SPACEDIM*2, 0, AMREX_SPACEDIM-1,Order::C> &r_qab,  Array2D<Real, 0, AMREX_SPACEDIM*2, 0, AMREX_SPACEDIM-1,Order::C> &dadq,  Array1D<Real, 0, AMREX_SPACEDIM*2> &q2, int a, int b){
    
    for(int p=0; p<r_qab.xlen(); p++){

		q2(p) = r_qab(p,X)*r_qab(p,X)+r_qab(p,Y)*r_qab(p,Y)
    #if(AMREX_SPACEDIM>2)
                +r_qab(p,Z)*r_qab(p,Z);
    #else
        ;
    #endif
		Real qx3 = r_qab(p,X)*r_qab(p,X)*r_qab(p,X);
        Real qy3 = r_qab(p,Y)*r_qab(p,Y)*r_qab(p,Y);
    #if(AMREX_SPACEDIM>2)
        Real qz3 = r_qab(p,Z)*r_qab(p,Z)*r_qab(p,Z);
    #endif

        Real q22 = q2(p)*q2(p);
        Real q23 = q2(p)*q2(p)*q2(p);
        Real q4 = qx3*r_qab(p,X)+qy3*r_qab(p,Y)
    #if(AMREX_SPACEDIM>2)
                +qz3*r_qab(p,Z);
    #else
        ;
    #endif

        if(fabs(q2(p))>1e-15){

            dadq(p,X) = 16.0*daby(a,b)*(qx3/q22-r_qab(p,X)*q4/q23);
            dadq(p,Y) = 16.0*daby(a,b)*(qy3/q22-r_qab(p,Y)*q4/q23);
        #if(AMREX_SPACEDIM>2)
            dadq(p,Z) = 16.0*daby(a,b)*(qz3/q22-r_qab(p,Z)*q4/q23);
        #endif

        }
        else{
            dadq(p,X) = 0.0;
            dadq(p,Y) = 0.0;
        #if(AMREX_SPACEDIM>2)
            dadq(p,Z) = 0.0;
        #endif
        }

    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void function_ac(int i, int j , int k, Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> daby, Array2D<Real, 0, AMREX_SPACEDIM*2, 0, AMREX_SPACEDIM-1,Order::C> &r_qab,  Array1D<Real, 0, AMREX_SPACEDIM*2> &ac,  Array1D<Real, 0, AMREX_SPACEDIM*2> &q2, int a, int b){
    
    for(int p=0; p<r_qab.xlen(); p++){
        
        Real qx2 = r_qab(p,X)*r_qab(p,X);
        Real qy2 = r_qab(p,Y)*r_qab(p,Y);
        Real qz2{0.0};
        Real qz4{0.0};
        if(AMREX_SPACEDIM>2){
            qz2 = r_qab(p,Z)*r_qab(p,Z);
            qz4 = qz2*qz2;
        }

        Real qx4 = qx2*qx2;
        Real qy4 = qy2*qy2;

        if(fabs(q2(p))>1e-15){
            ac(p) = 1.0-daby(a,b)*(3.0-4.0*(qx4+qy4+qz4)/(q2(p)*q2(p)));
        }
        else{
            ac(p) = 1.0;
        }

    }
}

#endif