#ifndef GP_UTILITY_H
#define GP_UTILITY_H

#include <AMReX_Utility.H>
using namespace amrex;


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Vec_rot(int i, int j, int k, int a, int b, Array2D<Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> &vec, Array2D<Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> &r_vec, Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_x,Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_y,Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_z){

    
#if(AMREX_SPACEDIM == 2)
    
        for(int p = 0; p < vec.xlen(); p++){
            r_vec(p,X) = vec(p,X)*cos(matrot_z(a,b))-vec(p,Y)*sin(matrot_z(a,b));
            r_vec(p,Y) = vec(p,X)*sin(matrot_z(a,b))+vec(p,Y)*cos(matrot_z(a,b));
        }

#endif

#if(AMREX_SPACEDIM>2)

    for(int p = 0; p < vec.xlen(); p++){
            r_vec(p,X) = vec(p,X)*(cos(matrot_y(a,b))*cos(matrot_z(a,b)))-vec(p,Y)*(cos(matrot_y(a,b))*sin(matrot_z(a,b)))-vec(p,Z)*sin(matrot_y(a,b));
            r_vec(p,Y) = vec(p,X)*(-cos(matrot_z(a,b))*sin(matrot_x(a,b))*sin(matrot_y(a,b))+cos(matrot_x(a,b))*sin(matrot_z(a,b)))+vec(p,Y)*(sin(matrot_x(a,b))*sin(matrot_y(a,b))*sin(matrot_z(a,b))+cos(matrot_x(a,b))*cos(matrot_z(a,b)))-vec(p,Z)*(sin(matrot_x(a,b))*cos(matrot_y(a,b)));
            r_vec(p,Z) = vec(p,X)*(cos(matrot_x(a,b))*sin(matrot_y(a,b))*cos(matrot_z(a,b))+sin(matrot_z(a,b))*sin(matrot_x(a,b)))+vec(p,Y)*(-cos(matrot_x(a,b))*sin(matrot_y(a,b))*sin(matrot_z(a,b))+cos(matrot_z(a,b))*sin(matrot_x(a,b)))+vec(p,Z)*(cos(matrot_x(a,b))*cos(matrot_y(a,b)));
        }

#endif

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Vec_rot(int i, int j, int k, int a, int b, Array1D <Real, 0, AMREX_SPACEDIM-1> &vec, Array1D <Real, 0, AMREX_SPACEDIM-1> &r_vec,Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_x,Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_y,Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_z){

#if(AMREX_SPACEDIM == 2)
    
            r_vec(X) = vec(X)*cos(matrot_z(a,b))-vec(Y)*sin(matrot_z(a,b));
            r_vec(Y) = vec(X)*sin(matrot_z(a,b))+vec(Y)*cos(matrot_z(a,b));


#endif

#if(AMREX_SPACEDIM>2)

            r_vec(X) = vec(X)*(cos(matrot_y(a,b))*cos(matrot_z(a,b)))-vec(Y)*(cos(matrot_y(a,b))*sin(matrot_z(a,b)))-vec(Z)*sin(matrot_y(a,b));
            r_vec(Y) = vec(X)*(-cos(matrot_z(a,b))*sin(matrot_x(a,b))*sin(matrot_y(a,b))+cos(matrot_x(a,b))*sin(matrot_z(a,b)))+vec(Y)*(sin(matrot_x(a,b))*sin(matrot_y(a,b))*sin(matrot_z(a,b))+cos(matrot_x(a,b))*cos(matrot_z(a,b)))-vec(Z)*(sin(matrot_x(a,b))*cos(matrot_y(a,b)));
            r_vec(Z) = vec(X)*(cos(matrot_x(a,b))*sin(matrot_y(a,b))*cos(matrot_z(a,b))+sin(matrot_z(a,b))*sin(matrot_x(a,b)))+vec(Y)*(-cos(matrot_x(a,b))*sin(matrot_y(a,b))*sin(matrot_z(a,b))+cos(matrot_z(a,b))*sin(matrot_x(a,b)))+vec(Z)*(cos(matrot_x(a,b))*cos(matrot_y(a,b)));
        

#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Inv_vec_rot(int i, int j, int k, int a, int b, Array2D <Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> &vec, Array2D <Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> &r_vec, Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_x, Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_y, Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_z){
    
#if(AMREX_SPACEDIM == 2)
    
        for(int p = 0; p < vec.xlen(); p++){
            r_vec(p,X) = vec(p,X)*cos(matrot_z(a,b))+vec(p,Y)*sin(matrot_z(a,b));
            r_vec(p,Y) = -vec(p,X)*sin(matrot_z(a,b))+vec(p,Y)*cos(matrot_z(a,b));
        }

#endif

#if(AMREX_SPACEDIM>2)

    for(int p = 0; p < vec.xlen(); p++){
            r_vec(p,X) = vec(p,X)*(cos(matrot_y(a,b))*cos(matrot_z(a,b)))+vec(p,Y)*(cos(matrot_y(a,b))*sin(matrot_z(a,b)))+vec(p,Z)*sin(matrot_y(a,b));
            r_vec(p,Y) = vec(p,X)*(-cos(matrot_z(a,b))*sin(matrot_x(a,b))*sin(matrot_y(a,b))-cos(matrot_x(a,b))*sin(matrot_z(a,b)))+vec(p,Y)*(-sin(matrot_x(a,b))*sin(matrot_y(a,b))*sin(matrot_z(a,b))+cos(matrot_x(a,b))*cos(matrot_z(a,b)))+vec(p,Z)*(sin(matrot_x(a,b))*cos(matrot_y(a,b)));
            r_vec(p,Z) = vec(p,X)*(-cos(matrot_x(a,b))*sin(matrot_y(a,b))*cos(matrot_z(a,b))+sin(matrot_z(a,b))*sin(matrot_x(a,b)))+vec(p,Y)*(-cos(matrot_x(a,b))*sin(matrot_y(a,b))*sin(matrot_z(a,b))-cos(matrot_z(a,b))*sin(matrot_x(a,b)))+vec(p,Z)*(cos(matrot_x(a,b))*cos(matrot_y(a,b)));
        }

#endif
}



#endif