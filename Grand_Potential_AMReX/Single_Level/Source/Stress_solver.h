#ifndef STRESS_SOLVER_H
#define STRESS_SOLVER_H

#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>
#include <Variables.h>
using namespace amrex;

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void calculate_eigen_strain(int i, int j,int k,int numphase,Array2D <Real, 0, phasecount-1,0,6,Order::C> eigenst, Array2D<Real,0,AMREX_SPACEDIM-1,0,5,Order::C> &eg_st, amrex::Array4<Real const> const& phi){
    for(int a=0; a<numphase; a++){
        
        eg_st(X,xx) += eigenst(a,1)*phi(i,j,k,a);  
        eg_st(X,yy) += eigenst(a,2)*phi(i,j,k,a);
        eg_st(X,xy) += eigenst(a,6)*phi(i,j,k,a);
    #if(AMREX_SPACEDIM>2)
        eg_st(X,zz) += eigenst(a,3)*phi(i,j,k,a);
        eg_st(X,yz) += eigenst(a,4)*phi(i,j,k,a);
        eg_st(X,xz) += eigenst(a,5)*phi(i,j,k,a);
    #endif

        eg_st(Y,xx) += eigenst(a,1)*phi(i,j,k,a);  
        eg_st(Y,yy) += eigenst(a,2)*phi(i,j,k,a);
        eg_st(Y,xy) += eigenst(a,6)*phi(i,j,k,a);
    #if(AMREX_SPACEDIM>2)
        eg_st(Y,zz) += eigenst(a,3)*phi(i,j,k,a);
        eg_st(Y,yz) += eigenst(a,4)*phi(i,j,k,a);
        eg_st(Y,xz) += eigenst(a,5)*phi(i,j,k,a);
    #endif

    #if(AMREX_SPACEDIM>2)
        eg_st(Z,xx) += eigenst(a,1)*phi(i,j,k,a);  
        eg_st(Z,yy) += eigenst(a,2)*phi(i,j,k,a);
        eg_st(Z,xy) += eigenst(a,6)*phi(i,j,k,a);
    
        eg_st(Z,zz) += eigenst(a,3)*phi(i,j,k,a);
        eg_st(Z,yz) += eigenst(a,4)*phi(i,j,k,a);
        eg_st(Z,xz) += eigenst(a,5)*phi(i,j,k,a);
    #endif

    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void calculate_stiffness(int i, int j, int k, int numphase,Array2D <Real, 0, phasecount-1,0,2,Order::C> stiff_phase, Array2D<Real,0,AMREX_SPACEDIM-1,0,2,Order::C> &stiffness_c, amrex::Array4<Real const> const& phi){
    for(int a=0; a<numphase; a++){
        stiffness_c(X,C11) += stiff_phase(a,1)*phi(i,j,k,a);
        stiffness_c(X,C12) += stiff_phase(a,2)*phi(i,j,k,a);
        stiffness_c(X,C44) += stiff_phase(a,3)*phi(i,j,k,a);

        stiffness_c(Y,C11) += stiff_phase(a,1)*phi(i,j,k,a);
        stiffness_c(Y,C12) += stiff_phase(a,2)*phi(i,j,k,a);
        stiffness_c(Y,C44) += stiff_phase(a,3)*phi(i,j,k,a); 
    
    #if(AMREX_SPACEDIM>2)
        stiffness_c(Z,C11) += stiff_phase(a,1)*phi(i,j,k,a);
        stiffness_c(Z,C12) += stiff_phase(a,2)*phi(i,j,k,a);
        stiffness_c(Z,C44) += stiff_phase(a,3)*phi(i,j,k,a);
    #endif
    }
}


void Stress_solver_2D(MultiFab& phi_old, MultiFab& disp_X, MultiFab& disp_Y, MultiFab& disp_Z, MultiFab& strain_X, MultiFab& strain_Y, MultiFab& strain_Z, Geometry const& geom){

    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();								//Defining the box for iteration space
		Array4<Real> const& phiOld = phi_old.array(mfi);				//Taking the Multifabs as array
		Array4<Real> const& dispx = disp_X.array(mfi);				//Taking the Multifabs as array
		Array4<Real> const& dispy = disp_Y.array(mfi);					//Taking the Multifabs as array
		Array4<Real> const& strainX = strain_X.array(mfi);
        Array4<Real> const& strainY = strain_Y.array(mfi);

        int numphase = nump;

        Array2D <Real, 0, phasecount-1,0,6,Order::C> eigenst{};
        for(int g=0; g<nump; g++){
            for(int h=0; h<7; h++){
                eigenst(g,h) = egstr[g][h];
            }
        }

        amrex::ParallelFor(vbx,
		[=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
		{   
            Array2D<Real,0,AMREX_SPACEDIM-1,0,5, Order::C> eg_st{};
                  //strain[X][xx] = 0,0 and so on

            calculate_eigen_strain(i,j,k,numphase,eigenst,eg_st,phiOld);

            strainX(i,j,k,xx) = 0.5*(dispx(i+1,j,k,2)-dispx(i-1,j,k,2)) -eg_st(X,xx);
            strainX(i,j,k,yy) = 0.5*(dispy(i,j+1,k,2)-dispy(i,j-1,k,2)) -eg_st(X,yy);
            strainX(i,j,k,xy) = 0.25*(dispx(i,j+1,k,2)-dispx(i,j-1,k,2)+dispy(i+1,j,k,2)-dispy(i-1,j,k,2));
            strainY(i,j,k,xx) = 0.5*(dispx(i+1,j,k,2)-dispx(i-1,j,k,2)) -eg_st(Y,xx);
            strainY(i,j,k,yy) = 0.5*(dispy(i,j+1,k,2)-dispy(i,j-1,k,2)) -eg_st(Y,yy);
            strainY(i,j,k,xy) = strainX(i,j,k,xy); 

        });
    }

    strain_X.FillBoundary(geom.periodicity());
	strain_Y.FillBoundary(geom.periodicity());


    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();								//Defining the box for iteration space
		Array4<Real> const& phiOld = phi_old.array(mfi);				//Taking the Multifabs as array
		Array4<Real> const& dispX = disp_X.array(mfi);				//Taking the Multifabs as array
		Array4<Real> const& dispY = disp_Y.array(mfi);
        Array4<Real> const& strainX = strain_X.array(mfi);
        Array4<Real> const& strainY = strain_Y.array(mfi);					//Taking the Multifabs as array
		
        int numphase = nump;

        Array2D <Real, 0, phasecount-1,0,2,Order::C> stiff_phase{};
        for(int g=0; g<nump; g++){
            for(int h=0; h<3; h++){
                stiff_phase(g,h) = stiffness_n[g][h];
            }
        }

        Real deltat_e = deltae;
        Real rho = ro; 
        Real damping_factor = dampfac;

        amrex::ParallelFor(vbx,
		[=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
		{   
            Array2D<Real,0,AMREX_SPACEDIM-1,0,5,Order::C> eg_st{};
            //Array2D<Real,0,2*AMREX_SPACEDIM-1,0,2,Order::C> str_X{};         //strain[X][xx] = 0,0 and so on
            //Array2D<Real,0,2*AMREX_SPACEDIM-1,0,2,Order::C> str_Y{};
            Array2D<Real,0,AMREX_SPACEDIM-1,0,2,Order::C> stiffness_c_ipo{};
            Array2D<Real,0,AMREX_SPACEDIM-1,0,2,Order::C> stiffness_c_imo{};
            Array2D<Real,0,AMREX_SPACEDIM-1,0,2,Order::C> stiffness_c_jpo{};
            Array2D<Real,0,AMREX_SPACEDIM-1,0,2,Order::C> stiffness_c_jmo{};
            Real trace_strain_ipo{0.0}, trace_strain_imo{0.0}, trace_strain_jpo{0.0},trace_strain_jmo{0.0}; 
            Real lambda_ipo{0.0}, lambda_imo{0.0}, lambda_jpo{0.0},lambda_jmo{0.0}; 
            Real mu1_ipo{0.0}, mu1_imo{0.0}, mu1_jpo{0.0},mu1_jmo{0.0}; 
            Real mu1_prime_ipo{0.0}, mu1_prime_imo{0.0}, mu1_prime_jpo{0.0},mu1_prime_jmo{0.0}; 
            Real sigma_xx_ipo{0.0}, sigma_xx_imo{0.0}, sigma_yx_jpo{0.0},sigma_yx_jmo{0.0}, sigma_xy_ipo{0.0}, sigma_xy_imo{0.0}, sigma_yy_jpo{0.0}, sigma_yy_jmo{0.0}; 
            Real ForceX{0.0}, ForceY{0.0};

            calculate_stiffness(i+1,j,k,numphase, stiff_phase,stiffness_c_ipo,phiOld);
            calculate_stiffness(i-1,j,k,numphase, stiff_phase,stiffness_c_imo,phiOld);
            calculate_stiffness(i,j+1,k,numphase, stiff_phase,stiffness_c_jpo,phiOld);
            calculate_stiffness(i,j-1,k,numphase, stiff_phase,stiffness_c_jmo,phiOld);

            trace_strain_ipo = strainX(i+1,j,k,xx)+strainX(i+1,j,k,yy);
            trace_strain_imo = strainX(i-1,j,k,xx)+strainX(i-1,j,k,yy);
            trace_strain_jpo = strainY(i,j+1,k,xx)+strainY(i,j+1,k,yy);
            trace_strain_jmo = strainY(i,j-1,k,xx)+strainY(i,j-1,k,yy);

            lambda_ipo = stiffness_c_ipo(X,C12);             //stiffness(,0)=c11, stiffness(,1)=c12, stiffness(,2)=c44
            lambda_imo = stiffness_c_imo(X,C12);
            lambda_jpo = stiffness_c_jpo(Y,C12);
            lambda_jmo = stiffness_c_jmo(Y,C12);

            mu1_ipo = stiffness_c_ipo(X,C44);
            mu1_imo = stiffness_c_imo(X,C44);
            mu1_jpo = stiffness_c_jpo(Y,C44);
            mu1_jmo = stiffness_c_jmo(Y,C44);

            mu1_prime_ipo = stiffness_c_ipo(X,C11)-stiffness_c_ipo(X,C12)-2.0*stiffness_c_ipo(X,C44);
            mu1_prime_imo = stiffness_c_imo(X,C11)-stiffness_c_imo(X,C12)-2.0*stiffness_c_imo(X,C44);
            mu1_prime_jpo = stiffness_c_jpo(Y,C11)-stiffness_c_jpo(Y,C12)-2.0*stiffness_c_jpo(Y,C44);
            mu1_prime_jmo = stiffness_c_jmo(Y,C11)-stiffness_c_jmo(Y,C12)-2.0*stiffness_c_jmo(Y,C44);

            sigma_xx_ipo = lambda_ipo*(trace_strain_ipo) + 2.0*mu1_ipo*strainX(i+1,j,k,xx) + mu1_prime_ipo*strainX(i+1,j,k,xx);
            sigma_xx_imo = lambda_imo*(trace_strain_imo) + 2.0*mu1_imo*strainX(i-1,j,k,xx) + mu1_prime_imo*strainX(i-1,j,k,xx);

            sigma_yx_jpo = 2.0*mu1_jpo*strainX(i,j+1,k,xy);
            sigma_yx_jmo = 2.0*mu1_jmo*strainX(i,j-1,k,xy);

            sigma_xy_ipo = 2.0*mu1_ipo*strainY(i+1,j,k,xy);
            sigma_xy_imo = 2.0*mu1_imo*strainY(i-1,j,k,xy);

            sigma_yy_jpo = lambda_jpo*(trace_strain_jpo) + 2.0*mu1_jpo*strainY(i,j+1,k,yy) + mu1_prime_jpo*strainY(i,j+1,k,yy);
            sigma_yy_jmo = lambda_jmo*(trace_strain_jmo) + 2.0*mu1_jmo*strainY(i,j-1,k,yy) + mu1_prime_jmo*strainY(i,j-1,k,yy);

            ForceX = 0.5*(sigma_xx_ipo-sigma_xx_imo+sigma_yx_jpo-sigma_yx_jmo);
            ForceY = 0.5*(sigma_xy_ipo-sigma_xy_imo+sigma_yy_jpo-sigma_yy_jmo);

            dispX(i,j,k,0) = dispX(i,j,k,1);
            dispY(i,j,k,0) = dispY(i,j,k,1);

            dispX(i,j,k,1) = dispX(i,j,k,2);
            dispY(i,j,k,1) = dispY(i,j,k,2);

            dispX(i,j,k,2) = (((deltat_e*deltat_e)/rho)*ForceX - (1 - damping_factor*deltat_e)*dispX(i,j,k,0) + 2*dispX(i,j,k,1))/(1.0 + damping_factor*deltat_e);
            dispY(i,j,k,2) = (((deltat_e*deltat_e)/rho)*ForceY - (1 - damping_factor*deltat_e)*dispY(i,j,k,0) + 2*dispY(i,j,k,1))/(1.0 + damping_factor*deltat_e);

        });

    }

    disp_X.FillBoundary(geom.periodicity());
	disp_Y.FillBoundary(geom.periodicity());

}

void Stress_solver_3D(MultiFab& phi_old, MultiFab& disp_X, MultiFab& disp_Y, MultiFab& disp_Z, MultiFab& strain_X, MultiFab& strain_Y, MultiFab& strain_Z, Geometry const& geom){
    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();								//Defining the box for iteration space
		Array4<Real> const& phiOld = phi_old.array(mfi);				//Taking the Multifabs as array
		Array4<Real> const& dispx = disp_X.array(mfi);				//Taking the Multifabs as array
		Array4<Real> const& dispy = disp_Y.array(mfi);					//Taking the Multifabs as array
		Array4<Real> const& dispz = disp_Z.array(mfi);
        Array4<Real> const& strainX = strain_X.array(mfi);
        Array4<Real> const& strainY = strain_Y.array(mfi);
        Array4<Real> const& strainZ = strain_Z.array(mfi);

        int numphase = nump;

        Array2D <Real, 0, phasecount-1,0,6,Order::C> eigenst{};
        for(int g=0; g<nump; g++){
            for(int h=0; h<7; h++){
                eigenst(g,h) = egstr[g][h];
            }
        }

        amrex::ParallelFor(vbx,
		[=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
		{   
            Array2D<Real,0,AMREX_SPACEDIM-1,0,5, Order::C> eg_st{};
                  //strain[X][xx] = 0,0 and so on

            calculate_eigen_strain(i,j,k,numphase,eigenst,eg_st,phiOld);

            strainX(i,j,k,xx) = 0.5*(dispx(i+1,j,k,2)-dispx(i-1,j,k,2)) -eg_st(X,xx);
            strainY(i,j,k,xx) = 0.5*(dispx(i+1,j,k,2)-dispx(i-1,j,k,2)) -eg_st(Y,xx);
            strainZ(i,j,k,xx) = 0.5*(dispx(i+1,j,k,2)-dispx(i-1,j,k,2)) -eg_st(Z,xx);

            strainX(i,j,k,yy) = 0.5*(dispy(i,j+1,k,2)-dispy(i,j-1,k,2)) -eg_st(X,yy);
            strainY(i,j,k,yy) = 0.5*(dispy(i,j+1,k,2)-dispy(i,j-1,k,2)) -eg_st(Y,yy);
            strainZ(i,j,k,yy) = 0.5*(dispy(i,j+1,k,2)-dispy(i,j-1,k,2)) -eg_st(Z,yy);

            strainX(i,j,k,zz) = 0.5*(dispz(i,j,k+1,2)-dispz(i,j,k-1,2)) -eg_st(X,zz);
            strainY(i,j,k,zz) = 0.5*(dispz(i,j,k+1,2)-dispz(i,j,k-1,2)) -eg_st(Y,zz);
            strainZ(i,j,k,zz) = 0.5*(dispz(i,j,k+1,2)-dispz(i,j,k-1,2)) -eg_st(Z,zz);
            
            strainX(i,j,k,xy) = 0.25*(dispx(i,j+1,k,2)-dispx(i,j-1,k,2)+dispy(i+1,j,k,2)-dispy(i-1,j,k,2));
            strainY(i,j,k,xy) = strainX(i,j,k,xy);
            strainZ(i,j,k,xy) = strainX(i,j,k,xy);

            strainX(i,j,k,xz) = 0.25*(dispx(i,j,k+1,2)-dispx(i,j,k-1,2)+dispz(i+1,j,k,2)-dispz(i-1,j,k,2));
            strainY(i,j,k,xz) = strainX(i,j,k,xz);
            strainZ(i,j,k,xz) = strainX(i,j,k,xz);

            strainX(i,j,k,yz) = 0.25*(dispy(i,j,k+1,2)-dispy(i,j,k-1,2)+dispz(i,j+1,k,2)-dispz(i,j-1,k,2));
            strainY(i,j,k,yz) = strainX(i,j,k,yz);
            strainZ(i,j,k,yz) = strainX(i,j,k,yz);

        });
    }

    strain_X.FillBoundary(geom.periodicity());
	strain_Y.FillBoundary(geom.periodicity());
    strain_Z.FillBoundary(geom.periodicity());


    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();								//Defining the box for iteration space
		Array4<Real> const& phiOld = phi_old.array(mfi);				//Taking the Multifabs as array
		Array4<Real> const& dispX = disp_X.array(mfi);				//Taking the Multifabs as array
		Array4<Real> const& dispY = disp_Y.array(mfi);
        Array4<Real> const& dispZ = disp_Z.array(mfi);
        Array4<Real> const& strainX = strain_X.array(mfi);
        Array4<Real> const& strainY = strain_Y.array(mfi);					//Taking the Multifabs as array
		Array4<Real> const& strainZ = strain_Z.array(mfi);

        int numphase = nump;

        Array2D <Real, 0, phasecount-1,0,2,Order::C> stiff_phase{};
        for(int g=0; g<nump; g++){
            for(int h=0; h<4; h++){
                stiff_phase(g,h) = stiffness_n[g][h];
            }
        }

        Real deltat_e = deltae;
        Real rho = ro; 
        Real damping_factor = dampfac;

        amrex::ParallelFor(vbx,
		[=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
		{   
            Array2D<Real,0,AMREX_SPACEDIM-1,0,5,Order::C> eg_st{};
            //Array2D<Real,0,2*AMREX_SPACEDIM-1,0,2,Order::C> str_X{};         //strain[X][xx] = 0,0 and so on
            //Array2D<Real,0,2*AMREX_SPACEDIM-1,0,2,Order::C> str_Y{};
            Array2D<Real,0,AMREX_SPACEDIM-1,0,2,Order::C> stiffness_c_ipo{};
            Array2D<Real,0,AMREX_SPACEDIM-1,0,2,Order::C> stiffness_c_imo{};
            Array2D<Real,0,AMREX_SPACEDIM-1,0,2,Order::C> stiffness_c_jpo{};
            Array2D<Real,0,AMREX_SPACEDIM-1,0,2,Order::C> stiffness_c_jmo{};
            Array2D<Real,0,AMREX_SPACEDIM-1,0,2,Order::C> stiffness_c_kpo{};
            Array2D<Real,0,AMREX_SPACEDIM-1,0,2,Order::C> stiffness_c_kmo{};
            Real trace_strain_ipo{0.0}, trace_strain_imo{0.0}, trace_strain_jpo{0.0},trace_strain_jmo{0.0}, trace_strain_kpo{0.0}, trace_strain_kmo{0.0}; 
            Real lambda_ipo{0.0}, lambda_imo{0.0}, lambda_jpo{0.0},lambda_jmo{0.0},lambda_kpo{0.0},lambda_kmo{0.0}; 
            Real mu1_ipo{0.0}, mu1_imo{0.0}, mu1_jpo{0.0},mu1_jmo{0.0},mu1_kpo{0.0},mu1_kmo{0.0}; 
            Real mu1_prime_ipo{0.0}, mu1_prime_imo{0.0}, mu1_prime_jpo{0.0},mu1_prime_jmo{0.0}, mu1_prime_kpo{0.0},mu1_prime_kmo{0.0}; 
            Real sigma_xx_ipo{0.0}, sigma_xx_imo{0.0}, sigma_xy_ipo{0.0}, sigma_xy_imo{0.0}, sigma_xz_ipo{0.0}, sigma_xz_imo{0.0};
            Real sigma_yy_jpo{0.0}, sigma_yy_jmo{0.0}, sigma_yx_jpo{0.0}, sigma_yx_jmo{0.0}, sigma_yz_jpo{0.0}, sigma_yz_jmo{0.0};
            Real sigma_zz_kpo{0.0}, sigma_zz_kmo{0.0}, sigma_zx_kpo{0.0}, sigma_zx_kmo{0.0}, sigma_zy_kpo{0.0}, sigma_zy_kmo{0.0}; 
            Real ForceX{0.0}, ForceY{0.0},  ForceZ{0.0};

            calculate_stiffness(i+1,j,k,numphase, stiff_phase,stiffness_c_ipo,phiOld);
            calculate_stiffness(i-1,j,k,numphase, stiff_phase,stiffness_c_imo,phiOld);
            calculate_stiffness(i,j+1,k,numphase, stiff_phase,stiffness_c_jpo,phiOld);
            calculate_stiffness(i,j-1,k,numphase, stiff_phase,stiffness_c_jmo,phiOld);
            calculate_stiffness(i,j,k+1,numphase, stiff_phase,stiffness_c_kpo,phiOld);
            calculate_stiffness(i,j,k-1,numphase, stiff_phase,stiffness_c_kmo,phiOld);

            trace_strain_ipo = strainX(i+1,j,k,xx)+strainX(i+1,j,k,yy)+strainX(i+1,j,k,zz);
            trace_strain_imo = strainX(i-1,j,k,xx)+strainX(i-1,j,k,yy)+strainX(i-1,j,k,zz);
            trace_strain_jpo = strainY(i,j+1,k,xx)+strainY(i,j+1,k,yy)+strainY(i,j+1,k,zz);
            trace_strain_jmo = strainY(i,j-1,k,xx)+strainY(i,j-1,k,yy)+strainY(i,j-1,k,zz);
            trace_strain_kpo = strainZ(i,j,k+1,xx)+strainZ(i,j,k+1,yy)+strainZ(i,j,k+1,zz);
            trace_strain_kmo = strainZ(i,j,k-1,xx)+strainZ(i,j,k-1,yy)+strainZ(i,j,k-1,zz);

            lambda_ipo = stiffness_c_ipo(X,C12);             //stiffness(,0)=c11, stiffness(,1)=c12, stiffness(,2)=c44
            lambda_imo = stiffness_c_imo(X,C12);
            lambda_jpo = stiffness_c_jpo(Y,C12);
            lambda_jmo = stiffness_c_jmo(Y,C12);
            lambda_kpo = stiffness_c_kpo(Z,C12);
            lambda_kmo = stiffness_c_kmo(Z,C12);

            mu1_ipo = stiffness_c_ipo(X,C44);
            mu1_imo = stiffness_c_imo(X,C44);
            mu1_jpo = stiffness_c_jpo(Y,C44);
            mu1_jmo = stiffness_c_jmo(Y,C44);
            mu1_kpo = stiffness_c_kpo(Z,C44);
            mu1_kmo = stiffness_c_kmo(Z,C44);

            mu1_prime_ipo = stiffness_c_ipo(X,C11)-stiffness_c_ipo(X,C12)-2.0*stiffness_c_ipo(X,C44);
            mu1_prime_imo = stiffness_c_imo(X,C11)-stiffness_c_imo(X,C12)-2.0*stiffness_c_imo(X,C44);
            mu1_prime_jpo = stiffness_c_jpo(Y,C11)-stiffness_c_jpo(Y,C12)-2.0*stiffness_c_jpo(Y,C44);
            mu1_prime_jmo = stiffness_c_jmo(Y,C11)-stiffness_c_jmo(Y,C12)-2.0*stiffness_c_jmo(Y,C44);
            mu1_prime_kpo = stiffness_c_kpo(Z,C11)-stiffness_c_kpo(Z,C12)-2.0*stiffness_c_kpo(Z,C44);
            mu1_prime_kmo = stiffness_c_kmo(Z,C11)-stiffness_c_kmo(Z,C12)-2.0*stiffness_c_kmo(Z,C44);

            sigma_xx_ipo = lambda_ipo*(trace_strain_ipo) + 2.0*mu1_ipo*strainX(i+1,j,k,xx) + mu1_prime_ipo*strainX(i+1,j,k,xx);
            sigma_xx_imo = lambda_imo*(trace_strain_imo) + 2.0*mu1_imo*strainX(i-1,j,k,xx) + mu1_prime_imo*strainX(i-1,j,k,xx);

            sigma_yx_jpo = 2.0*mu1_jpo*strainX(i,j+1,k,xy);
            sigma_yx_jmo = 2.0*mu1_jmo*strainX(i,j-1,k,xy);
            
            sigma_zx_kpo = 2.0*mu1_kpo*strainX(i,j,k+1,xz);
            sigma_zx_kmo = 2.0*mu1_kmo*strainX(i,j,k-1,xz);

            sigma_xz_ipo = 2.0*mu1_ipo*strainZ(i+1,j,k,xz);
            sigma_xz_imo = 2.0*mu1_imo*strainZ(i-1,j,k,xz);

            sigma_xy_ipo = 2.0*mu1_ipo*strainY(i+1,j,k,xy);
            sigma_xy_imo = 2.0*mu1_imo*strainY(i-1,j,k,xy);

            sigma_zy_kpo = 2.0*mu1_kpo*strainY(i,j,k+1,yz);
            sigma_zy_kmo = 2.0*mu1_kmo*strainY(i,j,k-1,yz);

            sigma_yz_jpo = 2.0*mu1_jpo*strainZ(i,j+1,k,yz);
            sigma_yz_jmo = 2.0*mu1_jmo*strainZ(i,j-1,k,yz); 

            sigma_yy_jpo = lambda_jpo*(trace_strain_jpo) + 2.0*mu1_jpo*strainY(i,j+1,k,yy) + mu1_prime_jpo*strainY(i,j+1,k,yy);
            sigma_yy_jmo = lambda_jmo*(trace_strain_jmo) + 2.0*mu1_jmo*strainY(i,j-1,k,yy) + mu1_prime_jmo*strainY(i,j-1,k,yy);

            sigma_zz_kpo = lambda_kpo*(trace_strain_kpo) + 2.0*mu1_kpo*strainZ(i,j,k+1,zz) + mu1_prime_kpo*strainZ(i,j,k+1,zz);
            sigma_zz_kmo = lambda_kmo*(trace_strain_kmo) + 2.0*mu1_kmo*strainZ(i,j,k-1,zz) + mu1_prime_kmo*strainZ(i,j,k-1,zz);

            ForceX = 0.5*((sigma_xx_ipo-sigma_xx_imo) + (sigma_yx_jpo-sigma_yx_jmo) + (sigma_zx_kpo-sigma_zx_kmo));
            ForceY = 0.5*((sigma_xy_ipo-sigma_xy_imo) + (sigma_yy_jpo-sigma_yy_jmo) + (sigma_zy_kpo-sigma_zy_kmo));
            ForceZ = 0.5*((sigma_xz_ipo-sigma_xz_imo) + (sigma_yz_jpo-sigma_yz_jmo) + (sigma_zz_kpo-sigma_zz_kmo));

            dispX(i,j,k,0) = dispX(i,j,k,1);
            dispY(i,j,k,0) = dispY(i,j,k,1);
            dispZ(i,j,k,0) = dispZ(i,j,k,1);

            dispX(i,j,k,1) = dispX(i,j,k,2);
            dispY(i,j,k,1) = dispY(i,j,k,2);
            dispZ(i,j,k,1) = dispZ(i,j,k,2);

            dispX(i,j,k,2) = (((deltat_e*deltat_e)/rho)*ForceX - (1 - damping_factor*deltat_e)*dispX(i,j,k,0) + 2*dispX(i,j,k,1))/(1.0 + damping_factor*deltat_e);
            dispY(i,j,k,2) = (((deltat_e*deltat_e)/rho)*ForceY - (1 - damping_factor*deltat_e)*dispY(i,j,k,0) + 2*dispY(i,j,k,1))/(1.0 + damping_factor*deltat_e);
            dispZ(i,j,k,2) = (((deltat_e*deltat_e)/rho)*ForceZ - (1 - damping_factor*deltat_e)*dispZ(i,j,k,0) + 2*dispZ(i,j,k,1))/(1.0 + damping_factor*deltat_e);

        });

    }

    disp_X.FillBoundary(geom.periodicity());
	disp_Y.FillBoundary(geom.periodicity());
    disp_Z.FillBoundary(geom.periodicity());
    
}


#endif