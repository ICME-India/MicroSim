#ifndef _FUNCTION_A_H_
#define _FUNCTION_A_H_

using namespace amrex;

#include "Function_Q.h"
#include "Anisotropy_01.h"
#include "GP_Utility.h"


//iph = i + 1/2
//imh = i - 1/2
//jph = j + 1/2
//jmh = j - 1/2
//kph = k + 1/2
//kmh = k - 1/2



//2-D isotropic function---------------------------------------------------------------------------------- 
void function_A_00_iso_2D(MultiFab& term1, MultiFab& phi_old, Geometry const& geom){
    
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif

	for(MFIter mfi(phi_old); mfi.isValid(); ++mfi)
	{
		const Box& bbx = mfi.validbox();									//Defining the box for iteration space
		auto const& phiOld = phi_old.const_array(mfi);						//Taking the MultiFab as arrays
		auto const& term1_val = term1.array(mfi);							//Taking the MultiFab as arrays

		//Turning the vector gamma to a GPU compatible array----------------------------------------------------------
		Array2D <Real,0,phasecount,0,phasecount,Order::C> gamm{};
		for(int a=0; a<nump; a++){
			for(int b=a; b<nump; b++){
				gamm(a,b) = gam[a][b];
				gamm(b,a) = gam[a][b];
			}
		}

		//Redifing variables in GPU space------------------------------------------------------------------------------
		Real numphase = nump;

		//delta stores dx,dy and dz------------------------------------------------------------------------------------
		GpuArray<Real,AMREX_SPACEDIM> delta = geom.CellSizeArray();

		//Iteration----------------------------------------------------------------------------------------------------
		amrex::ParallelFor(bbx,
		[=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
		{	
			
			//Defining variables---------------------------------------------------------------------------------------
			Array2D <Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> qab{};
			Array1D <Real, 0, AMREX_SPACEDIM-1> dqdphi{};
			amrex::Real dAdphi{0.0}, divdAdgradphi{0.0};
			//Array1D<Real,0,phasecount-1> div{};
			

			for(int a=0; a<numphase; a++){
				
					for(int b=0; b<numphase; b++){
							
							if(b!=a){
								
								//Calculating qab----------------------------------------------------------------------------------
								Function_Q(i,j,k, phiOld, delta, a, b, qab);
								
								//Calculating \frac{\partial \phi _{b}}{\partial x} ----------------------------------------------------------------------------------
								dqdphi(X) = (phiOld(i+1,j,k,b)-phiOld(i-1,j,k,b))/(2.0*delta[X]);

								//Calculating \frac{\partial \phi _{b}}{\partial y} ----------------------------------------------------------------------------------
								dqdphi(Y) = (phiOld(i,j+1,k,b)-phiOld(i,j-1,k,b))/(2.0*delta[Y]);

								//Calculating 2.0\cdot \gamma _{ab} \left (q_{ab}^{x}\cdot \frac{\partial \phi _{b}}{\partial x} + q_{ab}^{y}\cdot \frac{\partial \phi _{b}}{\partial y} \right ) ----------------------------------------------------------------------------------
								dAdphi += 2.0*gamm(a,b)*(qab(cent,X)*dqdphi(X)+qab(cent,Y)*dqdphi(Y));

								//Calculating \frac{\partial }{\partial x}\left ( 2.0 \cdot \gamma_{ab} \cdot \left ( - \phi_{b} \right )\cdot q_{ab}^{x} \right ) + \frac{\partial }{\partial y}\left ( 2.0 \cdot \gamma_{ab} \cdot \left ( - \phi_{b} \right )\cdot q_{ab}^{y} \right ) -----------------------------------
								divdAdgradphi += 2.0*gamm(a,b)*(qab(iph,X)*(-(phiOld(i+1,j,k,b)+phiOld(i,j,k,b))*0.5)-qab(imh,X)*(-(phiOld(i,j,k,b)+phiOld(i-1,j,k,b))*0.5))/delta[X] + 2.0*gamm(a,b)*(qab(jph,Y)*(-(phiOld(i,j+1,k,b)+phiOld(i,j,k,b))*0.5)-qab(jmh,Y)*(-(phiOld(i,j,k,b)+phiOld(i,j-1,k,b))*0.5))/delta[Y];									
							
						}
						
					}
					//Isotropic term for 2-D-------------------
					term1_val(i,j,k,a) = -dAdphi+divdAdgradphi;
					dAdphi = 0.0;
					divdAdgradphi = 0.0;
				
			}
		});
	}
}


//2-D isotropic function---------------------------------------------------------------------------------- 
void function_A_00_iso_3D(MultiFab& term1, MultiFab& phi_old, Geometry const& geom){

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif

	for(MFIter mfi(phi_old); mfi.isValid(); ++mfi)
	{
		const Box& bbx = mfi.validbox();										//Defining the box for iteration space
		auto const& phiOld = phi_old.const_array(mfi);							//Taking the MultiFab as arrays
		auto const& term1_val = term1.array(mfi);								//Taking the MultiFab as arrays

		//Turning the vector gamma to a GPU compatible array----------------------------------------------------------
		Array2D <Real,0,phasecount,0,phasecount,Order::C> gamm{};
		for(int a=0; a<nump; a++){
			for(int b=a; b<nump; b++){
				gamm(a,b) = gam[a][b];
				gamm(b,a) = gam[a][b];
			}
		}

		//Redifing variables in GPU space------------------------------------------------------------------------------
		Real numphase = nump;

		//delta stores dx,dy and dz------------------------------------------------------------------------------------
		GpuArray<Real,AMREX_SPACEDIM> delta = geom.CellSizeArray();
		

		//Iteration----------------------------------------------------------------------------------------------------
		amrex::ParallelFor(bbx,
		[=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
		{	
			
			//Defining variables---------------------------------------------------------------------------------------
			Array2D <Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> qab{};
			Array1D <Real, 0, AMREX_SPACEDIM-1> dqdphi{};
			amrex::Real dAdphi{0.0}, divdAdgradphi{0.0};
			
			
			for(int a=0; a<numphase; a++){
					
					for(int b=0; b<numphase; b++){
							
							if(b!=a){

								//Calculating qab----------------------------------------------------------------------------------	
								Function_Q(i,j,k, phiOld, delta, a, b, qab);

								//Calculating \frac{\partial \phi _{b}}{\partial x} ----------------------------------------------------------------------------------
								dqdphi(X) = (phiOld(i+1,j,k,b)-phiOld(i-1,j,k,b))/(2.0*delta[X]);
								
								//Calculating \frac{\partial \phi _{b}}{\partial x} ----------------------------------------------------------------------------------
								dqdphi(Y) = (phiOld(i,j+1,k,b)-phiOld(i,j-1,k,b))/(2.0*delta[Y]);
								
								//Calculating \frac{\partial \phi _{b}}{\partial z} ----------------------------------------------------------------------------------
								dqdphi(Z) = (phiOld(i,j,k+1,b)-phiOld(i,j,k-1,b))/(2.0*delta[Z]);

								//Calculating 2.0\cdot \gamma _{ab} \left (q_{ab}^{x}\cdot \frac{\partial \phi _{b}}{\partial x} + q_{ab}^{y}\cdot \frac{\partial \phi _{b}}{\partial y} + q_{ab}^{z}\cdot \frac{\partial \phi _{b}}{\partial z}\right ) ----------------------------------------------------------------------------------
								dAdphi += 2.0*gamm(a,b)*(qab(cent,X)*dqdphi(X) + qab(cent,Y)*dqdphi(Y) + qab(cent,Z)*dqdphi(Z));
								
								//Calculating \frac{\partial }{\partial x}\left ( 2.0 \cdot \gamma_{ab} \cdot \left ( - \phi_{b} \right )\cdot q_{ab}^{x} \right ) + \frac{\partial }{\partial y}\left ( 2.0 \cdot \gamma_{ab} \cdot \left ( - \phi_{b} \right )\cdot q_{ab}^{y} \right ) + \frac{\partial }{\partial z}\left ( 2.0 \cdot \gamma_{ab} \cdot \left ( - \phi_{b} \right )\cdot q_{ab}^{z} \right ) ----------------------------------------------------
								divdAdgradphi += 2.0*gamm(a,b)*(qab(iph,X)*(-(phiOld(i+1,j,k,b)+phiOld(i,j,k,b))*0.5)-qab(imh,X)*(-(phiOld(i,j,k,b)+phiOld(i-1,j,k,b))*0.5))/delta[X] 
											   + 2.0*gamm(a,b)*(qab(jph,Y)*(-(phiOld(i,j+1,k,b)+phiOld(i,j,k,b))*0.5)-qab(jmh,Y)*(-(phiOld(i,j,k,b)+phiOld(i,j-1,k,b))*0.5))/delta[Y] 
											   + 2.0*gamm(a,b)*(qab(kph,Z)*(-(phiOld(i,j,k+1,b)+phiOld(i,j,k,b))*0.5)-qab(kmh,Z)*(-(phiOld(i,j,k,b)+phiOld(i,j,k-1,b))*0.5))/delta[Z];										
							
						}
						
					}
					//Isotropic term for 3-D-------------------
					term1_val(i,j,k,a) = -dAdphi+divdAdgradphi;
					dAdphi = 0.0;
					divdAdgradphi = 0.0;
			}
		});
	}
}



//2-D anisotropic function---------------------------------------------------------------------------------- 
void function_A_01_ani_2D(MultiFab& term1, MultiFab& phi_old, Geometry const& geom){

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif

	for(MFIter mfi(phi_old); mfi.isValid(); ++mfi)
	{	
		const Box& bbx = mfi.validbox();										//Defining the box for iteration space
		auto const& phiOld = phi_old.const_array(mfi);							//Taking the MultiFab as arrays
		auto const& term1_val = term1.array(mfi);								//Taking the MultiFab as arrays

		//Turning the vector rotmat to a GPU compatible array----------------------------------------------------------
		Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_x{};
		Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_y{};
		Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_z{};
		for(int a=0; a<rotmat.size(); a++){
				int m = rotmat[a][0];
				int n = rotmat[a][1]; 
				matrot_x(m,n) = rotmat[a][2]*acos(-1)/180.0;
				matrot_x(n,m) = rotmat[a][2]*acos(-1)/180.0;
				matrot_y(m,n) = rotmat[a][3]*acos(-1)/180.0;
				matrot_y(n,m) = rotmat[a][3]*acos(-1)/180.0;
				matrot_z(m,n) = rotmat[a][4]*acos(-1)/180.0;
				matrot_z(n,m) = rotmat[a][4]*acos(-1)/180.0;
		
		}

		//Turning the vector gamma to a GPU compatible array----------------------------------------------------------
		Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> gamm{};
		for(int a=0; a<nump; a++){
			for(int b=a; b<nump; b++){
				gamm(a,b) = gam[a][b];
				gamm(b,a) = gam[a][b];
			}
		}

		//Turning the vector dabb to a GPU compatible array----------------------------------------------------------
		Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> dabbb{};
		for(int a=0; a<nump; a++){
			for(int b=a; b<nump; b++){
				dabbb(a,b) = dabb[a][b];
				dabbb(b,a) = dabb[a][b];
			}
		}

		//Redifing variables in GPU space------------------------------------------------------------------------------
		Real numphase = nump;
		int anisotropy_type = ANItype;

		//delta stores dx,dy and dz------------------------------------------------------------------------------------
		GpuArray<Real,AMREX_SPACEDIM> delta = geom.CellSizeArray();


		//Iteration----------------------------------------------------------------------------------------------------
		amrex::ParallelFor(bbx,
		[=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
		{	
			//Defining variables---------------------------------------------------------------------------------------
			Array2D <Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> qab{};
			Array2D <Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> r_qab{};
			// Array2D <Real,0,6,0,2,Order::C> qab{};
			// Array2D <Real,0,6,0,2,Order::C> r_qab{};
			Array2D <Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> dadq{};
			Array2D <Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> r_dadq{};
			Array1D <Real,0,AMREX_SPACEDIM*2> ac{};
			Array1D <Real,0,AMREX_SPACEDIM*2> q2{};
			Array1D <Real, 0, AMREX_SPACEDIM-1> dqdphi{};
			Array1D <Real, 0, AMREX_SPACEDIM-1> r_dqdphi{};
			// Array1D <Real, 0, 2> dqdphi{};
			// Array1D <Real, 0, 2> r_dqdphi{};
			amrex::Real dAdphi1{0.0}, dAdphi2{0.0}, divdAdgradphi1{0.0}, divdAdgradphi2{0.0};


			for(int a=0; a<numphase; a++){

					for(int b=0; b<numphase; b++){

							if(b!=a){
								
								//Calculating qab----------------------------------------------------------------------------------	
								Function_Q(i, j, k, phiOld, delta, a, b, qab);
						
								//Rotation of normal vector by the relevant angle--------------------------------------------------
								Vec_rot(i, j, k, a, b, qab, r_qab, matrot_x, matrot_y, matrot_z);

								if(anisotropy_type == 4){
									
									//Calculate dadq-------------------------------------------------------------------------------
									dAdq(i, j, k, dabbb, r_qab, dadq, q2, a, b);

									//Calculate ac---------------------------------------------------------------------------------
									function_ac(i, j, k, dabbb, r_qab, ac, q2, a, b);
								}

								//Calculating \frac{\partial \phi _{b}}{\partial x} ----------------------------------------------------------------------------------
								dqdphi(X) = (phiOld(i+1,j,k,b)-phiOld(i-1,j,k,b))/(2*delta[X]);

								//Calculating \frac{\partial \phi _{b}}{\partial y} ----------------------------------------------------------------------------------
								dqdphi(Y) = (phiOld(i,j+1,k,b)-phiOld(i,j-1,k,b))/(2*delta[Y]);
						

								//Rotating dqdphi by the desired angle ---------------------------------------------------------------
								Vec_rot(i, j, k, a, b, dqdphi, r_dqdphi, matrot_x, matrot_y, matrot_z);


								//Calculating 2.0\cdot \gamma _{ab} \cdot a_{c}^{2} \cdot \left (q_{ab}^{x}\cdot \frac{\partial \phi _{b}}{\partial x} + q_{ab}^{y}\cdot \frac{\partial \phi _{b}}{\partial y} \right ) ----------------------------------------------------------------------------------
								dAdphi1 += 2.0*gamm(a,b)*ac(cent)*ac(cent)*(r_qab(cent,X)*r_dqdphi(X)+r_qab(cent,Y)*r_dqdphi(Y));

								//Calculating 2.0 \gamma _{ab} a_{c} \left | q_{ab} \right |^{2} \left ( \left ( \frac{\partial a_{c}}{\partial q_{ab}}  \right ) \cdot \nabla\phi_{b} \right) -------------------------------------------------------
								dAdphi2 += 2.0*gamm(a,b)*ac(cent)*q2(cent)*(dadq(cent,X)*r_dqdphi(X)+dadq(cent,Y)*r_dqdphi(Y));
							

								//Inverse rotation of dadq by the desired angle ---------------------------------------------------------------------------
								Inv_vec_rot(i,j,k,a,b,dadq,r_dadq, matrot_x, matrot_y, matrot_z);
								
								//Calculating \frac{\partial }{\partial x}\left ( 2.0 \cdot \gamma_{ab} \cdot a_{c}^{2}\cdot \left ( - \phi_{b} \right )\cdot q_{ab}^{x} \right ) + \frac{\partial }{\partial y}\left ( 2.0 \cdot \gamma_{ab} \cdot a_{c}^{2} \cdot \left ( - \phi_{b} \right )\cdot q_{ab}^{y} \right ) \right ) ----------------------------------------------------------------------
								divdAdgradphi1 += 2*gamm(a,b)*(ac(iph)*ac(iph)*qab(iph,X)*(-(phiOld(i+1,j,k,b)+phiOld(i,j,k,b))*0.5)-ac(imh)*ac(imh)*qab(imh,X)*(-(phiOld(i,j,k,b)+phiOld(i-1,j,k,b))*0.5))/delta[X] + 2*gamm(a,b)*(ac(jph)*ac(jph)*qab(jph,Y)*(-(phiOld(i,j+1,k,b)+phiOld(i,j,k,b))*0.5)-ac(jmh)*ac(jmh)*qab(jmh,Y)*(-(phiOld(i,j,k,b)+phiOld(i,j-1,k,b))*0.5))/delta[Y];
								
								//Calculating \frac{\partial }{\partial x}\left ( 2.0 \cdot \gamma_{ab} \cdot a_{c} \cdot \left ( -\phi _{b} \right ) \cdot \left | q_{ab} \right |^{2} \cdot \left( \frac{\partial a_{c}}{\partial q_{ab}}\right)_{x} \right ) +  \frac{\partial }{\partial y}\left ( 2.0 \cdot \gamma_{ab} \cdot a_{c} \cdot \left ( -\phi _{b} \right ) \cdot \left | q_{ab} \right |^{2} \cdot \left( \frac{\partial a_{c}}{\partial q_{ab}}\right)_{y} \right ) ---------------
								divdAdgradphi2 += 2*gamm(a,b)*(ac(iph)*q2(iph)*(-(phiOld(i+1,j,k,b)+phiOld(i,j,k,b))*0.5)*r_dadq(iph,X)-ac(imh)*q2(imh)*(-(phiOld(i,j,k,b)+phiOld(i-1,j,k,b))*0.5)*r_dadq(imh,X))/delta[X] + 2*gamm(a,b)*(ac(jph)*q2(jph)*(-(phiOld(i,j+1,k,b)+phiOld(i,j,k,b))*0.5)*r_dadq(jph,Y)-ac(jmh)*q2(jmh)*(-(phiOld(i,j,k,b)+phiOld(i,j-1,k,b))*0.5)*r_dadq(jmh,Y))/delta[Y];


				 			}
					}

				//Final anisotropy term for 2-D
				term1_val(i,j,k,a) = (-dAdphi1-dAdphi2+divdAdgradphi1+divdAdgradphi2);

				dAdphi1=0.0;
				dAdphi2=0.0;
				divdAdgradphi1=0.0;
				divdAdgradphi2=0.0;

			}
		
	    });
	}
}


//3-D anisotropic function----------------------------------------------------------------------------------
void function_A_01_ani_3D(MultiFab& term1, MultiFab& phi_old, Geometry const& geom){

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif

	for(MFIter mfi(phi_old); mfi.isValid(); ++mfi)
	{
		const Box& bbx = mfi.validbox();										//Defining the box for iteration space
		auto const& phiOld = phi_old.const_array(mfi);							//Taking the MultiFabs as arrays
		auto const& term1_val = term1.array(mfi);								//Taking the MultiFabs as arrays

		//Turning the vector rotmat to a GPU compatible array----------------------------------------------------------
		Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_x{};
		Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_y{};
		Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_z{};
		for(int a=0; a<rotmat.size(); a++){
				int m = rotmat[a][0];
				int n = rotmat[a][1]; 
				matrot_x(m,n) = rotmat[a][2]*acos(-1)/180.0;
				matrot_x(n,m) = rotmat[a][2]*acos(-1)/180.0;
				matrot_y(m,n) = rotmat[a][3]*acos(-1)/180.0;
				matrot_y(n,m) = rotmat[a][3]*acos(-1)/180.0;
				matrot_z(m,n) = rotmat[a][4]*acos(-1)/180.0;
				matrot_z(n,m) = rotmat[a][4]*acos(-1)/180.0;
		
		}

		//Turning the vector gamma to a GPU compatible array----------------------------------------------------------
		Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> gamm{};
		for(int a=0; a<nump; a++){
			for(int b=a; b<nump; b++){
				gamm(a,b) = gam[a][b];
				gamm(b,a) = gam[a][b];
			}
		}

		//Turning the vector dabb to a GPU compatible array----------------------------------------------------------
		Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> dabbb{};
		for(int a=0; a<nump; a++){
			for(int b=a; b<nump; b++){
				dabbb(a,b) = dabb[a][b];
				dabbb(b,a) = dabb[a][b];
			}
		}

		//Redifing variables in GPU space------------------------------------------------------------------------------
		Real numphase = nump;
		int anisotropy_type = ANItype;

		//delta stores dx,dy and dz------------------------------------------------------------------------------------
		GpuArray<Real,AMREX_SPACEDIM> delta = geom.CellSizeArray();


		//Iteration ---------------------------------------------------------------------------------------------------
		amrex::ParallelFor(bbx,
		[=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
		{	
			//Defining variables---------------------------------------------------------------------------------------
			Array2D <Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> qab{};
			Array2D <Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> r_qab{};
			Array2D <Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> dadq{};
			Array2D <Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> r_dadq{};
			Array1D <Real,0,AMREX_SPACEDIM*2> ac{};
			Array1D <Real,0,AMREX_SPACEDIM*2> q2{};
			Array1D <Real, 0, AMREX_SPACEDIM-1> dqdphi{};
			Array1D <Real, 0, AMREX_SPACEDIM-1> r_dqdphi{};

			amrex::Real dAdphi1{0.0}, dAdphi2{0.0}, divdAdgradphi1{0.0}, divdAdgradphi2{0.0};

			for(int a=0; a<numphase; a++){

					for(int b=0; b<numphase; b++){

							if(b!=a){
								//Print()<<"in 3d loop\n";
								//Calculating qab----------------------------------------------------------------------------------	
								Function_Q(i, j, k, phiOld, delta, a, b, qab);
								
								//Rotation of normal vector by the relevant angle--------------------------------------------------
								Vec_rot(i, j, k, a, b, qab, r_qab, matrot_x, matrot_y, matrot_z);


								if(anisotropy_type == 6){
									//Calculate dadq ------------------------------------------------------------------------------
									dAdq(i, j, k, dabbb, r_qab, dadq, q2, a, b);

									//Calculate ac --------------------------------------------------------------------------------
									function_ac(i, j, k, dabbb, r_qab, ac, q2, a, b);
								}

								//Calculating \frac{\partial \phi _{b}}{\partial x} ----------------------------------------------------------------------------------
								dqdphi(X) = (phiOld(i+1,j,k,b)-phiOld(i-1,j,k,b))/(2*delta[X]);

								//Calculating \frac{\partial \phi _{b}}{\partial y} ----------------------------------------------------------------------------------
								dqdphi(Y) = (phiOld(i,j+1,k,b)-phiOld(i,j-1,k,b))/(2*delta[Y]);

								//Calculating \frac{\partial \phi _{b}}{\partial z} ----------------------------------------------------------------------------------
								dqdphi(Z) = (phiOld(i,j,k+1,b)-phiOld(i,j,k-1,b))/(2*delta[Z]);
							

								//Rotating dqdphi by the desired angle ---------------------------------------------------------------
								Vec_rot(i, j, k, a, b, dqdphi, r_dqdphi, matrot_x, matrot_y, matrot_z);


								//Calculating 2.0\cdot \gamma _{ab} \cdot a_{c}^{2} \cdot \left (q_{ab}^{x}\cdot \frac{\partial \phi _{b}}{\partial x} + q_{ab}^{y}\cdot \frac{\partial \phi _{b}}{\partial y} + q_{ab}^{z}\cdot \frac{\partial \phi _{b}}{\partial z} \right ) -------------------------------------------
								dAdphi1 += 2.0*gamm(a,b)*ac(cent)*ac(cent)*(r_qab(cent,X)*r_dqdphi(X)+r_qab(cent,Y)*r_dqdphi(Y))+2.0*gamm(a,b)*ac(cent)*ac(cent)*(r_qab(cent,Z)*r_dqdphi(Z));

								//Calculating 2.0 \gamma _{ab} a_{c} \left | q_{ab} \right |^{2} \left ( \left ( \frac{\partial a_{c}}{\partial q_{ab}}  \right ) \cdot \nabla\phi_{b} \right) -------------------------------------------------------
								dAdphi2 += 2.0*gamm(a,b)*ac(cent)*q2(cent)*(dadq(cent,X)*r_dqdphi(X)+dadq(cent,Y)*r_dqdphi(Y))+2.0*gamm(a,b)*ac(cent)*q2(cent)*(dadq(cent,Z)*r_dqdphi(Z));
							

								//Inverse rotation of dadq by the desired angle ---------------------------------------------------------------------------
								Inv_vec_rot(i,j,k,a,b,dadq,r_dadq, matrot_x, matrot_y, matrot_z);
								
								//Calculating \frac{\partial }{\partial x}\left ( 2.0 \cdot \gamma_{ab} \cdot a_{c}^{2}\cdot \left ( - \phi_{b} \right )\cdot q_{ab}^{x} \right ) + \frac{\partial }{\partial y}\left ( 2.0 \cdot \gamma_{ab} \cdot a_{c}^{2} \cdot \left ( - \phi_{b} \right )\cdot q_{ab}^{y} \right ) + \frac{\partial }{\partial z}\left ( 2.0 \cdot \gamma_{ab} \cdot a_{c}^{2} \cdot \left ( - \phi_{b} \right )\cdot q_{ab}^{z} \right ) ------------------------------------------------------------------------------------------------
								divdAdgradphi1 += 2*gamm(a,b)*(ac(iph)*ac(iph)*qab(iph,X)*(-(phiOld(i+1,j,k,b)+phiOld(i,j,k,b))*0.5)-ac(imh)*ac(imh)*qab(imh,X)*(-(phiOld(i,j,k,b)+phiOld(i-1,j,k,b))*0.5))/delta[X] + 2*gamm(a,b)*(ac(jph)*ac(jph)*qab(jph,Y)*(-(phiOld(i,j+1,k,b)+phiOld(i,j,k,b))*0.5)-ac(jmh)*ac(jmh)*qab(jmh,Y)*(-(phiOld(i,j,k,b)+phiOld(i,j-1,k,b))*0.5))/delta[Y]+ 2*gamm(a,b)*(ac(kph)*ac(kph)*qab(kph,Z)*(-(phiOld(i,j,k+1,b)+phiOld(i,j,k,b))*0.5)-ac(kmh)*ac(kmh)*qab(kmh,Z)*(-(phiOld(i,j,k,b)+phiOld(i,j,k-1,b))*0.5))/delta[Z];
							
								//Calculating \frac{\partial }{\partial x}\left ( 2.0 \cdot \gamma_{ab} \cdot a_{c} \cdot \left ( -\phi _{b} \right ) \cdot \left | q_{ab} \right |^{2} \cdot \left( \frac{\partial a_{c}}{\partial q_{ab}}\right)_{x} \right ) +  \frac{\partial }{\partial y}\left ( 2.0 \cdot \gamma_{ab} \cdot a_{c} \cdot \left ( -\phi _{b} \right ) \cdot \left | q_{ab} \right |^{2} \cdot \left( \frac{\partial a_{c}}{\partial q_{ab}}\right)_{y} \right ) + \frac{\partial }{\partial z}\left ( 2.0 \cdot \gamma_{ab} \cdot a_{c} \cdot \left ( -\phi _{b} \right ) \cdot \left | q_{ab} \right |^{2} \cdot \left( \frac{\partial a_{c}}{\partial q_{ab}}\right)_{z} \right ) ---------------
								divdAdgradphi2 += 2*gamm(a,b)*(ac(iph)*q2(iph)*(-(phiOld(i+1,j,k,b)+phiOld(i,j,k,b))*0.5)*r_dadq(iph,X)-ac(imh)*q2(imh)*(-(phiOld(i,j,k,b)+phiOld(i-1,j,k,b))*0.5)*r_dadq(imh,X))/delta[X] + 2*gamm(a,b)*(ac(jph)*q2(jph)*(-(phiOld(i,j+1,k,b)+phiOld(i,j,k,b))*0.5)*r_dadq(jph,Y)-ac(jmh)*q2(jmh)*(-(phiOld(i,j,k,b)+phiOld(i,j-1,k,b))*0.5)*r_dadq(jmh,Y))/delta[Y]+2*gamm(a,b)*(ac(kph)*q2(kph)*(-(phiOld(i,j,k,b)+phiOld(i,j,k+1,b))*0.5)*r_dadq(kph,Z)-ac(kmh)*q2(kmh)*(-(phiOld(i,j,k,b)+phiOld(i,j,k-1,b))*0.5)*r_dadq(kmh,Z))/delta[Z];

							}
					}

				//Final anisotropy term for 3-D --------------------------------------------------
				term1_val(i,j,k,a) = (-dAdphi1-dAdphi2+divdAdgradphi1+divdAdgradphi2);

				dAdphi1=0.0;
				dAdphi2=0.0;
				divdAdgradphi1=0.0;
				divdAdgradphi2=0.0;
			}
	    });
	}
}




#endif
