#ifndef _FILLING_FV_H_
#define _FILLING_FV_H_

using namespace amrex;

void init_phi_cyl(MultiFab& phi_new)
{
	for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
	{
		const Box& wbx = mfi.validbox();
		Array4<Real> const& phiNew = phi_new.array(mfi);
		
		for(int p=0; p<cylinder.size(); p++){
			Real cyl_comp = cylinder[p][0];
			Real cyl_X_cent = cylinder[p][1];
			Real cyl_Y_cent = cylinder[p][2];
			Real cyl_Z_strt = cylinder[p][3];
			Real cyl_Z_end = cylinder[p][4];
			Real cyl_rad = cylinder[p][5];

		amrex::ParallelFor( wbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
		{	

			#if (AMREX_SPACEDIM <= 2)
			if(((i-cyl_X_cent)*(i-cyl_X_cent) + (j-cyl_Y_cent)*(j-cyl_Y_cent)) < cyl_rad*cyl_rad)
			{
				phiNew(i,j,k,cyl_comp) = 1.0;
			}
			#endif
		

			#if (AMREX_SPACEDIM > 2)
			if(((i-cyl_X_cent)*(i-cyl_X_cent) + (j-cyl_Y_cent)*(j-cyl_Y_cent)) < cyl_rad*cyl_rad && cyl_Z_strt<=k<=cyl_Z_end)
			{
				phiNew(i,j,k,cyl_comp) = 1.0;
			}
			#endif 
		});
		}
	}
}	

void init_phi_sph(MultiFab& phi_new)
{
	for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
	{
		const Box& wbx = mfi.validbox();
		auto const& phiNew = phi_new.array(mfi);
		
		for(int p=0; p<sphere.size(); p++){
			Real sph_phase = sphere[p][0];
			Real sph_X_cent = sphere[p][1];
			Real sph_Y_cent = sphere[p][2];
			Real sph_Z_cent = sphere[p][3];
			Real sph_rad = sphere[p][4];

		amrex::ParallelFor( wbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
		{	
			#if (AMREX_SPACEDIM>2)
			if(((i-sph_X_cent)*(i-sph_X_cent) + (j-sph_Y_cent)*(j-sph_Y_cent)+(k-sph_Z_cent)*(k-sph_Z_cent)) < sph_rad*sph_rad)
			{
				phiNew(i,j,k,sph_phase) = 1.0;
			}
			#endif
			
		});
	}
	}
}

void init_phi_cube(MultiFab& phi_new)
{
	for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
	{
		const Box& wbx = mfi.validbox();
		auto const& phiNew = phi_new.array(mfi);
		
		for(int p=0; p<cube.size(); p++){
			Real cube_phase = cube[p][0];
			Real cube_X_strt = cube[p][1];
			Real cube_Y_strt = cube[p][2];
			Real cube_Z_strt = cube[p][3];
			Real cube_X_end = cube[p][4];
			Real cube_Y_end = cube[p][5];
			Real cube_Z_end = cube[p][6];

		amrex::ParallelFor( wbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
		{	
			#if (AMREX_SPACEDIM>2)
			if(cube_X_strt<=i<=cube_X_end && cube_Y_strt<=j<=cube_Y_end && cube_Z_strt<=k<=cube_Z_end)
			{
				phiNew(i,j,k,cube_phase) = 1.0;
			}
			#endif
			
		});
	}
}
}

void init_phi_ellip(MultiFab& phi_new)
{
	for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
	{
		const Box& wbx = mfi.validbox();
		Array4<Real> const& phiNew = phi_new.array(mfi);
		Real pi = acos(-1);
		for(int p=0; p<ellipse.size(); p++){
			Real rot_angle = ellipse[p][6]*pi/180.0;
			Real a = ellipse[p][4];
			Real b = sqrt(a*a*(1-ellipse[p][5]*ellipse[p][5]));
			Real ellip_X_cent = ellipse[p][1];
			Real ellip_Y_cent = ellipse[p][2];
			Real ellip_phase = ellipse[p][0]; 
			
		amrex::ParallelFor( wbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
		{	
			if((((i-ellip_X_cent)*cos(rot_angle)+(j-ellip_Y_cent)*sin(rot_angle))*((i-ellip_X_cent)*cos(rot_angle)+(j-ellip_Y_cent)*sin(rot_angle))/(a*a) + (-(i-ellip_X_cent)*sin(rot_angle)+(j-ellip_Y_cent)*cos(rot_angle))*(-(i-ellip_X_cent)*sin(rot_angle)+(j-ellip_Y_cent)*cos(rot_angle))/(b*b)) <= 1)
			{
				phiNew(i,j,k,ellip_phase) = 1.0;
			}
		});
	}
	}
}

void init_mu (MultiFab& mu_new)
{ 	
	for (MFIter mfi(mu_new); mfi.isValid(); ++mfi)
	{
		const Box& wbx = mfi.validbox();
		Array4<Real> const& muNew = mu_new.array(mfi);

		Array2D<Real,0,compcount-2,0,compcount-2> A_liq{};
		for(int l=0; l<numcom-1; l++){
			for(int m=0; m<numcom-1;m++){
				A_liq(l,m) = A[nump-1][l][m];
			}
		}
		
		Array1D<Real,0,compcount-2> ceq_liq{};
		for(int l=0; l<numcom-1; l++){
				ceq_liq(l) = conceq[nump-1][l];
		}
		
		Real numcomp = numcom;

		amrex::ParallelFor( wbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
		{	
			double sum =0.0;
			for(int l=0; l<numcomp-1; l++){
				for(int m=0; m<numcomp-1; m++){
					
					if(l==m){
						sum += 2.0*A_liq(l,m)*ceq_liq(m);
					}
					else{
						sum += A_liq(l,m)*ceq_liq(m);
					}
				}
				muNew(i,j,k,l) = sum;
				sum=0.0;
			}
		});
	}
}	




void Init_liq(MultiFab& phi_new){

		for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
	{
		const Box& wbx = mfi.validbox();
		Array4<Real> const& phiNew = phi_new.array(mfi);
		Real numphase = nump;
		
		amrex::ParallelFor( wbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
		{		
				Real sum{0.0};

				for(int a=0; a<numphase-1; a++){
					
					sum = sum + phiNew(i,j,k,a);
				
				}

				phiNew(i,j,k,numphase-1) = 1.0 - sum;
		});
		}
}



// void Init_comp(MultiFab& phi_new, MultiFab& comp_new){

// 		for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
// 	{	
// 		const Box& pbx = mfi.validbox();
// 		Array4<Real> const& phiNew = phi_new.array(mfi);
// 		Array4<Real> const& compNew = comp_new.array(mfi);
// 		int numphase = nump;
// 		int numcomp = numcom;

// 		Array2D<Real,0,phasecount-1,0,compcount-2,Order::C> co{};
// 		for(int a=0; a<nump; a++){
// 			for(int l=0; l<numcom-1; l++){
// 				co(a,l) = conc[a][l];
// 			}
// 		}

// 		Array2D<Real,0,phasecount-1,0,compcount-2,Order::C> coeq{};
// 		for(int a=0; a<nump; a++){
// 			for(int l=0; l<numcom-1; l++){
// 				coeq(a,l) = conceq[a][l];
// 			}
// 		} 
		
// 		amrex::ParallelFor( pbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
// 		{		
// 				Real sum{0.0};
// 				int val{0};

// 				for(int a=0; a<numphase-1; a++){
					
// 					if(phiNew(i,j,k,a)==1.0){
// 						for(int l=0; l<numcomp-1; l++){
// 						compNew(i,j,k,l) = coeq(a,l);
// 					}
// 					val=1;
// 					break;
// 					}
// 				}

// 				if(val==0){
// 					for(int l=0; l<numcomp-1; l++){
// 						compNew(i,j,k,l) = co(nump-1,l);;
// 					}
// 				}

// 		});
// 		}
// }

#endif
