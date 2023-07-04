#ifndef _INITIALIZE_H_
#define _INITIALIZE_H_

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
			Real cyl_rad = cylinder[p][5];

		amrex::ParallelFor( wbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
		{	

			#if (AMREX_SPACEDIM <= 2)
			if(((i-cyl_X_cent)*(i-cyl_X_cent) + (j-cyl_Y_cent)*(j-cyl_Y_cent)) < cyl_rad*cyl_rad)
			{
				phiNew(i,j,k,cyl_comp) = 1.0;
			}
			#endif
		

			// #if (AMREX_SPACEDIM > 2)
			// if(((i-cylinder[p][1])*(i-cylinder[p][1]) + (j-cylinder[p][2])*(j-cylinder[p][2])) < cylinder[p][5]*cylinder[p][5] && cylinder[p][3]<=k<=cylinder[p][4])
			// {
			// 	phiNew(i,j,k,cylinder[p][0]) = 1.0;
			// }
			// #endif 
		});
		}
	}
}	

void init_phi_sph(MultiFab& phi_new)
{
	// for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
	// {
	// 	const Box& wbx = mfi.validbox();
	// 	auto const& phiNew = phi_new.array(mfi);
		
	// 	for(int p=0; p<sphere.size(); p++){
	// 	amrex::ParallelFor( wbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
	// 	{	
	// 		#if (AMREX_SPACEDIM>2)
	// 		if(((i-sphere[p][1])*(i-sphere[p][1]) + (j-sphere[p][2])*(j-sphere[p][2])+(k-sphere[p][3])*(k-sphere[p][3])) < sphere[p][4]*sphere[p][4])
	// 		{
	// 			phiNew(i,j,k,sphere[p][0]) = 1.0;
	// 		}
	// 		#endif
			
	// 	});
	// }
	// }
}

void init_phi_cube(MultiFab& phi_new)
{
// 	for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
// 	{
// 		const Box& wbx = mfi.validbox();
// 		auto const& phiNew = phi_new.array(mfi);
		
// 		for(int p=0; p<cube.size(); p++){
// 		amrex::ParallelFor( wbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
// 		{	
// 			#if (AMREX_SPACEDIM>2)
// 			if(cube[p][1]<=i<=cube[p][4] && cube[p][2]<=j<=cube[p][5] && cube[p][3]<=k<=cube[p][6])
// 			{
// 				phiNew(i,j,k,cube[p][0]) = 1.0;
// 			}
// 			#endif
			
// 		});
// 	}
// }
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
		Real A_liq = A[nump-1][0][0];
		Real cequil = conceq[nump-1][0][0];

		amrex::ParallelFor( wbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
		{
			muNew(i,j,k) = 2*A_liq*cequil;
			//muNew(i,j,k) = 2*A_liq*ceq_liq;
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



#endif
