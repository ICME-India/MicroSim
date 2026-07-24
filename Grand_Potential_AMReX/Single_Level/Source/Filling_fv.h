#ifndef _FILLING_FV_H_
#define _FILLING_FV_H_

#include <gsl/gsl_rng.h>

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

void init_phi_cyl_rand(MultiFab& phi_new)
{	
	
	srand(time(0));
    int count{0};
	
    Vector<Vector<Vector<long>>> center(cylrand.size(),Vector<Vector<long>>(2,Vector<long>(0,0)));
	Vector<Vector<Real>> radi(cylrand.size(),Vector<Real>(0,0.0));

	for(int p =0; p<cylrand.size();p++){
			int cyl_comp = int(cylrand[p][0]);
			Real cyl_ppt_rad = cylrand[p][1];
			Real cyl_vol_frac = cylrand[p][2];
			Real cyl_shield = cylrand[p][3];
			Real cyl_spread = cylrand[p][4];

			Real mdev = cyl_spread*cyl_ppt_rad;
			Real volume_domain = ncellx*ncelly*ncellz;
			Real volume_per_particle = M_PI*cyl_ppt_rad*cyl_ppt_rad;


			int num_particles = ceil(volume_domain*cyl_vol_frac/volume_per_particle);

			Print()<<"num particle:"<<num_particles<<"\n";

			int part_id{0};

            center[cyl_comp][0].push_back(long(ncellx*((double)rand())/RAND_MAX));
			center[cyl_comp][1].push_back(long(ncelly*((double)rand())/RAND_MAX));
			radi[cyl_comp].push_back(cyl_ppt_rad + (1.0*(((double)rand())/RAND_MAX)-0.5)*mdev);

			while (part_id<num_particles-1){
                //Print()<<"part:"<<part_id<<"\n";
				long centx = long(ncellx*((double)rand())/RAND_MAX);
				long centy = long(ncelly*((double)rand())/RAND_MAX);
				Real rad = cyl_ppt_rad + (1.0*(((double)rand())/RAND_MAX)-0.5)*mdev;

               // Print()<<"center x:"<<centx<<", center y:"<<centy<<"rad:"<<rad<<"\n";
               // Print()<<center.size()<<","<<center[cyl_comp].size()<<","<<center[cyl_comp][0].size()<<"\n";

                int overlap=0;
                
                if(count!=0){
			
                for(int t=0; t<center[cyl_comp][0].size(); t++){
					if(((center[cyl_comp][0][t]-centx)*(center[cyl_comp][0][t]-centx) + 
					   (center[cyl_comp][1][t]-centy)*(center[cyl_comp][1][t]-centy)) <= (((cyl_shield+1)*cyl_ppt_rad)*((cyl_shield+1)*cyl_ppt_rad)) ){
							
                            overlap=1;
                            break;
					}

                    // if(cylrand.size()>0 && cyl_comp!=0){
                    //     for(int m=0; m<cylrand.size();m++){
                    //         if(m!=cyl_comp){
                    //         for(int t=0; t<center[m][0].size(); t++){
					//             if(((center[m][0][t]-centx)*(center[m][0][t]-centx) + 
					//                  (center[m][1][t]-centy)*(center[m][1][t]-centy)) <= (((cyl_shield+1)*cyl_ppt_rad)*((cyl_shield+1)*cyl_ppt_rad)) ){
							
                    //                     overlap=1;
                    //                     break;
					//             }
                    //         }
                    //         }
                    //    	}
					// }
                
                }
                }

                if(overlap==0){
                    count++;
                    part_id++;
                    center[cyl_comp][0].push_back(centx);
                    center[cyl_comp][1].push_back(centy);
					radi[cyl_comp].push_back(rad);
                }
                
			}

            
    }

	for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
	{
		const Box& wbx = mfi.validbox();
		Array4<Real> const& phiNew = phi_new.array(mfi);

		Array3D<long,0,0,0,1,0,30> centr{};		//fix this
		for(int a=0; a<nump-1; a++){
			for(int l=0; l<2; l++){
				for(int m=0; m<30; m++){
					centr(a,l,m) = center[a][l][m];
					//Print()<<"Cent("<<a<<","<<l<<","<<m<<"):"<<center[a][l][m];
				}
				Print()<<"\n";
			}
		}

		
		amrex::ParallelFor(wbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
		{	

			for(int m=0; m<cylrand.size();m++){						//fix this
                for(int t=0; t<center[m][0].size(); t++){
					if(((i-center[m][0][t])*(i-center[m][0][t]) + (j-center[m][1][t])*(j-center[m][1][t])) < radi[m][t]*radi[m][t])
				{
					phiNew(i,j,k,m) = 1.0;
				}
			}
			}
		});
		
	}

	Print()<<"Random cylinder filled\n";
}

void init_phi_sph_rand(MultiFab& phi_new)
{	
	
	srand(time(0));
    int count{0};
	
    Vector<Vector<Vector<Real>>> center(sphrand.size(),Vector<Vector<Real>>(3,Vector<Real>(0,0.0)));
	Vector<Vector<Real>> radi(sphrand.size(),Vector<Real>(0,0.0));

	for(int p =0; p<sphrand.size();p++){
			Real sph_comp = sphrand[p][0];
			Real sph_ppt_rad = sphrand[p][1];
			Real sph_vol_frac = sphrand[p][2];
			Real sph_shield = sphrand[p][3];
			Real sph_spread = sphrand[p][4];

			Real mdev = sph_spread*sph_ppt_rad;
			Real volume_domain = ncellx*ncelly*ncellz;
			Real volume_per_particle = (4.0/3.0)*M_PI*sph_ppt_rad*sph_ppt_rad*sph_ppt_rad;


			int num_particles = ceil(volume_domain*sph_vol_frac/volume_per_particle);

			Print()<<"num particle:"<<num_particles<<"\n";

			int part_id{0};

            center[sph_comp][0].push_back(ncellx*((double)rand())/RAND_MAX);
			center[sph_comp][1].push_back(ncelly*((double)rand())/RAND_MAX);
			center[sph_comp][2].push_back(ncellz*((double)rand())/RAND_MAX);
			radi[sph_comp].push_back(sph_ppt_rad + (1.0*(((double)rand())/RAND_MAX)-0.5)*mdev);

			while (part_id<num_particles-1){
                //Print()<<"part:"<<part_id<<"\n";
				Real centx = ncellx*((double)rand())/RAND_MAX;
				Real centy = ncelly*((double)rand())/RAND_MAX;
				Real centz = ncellz*((double)rand())/RAND_MAX;
				Real rad = sph_ppt_rad + (1.0*(((double)rand())/RAND_MAX)-0.5)*mdev;

               // Print()<<"center x:"<<centx<<", center y:"<<centy<<"rad:"<<rad<<"\n";
               // Print()<<center.size()<<","<<center[cyl_comp].size()<<","<<center[cyl_comp][0].size()<<"\n";

                int overlap=0;
                
                if(count!=0){
			
                for(int t=0; t<center[sph_comp][0].size(); t++){
					if(((center[sph_comp][0][t]-centx)*(center[sph_comp][0][t]-centx) + 
					   (center[sph_comp][1][t]-centy)*(center[sph_comp][1][t]-centy)+
					   (center[sph_comp][2][t]-centz)*(center[sph_comp][2][t]-centz)) <= (((sph_shield+1)*sph_ppt_rad)*((sph_shield+1)*sph_ppt_rad)) ){
							
                            overlap=1;
                            break;
					   }

                       if(sphrand.size()>0 && sph_comp!=0){
                        for(int m=0; m<sphrand.size();m++){
                            if(m!=sph_comp){
                            for(int t=0; t<center[m][0].size(); t++){
					            if(((center[m][0][t]-centx)*(center[m][0][t]-centx) + 
					                 (center[m][1][t]-centy)*(center[m][1][t]-centy) + 
									  (center[m][2][t]-centz)*(center[m][2][t]-centz)) <= (((sph_shield+1)*sph_ppt_rad)*((sph_shield+1)*sph_ppt_rad)) ){
							
                                        overlap=1;
                                        break;
					            }
                            }
                            }
                       }
				}
                
                }
                }

                if(overlap==0){
                    count++;
                    part_id++;
                    center[sph_comp][0].push_back(centx);
                    center[sph_comp][1].push_back(centy);
					center[sph_comp][2].push_back(centz);
					radi[sph_comp].push_back(rad);
                }
                
			}

            
    }

	for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
	{
		const Box& wbx = mfi.validbox();
		Array4<Real> const& phiNew = phi_new.array(mfi);

		Array3D<Real,0,0,0,1,0,30> centr{};
		for(int a=0; a<center.size(); a++){
			for(int l=0; l<center[0].size(); l++){
				for(int m=0; m<center[0][0].size(); m++){
					centr(a,l,m) = center[a][l][m];
				}
			}
		}

		
		amrex::ParallelFor(wbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
		{	

			for(int m=0; m<sphrand.size();m++){
                for(int t=0; t<center[m][0].size(); t++){
					if(((i-center[m][0][t])*(i-center[m][0][t]) + (j-center[m][1][t])*(j-center[m][1][t]) + (k-center[m][2][t])*(k-center[m][2][t])) < radi[m][t]*radi[m][t])
				{
					phiNew(i,j,k,m) = 1.0;
				}
			}
			}
		});
		
	}
}

void init_phi_cube_rand_variants(MultiFab& phi_new){

	//srand(time(0));
	gsl_rng *rng;

	gsl_rng_env_setup();
    rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, time(0));

	int var = cuberand[0];
	Real sx = cuberand[1];
	Real sy = cuberand[2];
	Real sz = cuberand[3];
	Real sfrac = cuberand[4];
	Real vf = cuberand[5];
	Real shield = cuberand[6];

	Real volume_domain = ncellx*ncelly*ncellz;
	Real volume_per_particle = sx*sy*sz;

	Real r{};

	Vector<long> data{};
	Vector<int> variant{};
	Vector<Vector<long>> cube_cent{};
	Vector<long> snow{};
	Vector<Vector<long>> cube_size{};

	// Print()<<"TV: "<<volume_domain<<"\n";
	// Print()<<"VPP: "<<volume_per_particle<<"\n";
	// Print()<<"VF: "<<vf<<"\n";


	int num_particles = ceil(volume_domain*vf/volume_per_particle);
	int count =0;
	int part_id{0};
	Print()<<"Number of particles = "<<num_particles<<"\n";

	while(part_id<num_particles){
		r = ((double)rand())/RAND_MAX;
		//r = gsl_rng_uniform(rng);

		snow.push_back(sx + (2*r-1)*sfrac*sx);
		snow.push_back(sy + (2*r-1)*sfrac*sy);
		snow.push_back(sz + (2*r-1)*sfrac*sz);

		// long xnow = ncellx*(gsl_rng_uniform(rng));
		// long ynow = ncelly*(gsl_rng_uniform(rng));
		// long znow = ncellz*(gsl_rng_uniform(rng));

		long xnow = ncellx*((double)rand())/RAND_MAX;
		long ynow = ncelly*((double)rand())/RAND_MAX;
		long znow = ncellz*((double)rand())/RAND_MAX;


		Print()<<"-------------------------------\n";
		Print()<<"Particle_id = "<<part_id<<"\n";
		// Print()<<"snow = "<<snow[0]<<" "<<snow[1]<<" "<<snow[2]<<"\n";
		// Print()<<"now = "<<xnow<<" "<<ynow<<" "<<znow<<"\n";


		data.push_back(xnow);
		data.push_back(ynow);
		data.push_back(znow);
		

		int ovlp=0;

		if(count!=0){
			if(dim == 2){
				
				for(int m=0; m<cube_cent.size(); m++){	

					if(std::abs(data[0]-(cube_cent[m][0]))<(0.5*snow[0]+shield+0.5*cube_size[m][0]) && std::abs(data[1]-(cube_cent[m][1]))<(0.5*snow[1]+shield+0.5*cube_size[m][1])){

							ovlp=1;
							
							break;
						}

				}
			}

			#if(AMREX_SPACEDIM == 3)
				for(int m=0; m<cube_cent.size(); m++){	

					if(std::abs(data[0]-(cube_cent[m][0]))<(0.5*snow[0]+shield+0.5*cube_size[m][0]) && std::abs(data[1]-(cube_cent[m][1]))<(0.5*snow[1]+shield+0.5*cube_size[m][1]) && std::abs(data[2]-(cube_cent[m][2]))<(0.5*snow[2]+shield+0.5*cube_size[m][2]) ){

							ovlp=1;
							
							break;
						}

				}
			#endif
			
		}

		//Print()<<"overlap: "<<ovlp<<"\n";

		if(ovlp==0){
			count++;
			part_id++;
			cube_cent.push_back(data);
			cube_size.push_back(snow);
			variant.push_back(int(r*var));
			data.clear();
			snow.clear();
		}

		else{
			data.clear();
			snow.clear();
		}
	}

	Vector<Vector<Real>> cube_data{};

	Vector<Real> val{};

	int index=0;

	while(index<cube_cent.size()){

		val.push_back(cube_cent[index][0]-cube_size[index][0]/2);
		val.push_back(cube_cent[index][0]+cube_size[index][0]/2);
		val.push_back(cube_cent[index][1]-cube_size[index][1]/2);
		val.push_back(cube_cent[index][1]+cube_size[index][1]/2);
	#if(AMREX_SPACEDIM>2)
		val.push_back(cube_cent[index][2]-cube_size[index][2]/2);
		val.push_back(cube_cent[index][2]+cube_size[index][2]/2);
	#endif

	cube_data.push_back(val);

	val.clear();
	index++;

	}

	//Print()<<"All 20 done\n";
	// Print()<<"ht cube_data: "<<cube_cent.size()<<"\n";
	// Print()<<"wd of cube_data: "<<cube_cent[0].size()<<"\n";
	// Print()<<"Variant: "<<variant.size()<<"\n";
	// Print()<<"variant = "<<variant[0]<<" "<<variant[1]<<" "<<variant[2]<<" "<<variant[3]<<"\n";

	for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
	{
		const Box& wbx = mfi.validbox();
		auto const& phiNew = phi_new.array(mfi);

		//Print()<<"here1\n";

		Array2D<Real,0,50,0,5,Order::C> cubes{};			//fix this
		for(int l=0; l<cube_data.size(); l++){
			for(int m=0; m<cube_data[0].size(); m++){
				cubes(l,m) = cube_data[l][m];
			}
		}
		//Print()<<"here2\n";
		

		Array1D<int, 0, 50> phas{};					//fix this
		for(int l=0; l<variant.size(); l++){
			phas(l) = variant[l];
		}
			
		//Print()<<"here3\n";
		amrex::ParallelFor( wbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
		{	
			for(int m=0; m< cube_cent.size();m++){
					
					#if (AMREX_SPACEDIM==2)
					if(i>cubes(m,0) && 
					   i<cubes(m,1) && 
					   j>cubes(m,2) && 
					   j<cubes(m,3))
					{
						phiNew(i,j,k,phas(m)) = 1.0;
					}
					#endif

					#if (AMREX_SPACEDIM==3)
					if(i>cubes(m,0) && 
					   i<cubes(m,1) && 
					   j>cubes(m,2) && 
					   j<cubes(m,3) &&
					   j<cubes(m,4) &&
					   j<cubes(m,5))
					{
						phiNew(i,j,k,phas(m)) = 1.0;
					}
					#endif

			}
			
		});
		}
	


}

void init_phi_cube_pattern(MultiFab& phi_new){
	
	long variant = cubepat[0];
	long sx = cubepat[1];
	long sy = cubepat[2];
	long sz = cubepat[3];
	Real sfrac = cubepat[4];
	long gap = cubepat[5];
	Real gfrac = cubepat[6];

	Vector<long> data{};
	Vector<long> phas{};
	Vector<Vector<long>> cube_data{};
	Real outof = 0;
	if(cubepat.size()>7){
		outof = cubepat[7];
	}

	gsl_rng *rng;
	gsl_rng_env_setup();
    rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, time(0));
	Real r{};

	ldiv_t resx, resy, resz;

	resx = ldiv(ncellx, sx+gap);
    resy = ldiv(ncelly, sy+gap);
    resz = ldiv(ncellz, sz+gap);

	if ( resx.quot == 0 )
        resx.quot = 1;
    if ( resy.quot == 0 )
        resy.quot = 1;
    if ( resz.quot == 0 )
        resz.quot = 1;
    long nparticles = resx.quot * resy.quot * resz.quot;
    
	Print()<<"number particles: "<<nparticles<<"\n";
	
	for (int m=0; m<resx.quot; m++ )
    {
        if ( ncellx > 1 )
        {
            r = gsl_rng_uniform(rng);
            long sgn = 2*lround(r) - 1;
            data.push_back(m*(gap+sx) + gap + sgn*gfrac*gap - (2*r-1)*sfrac*sx);
            data.push_back((m+1) * (gap+sx) + sgn*gfrac*gap + (2*r-1)*sfrac*sx);
        }
        for (int n=0; n<resy.quot; n++ )
        {
            if ( ncelly > 1 )
            {
                r = gsl_rng_uniform(rng);
                long sgn = 2*lround(r) - 1;
                data.push_back(n*(gap+sy) + gap + sgn*gfrac*gap - (2*r-1)*sfrac*sy);
                data.push_back((n+1) * (gap+sy) + sgn*gfrac*gap + (2*r-1)*sfrac*sy);
            }
            for (int l=0; l<resz.quot; l++ )
            {
                if ( outof > 0 )
                {
                    if ( !gsl_rng_uniform_int(rng, outof) )
                        continue;
                }
                if ( ncellz > 1 )
                {
                    r = gsl_rng_uniform(rng);
                    long sgn = 2*lround(r) - 1;
                    data.push_back(l*(gap+sz) + gap + sgn*gfrac*gap - (2*r-1)*sfrac*sz);
                    data.push_back((l+1) * (gap+sz) + sgn*gfrac*gap + (2*r-1)*sfrac*sz);
                }
                
				cube_data.push_back(data);

				#if(AMREX_SPACEDIM>2)
				data.pop_back();
				data.pop_back();
				#endif

				// Print()<<"cube_data size: "<<cube_data.size()<<"\n";
				// for(int a=0; a<cube_data[0].size();a++){
				// 	Print()<<cube_data[cube_data.size()-1][a]<<" ";
				// }
				// Print()<<"\n";
				
                
				if ( variant >= nump )
                    variant = nump - 1;
                r = gsl_rng_uniform(rng);
                phas.push_back(r * variant);

				// Print()<<"phase size: "<<phas.size()<<"\n";
				// Print()<<"phase: "<<phas[n]<<"\n";
                
            }

			data.pop_back();
			data.pop_back();
        }
		data.clear();
    }
    // Free GSL RNG.
    gsl_rng_free(rng);

	for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
	{
		const Box& wbx = mfi.validbox();
		auto const& phiNew = phi_new.array(mfi);

		//Print()<<"here1\n";

		Array2D<Real,0,50,0,5,Order::C> cubes{};			//fix this
		for(int l=0; l<cube_data.size(); l++){
			for(int m=0; m<cube_data[0].size(); m++){
				cubes(l,m) = cube_data[l][m];
			}
		}
		//Print()<<"here2\n";
		

		Array1D<int, 0, 50> ph{};					//fix this
		for(int l=0; l<phas.size(); l++){
			ph(l) = phas[l];
		}
			
		//Print()<<"here3\n";
		amrex::ParallelFor( wbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
		{	
			for(int m=0; m<cubes.size();m++){
					
					#if (AMREX_SPACEDIM==2)
					if(i>cubes(m,0) && 
					   i<cubes(m,1) && 
					   j>cubes(m,2) && 
					   j<cubes(m,3))
					{
						phiNew(i,j,k,ph(m)) = 1.0;
					}
					#endif

					#if (AMREX_SPACEDIM==3)
					if(i>cubes(m,0) && 
					   i<cubes(m,1) && 
					   j>cubes(m,2) && 
					   j<cubes(m,3) &&
					   j<cubes(m,4) &&
					   j<cubes(m,5))
					{
						phiNew(i,j,k,ph(m)) = 1.0;
					}
					#endif

			}
			
		});
		}

}

void init_phi_voronoi_2D(MultiFab& phi_new){

		int Vx_lo = voronoi2D[0];
		int Vx_hi = voronoi2D[1];
		int Vy_lo = voronoi2D[2];
		int Vy_hi = voronoi2D[3];
		int Vnp = voronoi2D[4];
		int Vsz = voronoi2D[5];
		int gidy{0};
    	int gidy1{0};
		long rand_x, rand_y;
    	int Phase_filled=0;
		Gpu::AsyncVector<Real> flag(ncellx*ncelly,0.0);
		Gpu::AsyncVector<Real> n(Vnp,0.0);
		Gpu::AsyncVector<Real> m(Vnp,0.0);
		Gpu::AsyncVector<Real> phase(Vnp,0.0);

		long limit_x = Vx_hi-Vx_lo;
		long limit_y = Vy_hi-Vy_lo;

		for(int i=0; i<Vnp; i++){
        while(Phase_filled!=1){
            rand_x = drand48()*limit_x;
            rand_y = drand48()*limit_y;
            gidy = rand_x*ncelly+rand_y;

            if(flag[gidy]!=1){
                n[i] = rand_x;
                m[i] = rand_y;

            for(int j=0; j<ncellx; j++){
                for(int k=0; k<ncelly; k++){
                    if ((n[i]-j)*(n[i]-j) +(m[i]-k)*(m[i]-k) <= Vsz*Vsz) {
                    gidy1 = j*ncelly + k;
                    flag[gidy1] = 1;
                }
            }
            }
            Phase_filled=1;
            }
        }
        Phase_filled=0;
    	}

		for(int i=0; i<Vnp; i++){
       		phase[i] = lrand48()%(nump-1); 
    	}

	for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
	{
		const Box& wbx = mfi.validbox();
		Array4<Real> const& phiNew = phi_new.array(mfi);

		amrex::ParallelFor( wbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
		{
			Real minimum{0.0};
    		int location{0};
			Gpu::AsyncVector<Real> l(Vnp,0.0);

			for(int j=0; j<ncellx; j++){
				for(int k=0; k<ncelly; k++){
					//gidy = j*ncelly + k;
					for(int r=0; r<Vnp; r++){
						l[r] = (n[r]-j)*(n[r]-j)+(m[r]-k)*(m[r]-k);
					}
					minimum = l[0];
					location=0;

					for(int t=0;t<Vnp;t++) {
						if(l[t]<minimum) {
							minimum  = l[t];
							location = phase[t];
						}
					}

				// Print()<<"location:"<<location<<"\n";

					for(int s=0;s<nump-1;s++) {
						if(location == s) {
							phiNew(i,j,k,s)= 1.0;
						}
						if(location != s) {
         					phiNew(i,j,k,s) = 0.0;
        				}
					}
				}
			}

		});
	}

		


}

void fill_voronoi_3D(){
	
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

					if(sum>1.0){
						sum=1.0;
					}
				
				}

				phiNew(i,j,k,numphase-1) = 1.0 - sum;
		});
		}
}

void Init_mu (MultiFab& mu_new, MultiFab& phi_new)
{ 	
	//Print()<<"fin1\n";

	for (MFIter mfi(mu_new); mfi.isValid(); ++mfi)
	{
		const Box& wbx = mfi.validbox();
		Array4<Real> const& muNew = mu_new.array(mfi);
		Array4<Real> const& phiNew = phi_new.array(mfi);
		//Print()<<"fin2\n";

		Array3D<Real,0,phasecount,0,compcount-2,0,compcount-2> A_fill{};
		for(int l=0; l<nump; l++){
			for(int m=0; m<numcom-1;m++){
				for(int n=0; n<numcom-1; n++){
					A_fill(l,m,n) = A[l][m][n];
				// 	Print()<<"A_fill ["<<l<<" , "<<m<<" , "<<n<<"]: "<<A_fill(l,m,n)<<"\n";
				// 	Print()<<"A ["<<l<<" , "<<m<<" , "<<n<<"]: "<<A[l][m][n]<<"\n";
				}
			}
		}

		//Print()<<"fin3\n";
		Array1D<Real,0,compcount-2> ceq_liq{};
		for(int l=0; l<numcom-1; l++){
				ceq_liq(l) = conceq[nump-1][l];
				//Print()<<"ceq_liq["<<l<<"]: "<<ceq_liq(l)<<"\n";
		}

		Array1D<Real,0,compcount-2> ceq_sol{};
		for(int l=0; l<numcom-1; l++){
				ceq_sol(l) = conceq[0][l];
				//Print()<<"ceq_sol["<<l<<"]: "<<ceq_sol(l)<<"\n";
		}

		Array2D <Real,0,phasecount-1,0,compcount-2, Order::C> BB{};
		for(int a=0; a<nump; a++){
			for(int l=0; l<numcom-1; l++){
				BB(a,l) = B[a][l];
			}
		}

		//Print()<<"fin4\n";
		Real numcomp = numcom;
		Real numphase = nump;
		//Print()<<numcomp<<"\n";
		//Real A_liq = A[nump-1][0][0];
		//Real cequil = conceq[nump-1][0][0];

		amrex::ParallelFor( wbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
		{	
			//Print()<<"fin5\n";
			double sum =0.0;
			int val = -1;
			Gpu::DeviceVector <Real> flag(numphase,0.0);

			for(int a=0; a<nump-1; a++){
				if(phiNew(i,j,k,a)==1.0){
					val=a;
					break;
				}
			}

			if(val>=0){
				for(int l=0; l<numcomp-1; l++){
					for(int m=0; m<numcomp-1; m++){
				
						if(l==m){
							sum += 2.0*A_fill(val,l,m)*ceq_sol(m);
						}
						else{
							sum += A_fill(val,l,m)*ceq_sol(m);
						}
					}
					muNew(i,j,k,l) = sum + BB(val,l);
					sum=0.0;
				}
			}

			else{
				for(int l=0; l<numcomp-1; l++){
					for(int m=0; m<numcomp-1; m++){
				
						if(l==m){
							sum += 2.0*A_fill(nump-1,l,m)*ceq_liq(m);
						}
						else{
							sum += A_fill(nump-1,l,m)*ceq_liq(m);
						}
					}
				muNew(i,j,k,l) = sum;
				sum=0.0;
				}
			}
		});
	}
}	

void Init_comp(MultiFab& phi_new, MultiFab& comp_new){

		for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
	{	
		const Box& pbx = mfi.validbox();
		Array4<Real> const& phiNew = phi_new.array(mfi);
		Array4<Real> const& compNew = comp_new.array(mfi);
		int numphase = nump;
		int numcomp = numcom;

		Array2D<Real,0,phasecount-1,0,compcount-2,Order::C> co{};
		for(int a=0; a<nump; a++){
			for(int l=0; l<numcom-1; l++){
				co(a,l) = conc[a][l];
			}
		}

		Array2D<Real,0,phasecount-1,0,compcount-2,Order::C> coeq{};
		for(int a=0; a<nump; a++){
			for(int l=0; l<numcom-1; l++){
				coeq(a,l) = conceq[a][l];
			}
		} 
		
		amrex::ParallelFor( pbx, [=] AMREX_GPU_DEVICE( int i, int j ,int k)
		{		
				Real sum{0.0};
				int val{0};

				for(int a=0; a<numphase-1; a++){
					
					if(phiNew(i,j,k,a)==1.0){
						for(int l=0; l<numcomp-1; l++){
						compNew(i,j,k,l) = coeq(a,l);
					}
					val=1;
					break;
					}
				}

				if(val==0){
					for(int l=0; l<numcomp-1; l++){
						compNew(i,j,k,l) = coeq(nump-1,l);;
					}
				}

		});
		}
}

#endif
