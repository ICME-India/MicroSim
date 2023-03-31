#ifndef _ANISOTROPY_H_
#define _ANISOTROPY_H_

using namespace amrex;



// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void compute_ac(Vector<amrex::Real>& ac, Vector<amrex::Real>& acdashx, Vector<amrex::Real>& acdashy, Vector<amrex::Real>& nx, Vector<amrex::Real>& ny, Real dabb)
// {	
// 	Vector<Real> modphisq(4,0);

// 	for(int i=0; i<ac.size(); i++){

// 	modphisq[i] = pow(nx[i],2)+pow(ny[i],2);

// 	if(modphisq[i]>1e-15){
	
// 		ac[i] = (1-3*dabb) + 4*dabb*(nx[i]*nx[i]*nx[i]*nx[i]+ny[i]*ny[i]*ny[i]*ny[i])/(modphisq[i]*modphisq[i]);

// 		acdashx[i] = 16*dabb*((nx[i]*nx[i]*nx[i])/(modphisq[i])-nx[i]*(nx[i]*nx[i]*nx[i]*nx[i]+ny[i]*ny[i]*ny[i]*ny[i])/(modphisq[i]*modphisq[i]));
	
// 		acdashy[i] = 16*dabb*((ny[i]*ny[i]*ny[i])/(modphisq[i])-ny[i]*(nx[i]*nx[i]*nx[i]*nx[i]+ny[i]*ny[i]*ny[i]*ny[i])/(modphisq[i]*modphisq[i]));

// 		//Print()<<"ac_"<<i<<" :"<<ac[i];
	
// 	}
	
// 	else {
// 		ac[i] = 1.0;

// 		acdashx[i] = 0.0;

// 		acdashy[i] = 0.0;
// 	}

	
// 	}

// }


void computeterm1(MultiFab& term1, MultiFab& phi_old, MultiFab& derivx, MultiFab& derivy, Geometry const& geom){
	
	BL_PROFILE("computeterm1()");	

    for (MFIter mfi(phi_old); mfi.isValid(); ++mfi)
	{
		const Box& bbx = mfi.validbox();
		auto const& phiOld = phi_old.const_array(mfi);
		auto const& term1_val = term1.array(mfi);
		auto const& derivxx = derivx.const_array(mfi);
        auto const& derivyy = derivy.const_array(mfi);
		Real GAM = gammaa;
		Real dabb = dab;
		Real rot_angle = rotmat.back();
		

		amrex::ParallelFor( bbx, 
		[=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
		{
		
		     GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

			Array1D <Real,0,3> nx;
			Array1D <Real,0,3> ny;
			Array1D <Real,0,3> ac;
			Array1D <Real,0,3> acdashx;
			Array1D <Real,0,3> acdashy;
			Array1D <Real,0,3> modphisq;
			Array1D <Real,0,3> acdash_rot;
			
			
			// Vector <amrex:: Real> nx(4,0);
			// Vector <amrex:: Real> ny(4,0);
			// Vector <amrex:: Real> ac(4,0);
			// Vector <amrex:: Real> acdashx(4,0);
			// Vector <amrex:: Real> acdashy(4,0);
			// Vector <amrex::Real> modphisq(4,0);
			// Vector <amrex::Real> acdash_rot(4,0);

			Real pi = acos(-1);
			//Rotation by theta	
	        Real thetaf = rot_angle*pi/180.0;

			nx(0) = (derivxx(i,j,k,0)*cos(thetaf)-derivyy(i,j,k,0)*sin(thetaf));
			nx(1) = (derivxx(i,j,k,1)*cos(thetaf)-derivyy(i,j,k,1)*sin(thetaf));
			nx(2) = (derivxx(i,j,k,2)*cos(thetaf)-derivyy(i,j,k,2)*sin(thetaf));
			nx(3) = (derivxx(i,j,k,3)*cos(thetaf)-derivyy(i,j,k,3)*sin(thetaf));
	
	
			ny(0) = (derivxx(i,j,k,0)*sin(thetaf)+derivyy(i,j,k,0)*cos(thetaf));
			ny(1) = (derivxx(i,j,k,1)*sin(thetaf)+derivyy(i,j,k,1)*cos(thetaf));
			ny(2) = (derivxx(i,j,k,2)*sin(thetaf)+derivyy(i,j,k,2)*cos(thetaf));
			ny(3) = (derivxx(i,j,k,3)*sin(thetaf)+derivyy(i,j,k,3)*cos(thetaf));

			//compute_ac(ac, acdashx, acdashy,  nx, ny,dabb);

			for(int p=0; p<ac.size(); p++){

				// nx[p] = (derivxx(i,j,k,p)*cos(thetaf)-derivyy(i,j,k,p)*sin(thetaf));
				// ny[p] = (derivxx(i,j,k,p)*sin(thetaf)+derivyy(i,j,k,p)*cos(thetaf));

				modphisq(p) = nx(p)*nx(p)+ny(p)*ny(p);

				if(modphisq(p)>1e-15){
	
					ac(p) = (1-3*dabb) + 4*dabb*(nx(p)*nx(p)*nx(p)*nx(p)+ny(p)*ny(p)*ny(p)*ny(p))/(modphisq(p)*modphisq(p));

					acdashx(p) = 16*dabb*((nx(p)*nx(p)*nx(p))/(modphisq(p))-nx(p)*(nx(p)*nx(p)*nx(p)*nx(p)+ny(p)*ny(p)*ny(p)*ny(p))/(modphisq(p)*modphisq(p)));
	
					acdashy(p) = 16*dabb*((ny(p)*ny(p)*ny(p))/(modphisq(p))-ny(p)*(nx(p)*nx(p)*nx(p)*nx(p)+ny(p)*ny(p)*ny(p)*ny(p))/(modphisq(p)*modphisq(p)));
	
				}
	
				else {
					ac(p) = 1.0;

					acdashx(p) = 0.0;

					acdashy(p) = 0.0;
				}

				//acdash_rot[p] = acdashx[p]*cos(thetaf)+acdashy[p]*sin(thetaf);

			}


			//Rotated derivative
			

			acdash_rot(0) = acdashx(0)*cos(thetaf)+acdashy(0)*sin(thetaf);
			acdash_rot(1) = acdashx(1)*cos(thetaf)+acdashy(1)*sin(thetaf);
			acdash_rot(2) = -acdashx(2)*sin(thetaf)+acdashy(2)*cos(thetaf);
			acdash_rot(3) = -acdashx(3)*sin(thetaf)+acdashy(3)*cos(thetaf);
			

			Real ani_term1 = 2.0*GAM*((ac(0)*ac(0)*derivxx(i,j,k,0) - ac(1)*ac(1)*derivxx(i,j,k,1))/dx[0] + (ac(2)*ac(2)*derivyy(i,j,k,2) - ac(3)*ac(3)*derivyy(i,j,k,3))/dx[1]);

			Real ani_term2 = 2.0*GAM*((ac(0)*acdash_rot(0)-ac(1)*acdash_rot(1))/(dx[0])+(ac(2)*acdash_rot(2)-ac(3)*acdash_rot(3))/(dx[1]));
			

			term1_val(i,j,k) = ani_term1 + ani_term2;

			// nx.clear();
			// ny.clear();
			// ac.clear();
			// acdashx.clear();
			// acdashy.clear();
			// acdash_rot.clear();
			// modphisq.clear();
		
		});
	}
	

}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////









#endif
