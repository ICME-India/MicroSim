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


// void computeterm1(MultiFab& term1, MultiFab& phi_old, MultiFab& derivx, MultiFab& derivy, Geometry const& geom){
	
// 	BL_PROFILE("computeterm1()");	

//     for (MFIter mfi(phi_old); mfi.isValid(); ++mfi)
// 	{
// 		const Box& bbx = mfi.validbox();
// 		auto const& phiOld = phi_old.const_array(mfi);
// 		auto const& term1_val = term1.array(mfi);
// 		auto const& derivxx = derivx.const_array(mfi);
//         auto const& derivyy = derivy.const_array(mfi);
// 		Real GAM = gammaa;
// 		Real dabb = dab;
// 		Real rot_angle = rotmat.back();
		

// 		amrex::ParallelFor( bbx, 
// 		[=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
// 		{
		
// 		     GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

// 			Array1D <Real,0,3> nx;
// 			Array1D <Real,0,3> ny;
// 			Array1D <Real,0,3> ac;
// 			Array1D <Real,0,3> acdashx;
// 			Array1D <Real,0,3> acdashy;
// 			Array1D <Real,0,3> modphisq;
// 			Array1D <Real,0,3> acdash_rot;
			
			
// 			// Vector <amrex:: Real> nx(4,0);
// 			// Vector <amrex:: Real> ny(4,0);
// 			// Vector <amrex:: Real> ac(4,0);
// 			// Vector <amrex:: Real> acdashx(4,0);
// 			// Vector <amrex:: Real> acdashy(4,0);
// 			// Vector <amrex::Real> modphisq(4,0);
// 			// Vector <amrex::Real> acdash_rot(4,0);

// 			Real pi = acos(-1);
// 			//Rotation by theta	
// 	        Real thetaf = rot_angle*pi/180.0;

// 			nx(0) = (derivxx(i,j,k,0)*cos(thetaf)-derivyy(i,j,k,0)*sin(thetaf));
// 			nx(1) = (derivxx(i,j,k,1)*cos(thetaf)-derivyy(i,j,k,1)*sin(thetaf));
// 			nx(2) = (derivxx(i,j,k,2)*cos(thetaf)-derivyy(i,j,k,2)*sin(thetaf));
// 			nx(3) = (derivxx(i,j,k,3)*cos(thetaf)-derivyy(i,j,k,3)*sin(thetaf));
	
	
// 			ny(0) = (derivxx(i,j,k,0)*sin(thetaf)+derivyy(i,j,k,0)*cos(thetaf));
// 			ny(1) = (derivxx(i,j,k,1)*sin(thetaf)+derivyy(i,j,k,1)*cos(thetaf));
// 			ny(2) = (derivxx(i,j,k,2)*sin(thetaf)+derivyy(i,j,k,2)*cos(thetaf));
// 			ny(3) = (derivxx(i,j,k,3)*sin(thetaf)+derivyy(i,j,k,3)*cos(thetaf));

// 			//compute_ac(ac, acdashx, acdashy,  nx, ny,dabb);

// 			for(int p=0; p<ac.size(); p++){

// 				// nx[p] = (derivxx(i,j,k,p)*cos(thetaf)-derivyy(i,j,k,p)*sin(thetaf));
// 				// ny[p] = (derivxx(i,j,k,p)*sin(thetaf)+derivyy(i,j,k,p)*cos(thetaf));

// 				modphisq(p) = nx(p)*nx(p)+ny(p)*ny(p);

// 				if(modphisq(p)>1e-15){
	
// 					ac(p) = (1-3*dabb) + 4*dabb*(nx(p)*nx(p)*nx(p)*nx(p)+ny(p)*ny(p)*ny(p)*ny(p))/(modphisq(p)*modphisq(p));

// 					acdashx(p) = 16*dabb*((nx(p)*nx(p)*nx(p))/(modphisq(p))-nx(p)*(nx(p)*nx(p)*nx(p)*nx(p)+ny(p)*ny(p)*ny(p)*ny(p))/(modphisq(p)*modphisq(p)));
	
// 					acdashy(p) = 16*dabb*((ny(p)*ny(p)*ny(p))/(modphisq(p))-ny(p)*(nx(p)*nx(p)*nx(p)*nx(p)+ny(p)*ny(p)*ny(p)*ny(p))/(modphisq(p)*modphisq(p)));
	
// 				}
	
// 				else {
// 					ac(p) = 1.0;

// 					acdashx(p) = 0.0;

// 					acdashy(p) = 0.0;
// 				}

// 				//acdash_rot[p] = acdashx[p]*cos(thetaf)+acdashy[p]*sin(thetaf);

// 			}


// 			//Rotated derivative
			

// 			acdash_rot(0) = acdashx(0)*cos(thetaf)+acdashy(0)*sin(thetaf);
// 			acdash_rot(1) = acdashx(1)*cos(thetaf)+acdashy(1)*sin(thetaf);
// 			acdash_rot(2) = -acdashx(2)*sin(thetaf)+acdashy(2)*cos(thetaf);
// 			acdash_rot(3) = -acdashx(3)*sin(thetaf)+acdashy(3)*cos(thetaf);
			

// 			Real ani_term1 = 2.0*GAM*((ac(0)*ac(0)*derivxx(i,j,k,0) - ac(1)*ac(1)*derivxx(i,j,k,1))/dx[0] + (ac(2)*ac(2)*derivyy(i,j,k,2) - ac(3)*ac(3)*derivyy(i,j,k,3))/dx[1]);

// 			Real ani_term2 = 2.0*GAM*((ac(0)*acdash_rot(0)-ac(1)*acdash_rot(1))/(dx[0])+(ac(2)*acdash_rot(2)-ac(3)*acdash_rot(3))/(dx[1]));
			

// 			term1_val(i,j,k) = ani_term1 + ani_term2;

// 			// nx.clear();
// 			// ny.clear();
// 			// ac.clear();
// 			// acdashx.clear();
// 			// acdashy.clear();
// 			// acdash_rot.clear();
// 			// modphisq.clear();
		
// 		});
// 	}
	

// }



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void computeterm1(MultiFab& term1, MultiFab& phi_old, MultiFab& derivx, MultiFab& derivy, Geometry const& geom){

	for(MFIter mfi(phi_old); mfi.isValid(); ++mfi)
	{
		const Box& bbx = mfi.validbox();
		auto const& phiOld = phi_old.const_array(mfi);
		auto const& term1_val = term1.array(mfi);
		//Real rot_angle = rotmat.back();
		Array2D <Real,0,3,0,3> rot_angle{0.0};
		//check this
		for(int a=0; a<rotmat.size(); a++){
				int m = rotmat[a][0];
				int n = rotmat[a][1]; 
				rot_angle(m,n) = rotmat[a][4];
				rot_angle(n,m) = rotmat[a][4];
				
		}
		int numphase = nump;
		// for(int a=0; a<nump; a++){
		// 	for(int b=0; b<nump; b++){
		// 		Print()<<rot_angle(a,b)<<" , ";
		// 	}
		// 	Print()<<"\n";
		// }
		int ani_func = funcANI;
		Real Gam = gammaa;
		Real dabb = dab;
		GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

		amrex::ParallelFor(bbx,
		[=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
		{	
			
			Array1D <Real,0,4> qabx{0.0};
			Array1D <Real,0,4> r_qabx{0.0};
			Array1D <Real,0,4> qaby{0.0};
			Array1D <Real,0,4> r_qaby{0.0};
			Array1D <Real,0,4> dadqx{0.0};
			Array1D <Real,0,4> dadqy{0.0};
			Array1D <Real,0,4> r_dadqx{0.0};
			Array1D <Real,0,4> r_dadqy{0.0};
			Array1D <Real,0,4> ac{0.0};
			Array1D <Real,0,4> q2{0.0};
			Real dqdphi_x{0.0};
			Real dqdphi_y{0.0};
			Real r_dqdphi_x{0.0};
			Real r_dqdphi_y{0.0};
			amrex::Real ani_term1{0.0}, ani_term2{0.0}, ani_term3{0.0}, ani_term4{0.0};

			for(int a=0; a<numphase; a++){

				//Print()<<"1\n";
				if(ani_func == 1){

				for(int b=0; b<numphase; b++){

					//Print()<<"2\n";

					if(b!=a){

						//Print()<<"3\n";

						qabx(0)=phiOld(i,j,k,a)*(phiOld(i+1,j,k,b)-phiOld(i-1,j,k,b))/(2*dx[0]) -phiOld(i,j,k,b)*(phiOld(i+1,j,k,a)-phiOld(i-1,j,k,a))/(2*dx[0]);
						qabx(1)=0.5*(phiOld(i,j,k,a)+phiOld(i+1,j,k,a))*(phiOld(i+1,j,k,b)-phiOld(i,j,k,b))/(dx[0]) -0.5*(phiOld(i,j,k,b)+phiOld(i+1,j,k,b))*(phiOld(i+1,j,k,a)-phiOld(i,j,k,a))/(dx[0]);
						qabx(2)=0.5*(phiOld(i,j,k,a)+phiOld(i-1,j,k,a))*(phiOld(i,j,k,b)-phiOld(i-1,j,k,b))/(dx[0]) -0.5*(phiOld(i,j,k,b)+phiOld(i-1,j,k,b))*(phiOld(i,j,k,a)-phiOld(i-1,j,k,a))/(dx[0]);
						qabx(3)=0.5*(phiOld(i,j+1,k,a)+phiOld(i,j,k,a))*(-phiOld(i-1,j+1,k,b)+phiOld(i+1,j+1,k,b)-phiOld(i-1,j,k,b)+phiOld(i+1,j,k,b))/(4*dx[0]) -0.5*(phiOld(i,j+1,k,b)+phiOld(i,j,k,b))*(-phiOld(i-1,j+1,k,a)+phiOld(i+1,j+1,k,a)-phiOld(i-1,j,k,a)+phiOld(i+1,j,k,a))/(4*dx[0]);
						qabx(4)=0.5*(phiOld(i,j,k,a)+phiOld(i,j-1,k,a))*(-phiOld(i-1,j,k,b)+phiOld(i+1,j,k,b)-phiOld(i-1,j-1,k,b)+phiOld(i+1,j-1,k,b))/(4*dx[0]) -0.5*(phiOld(i,j,k,b)+phiOld(i,j-1,k,b))*(-phiOld(i-1,j,k,a)+phiOld(i+1,j,k,a)-phiOld(i-1,j-1,k,a)+phiOld(i+1,j-1,k,a))/(4*dx[0]);

						qaby(0)=phiOld(i,j,k,a)*(phiOld(i,j+1,k,b)-phiOld(i,j-1,k,b))/(2*dx[1]) -phiOld(i,j,k,b)*(phiOld(i,j+1,k,a)-phiOld(i,j-1,k,a))/(2*dx[1]);
						qaby(1)=0.5*(phiOld(i,j,k,a)+phiOld(i+1,j,k,a))*(phiOld(i+1,j+1,k,b)-phiOld(i+1,j-1,k,b)+phiOld(i,j+1,k,b)-phiOld(i,j-1,k,b))/(4*dx[1]) -0.5*(phiOld(i,j,k,b)+phiOld(i+1,j,k,b))*(phiOld(i+1,j+1,k,a)-phiOld(i+1,j-1,k,a)+phiOld(i,j+1,k,a)-phiOld(i,j-1,k,a))/(4*dx[1]);
						qaby(2)=0.5*(phiOld(i,j,k,a)+phiOld(i-1,j,k,a))*(phiOld(i,j+1,k,b)-phiOld(i,j-1,k,b)+phiOld(i-1,j+1,k,b)-phiOld(i-1,j-1,k,b))/(4*dx[1]) -0.5*(phiOld(i,j,k,b)+phiOld(i-1,j,k,b))*(phiOld(i,j+1,k,a)-phiOld(i,j-1,k,a)+phiOld(i-1,j+1,k,a)-phiOld(i-1,j-1,k,a))/(4*dx[1]);
						qaby(3)=0.5*(phiOld(i,j+1,k,a)+phiOld(i,j,k,a))*(phiOld(i,j+1,k,b)-phiOld(i,j,k,b))/(dx[1]) -0.5*(phiOld(i,j+1,k,b)+phiOld(i,j,k,b))*(phiOld(i,j+1,k,a)-phiOld(i,j,k,a))/(dx[1]);
						qaby(4)=0.5*(phiOld(i,j,k,a)+phiOld(i,j-1,k,a))*(phiOld(i,j,k,b)-phiOld(i,j-1,k,b))/(dx[1]) -0.5*(phiOld(i,j,k,b)+phiOld(i,j-1,k,b))*(phiOld(i,j,k,a)-phiOld(i,j-1,k,a))/(dx[1]);

						//Print()<<"4\n";

						Real pi = acos(-1);	
						Real thetaf = rot_angle(a,b)*pi/180.0;
						

						r_qabx(0) = qabx(0)*cos(thetaf)-qaby(0)*sin(thetaf);
						r_qabx(1) = qabx(1)*cos(thetaf)-qaby(1)*sin(thetaf);
						r_qabx(2) = qabx(2)*cos(thetaf)-qaby(2)*sin(thetaf);
						r_qabx(3) = qabx(3)*cos(thetaf)-qaby(3)*sin(thetaf);
						r_qabx(4) = qabx(4)*cos(thetaf)-qaby(4)*sin(thetaf);

						r_qaby(0) = qabx(0)*sin(thetaf)+qaby(0)*cos(thetaf);
						r_qaby(1) = qabx(1)*sin(thetaf)+qaby(1)*cos(thetaf);
						r_qaby(2) = qabx(2)*sin(thetaf)+qaby(2)*cos(thetaf);
						r_qaby(3) = qabx(3)*sin(thetaf)+qaby(3)*cos(thetaf);
						r_qaby(4) = qabx(4)*sin(thetaf)+qaby(4)*cos(thetaf);

						for(int p=0; p<qabx.size(); p++){

							q2(p) = r_qabx(p)*r_qabx(p)+r_qaby(p)*r_qaby(p);
							Real qx2 = r_qabx(p)*r_qabx(p);
							Real qx3 = r_qabx(p)*r_qabx(p)*r_qabx(p);
							Real qy2 = r_qaby(p)*r_qaby(p);
							Real qy3 = r_qaby(p)*r_qaby(p)*r_qaby(p);
							Real q22 = q2(p)*q2(p);
							Real q23 = q2(p)*q2(p)*q2(p);
							Real qx4 = qx2*qx2;
							Real qy4 = qy2*qy2;
							Real q4 = qx3*r_qabx(p)+qy3*r_qaby(p);

							if(fabs(q2(p))>1e-15){

								ac(p) = 1.0-dabb*(3.0-4.0*(qx4+qy4)/(q2(p)*q2(p)));
								dadqx(p) = 16*dabb*(qx3/q22-r_qabx(p)*q4/q23);
								dadqy(p) = 16*dabb*(qy3/q22-r_qaby(p)*q4/q23);

							}

							else{
								ac(p) = 1.0;
								dadqx(p) = 0.0;
								dadqy(p) = 0.0;
							}

						}
	
						dqdphi_x = (phiOld(i+1,j,k,b)-phiOld(i-1,j,k,b))/(2*dx[0]);
						dqdphi_y = (phiOld(i,j+1,k,b)-phiOld(i,j-1,k,b))/(2*dx[1]);

						r_dqdphi_x = dqdphi_x*cos(thetaf)-dqdphi_y*sin(thetaf);
						r_dqdphi_y = dqdphi_x*sin(thetaf)+dqdphi_y*cos(thetaf);

						ani_term1 += 2.0*Gam*ac(0)*ac(0)*(r_qabx(0)*r_dqdphi_x+r_qaby(0)*r_dqdphi_y);
						ani_term2 += 2.0*Gam*ac(0)*q2(0)*(dadqx(0)*r_dqdphi_x+dadqy(0)*r_dqdphi_y);

						//Check this
						r_dadqx(0) = dadqx(0)*cos(thetaf)+dadqy(0)*sin(thetaf);
						r_dadqx(1) = dadqx(1)*cos(thetaf)+dadqy(1)*sin(thetaf);
						r_dadqx(2) = dadqx(2)*cos(thetaf)+dadqy(2)*sin(thetaf);
						r_dadqx(3) = dadqx(3)*cos(thetaf)+dadqy(3)*sin(thetaf);
						r_dadqx(4) = dadqx(4)*cos(thetaf)+dadqy(4)*sin(thetaf);

						r_dadqy(0) = -dadqx(0)*sin(thetaf)+dadqy(0)*cos(thetaf);
						r_dadqy(1) = -dadqx(1)*sin(thetaf)+dadqy(1)*cos(thetaf);
						r_dadqy(2) = -dadqx(2)*sin(thetaf)+dadqy(2)*cos(thetaf);
						r_dadqy(3) = -dadqx(3)*sin(thetaf)+dadqy(3)*cos(thetaf);
						r_dadqy(4) = -dadqx(4)*sin(thetaf)+dadqy(4)*cos(thetaf);
						//Check this

						//Print()<<"6\n";

						ani_term3 += 2*Gam*(ac(1)*ac(1)*qabx(1)*(-(phiOld(i+1,j,k,b)+phiOld(i,j,k,b))*0.5)-ac(2)*ac(2)*qabx(2)*(-(phiOld(i,j,k,b)+phiOld(i-1,j,k,b))*0.5))/dx[0] + 2*Gam*(ac(3)*ac(3)*qaby(3)*(-(phiOld(i,j+1,k,b)+phiOld(i,j,k,b))*0.5)-ac(4)*ac(4)*qaby(4)*(-(phiOld(i,j,k,b)+phiOld(i,j-1,k,b))*0.5))/dx[1];
						ani_term4 += 2*Gam*(ac(1)*q2(1)*(-(phiOld(i+1,j,k,b)+phiOld(i,j,k,b))*0.5)*r_dadqx(1)-ac(2)*q2(2)*(-(phiOld(i,j,k,b)+phiOld(i-1,j,k,b))*0.5)*r_dadqx(2))/dx[0] + 2*Gam*(ac(3)*q2(3)*(-(phiOld(i,j+1,k,b)+phiOld(i,j,k,b))*0.5)*r_dadqy(3)-ac(4)*q2(4)*(-(phiOld(i,j,k,b)+phiOld(i,j-1,k,b))*0.5)*r_dadqy(4))/dx[1];
						
				}
			}

			term1_val(i,j,k,a) = (-ani_term1-ani_term2+ani_term3+ani_term4);

			ani_term1=0.0;
			ani_term2=0.0;
			ani_term3=0.0;
			ani_term4=0.0;

			}
			 //Print()<<"7\n";

			else if(ani_func==0){
			term1_val(i,j,k,a) = 2.0*Gam*((phiOld(i+1,j,k,a)-2.0*phiOld(i,j,k,a)+phiOld(i-1,j,k,a))/(dx[0]*dx[0])+(phiOld(i,j+1,k,a)-2.0*phiOld(i,j,k,a)+phiOld(i,j-1,k,a))/(dx[1]*dx[1]));
			}
		}
	    });
	}
}






#endif
