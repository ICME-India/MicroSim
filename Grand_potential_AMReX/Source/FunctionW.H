#ifndef _FUNCTIONW_H_
#define _FUNCTIONW_H_

using namespace amrex;

void computeterm2(amrex::MultiFab& term2, amrex::MultiFab& phi_old, Geometry const& geom)
{
	BL_PROFILE("computeterm2()");
	
	for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();
		Array4<Real> const& phiOld = phi_old.array(mfi);
		Array4<Real> const& term = term2.array(mfi);
		Real GAM = gammaa;
		int numphase = nump;
	
		// amrex::ParallelFor(vbx,
        // [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        // {
		// 	term(i,j,k) = 9.0*GAM*2.0*phiold(i,j,k)*(1.0-phiold(i,j,k))*(1.0 - 2.0*phiold(i,j,k));
		// });
		GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

		amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {	
			Real sum{0.0};
			Real phi_phi{0.0};
			
			for(int a =0; a<numphase; a++){

				for(int b=0; b<numphase;b++){
					if(b!=a){
						sum = sum + 2.0*GAM*phiOld(i,j,k,a)*phiOld(i,j,k,b)*phiOld(i,j,k,b);
					}
				}
				sum*=9.0;
				
				for(int b =0; b<numphase; b++){
					for(int c = 0; c<numphase; c++){

						if(b!=a && c!=a && b<c){

							Real div_b = (phiOld(i+1,j,k,b)-2.0*phiOld(i,j,k,b)+phiOld(i-1,j,k,b))/(dx[0]*dx[0])+(phiOld(i,j+1,k,b)-2.0*phiOld(i,j,k,b)+phiOld(i,j-1,k,b))/(dx[1]*dx[1]);
							Real div_c = (phiOld(i+1,j,k,c)-2.0*phiOld(i,j,k,c)+phiOld(i-1,j,k,c))/(dx[0]*dx[0])+(phiOld(i,j+1,k,c)-2.0*phiOld(i,j,k,c)+phiOld(i,j-1,k,c))/(dx[1]*dx[1]);						

							if(fabs(div_b)>0.0 && fabs(div_c)>0.0){ 
							phi_phi = phiOld(i,j,k,b)*phiOld(i,j,k,c);
							sum += GAM*phi_phi;
							}
						}
					}
				}
				

			term(i,j,k,a) = sum;

			//term(i,j,k,a) = 18.0*GAM*phiOld(i,j,k,a)*(1.0-phiOld(i,j,k,a))*(1.0-2.0*phiOld(i,j,k,a));
			phi_phi=0.0;
			sum=0.0;

			}

			// term(i,j,k,0) = 18*GAM*phiOld(i,j,k,0)*(pow(phiOld(i,j,k,1),2)+pow(phiOld(i,j,k,2),2)) + GAM*(phiOld(i,j,k,1)*phiOld(i,j,k,2));
			// term(i,j,k,1) = 18*GAM*phiOld(i,j,k,1)*(pow(phiOld(i,j,k,0),2)+pow(phiOld(i,j,k,2),2)) + GAM*(phiOld(i,j,k,0)*phiOld(i,j,k,2));
			// term(i,j,k,2) = 18*GAM*phiOld(i,j,k,2)*(pow(phiOld(i,j,k,1),2)+pow(phiOld(i,j,k,0),2)) + GAM*(phiOld(i,j,k,1)*phiOld(i,j,k,0));
		});



	
	}
}



/*

void computeterm2(amrex::MultiFab& term2, amrex::MultiFab& phi_old)
{
	BL_PROFILE("computeterm2()");
	
	for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();
		Array4<Real> const& phiOld = phi_old.array(mfi);
		Array4<Real> const& term2val = term2.array(mfi);
		Real Gam = gammaa;
	
		amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {	
			Real db_well1, db_well2;

			for(int a=0; a<ntp; a++){
				for(int b=0; b<ntp; b++){

					if(a!=b && b>a){
					db_well1 += 18*Gam*phiOld(i,j,k,a)*phiOld(i,j,k,b)*phiOld(i,j,k,b);

					for(int c=0; c<ntp; c++){
						if(a!=b && b!=c && b<c){
							db_well2 += Gam*phiOld(i,j,k,b)*phiOld(i,j,k,c);
						}
					}
					}
				}
				term2val(i,j,k,a) = db_well1+db_well2;
			}

		});
	
	}
}

*/

#endif
