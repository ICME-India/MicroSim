#ifndef _CAL_JAT_H_
#define _CAL_JAT_H_

#include <AMReX_Utility.H>
#include <Variables.h>

using namespace amrex;


//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



//                                        THIS FUNCTION IS NOT IN USE



//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx




void Calc_jat(MultiFab& phi_new, MultiFab& phi_old, MultiFab& jat, MultiFab& mu_old,Geometry const& geom){


    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();
		Array4<Real> const& phiOld = phi_old.array(mfi);
		Array4<Real> const& phiNew = phi_new.array(mfi);
		Array4<Real> const& mu = mu_old.array(mfi);
		Array4<Real> const& jatt = jat.array(mfi);
		
		Array1D <Real,0,phasecount-1> BB;
		for(int a=0; a<nump; a++){
		BB(a) = B[a];
		}
		
		Array1D <Real,0,phasecount-1> dercmu;
		for(int a=0; a<nump; a++){
		dercmu(a) = dcdmu[a];
		}

		Array1D <Real,0,phasecount-1> diffs;
		for(int a=0; a<nump; a++){
		diffs(a) = diff[a];
		}

		Real time_step = dt;
		int numphase = nump;
		Print()<<"outside for jat\n";
		
		GpuArray<Real,AMREX_SPACEDIM> delta = geom.CellSizeArray();
	
		amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {	

        Real gradphix_l = (phiOld(i+1,j+1,k,nump-1)-phiOld(i-1,j+1,k,nump-1)+phiOld(i+1,j,k,nump-1)-phiOld(i-1,j,k,nump-1))/(4*delta[X]);
        Real normgradphiy_l = sqrt(gradphix_l*gradphix_l+pow(((phiOld(i,j+1,k,nump-1)-phiOld(i,j,k,nump-1))/delta[Y]),2));

        Real gradphiy_l = (phiOld(i+1,j+1,k,nump-1)-phiOld(i+1,j-1,k,nump-1)+phiOld(i,j+1,k,nump-1)-phiOld(i,j-1,k,nump-1))/(4*delta[Y]);
        Real normgradphix_l = sqrt(gradphiy_l*gradphiy_l+pow(((phiOld(i+1,j,k,nump-1)-phiOld(i,j,k,nump-1))/delta[Z]),2));

        Print()<<"grad_l calculated\n";

        Real scalprod{0.0};

        for (int a=0; a < nump-1; a++) {

          Real gradphix      = (phiOld(i+1,j+1,k,a)-phiOld(i-1,j+1,k,a)+phiOld(i+1,j,k,a)-phiOld(i-1,j,k,a))/(4*delta[X]);
          Real normgradphi   = sqrt(gradphix*gradphix+pow(((phiOld(i,j+1,k,a)-phiOld(i,j,k,a))/delta[Y]),2));

          Print()<<"grad_s calculated\n";

            if (normgradphi != 0.0) {
                
                Real cl_right  = (mu(i,j+1,k)-BB(nump-1))*dercmu(nump-1);
                Real cl_center = (mu(i,j,k)-BB(nump-1))*dercmu(nump-1);

                Real cs_right  = (mu(i,j+1,k)-BB(a))*dercmu(a);
                Real cs_center = (mu(i,j,k)-BB(a))*dercmu(a);

                Print()<<"cs & cl calculated\n";

                scalprod  = ((phiOld(i,j+1,k,a)-phiOld(i,j,k,a))/delta[Y])*((phiOld(i,j+1,k,nump-1)-phiOld(i,j,k,nump-1))/delta[Y]);
                scalprod += gradphix*gradphix_l;
                
                if (normgradphiy_l > 0.0) {
                  scalprod /= (normgradphiy_l*normgradphi);
                }

                jatt(i,j,k,Y) += (1.0-diffs(a)/diffs(nump-1))*fabs(scalprod)*(0.5/sqrt(2.0))*((cl_right-cs_right)*(phiNew(i,j+1,k,a)-phiOld(i,j+1,k,a)) + (cl_center-cs_center)*(phiNew(i,j,k,a)-phiOld(i,j,k,a)))*((phiOld(i,j+1,k,a)-phiOld(i,j,k,a))/delta[Y])/normgradphi;

                Print()<<"jat_y calculated\n";
            
          }
        
        Real gradphiy = (phiOld(i+1,j+1,k,a)-phiOld(i-1,j+1,k,a)+phiOld(i+1,j,k,a)-phiOld(i-1,j,k,a))/(4*delta[X]);
        normgradphi = sqrt(gradphiy*gradphiy + pow(((phiOld(i+1,j,k,a)-phiOld(i,j,k,a))/delta[X]),2));

        Print()<<"grad_s calculated\n";

            if (normgradphi != 0.0) {
              Real cl_front  = (mu(i+1,j,k)-BB(nump-1))*dercmu(nump-1);
              Real cl_center = (mu(i,j,k)-BB(nump-1))*dercmu(nump-1);
              
              Real cs_front  = (mu(i+1,j,k)-BB(a))*dercmu(a);
              Real cs_center = (mu(i,j,k)-BB(a))*dercmu(a);

              Print()<<"cs & cl calculated\n";

              scalprod  = ((phiOld(i+1,j,k,a)-phiOld(i,j,k,a))/delta[X])*((phiOld(i+1,j,k,nump-1)-phiOld(i,j,k,nump-1))/delta[X]);
              scalprod += gradphiy*gradphiy_l;

              if (normgradphix_l > 0.0) {
                scalprod /= (normgradphix_l*normgradphi);
              }

              jatt(i,j,k,X) += (1.0-diffs(a)/diffs(nump-1))*fabs(scalprod)*(0.5/sqrt(2.0))*((cl_front-cs_front)*(phiNew(i+1,j,k,a)-phiOld(i+1,j,k,a)) + (cl_center-cs_center)*(phiNew(i,j,k,a)-phiOld(i,j,k,a)))*((phiOld(i+1,j,k,a)-phiOld(i,j,k,a))/delta[X])/normgradphi;

              Print()<<"jatx calculated\n";
            }
        
        }

      });
    }
}

#endif