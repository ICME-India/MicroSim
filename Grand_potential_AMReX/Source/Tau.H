#ifndef _TAU_H_
#define _TAU_H_

using namespace amrex; 

double Function_tau(MultiFab& phi_old)
{   
    Vector<Vector<Real>> tau_ab(nump,Vector<Real>(nump,0));

	for(int a=0; a<nump; a++)
    {
        for(int b=a+1; b<nump; b++)
            tau_ab[a][b] = 0.182223*eps*(ceq[0][2]-ceq[1][2])*(ceq[0][2]-ceq[1][2])*2*Aeq[1][0][0]/(diff[1][2]*Vm);
            
            //Print()<<"A_liq: "<<Aeq[1][0][0]<<"\n";
    }

    // if(nump>2)
    // {
    //     for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    //     {
    //         const Box& vbx = mfi.validbox();
    //         auto const& phiold = phi_old.array(mfi);
            
    //         for(int a=0; a<nump; a++){
    //             for(int b=a+1; b<nump; b++){

    //                 amrex::ParallelFor(vbx,[=] AMREX_GPU_DEVICE (int i, int j, int k)
    //                 {
    //                     sum = sum + tau_ab[a][b]*phiold(i,j,k)*(1-phiold(i,j,k));
    //                     sum1 = sum1 + phiold(i,j,k)*(1-phiold(i,j,k));
    //                 });
    //             }
    //         }
            
      
    //     }
    // } 
    // else
    // {
    //     sum1=0.0;
    // } 

    // if(sum1!=0)
    //     {
    //         return sum/sum1;
    //     }
    // else{
            return tau_ab[0][1];
    //    }

}
 

#endif
