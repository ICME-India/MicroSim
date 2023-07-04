#ifndef _TAU_H_
#define _TAU_H_

using namespace amrex; 

double Function_tau(MultiFab& phi_old)
{   
    Vector<Vector<Real>> tau_ab(nump-1,Vector<Real>(nump-1,0));
    Real min_tau{0};

	// for(int a=0; a<nump; a++)
    // {
    //     for(int b=a+1; b<nump; b++)
    //         tau_ab[a][b] = 0.182223*eps*(conceq[b][0][0]-conceq[a][0][0])*(conceq[b][0][0]-conceq[a][0][0])*2*Aeq[1][0][0]/(diff[1][2]*Vm);
            
    //         //Print()<<"A_liq: "<<Aeq[1][0][0]<<"\n";
    // }


    for(int a=0; a<nump-1; a++){

        tau_ab[a][nump-1] = 0.182223*eps*(conceq[nump-1][0][0]-conceq[a][0][0])*(conceq[nump-1][0][0]-conceq[a][0][0])/(dcdmu_eq[nump-1]*diff[nump-1]*Vm);
    
        if(a==0){
            min_tau = tau_ab[a][nump-1];
        }
        if(a!=0 && (tau_ab[a][nump-1]<min_tau)){
            min_tau = tau_ab[a][nump-1];
        }
        Print()<<"Tau "<<a<<" : "<<tau_ab[a][nump-1]<<"\n";
        Print()<<"conc_liq : "<<conceq[nump-1][0][0]<<"\n";
        Print()<<"conc_sol "<<a<<" : "<<conceq[a][0][0]<<"\n";
        Print()<<"dcdmu_eq : "<<dcdmu_eq[nump-1]<<"\n";
        Print()<<"diff_liq : "<<diff[0]<<" , "<<diff[nump-1]<<"\n";
        Print()<<"Vm : "<< Vm<<"\n";
        Print()<<"eps : "<< eps<<"\n";
       
    }
 
    return min_tau*0.5;
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
            //return tau_ab[0][1];
    //    }

}
 

#endif
