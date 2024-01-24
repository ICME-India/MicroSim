#ifndef _TAU_H_
#define _TAU_H_

#include <AMReX_Utility.H>
#include <matrix.h>

using namespace amrex; 

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
double FunctionTau(int i, int j, int k, int numphase, Real Tauu, amrex::Array4<Real const> const& phi){
    double sum=0.0, sum1=0.0;
    long a, b;
    for (a=0; a<numphase; a++) {
        for (b=0; b<numphase; b++) {
            if (a<b) {
                sum  += Tauu*phi(i,j,k,a)*phi(i,j,k,b);
                sum1 += phi(i,j,k,a)*phi(i,j,k,b);
            }
         }
    }
    if (sum1) {
        return sum/sum1;
    } else {
        return Tauu;
    }
}



void Calculate_Tau(MultiFab& phi_old)
{   
    Real min_tau{0.0};

        Vector<Real> deltac(numcom-1,0.0);
        Vector<Real> deltamu(numcom-1,0.0);
        Vector<Vector<Real>>prod(numcom-1,Vector<Real>(numcom-1,0.0));
        Vector<Vector<Real>>inv_dcdmu(numcom-1,Vector<Real>(numcom-1,0.0));
        Vector<Vector<Real>> tau_ab(nump,Vector<Real>(nump,0.0));

    for(int a=0; a<nump-1; a++){
    
        for(int k=0; k<numcom-1; k++){
            deltac[k] = conceq[nump-1][k]-conceq[a][k];
        }

        Print()<<"deltac : "<<deltac[0]<<"\n";

        if(numcom ==3){
           
        prod[0][0] = dcdmu[nump-1][0][0]*diff[nump-1][0][0] + dcdmu[nump-1][0][1]*diff[nump-1][1][0];
        prod[0][1] = dcdmu[nump-1][0][1]*diff[nump-1][0][0] + dcdmu[nump-1][1][1]*diff[nump-1][0][1];
        prod[1][0] = dcdmu[nump-1][0][0]*diff[nump-1][1][0] + dcdmu[nump-1][1][0]*diff[nump-1][1][1];
        prod[1][1] = dcdmu[nump-1][1][0]*diff[nump-1][0][1] + dcdmu[nump-1][1][1]*diff[nump-1][1][1];

        Real det = prod[0][0]*prod[1][1] - prod[0][1]*prod[1][0];
        inv_dcdmu[0][0] = prod[1][1]/det;
        inv_dcdmu[0][1] = -prod[0][1]/det;
        inv_dcdmu[1][0] = -prod[1][0]/det;
        inv_dcdmu[1][1] = prod[0][0]/det; 

        deltamu[0] = inv_dcdmu[0][0]*deltac[0] + inv_dcdmu[0][1]*deltac[1];
        deltamu[1] = inv_dcdmu[1][0]*deltac[0] + inv_dcdmu[1][1]*deltac[1]; 

        }

        else if(numcom ==2){
            prod[0][0] = dcdmu_eq[nump-1][0][0]*diff[nump-1][0][0];

            inv_dcdmu[0][0] = 1/prod[0][0];

            deltamu[0] = inv_dcdmu[0][0]*deltac[0];
        }

        // Print()<<"prod: "<<prod[0][0]<<"\n";

        // Print()<<"inv_dcdmu"<<inv_dcdmu[0][0]<<"\n";

        // Print()<<"deltamu"<<deltamu[0]<<"\n";

        double sum=0.0;
        for (int k=0; k<numcom-1; k++) {
        sum += deltamu[k]*deltac[k];
        }

        //Print()<<"sum: "<<sum<<"\n";

        tau_ab[a][nump-1] = sum*eps*(0.182223)/Vm;
        
        //Print()<<"Tau["<<a<<",nump-1]"<<tau_ab[a][nump-1]<<"\n";
       
        tau_ab[nump-1][a] = tau_ab[a][nump-1];
        

        if (a==0) {
        min_tau = tau_ab[a][nump-1];
        }

        if (tau_ab[a][nump-1] < min_tau) {
        min_tau = tau_ab[a][nump-1];
        }

        deltac.clear();
        deltamu.clear();
        prod.clear();
        inv_dcdmu.clear();

        deltac = Vector<Real>(numcom-1,0.0);
        deltamu = Vector<Real>(numcom-1,0.0);
        prod = Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0.0));
        inv_dcdmu = Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0.0));
    }
        for (int a=0; a<nump; a++) {
            for (int b=0; b<nump; b++) {
                tau_ab[a][b] = 0.5*min_tau;
            }
        }
        
        tau = 0.5*min_tau;
        
        Print()<<"Tau : "<<tau<<"\n";
       
}

 

#endif
