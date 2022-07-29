#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#include "solverloop/defines.h"
#include "solverloop/defines1.h"
#include "solverloop/ThermoCL_2.h"
#include "solverloop/ThermoCL_2.c"
//#include "solverloop/GibbsEnergyData1.h"
#include "solverloop/GibbsEnergyData2.h"
#include "solverloop/FunctionF_2_c_mu_NR_GEData2.h"
//#include "solverloop/GibbsEnergyData_1.h"

/*
      0          ---- j ----->           nx (X)
    
 0    #  #  #  #  #  #  #  #  #  #  #  #  # 
               
      #  #  #  #  #  #  #  #  #  #  #  #  #
               
 |    #  #  #  #  #  #  #  #  #  #  #  #  # 
 |             
 |    #  #  #  #  0  1  2  #  #  #  #  #  # 
               
 i    #  #  #  #  3  4  5  #  #  #  #  #  # 
               
 |    #  #  #  #  6  7  8  #  #  #  #  #  # 
 |             
 V    #  #  #  #  #  #  #  #  #  #  #  #  # 
               
      #  #  #  #  #  #  #  #  #  #  #  #  # 
               
 ny   #  #  #  #  #  #  #  #  #  #  #  #  # 
 
 (Y)

*/

struct grid {
    double phi;
    double c1;
    double mu[1];
};
struct csle {
    double c1l;
    double c1s;
};
struct pfmpar {
    double E0;
    double surfTen;
    double ee;
    double w;
    double deltax;
    double deltay;
    double deltat;
    //double dx;
    //double dy;
    //double dt;
    double eesqrt;
    double Er;
    double IntMob;
    double IntMobInv;
};
struct pfmval {
  double Rotation_matrix[2][2][3][3];
  double Inv_Rotation_matrix[2][2][3][3];
    double Tr;
    double sigma;
    double Vm;
    double D11l;
    double D11s;
    double phisolid;
    double philiquid;
    double Rg;
    double T0;
    double Teq;
    double Tfill;
    double lrep;
    double c1l_Initial;
    double c1s_Initial;
    double c1l_1stguess;
    double c1s_1stguess;
    double a2;
    double rad;
    double epsc;
    double intwidth;
    double dxIntfcPoints;
    double dtParam;
    double epsm;
    double InterfaceMobility;
    double RefD;
    double angle;
  double TLiquidus;
  double Toffset;
  //double PosOffset;
  //double TG;
  //double Vp;
  double TPosOffset; 
  double TGRADIENT;
  double velocity;
  double NoiseFac;
  double interfaceUplimit;
  double interfaceDownlimit;
  long shift_OFFSET;
  int thermophase[npha];
    int   nproc;
    int   jNx;
    int   iNy;
    int   jDimX;
    int   iDimY;
    int   ntimesteps;
    int   savetime;
    int   myrank;
  int   ISOTHERMAL;
};
struct propmatf3 {
  double ceq[npha][npha][nsol];
  double cfill[npha][npha][nsol];
  double slopes[npha][npha][nsol];
  double dcbdT[npha][npha][nsol];
  double A[npha][nsol][nsol]; 
  double DELTA_T[npha][npha]; 
  double DELTA_C[npha][nsol]; 
  double dcbdT_phase[npha][nsol];
  double B[npha][nsol];
  double Beq[npha][nsol];
  double dBbdT[npha][nsol];
  double C[npha];
  double cmu[npha][nsol][nsol];
  double muc[npha][nsol][nsol];
};
struct propmatf4 {
  double ceq[npha][npha][nsol];
  double cfill[npha][npha][nsol];
  double slopes[npha][npha][nsol];
  double dcbdT[npha][npha][nsol];
  double A[npha][nsol][nsol]; 
  double DELTA_T[npha][npha]; 
  double DELTA_C[npha][nsol]; 
  double dcbdT_phase[npha][nsol];
  double B[npha][nsol];
  double Beq[npha][nsol];
  double dBbdT[npha][nsol];
  double C[npha];
  double cmu[npha][nsol][nsol];
  double muc[npha][nsol][nsol];
};
struct propmatf4spline {
  double A[npha][nsol][nsol]; 
  double B[npha][nsol];
  double C[npha];
};


__kernel void SolverCsClEq_2(__global struct grid *gridOld, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf3 *propf3) {
    
    int i;
    int j;
    int nx;
    int ny;
    int index;
    int interface, phaseliquid, phasesolid, bulkphase;
    double cg[nsol], mu[nsol], cphas[nsol];
    
    double Ti, tmp0, fi;

    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;

    index = (i*nx + j);
    
    Ti = temp[j];
    fi = gridOld[index].phi;
    mu[0] = gridOld[index].mu[0];
    
    interface = 1; 
    phaseliquid = 1; 
    phasesolid = 0; 
    
    if ( fi >= pfmdat->interfaceUplimit ) { 
      interface = 0;
    }
    else if ( fi <= pfmdat->interfaceDownlimit ) { 
      interface = 0; 
    }
    
    mu[0] = gridOld[index].mu[0];
    
    if ( interface ) {
      
      cg[0] = pfmdat->c1s_1stguess; 
      
      c_mu(mu, cphas, Ti, phasesolid, cg, tstep[0], i, j, pfmdat->thermophase[phasesolid]);
      
      cscl[index].c1s = cphas[0];
      
      cg[0] = pfmdat->c1l_1stguess; 
      
      c_mu(mu, cphas, Ti, phaseliquid, cg, tstep[0], i, j, pfmdat->thermophase[phaseliquid]);
      
      cscl[index].c1l = cphas[0];
      
    }
 
}

__kernel void SolverCsClEq_3(__global struct grid *gridOld, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf3 *propf3) {
    
    int i, j;
    int nx, ny;
    int index;
    
    int interface, phaseliquid, phasesolid;
    int bulkphase; 
    double cg[nsol], ci[nsol];
    double Ti;
    double tmp0;
    double fi, mu[1];

    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;

    index = (i*nx + j);
    
    Ti = temp[j];
 
    fi = gridOld[index].phi;
    
    mu[0] = gridOld[index].mu[0];
    
    interface = 1; 
    phaseliquid = 1; 
    phasesolid = 0; 
    
    if ( fi >= pfmdat->interfaceUplimit ) { 
      interface = 0;
    }
    else if ( fi <= pfmdat->interfaceDownlimit ) { 
      interface = 0; 
    }
    
    if ( interface ) {

      cscl[index].c1s = propf3->cmu[0][0][0] * ( gridOld[index].mu[0] - ( propf3->Beq[0][0] + propf3->dBbdT[0][0] * (Ti-pfmdat->Teq) ) );

      cscl[index].c1l = propf3->cmu[1][0][0] * ( gridOld[index].mu[0] - ( propf3->Beq[1][0] + propf3->dBbdT[1][0] * (Ti-pfmdat->Teq) ) );
      
    }
    
 
}

__kernel void SolverCsClEq_4(__global struct grid *gridOld, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf4 *propf4, __constant struct propmatf4spline *propf4spline) {
    
    int i, j;
    int nx, ny;
    int index;
    
    int interface, phaseliquid, phasesolid;
    int bulkphase; 
    double cg[nsol], ci[nsol];
    double Ti;
    double tmp0;
    double fi, mu[1];
    double muc[1], cmu[1];

    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;

    index = (i*nx + j);
    
    Ti = temp[j];
    
    fi = gridOld[index].phi;
    
    mu[0] = gridOld[index].mu[0];
    
    interface = 1; 
    phaseliquid = 1; 
    phasesolid = 0; 
    
    if ( fi >= pfmdat->interfaceUplimit ) { 
      interface = 0;
    }
    else if ( fi <= pfmdat->interfaceDownlimit ) { 
      interface = 0; 
    }
    
    
    if ( interface ) {

        if ( pfmdat->ISOTHERMAL ) { 
          cscl[index].c1s = propf4->cmu[0][0][0] * ( gridOld[index].mu[0] - propf4->B[0][0] );
          
          cscl[index].c1l = propf4->cmu[1][0][0] * ( gridOld[index].mu[0] - propf4->B[1][0] );
          
        }
        else { 
          
          muc[0] = 2.0 * propf4spline[j].A[0][0][0]; 
          
          cmu[0] = 1.0 / muc[0]; 
          
          cscl[index].c1s = cmu[0] * ( gridOld[index].mu[0] - propf4spline[j].B[0][0] );
          
          
          muc[0] = 2.0 * propf4spline[j].A[1][0][0]; 
          
          cmu[0] = 1.0 / muc[0]; 
          
          cscl[index].c1l = cmu[0] * ( gridOld[index].mu[0] - propf4spline[j].B[1][0] );
          
        }
        
    }
    
 
}

__kernel void SolverCWoatr_2(__global struct grid *gridOld, __global struct grid *gridNew, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf3 *propf3) { 

    int i;
    int j;
    int k;
    int ii;    
    int nx;
    int ny;
    int nz;
    int index;
    int interface, phaseliquid, phasesolid, bulkphase;

    double hphi[5], fhi[5], hphid[5];
    double term1lx;
    double term1sx;
    double term1ly;
    double term1sy;
    double term1l;
    double term1s;
    double c1dot;

    double gradmux1[1], gradmux2[1], gradmuy1[1], gradmuy2[1], mu[1]; 
    double Ti, ddgldx1ldx1l, ddgsdx1sdx1s, dcdmu[2][1][1], dc_dmu[5][2][1][1];
    double Da[5], Damidx1, Damidx2, Damidy1, Damidy2, divflux, deltamu; 
    double suma, deltac, dcdmudenom, y[2], retmu[1], Tij[5]; 
    double retdmuphase[1], DELTAT, cg[1], dcbdT_phase[2][1], c_tdt[1], sum_dcbdT; 

    struct grid stgridO[9];
    struct csle stcscl[5];

    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    Ti = temp[j];
    
    if ( i != 0 && i != ny-1 && j != 0 && j != nx-1 ) {

        index = (i*nx + j);

        stgridO[0] = gridOld[(i-1)*nx + (j-1)];
        stgridO[1] = gridOld[(i-1)*nx + ( j )];
        stgridO[2] = gridOld[(i-1)*nx + (j+1)];
        stgridO[3] = gridOld[( i )*nx + (j-1)];
        stgridO[4] = gridOld[( i )*nx + ( j )];
        stgridO[5] = gridOld[( i )*nx + (j+1)];
        stgridO[6] = gridOld[(i+1)*nx + (j-1)];
        stgridO[7] = gridOld[(i+1)*nx + ( j )];
        stgridO[8] = gridOld[(i+1)*nx + (j+1)];

        stcscl[0] = cscl[(i-1)*nx + ( j )];
        stcscl[1] = cscl[( i )*nx + (j-1)];
        stcscl[4] = cscl[( i )*nx + ( j )];
        stcscl[2] = cscl[( i )*nx + (j+1)];
        stcscl[3] = cscl[(i+1)*nx + ( j )];

        Tij[0] = temp[( j )];
        Tij[1] = temp[(j-1)];
        Tij[4] = temp[( j )];
        Tij[2] = temp[(j+1)];
        Tij[3] = temp[( j )];

        fhi[0] = stgridO[1].phi;
        fhi[1] = stgridO[3].phi;
        fhi[4] = stgridO[4].phi;
        fhi[2] = stgridO[5].phi;
        fhi[3] = stgridO[7].phi;

        hphid[0] = fhi[0]; //stgridO[1].phi;
        hphid[1] = fhi[1]; //stgridO[3].phi;
        hphid[4] = fhi[4]; //stgridO[4].phi;
        hphid[2] = fhi[2]; //stgridO[5].phi;
        hphid[3] = fhi[3]; //stgridO[7].phi;
        

        for ( ii = 0; ii < 5; ii++ ) { 
          //printf("%d, %d, %d, %d, %le\n", tstep[0], i, j, ii, fhi[ii]);

          interface = 1;
          phaseliquid = 1; 
          phasesolid = 0; 
          
          if (fhi[ii] >= pfmdat->interfaceUplimit) { 
            bulkphase = 0;
            interface = 0;
          }
          else if (fhi[ii] <= pfmdat->interfaceDownlimit) { 
            bulkphase = 1;
            interface = 0;
          }

          if ( interface ) { 
            
            y[0] = stcscl[ii].c1s; 
            y[1] = 1.0-y[0];
            
            dMudc(Tij[ii], y, retdmuphase, pfmdat->thermophase[phasesolid]);
            dc_dmu[ii][0][0][0] = 1.0 / retdmuphase[0];
            
            y[0] = stcscl[ii].c1l; 
            y[1] = 1.0-y[0];
            
            dMudc(Tij[ii], y, retdmuphase, pfmdat->thermophase[phaseliquid]);
            dc_dmu[ii][1][0][0] = 1.0 / retdmuphase[0];

            Da[ii] = pfmdat->D11l * (1.0-hphid[ii]) * dc_dmu[ii][1][0][0] + pfmdat->D11s * (hphid[ii]) * dc_dmu[ii][0][0][0];
          }
          else { 
              
              y[0] = stgridO[ii].c1;
              y[1] = 1.0 -y[0];
              
            if ( bulkphase == 0 ) { 
              
              dMudc(Tij[ii], y, retdmuphase, pfmdat->thermophase[0]);
              
              Da[ii] = pfmdat->D11s * ( 1.0 / retdmuphase[0] );
              
            }
            else if ( bulkphase == 1 ) { 
              
              dMudc(Tij[ii], y, retdmuphase, pfmdat->thermophase[1]);
              
              Da[ii] = pfmdat->D11l * ( 1.0 / retdmuphase[0] );
              
            }
          }
        }
        
        Damidx1 = ( Da[2] + Da[4] ) / 2.0;
        Damidx2 = ( Da[4] + Da[1] ) / 2.0;
        Damidy1 = ( Da[3] + Da[4] ) / 2.0;
        Damidy2 = ( Da[4] + Da[0] ) / 2.0;
        
        gradmux1[0] = ( stgridO[5].mu[0] - stgridO[4].mu[0] ) / pfmvar->deltax;
        gradmux2[0] = ( stgridO[4].mu[0] - stgridO[3].mu[0] ) / pfmvar->deltax;
        gradmuy1[0] = ( stgridO[7].mu[0] - stgridO[4].mu[0] ) / pfmvar->deltay;
        gradmuy2[0] = ( stgridO[4].mu[0] - stgridO[1].mu[0] ) / pfmvar->deltay;

        divflux = ( Damidx1*gradmux1[0] - Damidx2*gradmux2[0] ) / pfmvar->deltax + ( Damidy1*gradmuy1[0] - Damidy2*gradmuy2[0] ) / pfmvar->deltay;
        
        interface = 1;
        phaseliquid = 1; 
        phasesolid = 0; 
        
        if (stgridO[4].phi >= pfmdat->interfaceUplimit) {
          bulkphase=0;
          interface = 0;
        }
        else if (stgridO[4].phi <= pfmdat->interfaceDownlimit) {
          bulkphase = 1;
          interface = 0; 
        }

        if ( interface == 0) { 
          
          gridNew[index].c1 = stgridO[4].c1 + pfmvar->deltat * divflux;
          
          y[0] = gridNew[index].c1; 
          y[1] = 1.0 - y[0];
          
          Mu(Ti, y, retmu, pfmdat->thermophase[bulkphase]);
          
          // if ( phasesolid ) { 
          //   //Mu_0(Ti, y, retmu);
          //   Mu(Ti, y, retmu, pfmdat->thermophase[0]);
          // }
          // else if ( phaseliquid ) { 
          //   //Mu_1(Ti, y, retmu);
          //   Mu(Ti, y, retmu, pfmdat->thermophase[1]);
          // }
          
          deltamu = retmu[0] - stgridO[4].mu[0];
          gridNew[index].mu[0] = retmu[0];
          
        }
        else if ( interface ) { 
        
        if ( !pfmdat->ISOTHERMAL ) {
          DELTAT = pfmvar->deltat * ( -pfmdat->TGRADIENT * pfmdat->velocity ); 
          
          cg[0] = pfmdat->c1s_1stguess; 
          
          mu[0] = gridOld[index].mu[0];
          
          c_mu(mu, c_tdt, Ti+DELTAT, phasesolid, cg, tstep[0], i, j, pfmdat->thermophase[phasesolid]);
          
          dcbdT_phase[0][0] = c_tdt[0] - stcscl[4].c1l;
          
          cg[0] = pfmdat->c1l_1stguess; 
          
          c_mu(mu, c_tdt, Ti+DELTAT, phaseliquid, cg, tstep[0], i, j, pfmdat->thermophase[phaseliquid]);
          
          dcbdT_phase[1][0] = c_tdt[0] - stcscl[4].c1s;
          
        } 

          suma = ( stcscl[4].c1s - stcscl[4].c1l ) * (1.0) * ( gridNew[index].phi - stgridO[4].phi );
          
          sum_dcbdT = fhi[4] * dcbdT_phase[0][0] + (1.0-fhi[4]) * dcbdT_phase[1][0];

          deltac = pfmvar->deltat * divflux; 

          gridNew[index].c1 = stgridO[4].c1 + deltac;
          
          if ( pfmdat->ISOTHERMAL ) { 
            deltac = deltac - suma;
          }
          else { 
            deltac = deltac - suma  - sum_dcbdT;
          }
          
          hphi[4] = fhi[4];
          dcdmudenom = dc_dmu[4][0][0][0]*hphi[4] + dc_dmu[4][1][0][0]*(1.0-hphi[4]);
          
          gridNew[index].mu[0] = stgridO[4].mu[0] + deltac / dcdmudenom;
          

        }
        
    }
}

__kernel void SolverCWoatr_3(__global struct grid *gridOld, __global struct grid *gridNew, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf3 *propf3) { 

    int i;
    int j;
    int k;
    int ii;    
    int nx;
    int ny;
    int nz;
    int index;
    int interface, phaseliquid, phasesolid, bulkphase;

    double hphi[5], fhi[5], hphid[5];
    double term1lx;
    double term1sx;
    double term1ly;
    double term1sy;
    double term1l;
    double term1s;
    double c1dot;

    double gradmux1[1], gradmux2[1], gradmuy1[1], gradmuy2[1], mu[1]; 
    double Ti, ddgldx1ldx1l, ddgsdx1sdx1s, dcdmu[2][1][1], dc_dmu[5][2][1][1];
    double Da[5], Damidx1, Damidx2, Damidy1, Damidy2, divflux, deltamu; 
    double suma, deltac, dcdmudenom, y[2], retmu[1]; 
    double retdmuphase[1], DELTAT, cgs, cgl, dcbdT_phase[2][1], c_tdt[1], sum_dcbdT; 

    struct grid stgridO[9];
    struct csle stcscl[5];

    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    Ti = temp[j];
    
    if ( i != 0 && i != ny-1 && j != 0 && j != nx-1 ) {

        index = (i*nx + j);

        stgridO[0] = gridOld[(i-1)*nx + (j-1)];
        stgridO[1] = gridOld[(i-1)*nx + ( j )];
        stgridO[2] = gridOld[(i-1)*nx + (j+1)];
        stgridO[3] = gridOld[( i )*nx + (j-1)];
        stgridO[4] = gridOld[( i )*nx + ( j )];
        stgridO[5] = gridOld[( i )*nx + (j+1)];
        stgridO[6] = gridOld[(i+1)*nx + (j-1)];
        stgridO[7] = gridOld[(i+1)*nx + ( j )];
        stgridO[8] = gridOld[(i+1)*nx + (j+1)];

        stcscl[0] = cscl[(i-1)*nx + ( j )];
        stcscl[1] = cscl[( i )*nx + (j-1)];
        stcscl[4] = cscl[( i )*nx + ( j )];
        stcscl[2] = cscl[( i )*nx + (j+1)];
        stcscl[3] = cscl[(i+1)*nx + ( j )];

        fhi[0] = stgridO[1].phi;
        fhi[1] = stgridO[3].phi;
        fhi[4] = stgridO[4].phi;
        fhi[2] = stgridO[5].phi;
        fhi[3] = stgridO[7].phi;

        hphid[0] = fhi[0]; //stgridO[1].phi;
        hphid[1] = fhi[1]; //stgridO[3].phi;
        hphid[4] = fhi[4]; //stgridO[4].phi;
        hphid[2] = fhi[2]; //stgridO[5].phi;
        hphid[3] = fhi[3]; //stgridO[7].phi;
        

        for ( ii = 0; ii < 5; ii++ ) { 
          //printf("%d, %d, %d, %d, %le\n", tstep[0], i, j, ii, fhi[ii]);

          interface = 1;
          phaseliquid = 1; 
          phasesolid = 0; 
          
          if (fhi[ii] >= pfmdat->interfaceUplimit) { 
            bulkphase = 0;
            interface = 0;
          }
          else if (fhi[ii] <= pfmdat->interfaceDownlimit) { 
            bulkphase = 1;
            interface = 0;
          }

          if ( interface ) { 

            //dc_dmu[ii][0][0][0] = 1.0 / (2.0*propf3->A[0][0][0]);
            //dc_dmu[ii][1][0][0] = 1.0 / (2.0*propf3->A[1][0][0]);

            dc_dmu[ii][0][0][0] = propf3->cmu[0][0][0];
            dc_dmu[ii][1][0][0] = propf3->cmu[1][0][0];

            Da[ii] = pfmdat->D11l * (1.0-hphid[ii]) * dc_dmu[ii][1][0][0] + pfmdat->D11s * (hphid[ii]) * dc_dmu[ii][0][0][0];
            
          }
          else { 
            if ( bulkphase == 0 ) { 
              
              //Da[ii] = pfmdat->D11s * ( 1.0 / (2.0*propf3->A[0][0][0]) );
              
              Da[ii] = pfmdat->D11s * propf3->cmu[0][0][0];
              
            }
            else if ( bulkphase == 1 ) { 
              
              Da[ii] = pfmdat->D11l * propf3->cmu[1][0][0];
              
            }
          }
        }
        
        Damidx1 = ( Da[2] + Da[4] ) / 2.0;
        Damidx2 = ( Da[4] + Da[1] ) / 2.0;
        Damidy1 = ( Da[3] + Da[4] ) / 2.0;
        Damidy2 = ( Da[4] + Da[0] ) / 2.0;
        
        
        gradmux1[0] = ( stgridO[5].mu[0] - stgridO[4].mu[0] ) / pfmvar->deltax;
        gradmux2[0] = ( stgridO[4].mu[0] - stgridO[3].mu[0] ) / pfmvar->deltax;
        gradmuy1[0] = ( stgridO[7].mu[0] - stgridO[4].mu[0] ) / pfmvar->deltay;
        gradmuy2[0] = ( stgridO[4].mu[0] - stgridO[1].mu[0] ) / pfmvar->deltay;

        divflux = ( Damidx1*gradmux1[0] - Damidx2*gradmux2[0] ) / pfmvar->deltax + ( Damidy1*gradmuy1[0] - Damidy2*gradmuy2[0] ) / pfmvar->deltay;
        
        
        interface = 1;
        phaseliquid = 1; 
        phasesolid = 0;
        
        if (stgridO[4].phi >= pfmdat->interfaceUplimit) {
          bulkphase = 0;
          interface = 0;
        }
        else if (stgridO[4].phi <= pfmdat->interfaceDownlimit) {
          bulkphase = 1;
          interface = 0; 
        }

        if ( interface == 0) { 

          gridNew[index].c1 = stgridO[4].c1 + pfmvar->deltat * divflux;
          
          if ( bulkphase == 0 ) { 
            retmu[0] = 2.0*propf3->A[0][0][0]*gridNew[index].c1 + (propf3->Beq[0][0] + propf3->dBbdT[0][0]*(Ti-pfmdat->Teq));
          }
          else if ( bulkphase == 1 ) { 
            retmu[0] = 2.0*propf3->A[1][0][0]*gridNew[index].c1 + (propf3->Beq[1][0] + propf3->dBbdT[1][0]*(Ti-pfmdat->Teq));
          }
          deltamu = retmu[0] - stgridO[4].mu[0];
          gridNew[index].mu[0] = retmu[0];
          
          
        }
        else if ( interface ) { 
          
          if ( !pfmdat->ISOTHERMAL ) {
            
            DELTAT = pfmvar->deltat * ( -pfmdat->TGRADIENT * pfmdat->velocity ); 
            
            c_tdt[0] = propf3->cmu[0][0][0] * ( stgridO[4].mu[0] - ( propf3->Beq[0][0] + propf3->dBbdT[0][0] * (Ti+DELTAT-pfmdat->Teq) ) );
            
            dcbdT_phase[0][0] = c_tdt[0] - stcscl[4].c1s;
            
            c_tdt[0] = propf3->cmu[1][0][0] * ( stgridO[4].mu[0] - ( propf3->Beq[1][0] + propf3->dBbdT[1][0] * (Ti+DELTAT-pfmdat->Teq) ) );
            
            dcbdT_phase[1][0] = c_tdt[0] - stcscl[4].c1l;
          
          }
        
          
          suma = ( stcscl[4].c1s - stcscl[4].c1l ) * (1.0) * ( gridNew[index].phi - stgridO[4].phi );
          
          sum_dcbdT = fhi[4] * dcbdT_phase[0][0] + (1.0-fhi[4]) * dcbdT_phase[1][0];

          deltac = pfmvar->deltat * divflux; 

          gridNew[index].c1 = stgridO[4].c1 + deltac;
          
          if ( pfmdat->ISOTHERMAL ) { 
            deltac = deltac - suma;
          }
          else { 
            deltac = deltac - suma  - sum_dcbdT;
          }
          
          hphi[4] = fhi[4];
          dcdmudenom = dc_dmu[4][0][0][0]*hphi[4] + dc_dmu[4][1][0][0]*(1.0-hphi[4]);
          
          gridNew[index].mu[0] = stgridO[4].mu[0] + deltac / dcdmudenom;
        }
        
    }
}

__kernel void SolverCWoatr_4(__global struct grid *gridOld, __global struct grid *gridNew, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf4 *propf4, __constant struct propmatf4spline *propf4spline, __constant struct propmatf4spline *propf4spline1) { 

    int i;
    int j;
    int k;
    int ii;    
    int nx;
    int ny;
    int nz;
    int index;
    int interface, phaseliquid, phasesolid, bulkphase, indij[5];

    double hphi[5], fhi[5], hphid[5], Tij[5];
    double term1lx;
    double term1sx;
    double term1ly;
    double term1sy;
    double term1l;
    double term1s;
    double c1dot;

    double gradmux1[1], gradmux2[1], gradmuy1[1], gradmuy2[1], mu[1]; 
    double Ti, ddgldx1ldx1l, ddgsdx1sdx1s, dcdmu[2][1][1], dc_dmu[5][2][1][1];
    double Da[5], Damidx1, Damidx2, Damidy1, Damidy2, divflux, deltamu; 
    double suma, deltac, dcdmudenom, y[2], retmu[1]; 
    double retdmuphase[1], DELTAT, cgs, cgl, dcbdT_phase[2][1], c_tdt[1], sum_dcbdT, muc1[1], cmu1[1]; 

    struct grid stgridO[9];
    struct csle stcscl[5];

    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    Ti = temp[j];
    
    if ( i != 0 && i != ny-1 && j != 0 && j != nx-1 ) {

        index = (i*nx + j);

        stgridO[0] = gridOld[(i-1)*nx + (j-1)];
        stgridO[1] = gridOld[(i-1)*nx + ( j )];
        stgridO[2] = gridOld[(i-1)*nx + (j+1)];
        stgridO[3] = gridOld[( i )*nx + (j-1)];
        stgridO[4] = gridOld[( i )*nx + ( j )];
        stgridO[5] = gridOld[( i )*nx + (j+1)];
        stgridO[6] = gridOld[(i+1)*nx + (j-1)];
        stgridO[7] = gridOld[(i+1)*nx + ( j )];
        stgridO[8] = gridOld[(i+1)*nx + (j+1)];

        stcscl[0] = cscl[(i-1)*nx + ( j )];
        stcscl[1] = cscl[( i )*nx + (j-1)];
        stcscl[4] = cscl[( i )*nx + ( j )];
        stcscl[2] = cscl[( i )*nx + (j+1)];
        stcscl[3] = cscl[(i+1)*nx + ( j )];

        fhi[0] = stgridO[1].phi;
        fhi[1] = stgridO[3].phi;
        fhi[4] = stgridO[4].phi;
        fhi[2] = stgridO[5].phi;
        fhi[3] = stgridO[7].phi;

        hphid[0] = fhi[0]; //stgridO[1].phi;
        hphid[1] = fhi[1]; //stgridO[3].phi;
        hphid[4] = fhi[4]; //stgridO[4].phi;
        hphid[2] = fhi[2]; //stgridO[5].phi;
        hphid[3] = fhi[3]; //stgridO[7].phi;

        Tij[0] = temp[( j )];
        Tij[1] = temp[(j-1)];
        Tij[4] = temp[( j )];
        Tij[2] = temp[(j+1)];
        Tij[3] = temp[( j )];
        
        indij[0] = ( j );
        indij[1] = (j-1);
        indij[4] = ( j );
        indij[2] = (j+1);
        indij[3] = ( j );
        

        for ( ii = 0; ii < 5; ii++ ) { 
          //printf("%d, %d, %d, %d, %le\n", tstep[0], i, j, ii, fhi[ii]);

          interface = 1;
          phaseliquid = 1; 
          phasesolid = 0; 
          
          if (fhi[ii] >= pfmdat->interfaceUplimit) { 
            bulkphase = 0;
            interface = 0;
          }
          else if (fhi[ii] <= pfmdat->interfaceDownlimit) { 
            bulkphase = 1;
            interface = 0;
          }

          if ( interface ) { 
             if ( pfmdat->ISOTHERMAL ) { 
               dc_dmu[ii][0][0][0] = propf4->cmu[0][0][0];
               dc_dmu[ii][1][0][0] = propf4->cmu[1][0][0];
             }
             else { 
               dc_dmu[ii][0][0][0] = 1 / (2.0 * propf4spline[indij[ii]].A[0][0][0]);
               dc_dmu[ii][1][0][0] = 1 / (2.0 * propf4spline[indij[ii]].A[1][0][0]);
             }

            Da[ii] = pfmdat->D11l * (1.0-hphid[ii]) * dc_dmu[ii][1][0][0] + pfmdat->D11s * (hphid[ii]) * dc_dmu[ii][0][0][0];
            
          }
          else {
           if ( pfmdat->ISOTHERMAL ) { 
            if ( bulkphase == 0 ) { 
              
              dc_dmu[ii][0][0][0] = propf4->cmu[bulkphase][0][0];
              
              
              Da[ii] = pfmdat->D11s * dc_dmu[ii][0][0][0];
              
            }
            else if ( bulkphase == 1 ) { 
              
              dc_dmu[ii][1][0][0] = propf4->cmu[bulkphase][0][0];
              
              Da[ii] = pfmdat->D11l * dc_dmu[ii][1][0][0];
              
            }
           }
           else {
              
            if ( bulkphase == 0 ) { 
              
              dc_dmu[ii][0][0][0] = 1.0 / (2.0 * propf4spline[indij[ii]].A[bulkphase][0][0]);
              
              
              Da[ii] = pfmdat->D11s * dc_dmu[ii][0][0][0];
              
            }
            else if ( bulkphase == 1 ) { 
              
              dc_dmu[ii][1][0][0] = 1.0 / (2.0 * propf4spline[indij[ii]].A[bulkphase][0][0]);
              
              Da[ii] = pfmdat->D11l * dc_dmu[ii][1][0][0];
              
            }
           }
          }
        }
        
        Damidx1 = ( Da[2] + Da[4] ) / 2.0;
        Damidx2 = ( Da[4] + Da[1] ) / 2.0;
        Damidy1 = ( Da[3] + Da[4] ) / 2.0;
        Damidy2 = ( Da[4] + Da[0] ) / 2.0;
        
        
        gradmux1[0] = ( stgridO[5].mu[0] - stgridO[4].mu[0] ) / pfmvar->deltax;
        gradmux2[0] = ( stgridO[4].mu[0] - stgridO[3].mu[0] ) / pfmvar->deltax;
        gradmuy1[0] = ( stgridO[7].mu[0] - stgridO[4].mu[0] ) / pfmvar->deltay;
        gradmuy2[0] = ( stgridO[4].mu[0] - stgridO[1].mu[0] ) / pfmvar->deltay;

        divflux = ( Damidx1*gradmux1[0] - Damidx2*gradmux2[0] ) / pfmvar->deltax + ( Damidy1*gradmuy1[0] - Damidy2*gradmuy2[0] ) / pfmvar->deltay;
        
        
        interface = 1;
        phaseliquid = 1; 
        phasesolid = 0;
        
        if (stgridO[4].phi >= pfmdat->interfaceUplimit) {
          bulkphase = 0;
          interface = 0;
        }
        else if (stgridO[4].phi <= pfmdat->interfaceDownlimit) {
          bulkphase = 1;
          interface = 0; 
        }

        if ( interface == 0) { 

          gridNew[index].c1 = stgridO[4].c1 + pfmvar->deltat * divflux;
          
         if ( pfmdat->ISOTHERMAL ) { 
          
          if ( bulkphase == 0 ) { 
            retmu[0] = 2.0*propf4->A[bulkphase][0][0]*gridNew[index].c1 + propf4->B[bulkphase][0];
          }
          else if ( bulkphase == 1 ) { 
            retmu[0] = 2.0*propf4->A[bulkphase][0][0]*gridNew[index].c1 + propf4->B[bulkphase][0];
          }
          
         }
         else { 
           
           if ( bulkphase == 0 ) { 
            retmu[0] = 2.0*propf4spline[j].A[bulkphase][0][0]*gridNew[index].c1 + propf4spline[j].B[bulkphase][0];
          }
          else if ( bulkphase == 1 ) { 
            retmu[0] = 2.0*propf4spline[j].A[bulkphase][0][0]*gridNew[index].c1 + propf4spline[j].B[bulkphase][0];
          }
           
         }
          
          
          deltamu = retmu[0] - stgridO[4].mu[0];
          gridNew[index].mu[0] = retmu[0];
          
        }
        else if ( interface ) {
          
          if ( pfmdat->ISOTHERMAL ) { 
               dcbdT_phase[0][0] = 0.0;
               dcbdT_phase[1][0] = 0.0;
          }
          else { 
            
            DELTAT = ( pfmvar->deltat ) * ( -pfmdat->TGRADIENT * pfmdat->velocity ); 
            
            muc1[0] = 2.0 * propf4spline1[j].A[0][0][0]; 
            
            cmu1[0] = 1.0 / muc1[0];
            
            c_tdt[0] = cmu1[0] * ( gridOld[index].mu[0] - propf4spline1[j].B[0][0] );
            
            dcbdT_phase[0][0] = c_tdt[0] - stcscl[4].c1s;
            
            muc1[0] = 2.0 * propf4spline1[j].A[1][0][0]; 
            
            cmu1[0] = 1.0 / muc1[0];
            
            c_tdt[0] = cmu1[0] * ( gridOld[index].mu[0] - propf4spline1[j].B[1][0] );
            
            dcbdT_phase[1][0] = c_tdt[0] - stcscl[4].c1l;
            
          }
        
          
          suma = ( stcscl[4].c1s - stcscl[4].c1l ) * (1.0) * ( gridNew[index].phi - stgridO[4].phi );
          
          sum_dcbdT = fhi[4] * dcbdT_phase[0][0] + (1.0-fhi[4]) * dcbdT_phase[1][0];

          deltac = pfmvar->deltat * divflux; 

          gridNew[index].c1 = stgridO[4].c1 + deltac;
          
          if ( pfmdat->ISOTHERMAL ) { 
            deltac = deltac - suma;
          }
          else { 
            deltac = deltac - suma  - sum_dcbdT;
          }
          
          hphi[4] = fhi[4];
          dcdmudenom = dc_dmu[4][0][0][0]*hphi[4] + dc_dmu[4][1][0][0]*(1.0-hphi[4]);
          
          gridNew[index].mu[0] = stgridO[4].mu[0] + deltac / dcdmudenom;
        }
        
    }
}


__kernel void SolverPhiIso_2(__global struct grid *gridOld, __global struct grid *gridNew, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf3 *propf3) {   
  
    int i;
    int j;
    int k;
    int ii;    
    int nx;
    int ny;
    int nz;
    int index;
    int interface, phaseliquid, phasesolid, bulkphase;

    double gphi;
    double dgphidphi;
    double dhphidphi; 
    double DetDL;
    double D11invL;
    double gl;
    double gs;
    double dgldc1l;
    double ddgldc1ldc1l;
    double zeta;
    double div_phi;
    double Mphi;
    double Mphi_aniso;
    double Mphi_FiniteMobility;
    double Mphi_FiniteMobilityAniso;
    double phidot_iso;
    
    double cliq1;
    double csol1;
    double cliq2;
    double csol2;

    double Bf3[npha][nsol], Cf3[npha];

    double Ti;

    double tmp1;
    
    double y[2], retGphase[1], retdmuphase[1], c_sol[1], c_liq[1];
    
    struct grid stgridO[9];

    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    if ( i != 0 && i != ny-1 && j != 0 && j != nx-1 ) {

        Ti = temp[j];

        index = (i)*nx + (j);

        stgridO[0] = gridOld[(i-1)*nx + (j-1)];
        stgridO[1] = gridOld[(i-1)*nx + ( j )];
        stgridO[2] = gridOld[(i-1)*nx + (j+1)];
        stgridO[3] = gridOld[( i )*nx + (j-1)];
        stgridO[4] = gridOld[( i )*nx + ( j )];
        stgridO[5] = gridOld[( i )*nx + (j+1)];
        stgridO[6] = gridOld[(i+1)*nx + (j-1)];
        stgridO[7] = gridOld[(i+1)*nx + ( j )];
        stgridO[8] = gridOld[(i+1)*nx + (j+1)];

        cliq1 = cscl[( i )*nx + ( j )].c1l;
        csol1 = cscl[( i )*nx + ( j )].c1s;

        div_phi = ( 4.0*(stgridO[1].phi+stgridO[3].phi+stgridO[7].phi+stgridO[5].phi ) + stgridO[0].phi+stgridO[6].phi+stgridO[8].phi+stgridO[2].phi - 20.0*stgridO[4].phi ) / ( 6.0*(pfmvar->deltax)*(pfmvar->deltax) );

          D11invL = 1.0/pfmdat->D11l;
          
          y[0] = csol1;
          y[1] = 1.0 - y[0];
          
          Ge(Ti, y, retGphase, pfmdat->thermophase[0]);
          gs = retGphase[0]/(pfmdat->Vm);
          
          y[0] = cliq1;
          y[1] = 1.0 - y[0];
          
          Ge(Ti, y, retGphase, pfmdat->thermophase[1]); 
          
          gl = retGphase[0]/(pfmdat->Vm);
          
          dgldc1l = stgridO[4].mu[0] / (pfmdat->Vm);
          
          //ddgldc1ldc1l = (ddGLIQdX1LdX1L(Ti, cliq1))/(pfmdat->Vm);
          
          dMudc(Ti, y, retdmuphase, pfmdat->thermophase[1]);
          ddgldc1ldc1l = retdmuphase[0] / pfmdat->Vm;
          
          zeta = ((cliq1-csol1)*(cliq1-csol1)*ddgldc1ldc1l)*D11invL;
          Mphi = (pfmvar->w)/(3.0*(pfmvar->ee)*(pfmdat->a2)*zeta);
          
          gphi = stgridO[4].phi*stgridO[4].phi*(1.0-stgridO[4].phi)*(1.0-stgridO[4].phi);
          dgphidphi = 2.0*stgridO[4].phi*(1.0-stgridO[4].phi)*(1.0-2.0*stgridO[4].phi);
          dhphidphi = 30.0*gphi;
          

        if ( fabs(div_phi) > 0.0 )  { 
          
          phidot_iso =  Mphi*((pfmvar->ee)*div_phi - (pfmvar->w)*dgphidphi - dhphidphi*(gs-gl-(csol1-cliq1)*dgldc1l));

          //printf("%d, %d, %d, %le\n", tstep[0], i, j, phidot_iso*pfmvar->deltat);
          //printf("%d, %d, %d, %le, %le, %le, %le\n", tstep[0], i, j, csol1, gs*pfmdat->Vm, cliq1, gl*pfmdat->Vm);
          
          gridNew[ ( ( i )*nx + ( j ) ) ].phi = stgridO[4].phi + (pfmvar->deltat)*phidot_iso;

        }
        // else { 
        //   phidot_iso = 0.0;
        //   gridNew[ ( ( i )*nx + ( j ) ) ].phi = stgridO[4].phi + (pfmvar->deltat)*phidot_iso;
        //   //printf("%d, %d, %d, %le\n", tstep[0], i, j, phidot_iso*pfmvar->deltat);
        //   //printf("%d, %d, %d, %le, %le\n", tstep[0], i, j, phidot_iso, phidot_iso);
        // }
        
        //printf(">> %d, %d, %d, %le, %le, %le, %le, %le, %le, %le, %le, %le\n", tstep[0], i, j, Mphi, div_phi, dgphidphi, dhphidphi, gs, gl,csol1, cliq1,dgldc1l);

        
        //printf("%le,\t%le,\t%le,\t%le,\t%le,\t%le,\t%le,\t%le\n", Mphi, pfmvar->w, pfmvar->ee, pfmdat->a2, zeta, cliq1, csol1, ddgldc1ldc1l);
    }
}

__kernel void SolverPhiIso_3(__global struct grid *gridOld, __global struct grid *gridNew, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf3 *propf3) {   
  
    int i;
    int j;
    int k;
    int ii;    
    int nx;
    int ny;
    int nz;
    int index;

    double gphi;
    double dgphidphi;
    double dhphidphi; 
    double DetDL;
    double D11invL;
    double gl;
    double gs;
    double dgldc1l;
    double ddgldc1ldc1l;
    double zeta;
    double div_phi;
    double Mphi;
    double Mphi_aniso;
    double Mphi_FiniteMobility;
    double Mphi_FiniteMobilityAniso;
    double phidot_iso;
    
    double cliq1;
    double csol1;
    double cliq2;
    double csol2;

    double Bf3[npha][nsol], Cf3[npha];
    
    double y[2], retGphase[1], retdmuphase[1], c_sol[1], c_liq[1];

    double Ti;

    double tmp1;
    
    struct grid stgridO[9];

    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    if ( i != 0 && i != ny-1 && j != 0 && j != nx-1 ) {

        Ti = temp[j];

        index = (i)*nx + (j);

        stgridO[0] = gridOld[(i-1)*nx + (j-1)];
        stgridO[1] = gridOld[(i-1)*nx + ( j )];
        stgridO[2] = gridOld[(i-1)*nx + (j+1)];
        stgridO[3] = gridOld[( i )*nx + (j-1)];
        stgridO[4] = gridOld[( i )*nx + ( j )];
        stgridO[5] = gridOld[( i )*nx + (j+1)];
        stgridO[6] = gridOld[(i+1)*nx + (j-1)];
        stgridO[7] = gridOld[(i+1)*nx + ( j )];
        stgridO[8] = gridOld[(i+1)*nx + (j+1)];

        cliq1 = cscl[( i )*nx + ( j )].c1l;
        csol1 = cscl[( i )*nx + ( j )].c1s;

        div_phi = ( 4.0*(stgridO[1].phi+stgridO[3].phi+stgridO[7].phi+stgridO[5].phi ) + stgridO[0].phi+stgridO[6].phi+stgridO[8].phi+stgridO[2].phi - 20.0*stgridO[4].phi ) / ( 6.0*(pfmvar->deltax)*(pfmvar->deltax) );

        D11invL = 1.0/pfmdat->D11l;
        
        if ( !pfmdat->ISOTHERMAL ) { 
          Bf3[0][0] = propf3->Beq[0][0] + propf3->dBbdT[0][0]*(Ti-pfmdat->Teq);
          Bf3[1][0] = propf3->Beq[1][0] + propf3->dBbdT[1][0]*(Ti-pfmdat->Teq);
        }
        else {
          Bf3[0][0] = propf3->B[0][0];
          Bf3[1][0] = propf3->B[1][0];
        }
        
        if ( !pfmdat->ISOTHERMAL ) { 
          c_liq[0] = propf3->ceq[0][1][0] - propf3->DELTA_C[0][0]*(pfmdat->Teq-Ti) / propf3->DELTA_T[0][1];
          
          c_sol[0] = propf3->ceq[0][0][0] - propf3->DELTA_C[0][0]*(pfmdat->Teq-Ti) / propf3->DELTA_T[0][0]; 
          
          Cf3[0] = propf3->A[0][0][0] * c_sol[0]  * c_sol[0] - propf3->A[1][0][0] * c_liq[0] * c_liq[0]; 
          
          Cf3[1] = propf3->C[1];
        }
        else { 
          Cf3[0] = propf3->C[0]; 
          Cf3[1] = propf3->C[1];
        }
          
        
        
        gs = ( propf3->A[0][0][0]*csol1*csol1 + Bf3[0][0]*csol1 + Cf3[0] ) / pfmdat->Vm; 
        gl = ( propf3->A[1][0][0]*cliq1*cliq1 + Bf3[1][0]*cliq1 + Cf3[1] ) / pfmdat->Vm; 
        //printf("%d, %d, %d, %le, %le, %le, %le, %le, %le\n", tstep[0], i, j, gs, propf3->A[0][0][0],  propf3->Beq[0][0], Bf3[0][0], csol1,  Cf3[0]);
        
         dgldc1l = stgridO[4].mu[0] / (pfmdat->Vm);
         
        ddgldc1ldc1l = ( 2.0*propf3->A[1][0][0] )/(pfmdat->Vm);
        
        zeta = ((cliq1-csol1)*(cliq1-csol1)*ddgldc1ldc1l)*D11invL;

        gphi = stgridO[4].phi*stgridO[4].phi*(1.0-stgridO[4].phi)*(1.0-stgridO[4].phi);
        dgphidphi = 2.0*stgridO[4].phi*(1.0-stgridO[4].phi)*(1.0-2.0*stgridO[4].phi);
        dhphidphi = 30.0*gphi;

        Mphi = (pfmvar->w)/(3.0*(pfmvar->ee)*(pfmdat->a2)*zeta);

        if ( fabs(div_phi) > 0.0 )  { 

        phidot_iso =  Mphi*((pfmvar->ee)*div_phi - (pfmvar->w)*dgphidphi - dhphidphi*(gs-gl-(csol1-cliq1)*dgldc1l));

        gridNew[ ( ( i )*nx + ( j ) ) ].phi = stgridO[4].phi + (pfmvar->deltat)*phidot_iso;
        
        //printf("%d, %d, %d, %le, %le, %le, %le, %le, %le, %le\n", tstep[0], i, j, stgridO[4].phi, pfmvar->deltat, phidot_iso/Mphi, Mphi, gridNew[ ( ( i )*nx + ( j ) ) ].phi, stgridO[4].c1, stgridO[4].mu[0]);
        
        //printf("%d, %d, %d, %le, %le, %le, %le, %le, %le, %le, %le, %le, %le\n", tstep[0], i, j, pfmvar->ee, div_phi, pfmvar->w, dgphidphi, dhphidphi, gs, gl, csol1, cliq1, dgldc1l);

        }

    }
}

__kernel void SolverPhiIso_4(__global struct grid *gridOld, __global struct grid *gridNew, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf4 *propf4, __constant struct propmatf4spline *propf4spline) {   
  
    int i;
    int j;
    int k;
    int ii;    
    int nx;
    int ny;
    int nz;
    int index;

    double gphi;
    double dgphidphi;
    double dhphidphi; 
    double DetDL;
    double D11invL;
    double gl;
    double gs;
    double dgldc1l;
    double ddgldc1ldc1l;
    double zeta;
    double div_phi;
    double Mphi;
    double Mphi_aniso;
    double Mphi_FiniteMobility;
    double Mphi_FiniteMobilityAniso;
    double phidot_iso;
    
    double cliq1;
    double csol1;
    double cliq2;
    double csol2;

    double Bf3[npha][nsol], Cf3[npha];
    
    double y[2], retGphase[1], retdmuphase[1], c_sol[1], c_liq[1];

    double Ti;

    double tmp1, tmp0;
    
    struct grid stgridO[9];

    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    if ( i != 0 && i != ny-1 && j != 0 && j != nx-1 ) {

        Ti = temp[j];

        index = (i)*nx + (j);

        stgridO[0] = gridOld[(i-1)*nx + (j-1)];
        stgridO[1] = gridOld[(i-1)*nx + ( j )];
        stgridO[2] = gridOld[(i-1)*nx + (j+1)];
        stgridO[3] = gridOld[( i )*nx + (j-1)];
        stgridO[4] = gridOld[( i )*nx + ( j )];
        stgridO[5] = gridOld[( i )*nx + (j+1)];
        stgridO[6] = gridOld[(i+1)*nx + (j-1)];
        stgridO[7] = gridOld[(i+1)*nx + ( j )];
        stgridO[8] = gridOld[(i+1)*nx + (j+1)];

        cliq1 = cscl[( i )*nx + ( j )].c1l;
        csol1 = cscl[( i )*nx + ( j )].c1s;

        div_phi = ( 4.0*(stgridO[1].phi+stgridO[3].phi+stgridO[7].phi+stgridO[5].phi ) + stgridO[0].phi+stgridO[6].phi+stgridO[8].phi+stgridO[2].phi - 20.0*stgridO[4].phi ) / ( 6.0*(pfmvar->deltax)*(pfmvar->deltax) );

        D11invL = 1.0/pfmdat->D11l;
        
        if ( pfmdat->ISOTHERMAL ) { 
          
          gs = ( propf4->A[0][0][0]*csol1*csol1 + propf4->B[0][0]*csol1 + propf4->C[0] ) / pfmdat->Vm;
          
          gl = ( propf4->A[1][0][0]*cliq1*cliq1 + propf4->B[1][0]*cliq1 + propf4->C[1] ) / pfmdat->Vm;
          
        }
        else { 
          tmp0 = propf4spline[j].A[0][0][0]*csol1*csol1;
          
          tmp1 = propf4spline[j].B[0][0]*csol1; 
        
          gs = ( tmp0 + tmp1 + propf4spline[j].C[0] ) / pfmdat->Vm;
        
          tmp0 = propf4spline[j].A[1][0][0]*cliq1*cliq1;
          
          tmp1 = propf4spline[j].B[1][0]*cliq1; 
        
          gl = ( tmp0 + tmp1 + propf4spline[j].C[1] ) / pfmdat->Vm;
        }
        
        if ( pfmdat->ISOTHERMAL ) { 
          ddgldc1ldc1l = 2.0*propf4->A[1][0][0] / pfmdat->Vm;
        }
        else { 
          ddgldc1ldc1l = 2.0*propf4spline[j].A[1][0][0] / pfmdat->Vm;
        }
        
        dgldc1l = gridOld[index].mu[0] / pfmdat->Vm;
        
        zeta = ((cliq1-csol1)*(cliq1-csol1)*ddgldc1ldc1l)*D11invL;

        gphi = stgridO[4].phi*stgridO[4].phi*(1.0-stgridO[4].phi)*(1.0-stgridO[4].phi);
        dgphidphi = 2.0*stgridO[4].phi*(1.0-stgridO[4].phi)*(1.0-2.0*stgridO[4].phi);
        dhphidphi = 30.0*gphi;

        Mphi = (pfmvar->w)/(3.0*(pfmvar->ee)*(pfmdat->a2)*zeta);

        if ( fabs(div_phi) > 0.0 )  { 

        phidot_iso =  Mphi*((pfmvar->ee)*div_phi - (pfmvar->w)*dgphidphi - dhphidphi*(gs-gl-(csol1-cliq1)*dgldc1l));

        gridNew[ ( ( i )*nx + ( j ) ) ].phi = stgridO[4].phi + (pfmvar->deltat)*phidot_iso;

        }

    }
}

__kernel void SolverPhi_2(__global struct grid *gridOld, __global struct grid *gridNew, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf3 *propf3) {   
  
    int i;
    int j;
    int k;
    int ii, i1, ii1, jj1;  
    int nx;
    int ny;
    int nz;
    int index;

    double gphi;
    double dgphidphi;
    double dhphidphi; 
    double DetDL;
    double D11invL;
    double gl;
    double gs;
    double dgldc1l;
    double ddgldc1ldc1l;
    double zeta;
    double gradx_phi[3];
    double grady_phi[3];
    double phix[5];
    double phiy[5];
    double div_phi;
    double angleig;
    double Nnx;
    double Nny;
    double neta[5];
    double etam[5];
    double B[5];
    double neta_x;
    double neta_y;
    double alp1;
    double eep_x;
    double eep_y;
    double alp2;
    double alp3;
    double alp;
    double Mphi;
    double Mphi_aniso;
    double Mphi_FiniteMobility;
    double Mphi_FiniteMobilityAniso;
    double phidot_aniso;
    
    double cliq1;
    double csol1;
    double cliq2;
    double csol2;

    double Bf3[npha][nsol];
    
     double y[2], retGphase[1], retdmuphase[1], c_sol[1], c_liq[1];

    double Ti;

    double tmp1;

    double phiz[5], phia[5][3], tmp0, phir[5][3], phixr[5], phiyr[5], phizr[5], Bx[5], By[5], Bz[5], modgradphi1, modgradphi4, Bxir[5], Byir[5], Bzir[5], Ba[5][3], Bir[5][3];
    
    struct grid stgridO[9];

    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    if ( i != 0 && i != ny-1 && j != 0 && j != nx-1 ) {

        Ti = temp[j];

        index = (i)*nx + (j);

        stgridO[0] = gridOld[(i-1)*nx + (j-1)];
        stgridO[1] = gridOld[(i-1)*nx + ( j )];
        stgridO[2] = gridOld[(i-1)*nx + (j+1)];
        stgridO[3] = gridOld[( i )*nx + (j-1)];
        stgridO[4] = gridOld[( i )*nx + ( j )];
        stgridO[5] = gridOld[( i )*nx + (j+1)];
        stgridO[6] = gridOld[(i+1)*nx + (j-1)];
        stgridO[7] = gridOld[(i+1)*nx + ( j )];
        stgridO[8] = gridOld[(i+1)*nx + (j+1)];

        cliq1 = cscl[( i )*nx + ( j )].c1l;
        csol1 = cscl[( i )*nx + ( j )].c1s;

        gradx_phi[0] = ( stgridO[5].phi - stgridO[3].phi ) / ( 2.0*(pfmvar->deltax) );
        gradx_phi[1] = ( stgridO[8].phi - stgridO[6].phi ) / ( 2.0*(pfmvar->deltax) );
        gradx_phi[2] = ( stgridO[2].phi - stgridO[0].phi ) / ( 2.0*(pfmvar->deltax) );

        grady_phi[0] = ( stgridO[7].phi - stgridO[1].phi ) / ( 2.0*(pfmvar->deltay) );
        grady_phi[1] = ( stgridO[8].phi - stgridO[2].phi ) / ( 2.0*(pfmvar->deltay) );
        grady_phi[2] = ( stgridO[6].phi - stgridO[0].phi ) / ( 2.0*(pfmvar->deltay) );

        phix[0] = gradx_phi[0];
        phix[1] = ( stgridO[5].phi - stgridO[4].phi ) / ( (pfmvar->deltax) );
        phix[2] = ( stgridO[4].phi - stgridO[3].phi ) / ( (pfmvar->deltax) );
        phix[3] = ( gradx_phi[0] + gradx_phi[1] ) / ( 2.0 );
        phix[4] = ( gradx_phi[0] + gradx_phi[2] ) / ( 2.0 );

        phiy[0] = grady_phi[0];
        phiy[1] = ( grady_phi[0] + grady_phi[1] ) / ( 2.0 );
        phiy[2] = ( grady_phi[0] + grady_phi[2] ) / ( 2.0 );
        phiy[3] = ( stgridO[7].phi - stgridO[4].phi ) / ( (pfmvar->deltay) );
        phiy[4] = ( stgridO[4].phi - stgridO[1].phi ) / ( (pfmvar->deltay) );

        div_phi = ( 4.0*(stgridO[1].phi+stgridO[3].phi+stgridO[7].phi+stgridO[5].phi ) + stgridO[0].phi+stgridO[6].phi+stgridO[8].phi+stgridO[2].phi - 20.0*stgridO[4].phi ) / ( 6.0*(pfmvar->deltax)*(pfmvar->deltax) );
        
          
        D11invL = 1.0/pfmdat->D11l;

        y[0] = csol1;
        y[1] = 1.0 - y[0];
        
        Ge(Ti, y, retGphase, pfmdat->thermophase[0]);
        gs = retGphase[0]/(pfmdat->Vm);
        
        y[0] = cliq1;
        y[1] = 1.0 - y[0];
        
        Ge(Ti, y, retGphase, pfmdat->thermophase[1]); 
        
        gl = retGphase[0]/(pfmdat->Vm);
        
        dgldc1l = stgridO[4].mu[0] / (pfmdat->Vm);
        
        //ddgldc1ldc1l = (ddGLIQdX1LdX1L(Ti, cliq1))/(pfmdat->Vm);
        
        dMudc(Ti, y, retdmuphase, pfmdat->thermophase[1]);
        ddgldc1ldc1l = retdmuphase[0] / pfmdat->Vm;
          
        zeta = ((cliq1-csol1)*(cliq1-csol1)*ddgldc1ldc1l)*D11invL;
        Mphi = (pfmvar->w)/(3.0*(pfmvar->ee)*(pfmdat->a2)*zeta);
        
        // ddgldc1ldc1l = (ddGLIQdX1LdX1L(Ti, pfmdat->c1l_Initial))/(pfmdat->Vm);
        // zeta = ((pfmdat->c1l_Initial-pfmdat->c1s_Initial)*(pfmdat->c1l_Initial-pfmdat->c1s_Initial)*ddgldc1ldc1l)*D11invL;
        // Mphi = (pfmvar->w)/(3.0*(pfmvar->ee)*(pfmdat->a2)*zeta);
        
        //printf("%d, %d, %d, %le %le\n", tstep[0], i, j, zeta, Mphi)

        gphi = stgridO[4].phi*stgridO[4].phi*(1.0-stgridO[4].phi)*(1.0-stgridO[4].phi);
        dgphidphi = 2.0*stgridO[4].phi*(1.0-stgridO[4].phi)*(1.0-2.0*stgridO[4].phi);
        dhphidphi = 30.0*gphi;
        
        // gphi = stgridO[4].phi*stgridO[4].phi*(1.0-stgridO[4].phi)*(1.0-stgridO[4].phi);
        // dgphidphi = 2.0*stgridO[4].phi*(1.0-stgridO[4].phi)*(1.0-2.0*stgridO[4].phi);
        // dhphidphi = 6.0*stgridO[4].phi*(1.0-stgridO[4].phi);
        
        phiz[0] = 0.0; 
        phiz[1] = 0.0; 
        phiz[2] = 0.0; 
        phiz[3] = 0.0; 
        phiz[4] = 0.0; 
        
        for ( i1 = 0; i1 < 5; i1++ ) { 
          phia[i1][0] = phix[i1];
          phia[i1][1] = phiy[i1];
          phia[i1][2] = phiz[i1];
        }
        
        // Rotation Matrix multiplication with gradients vector
        for ( i1 = 0; i1 < 5; i1++ ) { 
          for ( ii1 = 0; ii1 < 3; ii1++ ) { 
            tmp0 = 0.0;
            for ( jj1 = 0; jj1 < 3; jj1++ ) { 
              tmp0 += pfmdat->Rotation_matrix[0][1][ii1][jj1]*phia[i1][jj1];
              // if ( i == 1 && j == 1 && i1 ==0 && ip1 == 0 && ip2 == 1 ) { 
              //   printf("%le\t", pfmdat->Rotation_matrix[ii1][jj1]);
              // }
            }
            // if ( i == 1 && j == 1 && i1 ==0 && ip1 == 0 && ip2 == 1 ) { 
            //   printf("\n");
            // }
            phir[i1][ii1] = tmp0;
          }
        }
        
        for ( i1 = 0; i1 < 5; i1++ ) { 
          phixr[i1] = phir[i1][0];
          phiyr[i1] = phir[i1][1];
          phizr[i1] = phir[i1][2];
        }
        
        for (ii = 0; ii < 5; ii++) { 
          neta[ii] = 0.0;
          Bx[ii] = 0.0; 
          By[ii] = 0.0; 
          Bz[ii] = 0.0; 
          B[ii] = 0.0; 
        }
        
        for (ii = 0; ii <5; ii++) { 
          
          modgradphi1 = sqrt( phixr[ii]*phixr[ii] + phiyr[ii]*phiyr[ii] );
          modgradphi4 = modgradphi1*modgradphi1*modgradphi1*modgradphi1;
          
          if ( modgradphi4 != 0.0 ) { 
            
            Nnx = phixr[ii] / modgradphi1;
            Nny = phiyr[ii] / modgradphi1;
            
            tmp1 = ( Nnx*Nnx*Nnx*Nnx + Nny*Nny*Nny*Nny );
            
            //tmp2 = ( phixr[ii]*phixr[ii] + phiyr[ii]*phiyr[ii] ); 
            
            neta[ii] = ( 1.0 - 3.0*(pfmdat->epsc) ) * ( 1.0 + ( (4.0*(pfmdat->epsc))/(1.0-3.0*(pfmdat->epsc)) )*tmp1 );
            
            //Bx[ii] = neta[ii] *  16.0 * pfmdat->epsc * ( Nnx*Nnx - tmp1 ) * tmp2;
            //By[ii] = neta[ii] *  16.0 * pfmdat->epsc * ( Nny*Nny - tmp1 ) * tmp2;
            
            Bx[ii] = 16.0 * pfmdat->epsc * ( Nnx*Nnx - tmp1 ) * phixr[ii];
            By[ii] = 16.0 * pfmdat->epsc * ( Nny*Nny - tmp1 ) * phiyr[ii];
            
            //etam[ii] = ( 1.0 - 3.0*(pfmdat->epsm) ) * ( 1.0 + ( (4.0*(pfmdat->epsm))/(1.0-3.0*(pfmdat->epsm)) )*tmp1 );
             
          }
          else {
            
            neta[ii] = ( 1.0 - 3.0*(pfmdat->epsc) );
            
            Bx[ii] = 0.0;
            By[ii] = 0.0;
            
            //etam[ii] = ( 1.0 - 3.0*(pfmdat->epsm) );

          }
        }
        
        for ( i1 = 0; i1 < 5; i1++ ) { 
          Ba[i1][0] = Bx[i1];
          Ba[i1][1] = By[i1];
          Ba[i1][2] = Bz[i1];
        }
        
        // Inverse Rotation Matrix multiplication with depsilon/dphix, depsilon/dphiy
        for ( i1 = 0; i1 < 5; i1++ ) { 
          for ( ii1 = 0; ii1 < 3; ii1++ ) { 
            tmp0 = 0.0;
            for ( jj1 = 0; jj1 < 3; jj1++ ) { 
              tmp0 += pfmdat->Inv_Rotation_matrix[0][1][ii1][jj1]*Ba[i1][jj1];
            }
            Bir[i1][ii1] = tmp0;
          }
        }
        
        for ( i1 = 0; i1 < 5; i1++ ) { 
          Bxir[i1] = Bir[i1][0];
          Byir[i1] = Bir[i1][1];
          Bzir[i1] = Bir[i1][2];
        }
        
        neta_x = ( neta[1] - neta[2] ) / ( (pfmvar->deltax) ); 
        neta_y = ( neta[3] - neta[4] ) / ( (pfmvar->deltay) );

        alp1 = neta[0]*neta[0]*div_phi + 2.0*neta[0] * ( neta_x * phix[0] + neta_y * phiy[0] );
        
        alp2 = ( neta[1] * Bxir[1] - neta[2] * Bxir[2] ) / ( (pfmvar->deltax) );
        alp3 = ( neta[3] * Byir[3] - neta[4] * Byir[4] ) / ( (pfmvar->deltay) );

        alp = alp1 + alp2 + alp3;
        
        if ( fabs(div_phi) > 0.0 )  { 
        phidot_aniso =  Mphi*((pfmvar->ee)*alp - (pfmvar->w)*dgphidphi - dhphidphi*(gs-gl-(csol1-cliq1)*dgldc1l));
        gridNew[ ( ( i )*nx + ( j ) ) ].phi = stgridO[4].phi + (pfmvar->deltat)*phidot_aniso;

      }
      // else { 
      //   phidot_aniso = 0.0;
      //   gridNew[ ( ( i )*nx + ( j ) ) ].phi = stgridO[4].phi + (pfmvar->deltat)*phidot_aniso;
      //   //printf("%d, %d, %d, %le\n", tstep[0], i, j, phidot_iso*pfmvar->deltat);
      //   //printf("%d, %d, %d, %le, %le\n", tstep[0], i, j, phidot_iso, phidot_iso);
      // }
      
      
      //printf("%le,\t%le,\t%le,\t%le,\t%le,\t%le,\t%le,\t%le\n", Mphi, pfmvar->w, pfmvar->ee, pfmdat->a2, zeta, cliq1, csol1, ddgldc1ldc1l);
		

    }
}

__kernel void SolverPhi_3(__global struct grid *gridOld, __global struct grid *gridNew, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf3 *propf3) {   
  
    int i;
    int j;
    int k;
    int ii, i1, ii1, jj1;  
    int nx;
    int ny;
    int nz;
    int index;

    double gphi;
    double dgphidphi;
    double dhphidphi; 
    double DetDL;
    double D11invL;
    double gl;
    double gs;
    double dgldc1l;
    double ddgldc1ldc1l;
    double zeta;
    double gradx_phi[3];
    double grady_phi[3];
    double phix[5];
    double phiy[5];
    double div_phi;
    double angleig;
    double Nnx;
    double Nny;
    double neta[5];
    double etam[5];
    double B[5];
    double neta_x;
    double neta_y;
    double alp1;
    double eep_x;
    double eep_y;
    double alp2;
    double alp3;
    double alp;
    double Mphi;
    double Mphi_aniso;
    double Mphi_FiniteMobility;
    double Mphi_FiniteMobilityAniso;
    double phidot_aniso;
    
    double cliq1;
    double csol1;
    double cliq2;
    double csol2;

    double Bf3[npha][nsol], Cf3[npha];
    
     double y[2], retGphase[1], retdmuphase[1], c_sol[1], c_liq[1];

    double Ti;

    double tmp1;

    double phiz[5], phia[5][3], tmp0, phir[5][3], phixr[5], phiyr[5], phizr[5], Bx[5], By[5], Bz[5], modgradphi1, modgradphi4, Bxir[5], Byir[5], Bzir[5], Ba[5][3], Bir[5][3];
    
    struct grid stgridO[9];

    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    if ( i != 0 && i != ny-1 && j != 0 && j != nx-1 ) {

        Ti = temp[j];

        index = (i)*nx + (j);

        stgridO[0] = gridOld[(i-1)*nx + (j-1)];
        stgridO[1] = gridOld[(i-1)*nx + ( j )];
        stgridO[2] = gridOld[(i-1)*nx + (j+1)];
        stgridO[3] = gridOld[( i )*nx + (j-1)];
        stgridO[4] = gridOld[( i )*nx + ( j )];
        stgridO[5] = gridOld[( i )*nx + (j+1)];
        stgridO[6] = gridOld[(i+1)*nx + (j-1)];
        stgridO[7] = gridOld[(i+1)*nx + ( j )];
        stgridO[8] = gridOld[(i+1)*nx + (j+1)];

        cliq1 = cscl[( i )*nx + ( j )].c1l;
        csol1 = cscl[( i )*nx + ( j )].c1s;

        
        div_phi = ( 4.0*(stgridO[1].phi+stgridO[3].phi+stgridO[7].phi+stgridO[5].phi ) + stgridO[0].phi+stgridO[6].phi+stgridO[8].phi+stgridO[2].phi - 20.0*stgridO[4].phi ) / ( 6.0*(pfmvar->deltax)*(pfmvar->deltax) );

        gradx_phi[0] = ( stgridO[5].phi - stgridO[3].phi ) / ( 2.0*(pfmvar->deltax) );
        gradx_phi[1] = ( stgridO[8].phi - stgridO[6].phi ) / ( 2.0*(pfmvar->deltax) );
        gradx_phi[2] = ( stgridO[2].phi - stgridO[0].phi ) / ( 2.0*(pfmvar->deltax) );

        grady_phi[0] = ( stgridO[7].phi - stgridO[1].phi ) / ( 2.0*(pfmvar->deltay) );
        grady_phi[1] = ( stgridO[8].phi - stgridO[2].phi ) / ( 2.0*(pfmvar->deltay) );
        grady_phi[2] = ( stgridO[6].phi - stgridO[0].phi ) / ( 2.0*(pfmvar->deltay) );

        phix[0] = gradx_phi[0];
        phix[1] = ( stgridO[5].phi - stgridO[4].phi ) / ( (pfmvar->deltax) );
        phix[2] = ( stgridO[4].phi - stgridO[3].phi ) / ( (pfmvar->deltax) );
        phix[3] = ( gradx_phi[0] + gradx_phi[1] ) / ( 2.0 );
        phix[4] = ( gradx_phi[0] + gradx_phi[2] ) / ( 2.0 );

        phiy[0] = grady_phi[0];
        phiy[1] = ( grady_phi[0] + grady_phi[1] ) / ( 2.0 );
        phiy[2] = ( grady_phi[0] + grady_phi[2] ) / ( 2.0 );
        phiy[3] = ( stgridO[7].phi - stgridO[4].phi ) / ( (pfmvar->deltay) );
        phiy[4] = ( stgridO[4].phi - stgridO[1].phi ) / ( (pfmvar->deltay) );

        

        D11invL = 1.0/pfmdat->D11l;
        
        if ( !pfmdat->ISOTHERMAL ) { 
          Bf3[0][0] = propf3->Beq[0][0] + propf3->dBbdT[0][0]*(Ti-pfmdat->Teq);
          Bf3[1][0] = propf3->Beq[1][0] + propf3->dBbdT[1][0]*(Ti-pfmdat->Teq);
        }
        else {
          Bf3[0][0] = propf3->B[0][0];
          Bf3[1][0] = propf3->B[1][0];
        }
        
        if ( !pfmdat->ISOTHERMAL ) { 
          c_liq[0] = propf3->ceq[0][1][0] - propf3->DELTA_C[0][0]*(pfmdat->Teq-Ti) / propf3->DELTA_T[0][1];
          
          c_sol[0] = propf3->ceq[0][0][0] - propf3->DELTA_C[0][0]*(pfmdat->Teq-Ti) / propf3->DELTA_T[0][0]; 
          
          Cf3[0] = propf3->A[0][0][0] * c_sol[0] * c_sol[0] - propf3->A[1][0][0] * c_liq[0] * c_liq[0]; 
          
          Cf3[1] = propf3->C[1];
        }
        else { 
          Cf3[0] = propf3->C[0]; 
          Cf3[1] = propf3->C[1];
        }
        
        gs = ( propf3->A[0][0][0]*csol1*csol1 + Bf3[0][0]*csol1 + Cf3[0] ) / pfmdat->Vm; 
        gl = ( propf3->A[1][0][0]*cliq1*cliq1 + Bf3[1][0]*cliq1 + Cf3[1] ) / pfmdat->Vm; 
        
         dgldc1l = stgridO[4].mu[0] / (pfmdat->Vm);
         
        ddgldc1ldc1l = ( 2.0*propf3->A[1][0][0] )/(pfmdat->Vm);
        
        zeta = ((cliq1-csol1)*(cliq1-csol1)*ddgldc1ldc1l*D11invL);

        gphi = stgridO[4].phi*stgridO[4].phi*(1.0-stgridO[4].phi)*(1.0-stgridO[4].phi);
        dgphidphi = 2.0*stgridO[4].phi*(1.0-stgridO[4].phi)*(1.0-2.0*stgridO[4].phi);
        dhphidphi = 30.0*gphi;


        phiz[0] = 0.0; 
        phiz[1] = 0.0; 
        phiz[2] = 0.0; 
        phiz[3] = 0.0; 
        phiz[4] = 0.0; 
        
        for ( i1 = 0; i1 < 5; i1++ ) { 
          phia[i1][0] = phix[i1];
          phia[i1][1] = phiy[i1];
          phia[i1][2] = phiz[i1];
        }
        
        // Rotation Matrix multiplication with gradients vector
        for ( i1 = 0; i1 < 5; i1++ ) { 
          for ( ii1 = 0; ii1 < 3; ii1++ ) { 
            tmp0 = 0.0;
            for ( jj1 = 0; jj1 < 3; jj1++ ) { 
              tmp0 += pfmdat->Rotation_matrix[0][1][ii1][jj1]*phia[i1][jj1];
              // if ( i == 1 && j == 1 && i1 ==0 && ip1 == 0 && ip2 == 1 ) { 
              //   printf("%le\t", pfmdat->Rotation_matrix[ii1][jj1]);
              // }
            }
            // if ( i == 1 && j == 1 && i1 ==0 && ip1 == 0 && ip2 == 1 ) { 
            //   printf("\n");
            // }
            phir[i1][ii1] = tmp0;
          }
        }
        
        for ( i1 = 0; i1 < 5; i1++ ) { 
          phixr[i1] = phir[i1][0];
          phiyr[i1] = phir[i1][1];
          phizr[i1] = phir[i1][2];
        }
        
        for (ii = 0; ii < 5; ii++) { 
          neta[ii] = 0.0;
          Bx[ii] = 0.0; 
          By[ii] = 0.0; 
          Bz[ii] = 0.0; 
          B[ii] = 0.0; 
        }
        
        for (ii = 0; ii <5; ii++) { 
          
          modgradphi1 = sqrt( phixr[ii]*phixr[ii] + phiyr[ii]*phiyr[ii] );
          modgradphi4 = modgradphi1*modgradphi1*modgradphi1*modgradphi1;
          
          if ( modgradphi4 != 0.0 ) { 
            
            Nnx = phixr[ii] / modgradphi1;
            Nny = phiyr[ii] / modgradphi1;
            
            tmp1 = ( Nnx*Nnx*Nnx*Nnx + Nny*Nny*Nny*Nny );
            
            //tmp2 = ( phixr[ii]*phixr[ii] + phiyr[ii]*phiyr[ii] ); 
            
            neta[ii] = ( 1.0 - 3.0*(pfmdat->epsc) ) * ( 1.0 + ( (4.0*(pfmdat->epsc))/(1.0-3.0*(pfmdat->epsc)) )*tmp1 );
            
            //Bx[ii] = neta[ii] *  16.0 * pfmdat->epsc * ( Nnx*Nnx - tmp1 ) * tmp2;
            //By[ii] = neta[ii] *  16.0 * pfmdat->epsc * ( Nny*Nny - tmp1 ) * tmp2;
            
            Bx[ii] = 16.0 * pfmdat->epsc * ( Nnx*Nnx - tmp1 ) * phixr[ii];
            By[ii] = 16.0 * pfmdat->epsc * ( Nny*Nny - tmp1 ) * phiyr[ii];
            
            //etam[ii] = ( 1.0 - 3.0*(pfmdat->epsm) ) * ( 1.0 + ( (4.0*(pfmdat->epsm))/(1.0-3.0*(pfmdat->epsm)) )*tmp1 );
             
          }
          else {
            
            neta[ii] = ( 1.0 - 3.0*(pfmdat->epsc) );
            
            Bx[ii] = 0.0;
            By[ii] = 0.0;
            
            //etam[ii] = ( 1.0 - 3.0*(pfmdat->epsm) );

          }
        }
        
        for ( i1 = 0; i1 < 5; i1++ ) { 
          Ba[i1][0] = Bx[i1];
          Ba[i1][1] = By[i1];
          Ba[i1][2] = Bz[i1];
        }
        
        // Inverse Rotation Matrix multiplication with depsilon/dphix, depsilon/dphiy
        for ( i1 = 0; i1 < 5; i1++ ) { 
          for ( ii1 = 0; ii1 < 3; ii1++ ) { 
            tmp0 = 0.0;
            for ( jj1 = 0; jj1 < 3; jj1++ ) { 
              tmp0 += pfmdat->Inv_Rotation_matrix[0][1][ii1][jj1]*Ba[i1][jj1];
            }
            Bir[i1][ii1] = tmp0;
          }
        }
        
        for ( i1 = 0; i1 < 5; i1++ ) { 
          Bxir[i1] = Bir[i1][0];
          Byir[i1] = Bir[i1][1];
          Bzir[i1] = Bir[i1][2];
        }
        
        neta_x = ( neta[1] - neta[2] ) / ( (pfmvar->deltax) ); 
        neta_y = ( neta[3] - neta[4] ) / ( (pfmvar->deltay) );

        alp1 = neta[0]*neta[0]*div_phi + 2.0*neta[0] * ( neta_x * phix[0] + neta_y * phiy[0] );
        
        alp2 = ( neta[1] * Bxir[1] - neta[2] * Bxir[2] ) / ( (pfmvar->deltax) );
        alp3 = ( neta[3] * Byir[3] - neta[4] * Byir[4] ) / ( (pfmvar->deltay) );

        alp = alp1 + alp2 + alp3;
        
        Mphi = (pfmvar->w)/(3.0*(pfmvar->ee)*(pfmdat->a2)*zeta);
        
        if ( fabs(div_phi) > 0.0 )  { 

        phidot_aniso =  Mphi*((pfmvar->ee)*alp - (pfmvar->w)*dgphidphi - dhphidphi*(gs-gl-(csol1-cliq1)*dgldc1l));

        gridNew[ ( ( i )*nx + ( j ) ) ].phi = stgridO[4].phi + (pfmvar->deltat)*phidot_aniso;
        }
		//printf("%le,\t%le,\t%le,\t%le,\t%le,\t%le,\t%le,\t%le\n", Mphi, pfmvar->w, pfmvar->ee, pfmdat->a2, zeta, cliq1, csol1, ddgldc1ldc1l);
		

    }
}

__kernel void SolverPhi_4(__global struct grid *gridOld, __global struct grid *gridNew, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf4 *propf4, __constant struct propmatf4spline *propf4spline) {   
  
    int i;
    int j;
    int k;
    int ii, i1, ii1, jj1;  
    int nx;
    int ny;
    int nz;
    int index;

    double gphi;
    double dgphidphi;
    double dhphidphi; 
    double DetDL;
    double D11invL;
    double gl;
    double gs;
    double dgldc1l;
    double ddgldc1ldc1l;
    double zeta;
    double gradx_phi[3];
    double grady_phi[3];
    double phix[5];
    double phiy[5];
    double div_phi;
    double angleig;
    double Nnx;
    double Nny;
    double neta[5];
    double etam[5];
    double B[5];
    double neta_x;
    double neta_y;
    double alp1;
    double eep_x;
    double eep_y;
    double alp2;
    double alp3;
    double alp;
    double Mphi;
    double Mphi_aniso;
    double Mphi_FiniteMobility;
    double Mphi_FiniteMobilityAniso;
    double phidot_aniso;
    
    double cliq1;
    double csol1;
    double cliq2;
    double csol2;

    double Bf3[npha][nsol], Cf3[npha];
    
     double y[2], retGphase[1], retdmuphase[1], c_sol[1], c_liq[1];

    double Ti;

    double tmp1;

    double phiz[5], phia[5][3], tmp0, phir[5][3], phixr[5], phiyr[5], phizr[5], Bx[5], By[5], Bz[5], modgradphi1, modgradphi4, Bxir[5], Byir[5], Bzir[5], Ba[5][3], Bir[5][3];
    
    struct grid stgridO[9];

    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    if ( i != 0 && i != ny-1 && j != 0 && j != nx-1 ) {

        Ti = temp[j];

        index = (i)*nx + (j);

        stgridO[0] = gridOld[(i-1)*nx + (j-1)];
        stgridO[1] = gridOld[(i-1)*nx + ( j )];
        stgridO[2] = gridOld[(i-1)*nx + (j+1)];
        stgridO[3] = gridOld[( i )*nx + (j-1)];
        stgridO[4] = gridOld[( i )*nx + ( j )];
        stgridO[5] = gridOld[( i )*nx + (j+1)];
        stgridO[6] = gridOld[(i+1)*nx + (j-1)];
        stgridO[7] = gridOld[(i+1)*nx + ( j )];
        stgridO[8] = gridOld[(i+1)*nx + (j+1)];

        cliq1 = cscl[( i )*nx + ( j )].c1l;
        csol1 = cscl[( i )*nx + ( j )].c1s;

        
        div_phi = ( 4.0*(stgridO[1].phi+stgridO[3].phi+stgridO[7].phi+stgridO[5].phi ) + stgridO[0].phi+stgridO[6].phi+stgridO[8].phi+stgridO[2].phi - 20.0*stgridO[4].phi ) / ( 6.0*(pfmvar->deltax)*(pfmvar->deltax) );

        gradx_phi[0] = ( stgridO[5].phi - stgridO[3].phi ) / ( 2.0*(pfmvar->deltax) );
        gradx_phi[1] = ( stgridO[8].phi - stgridO[6].phi ) / ( 2.0*(pfmvar->deltax) );
        gradx_phi[2] = ( stgridO[2].phi - stgridO[0].phi ) / ( 2.0*(pfmvar->deltax) );

        grady_phi[0] = ( stgridO[7].phi - stgridO[1].phi ) / ( 2.0*(pfmvar->deltay) );
        grady_phi[1] = ( stgridO[8].phi - stgridO[2].phi ) / ( 2.0*(pfmvar->deltay) );
        grady_phi[2] = ( stgridO[6].phi - stgridO[0].phi ) / ( 2.0*(pfmvar->deltay) );

        phix[0] = gradx_phi[0];
        phix[1] = ( stgridO[5].phi - stgridO[4].phi ) / ( (pfmvar->deltax) );
        phix[2] = ( stgridO[4].phi - stgridO[3].phi ) / ( (pfmvar->deltax) );
        phix[3] = ( gradx_phi[0] + gradx_phi[1] ) / ( 2.0 );
        phix[4] = ( gradx_phi[0] + gradx_phi[2] ) / ( 2.0 );

        phiy[0] = grady_phi[0];
        phiy[1] = ( grady_phi[0] + grady_phi[1] ) / ( 2.0 );
        phiy[2] = ( grady_phi[0] + grady_phi[2] ) / ( 2.0 );
        phiy[3] = ( stgridO[7].phi - stgridO[4].phi ) / ( (pfmvar->deltay) );
        phiy[4] = ( stgridO[4].phi - stgridO[1].phi ) / ( (pfmvar->deltay) );

        

        D11invL = 1.0/pfmdat->D11l;
        
        if ( pfmdat->ISOTHERMAL ) { 
          
          gs = ( propf4->A[0][0][0]*csol1*csol1 + propf4->B[0][0]*csol1 + propf4->C[0] ) / pfmdat->Vm;
          
          gl = ( propf4->A[1][0][0]*cliq1*cliq1 + propf4->B[1][0]*cliq1 + propf4->C[1] ) / pfmdat->Vm;
          
        }
        else { 
          tmp0 = propf4spline[j].A[0][0][0]*csol1*csol1;
          
          tmp1 = propf4spline[j].B[0][0]*csol1; 
        
          gs = ( tmp0 + tmp1 + propf4spline[j].C[0] ) / pfmdat->Vm;
        
          tmp0 = propf4spline[j].A[1][0][0]*cliq1*cliq1;
          
          tmp1 = propf4spline[j].B[1][0]*cliq1; 
        
          gl = ( tmp0 + tmp1 + propf4spline[j].C[1] ) / pfmdat->Vm;
        }
        
        
        if ( pfmdat->ISOTHERMAL ) { 
          ddgldc1ldc1l = 2.0*propf4->A[1][0][0] / pfmdat->Vm;
        }
        else { 
          ddgldc1ldc1l = 2.0*propf4spline[j].A[1][0][0] / pfmdat->Vm;
        }
        
         dgldc1l = stgridO[4].mu[0] / (pfmdat->Vm);
        
        zeta = ((cliq1-csol1)*(cliq1-csol1)*ddgldc1ldc1l*D11invL);

        gphi = stgridO[4].phi*stgridO[4].phi*(1.0-stgridO[4].phi)*(1.0-stgridO[4].phi);
        dgphidphi = 2.0*stgridO[4].phi*(1.0-stgridO[4].phi)*(1.0-2.0*stgridO[4].phi);
        dhphidphi = 30.0*gphi;


        phiz[0] = 0.0; 
        phiz[1] = 0.0; 
        phiz[2] = 0.0; 
        phiz[3] = 0.0; 
        phiz[4] = 0.0; 
        
        for ( i1 = 0; i1 < 5; i1++ ) { 
          phia[i1][0] = phix[i1];
          phia[i1][1] = phiy[i1];
          phia[i1][2] = phiz[i1];
        }
        
        // Rotation Matrix multiplication with gradients vector
        for ( i1 = 0; i1 < 5; i1++ ) { 
          for ( ii1 = 0; ii1 < 3; ii1++ ) { 
            tmp0 = 0.0;
            for ( jj1 = 0; jj1 < 3; jj1++ ) { 
              tmp0 += pfmdat->Rotation_matrix[0][1][ii1][jj1]*phia[i1][jj1];
              // if ( i == 1 && j == 1 && i1 ==0 && ip1 == 0 && ip2 == 1 ) { 
              //   printf("%le\t", pfmdat->Rotation_matrix[ii1][jj1]);
              // }
            }
            // if ( i == 1 && j == 1 && i1 ==0 && ip1 == 0 && ip2 == 1 ) { 
            //   printf("\n");
            // }
            phir[i1][ii1] = tmp0;
          }
        }
        
        for ( i1 = 0; i1 < 5; i1++ ) { 
          phixr[i1] = phir[i1][0];
          phiyr[i1] = phir[i1][1];
          phizr[i1] = phir[i1][2];
        }
        
        for (ii = 0; ii < 5; ii++) { 
          neta[ii] = 0.0;
          Bx[ii] = 0.0; 
          By[ii] = 0.0; 
          Bz[ii] = 0.0; 
          B[ii] = 0.0; 
        }
        
        for (ii = 0; ii <5; ii++) { 
          
          modgradphi1 = sqrt( phixr[ii]*phixr[ii] + phiyr[ii]*phiyr[ii] );
          modgradphi4 = modgradphi1*modgradphi1*modgradphi1*modgradphi1;
          
          if ( modgradphi4 != 0.0 ) { 
            
            Nnx = phixr[ii] / modgradphi1;
            Nny = phiyr[ii] / modgradphi1;
            
            tmp1 = ( Nnx*Nnx*Nnx*Nnx + Nny*Nny*Nny*Nny );
            
            //tmp2 = ( phixr[ii]*phixr[ii] + phiyr[ii]*phiyr[ii] ); 
            
            neta[ii] = ( 1.0 - 3.0*(pfmdat->epsc) ) * ( 1.0 + ( (4.0*(pfmdat->epsc))/(1.0-3.0*(pfmdat->epsc)) )*tmp1 );
            
            //Bx[ii] = neta[ii] *  16.0 * pfmdat->epsc * ( Nnx*Nnx - tmp1 ) * tmp2;
            //By[ii] = neta[ii] *  16.0 * pfmdat->epsc * ( Nny*Nny - tmp1 ) * tmp2;
            
            Bx[ii] = 16.0 * pfmdat->epsc * ( Nnx*Nnx - tmp1 ) * phixr[ii];
            By[ii] = 16.0 * pfmdat->epsc * ( Nny*Nny - tmp1 ) * phiyr[ii];
            
            //etam[ii] = ( 1.0 - 3.0*(pfmdat->epsm) ) * ( 1.0 + ( (4.0*(pfmdat->epsm))/(1.0-3.0*(pfmdat->epsm)) )*tmp1 );
             
          }
          else {
            
            neta[ii] = ( 1.0 - 3.0*(pfmdat->epsc) );
            
            Bx[ii] = 0.0;
            By[ii] = 0.0;
            
            //etam[ii] = ( 1.0 - 3.0*(pfmdat->epsm) );

          }
        }
        
        for ( i1 = 0; i1 < 5; i1++ ) { 
          Ba[i1][0] = Bx[i1];
          Ba[i1][1] = By[i1];
          Ba[i1][2] = Bz[i1];
        }
        
        // Inverse Rotation Matrix multiplication with depsilon/dphix, depsilon/dphiy
        for ( i1 = 0; i1 < 5; i1++ ) { 
          for ( ii1 = 0; ii1 < 3; ii1++ ) { 
            tmp0 = 0.0;
            for ( jj1 = 0; jj1 < 3; jj1++ ) { 
              tmp0 += pfmdat->Inv_Rotation_matrix[0][1][ii1][jj1]*Ba[i1][jj1];
            }
            Bir[i1][ii1] = tmp0;
          }
        }
        
        for ( i1 = 0; i1 < 5; i1++ ) { 
          Bxir[i1] = Bir[i1][0];
          Byir[i1] = Bir[i1][1];
          Bzir[i1] = Bir[i1][2];
        }
        
        neta_x = ( neta[1] - neta[2] ) / ( (pfmvar->deltax) ); 
        neta_y = ( neta[3] - neta[4] ) / ( (pfmvar->deltay) );

        alp1 = neta[0]*neta[0]*div_phi + 2.0*neta[0] * ( neta_x * phix[0] + neta_y * phiy[0] );
        
        alp2 = ( neta[1] * Bxir[1] - neta[2] * Bxir[2] ) / ( (pfmvar->deltax) );
        alp3 = ( neta[3] * Byir[3] - neta[4] * Byir[4] ) / ( (pfmvar->deltay) );

        alp = alp1 + alp2 + alp3;
        
        Mphi = (pfmvar->w)/(3.0*(pfmvar->ee)*(pfmdat->a2)*zeta);
        
        if ( fabs(div_phi) > 0.0 )  { 

        phidot_aniso =  Mphi*((pfmvar->ee)*alp - (pfmvar->w)*dgphidphi - dhphidphi*(gs-gl-(csol1-cliq1)*dgldc1l));

        gridNew[ ( ( i )*nx + ( j ) ) ].phi = stgridO[4].phi + (pfmvar->deltat)*phidot_aniso;
        }
		//printf("%le,\t%le,\t%le,\t%le,\t%le,\t%le,\t%le,\t%le\n", Mphi, pfmvar->w, pfmvar->ee, pfmdat->a2, zeta, cliq1, csol1, ddgldc1ldc1l);
		

    }
}


__kernel void SolverCatr_2(__global struct grid *gridOld, __global struct grid *gridNew, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf3 *propf3) { 

    int i;
    int j;
    int k;
    int jj, ii, i1, i2, j1, j2;
    int nx;
    int ny;
    int nz;
    int index;

    double hphi[5];
    double A1;
    double alpha1[5];
    double dphidt[5];
    double gradx_phi[3];
    double grady_phi[3];
    double phix[5];
    double phiy[5];
    double modgradphi[5];
    double alphidot1[4];
    double jat1[4];
    double c1jatx;
    double c1jaty;
    double c1jat;
    double term1lx;
    double term1ly;
    double term1l;
    double c1dot;
    
    double fhi[5], hphid[5];
    double gradmux1[1], gradmux2[1], gradmuy1[1], gradmuy2[1], mu[1], Tij[5]; 
    double Ti, ddgldx1ldx1l, ddgsdx1sdx1s, dcdmu[2][1][1], dc_dmu[5][2][1][1];
    double Da[5], Damidx1, Damidx2, Damidy1, Damidy2, divflux, deltamu; 
    double suma, deltac, dcdmudenom, y[2], retmu[1]; 
    int interface, phaseliquid, phasesolid, bulkphase;
    double retdmuphase[1], DELTAT, cg[2], dcbdT_phase[2][1], c_tdt[1], sum_dcbdT; 
    
    struct grid stgridO[9];
    struct grid stgridN[5];
    struct csle stcscl[5];

    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;

    if ( i != 0 && i != ny-1 && j != 0 && j != nx-1 ) {

      Ti = temp[j];

        index = (i*nx + j);

        stgridO[0] = gridOld[(i-1)*nx + (j-1)];
        stgridO[1] = gridOld[(i-1)*nx + ( j )];
        stgridO[2] = gridOld[(i-1)*nx + (j+1)];
        stgridO[3] = gridOld[( i )*nx + (j-1)];
        stgridO[4] = gridOld[( i )*nx + ( j )];
        stgridO[5] = gridOld[( i )*nx + (j+1)];
        stgridO[6] = gridOld[(i+1)*nx + (j-1)];
        stgridO[7] = gridOld[(i+1)*nx + ( j )];
        stgridO[8] = gridOld[(i+1)*nx + (j+1)];

        fhi[0] = stgridO[1].phi;
        fhi[1] = stgridO[3].phi;
        fhi[4] = stgridO[4].phi;
        fhi[2] = stgridO[5].phi;
        fhi[3] = stgridO[7].phi;

        hphi[0] = stgridO[1].phi;
        hphi[1] = stgridO[3].phi;
        hphi[4] = stgridO[4].phi;
        hphi[2] = stgridO[5].phi;
        hphi[3] = stgridO[7].phi;

        hphid[0] = fhi[0]; //stgridO[1].phi;
        hphid[1] = fhi[1]; //stgridO[3].phi;
        hphid[4] = fhi[4]; //stgridO[4].phi;
        hphid[2] = fhi[2]; //stgridO[5].phi;
        hphid[3] = fhi[3]; //stgridO[7].phi;

        stcscl[0] = cscl[(i-1)*nx + ( j )];
        stcscl[1] = cscl[( i )*nx + (j-1)];
        stcscl[4] = cscl[( i )*nx + ( j )];
        stcscl[2] = cscl[( i )*nx + (j+1)];
        stcscl[3] = cscl[(i+1)*nx + ( j )];

        Tij[0] = temp[(i-1)];
        Tij[1] = temp[( i )];
        Tij[4] = temp[( i )];
        Tij[2] = temp[( i )];
        Tij[3] = temp[(i+1)];

        stgridN[0].phi = gridNew[(i-1)*nx + ( j )].phi;
	    stgridN[1].phi = gridNew[( i )*nx + (j-1)].phi;
	    stgridN[4].phi = gridNew[( i )*nx + ( j )].phi;
	    stgridN[2].phi = gridNew[( i )*nx + (j+1)].phi;
	    stgridN[3].phi = gridNew[(i+1)*nx + ( j )].phi;

        dphidt[0] = ( stgridN[4].phi - stgridO[4].phi ) / ( (pfmvar->deltat) );
	    dphidt[1] = ( stgridN[2].phi - stgridO[5].phi ) / ( (pfmvar->deltat) );
	    dphidt[2] = ( stgridN[1].phi - stgridO[3].phi ) / ( (pfmvar->deltat) );
	    dphidt[3] = ( stgridN[3].phi - stgridO[7].phi ) / ( (pfmvar->deltat) );
	    dphidt[4] = ( stgridN[0].phi - stgridO[1].phi ) / ( (pfmvar->deltat) );

        A1 = (sqrt((pfmvar->ee)))/(sqrt(2.0*(pfmvar->w)));

        alpha1[0] = (A1)*( stcscl[4].c1l - stcscl[4].c1s );
        alpha1[1] = (A1)*( stcscl[2].c1l - stcscl[2].c1s );
        alpha1[2] = (A1)*( stcscl[1].c1l - stcscl[1].c1s );
        alpha1[3] = (A1)*( stcscl[3].c1l - stcscl[3].c1s );
        alpha1[4] = (A1)*( stcscl[0].c1l - stcscl[0].c1s );

        gradx_phi[0] = ( stgridO[5].phi - stgridO[3].phi ) / ( 2.0*(pfmvar->deltax) );
        gradx_phi[1] = ( stgridO[8].phi - stgridO[6].phi ) / ( 2.0*(pfmvar->deltax) );
        gradx_phi[2] = ( stgridO[2].phi - stgridO[0].phi ) / ( 2.0*(pfmvar->deltax) );

        grady_phi[0] = ( stgridO[7].phi - stgridO[1].phi ) / ( 2.0*(pfmvar->deltay) );
        grady_phi[1] = ( stgridO[8].phi - stgridO[2].phi ) / ( 2.0*(pfmvar->deltay) );
        grady_phi[2] = ( stgridO[6].phi - stgridO[0].phi ) / ( 2.0*(pfmvar->deltay) );

        phix[0] = gradx_phi[0];
        phix[1] = ( stgridO[5].phi - stgridO[4].phi ) / ( (pfmvar->deltax) );
        phix[2] = ( stgridO[4].phi - stgridO[3].phi ) / ( (pfmvar->deltax) );
        phix[3] = ( gradx_phi[0] + gradx_phi[1] ) / ( 2.0 );
        phix[4] = ( gradx_phi[0] + gradx_phi[2] ) / ( 2.0 );

        phiy[0] = grady_phi[0];
        phiy[1] = ( grady_phi[0] + grady_phi[1] ) / ( 2.0 );
        phiy[2] = ( grady_phi[0] + grady_phi[2] ) / ( 2.0 );
        phiy[3] = ( stgridO[7].phi - stgridO[4].phi ) / ( (pfmvar->deltay) );
        phiy[4] = ( stgridO[4].phi - stgridO[1].phi ) / ( (pfmvar->deltay) );

        alphidot1[0] = ( ( alpha1[0] * dphidt[0] ) + ( alpha1[1] * dphidt[1] ) ) / ( 2.0 );
        alphidot1[1] = ( ( alpha1[0] * dphidt[0] ) + ( alpha1[2] * dphidt[2] ) ) / ( 2.0 );
        alphidot1[2] = ( ( alpha1[0] * dphidt[0] ) + ( alpha1[3] * dphidt[3] ) ) / ( 2.0 );
        alphidot1[3] = ( ( alpha1[0] * dphidt[0] ) + ( alpha1[4] * dphidt[4] ) ) / ( 2.0 );

        for (jj=0; jj<5; jj++) {
            modgradphi[jj] = sqrt( phix[jj]*phix[jj] + phiy[jj]*phiy[jj] );
        }

        jat1[0] = ( alphidot1[0] * phix[1] ) / ( modgradphi[1] );
        jat1[1] = ( alphidot1[1] * phix[2] ) / ( modgradphi[2] );
        jat1[2] = ( alphidot1[2] * phiy[3] ) / ( modgradphi[3] );
        jat1[3] = ( alphidot1[3] * phiy[4] ) / ( modgradphi[4] );

        for (jj=0; jj<4; jj++) {
            if ( modgradphi[jj+1] == 0.0 ) {
                jat1[jj] = 0.0;
            }
        }

        c1jatx = ( jat1[0] - jat1[1] ) / ( (pfmvar->deltax) );
        c1jaty = ( jat1[2] - jat1[3] ) / ( (pfmvar->deltay) );

        c1jat = c1jatx + c1jaty;
        
        for ( ii = 0; ii < 5; ii++ ) { 

          //printf("%d, %d, %d, %d, %le\n", tstep[0], i, j, ii, fhi[ii]);

          interface = 1;
          phaseliquid = 1; 
          phasesolid = 0;
          
          if (fhi[ii] >= pfmdat->interfaceUplimit) { 
            bulkphase = 0;
            interface = 0;
          }
          else if (fhi[ii] <= pfmdat->interfaceDownlimit) { 
            bulkphase = 1;
            interface = 0;
          }

          if ( interface ) { 
            
            y[0] = stcscl[ii].c1s; 
            y[1] = 1.0-y[0];
            
            dMudc(Tij[ii], y, retdmuphase, pfmdat->thermophase[phasesolid]);
            dc_dmu[ii][0][0][0] = 1.0 / retdmuphase[0];
            
            y[0] = stcscl[ii].c1l; 
            y[1] = 1.0-y[0];
            
            dMudc(Tij[ii], y, retdmuphase, pfmdat->thermophase[phaseliquid]);
            dc_dmu[ii][1][0][0] = 1.0 / retdmuphase[0];

            Da[ii] = pfmdat->D11l * (1.0-hphid[ii]) * dc_dmu[ii][1][0][0] + pfmdat->D11s * (hphid[ii]) * dc_dmu[ii][0][0][0];

          }
          else { 
              
              y[0] = stgridO[ii].c1;
              y[1] = 1.0 -y[0];
              
            if ( bulkphase == 0 ) { 
              
              dMudc(Tij[ii], y, retdmuphase, pfmdat->thermophase[0]);
              
              Da[ii] = pfmdat->D11s * ( 1.0 / retdmuphase[0] );
              
            }
            else if ( bulkphase == 1 ) { 
              
              dMudc(Tij[ii], y, retdmuphase, pfmdat->thermophase[1]);
              
              Da[ii] = pfmdat->D11l * ( 1.0 / retdmuphase[0] );
              
            }
          }
        }

        Damidx1 = ( Da[2] + Da[4] ) / 2.0;
        Damidx2 = ( Da[4] + Da[1] ) / 2.0;
        Damidy1 = ( Da[3] + Da[4] ) / 2.0;
        Damidy2 = ( Da[4] + Da[0] ) / 2.0;
        
        gradmux1[0] = ( stgridO[5].mu[0] - stgridO[4].mu[0] ) / pfmvar->deltax;
        gradmux2[0] = ( stgridO[4].mu[0] - stgridO[3].mu[0] ) / pfmvar->deltax;
        gradmuy1[0] = ( stgridO[7].mu[0] - stgridO[4].mu[0] ) / pfmvar->deltay;
        gradmuy2[0] = ( stgridO[4].mu[0] - stgridO[1].mu[0] ) / pfmvar->deltay;

        divflux = ( Damidx1*gradmux1[0] - Damidx2*gradmux2[0] ) / pfmvar->deltax + ( Damidy1*gradmuy1[0] - Damidy2*gradmuy2[0] ) / pfmvar->deltay;
        
        interface = 1;
        phaseliquid = 1; 
        phasesolid = 0; 
        
        if (stgridO[4].phi >= pfmdat->interfaceUplimit) {
          bulkphase=0;
          interface = 0;
        }
        else if (stgridO[4].phi <= pfmdat->interfaceDownlimit) {
          bulkphase = 1;
          interface = 0; 
        }

        if ( interface == 0) { 

          gridNew[index].c1 = stgridO[4].c1 + pfmvar->deltat * (divflux + c1jat);
          //printf("%d, %d, %d, %d, %le %le, %le\n", tstep[0], i, j, interface, stgridO[4].c1, gridNew[index].c1, divflux);
          
          y[0] = gridNew[index].c1; 
          y[1] = 1.0 - y[0];
          
          Mu(Ti, y, retmu, pfmdat->thermophase[bulkphase]);
          
          // if ( phasesolid ) { 
          //   //Mu_0(Ti, y, retmu);
          //   Mu(Ti, y, retmu, pfmdat->thermophase[0]);
          // }
          // else if ( phaseliquid ) { 
          //   //Mu_1(Ti, y, retmu);
          //   Mu(Ti, y, retmu, pfmdat->thermophase[1]);
          // }
          
          deltamu = retmu[0] - stgridO[4].mu[0];
          gridNew[index].mu[0] = retmu[0];
          //printf("0# %d, %d, %d, %d, %le, %le, %le, %le\n", tstep[0], i, j, interface, stgridO[4].phi, stgridO[4].mu[0], gridNew[index].mu[0], deltamu);

        }
        else if ( interface ) { 
        
        if ( !pfmdat->ISOTHERMAL ) {
          DELTAT = pfmvar->deltat * ( -pfmdat->TGRADIENT * pfmdat->velocity ); 
          
          cg[0] = pfmdat->c1s_1stguess; 
          
          mu[0] = gridOld[index].mu[0];
          
          c_mu(mu, c_tdt, Ti+DELTAT, phasesolid, cg, tstep[0], i, j, pfmdat->thermophase[phasesolid]);
          
          dcbdT_phase[0][0] = c_tdt[0] - stcscl[4].c1l;
          
          cg[0] = pfmdat->c1l_1stguess; 
          
          c_mu(mu, c_tdt, Ti+DELTAT, phaseliquid, cg, tstep[0], i, j, pfmdat->thermophase[phaseliquid]);
          
          dcbdT_phase[1][0] = c_tdt[0] - stcscl[4].c1s;
          
        } 
          
          suma = ( stcscl[4].c1s - stcscl[4].c1l ) * (1.0) * ( gridNew[index].phi - stgridO[4].phi );
          
          sum_dcbdT = fhi[4] * dcbdT_phase[0][0] + (1.0-fhi[4]) * dcbdT_phase[1][0];

          deltac = pfmvar->deltat * (divflux + c1jat); 

          gridNew[index].c1 = stgridO[4].c1 + deltac;
          //printf("%d, %d, %d, %d, %le %le, %le\n", tstep[0], i, j, interface, stgridO[4].c1, gridNew[index].c1, divflux);
          
          if ( pfmdat->ISOTHERMAL ) { 
            deltac = deltac - suma;
          }
          else { 
            deltac = deltac - suma  - sum_dcbdT;
          }

          //hphi[4] = 3.0*fhi[4]*fhi[4] - 2.0*fhi[4]*fhi[4]*fhi[4];
          hphi[4] = fhi[4];
          dcdmudenom = dc_dmu[4][0][0][0]*hphi[4] + dc_dmu[4][1][0][0]*(1.0-hphi[4]);
          
          gridNew[index].mu[0] = stgridO[4].mu[0] + deltac / dcdmudenom;
          //printf("1# %d, %d, %d, %d, %le, %le, %le, %le, %le, %le, %le, %le\n", tstep[0], i, j, interface, stgridO[4].phi, stgridO[4].mu[0], dcdmu[0][0][0], hphi[4], dcdmu[1][0][0], (1.0-hphi[4]), gridNew[index].mu[0], dcdmudenom);
          //printf("1# %d, %d, %d, %d, %le, %le, %le, %le, %le, %le, %le\n", tstep[0], i, j, interface, stgridO[4].phi, stgridO[4].mu[0],  gridNew[index].mu[0], deltac / dcdmudenom, suma, ( gridNew[index].phi - stgridO[4].phi ), (6.0*fhi[4]*(1.0-fhi[4])) );

        }
    }
}

__kernel void SolverCatr_3(__global struct grid *gridOld, __global struct grid *gridNew, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf3 *propf3) { 

    int i;
    int j;
    int k;
    int jj, ii, i1, i2, j1, j2;
    int nx;
    int ny;
    int nz;
    int index;

    double hphi[5];
    double A1;
    double alpha1[5];
    double dphidt[5];
    double gradx_phi[3];
    double grady_phi[3];
    double phix[5];
    double phiy[5];
    double modgradphi[5];
    double alphidot1[4];
    double jat1[4];
    double c1jatx;
    double c1jaty;
    double c1jat;
    double term1lx;
    double term1ly;
    double term1l;
    double c1dot;
    
    double fhi[5], hphid[5];
    double gradmux1[1], gradmux2[1], gradmuy1[1], gradmuy2[1], mu[1]; 
    double Ti, ddgldx1ldx1l, ddgsdx1sdx1s, dcdmu[2][1][1], dc_dmu[5][2][1][1];
    double Da[5], Damidx1, Damidx2, Damidy1, Damidy2, divflux, deltamu; 
    double suma, deltac, dcdmudenom, call[2], retmu[1], Tij[5]; 
    int interface, phaseliquid, phasesolid, bulkphase;
    double retdmuphase[1], DELTAT, cg[2], dcbdT_phase[2][1], c_tdt[1], sum_dcbdT; 
    
    struct grid stgridO[9];
    struct grid stgridN[5];
    struct csle stcscl[5];

    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;

    if ( i != 0 && i != ny-1 && j != 0 && j != nx-1 ) {

      Ti = temp[j];

        index = (i*nx + j);

        stgridO[0] = gridOld[(i-1)*nx + (j-1)];
        stgridO[1] = gridOld[(i-1)*nx + ( j )];
        stgridO[2] = gridOld[(i-1)*nx + (j+1)];
        stgridO[3] = gridOld[( i )*nx + (j-1)];
        stgridO[4] = gridOld[( i )*nx + ( j )];
        stgridO[5] = gridOld[( i )*nx + (j+1)];
        stgridO[6] = gridOld[(i+1)*nx + (j-1)];
        stgridO[7] = gridOld[(i+1)*nx + ( j )];
        stgridO[8] = gridOld[(i+1)*nx + (j+1)];

        fhi[0] = stgridO[1].phi;
        fhi[1] = stgridO[3].phi;
        fhi[4] = stgridO[4].phi;
        fhi[2] = stgridO[5].phi;
        fhi[3] = stgridO[7].phi;

        hphi[0] = stgridO[1].phi;
        hphi[1] = stgridO[3].phi;
        hphi[4] = stgridO[4].phi;
        hphi[2] = stgridO[5].phi;
        hphi[3] = stgridO[7].phi;

        hphid[0] = fhi[0]; //stgridO[1].phi;
        hphid[1] = fhi[1]; //stgridO[3].phi;
        hphid[4] = fhi[4]; //stgridO[4].phi;
        hphid[2] = fhi[2]; //stgridO[5].phi;
        hphid[3] = fhi[3]; //stgridO[7].phi;

        stcscl[0] = cscl[(i-1)*nx + ( j )];
        stcscl[1] = cscl[( i )*nx + (j-1)];
        stcscl[4] = cscl[( i )*nx + ( j )];
        stcscl[2] = cscl[( i )*nx + (j+1)];
        stcscl[3] = cscl[(i+1)*nx + ( j )];

        Tij[0] = temp[(i-1)];
        Tij[1] = temp[( i )];
        Tij[4] = temp[( i )];
        Tij[2] = temp[( i )];
        Tij[3] = temp[(i+1)];

        stgridN[0].phi = gridNew[(i-1)*nx + ( j )].phi;
	    stgridN[1].phi = gridNew[( i )*nx + (j-1)].phi;
	    stgridN[4].phi = gridNew[( i )*nx + ( j )].phi;
	    stgridN[2].phi = gridNew[( i )*nx + (j+1)].phi;
	    stgridN[3].phi = gridNew[(i+1)*nx + ( j )].phi;

        dphidt[0] = ( stgridN[4].phi - stgridO[4].phi ) / ( (pfmvar->deltat) );
	    dphidt[1] = ( stgridN[2].phi - stgridO[5].phi ) / ( (pfmvar->deltat) );
	    dphidt[2] = ( stgridN[1].phi - stgridO[3].phi ) / ( (pfmvar->deltat) );
	    dphidt[3] = ( stgridN[3].phi - stgridO[7].phi ) / ( (pfmvar->deltat) );
	    dphidt[4] = ( stgridN[0].phi - stgridO[1].phi ) / ( (pfmvar->deltat) );

        A1 = (sqrt((pfmvar->ee)))/(sqrt(2.0*(pfmvar->w)));

        alpha1[0] = (A1)*( stcscl[4].c1l - stcscl[4].c1s );
        alpha1[1] = (A1)*( stcscl[2].c1l - stcscl[2].c1s );
        alpha1[2] = (A1)*( stcscl[1].c1l - stcscl[1].c1s );
        alpha1[3] = (A1)*( stcscl[3].c1l - stcscl[3].c1s );
        alpha1[4] = (A1)*( stcscl[0].c1l - stcscl[0].c1s );

        gradx_phi[0] = ( stgridO[5].phi - stgridO[3].phi ) / ( 2.0*(pfmvar->deltax) );
        gradx_phi[1] = ( stgridO[8].phi - stgridO[6].phi ) / ( 2.0*(pfmvar->deltax) );
        gradx_phi[2] = ( stgridO[2].phi - stgridO[0].phi ) / ( 2.0*(pfmvar->deltax) );

        grady_phi[0] = ( stgridO[7].phi - stgridO[1].phi ) / ( 2.0*(pfmvar->deltay) );
        grady_phi[1] = ( stgridO[8].phi - stgridO[2].phi ) / ( 2.0*(pfmvar->deltay) );
        grady_phi[2] = ( stgridO[6].phi - stgridO[0].phi ) / ( 2.0*(pfmvar->deltay) );

        phix[0] = gradx_phi[0];
        phix[1] = ( stgridO[5].phi - stgridO[4].phi ) / ( (pfmvar->deltax) );
        phix[2] = ( stgridO[4].phi - stgridO[3].phi ) / ( (pfmvar->deltax) );
        phix[3] = ( gradx_phi[0] + gradx_phi[1] ) / ( 2.0 );
        phix[4] = ( gradx_phi[0] + gradx_phi[2] ) / ( 2.0 );

        phiy[0] = grady_phi[0];
        phiy[1] = ( grady_phi[0] + grady_phi[1] ) / ( 2.0 );
        phiy[2] = ( grady_phi[0] + grady_phi[2] ) / ( 2.0 );
        phiy[3] = ( stgridO[7].phi - stgridO[4].phi ) / ( (pfmvar->deltay) );
        phiy[4] = ( stgridO[4].phi - stgridO[1].phi ) / ( (pfmvar->deltay) );

        alphidot1[0] = ( ( alpha1[0] * dphidt[0] ) + ( alpha1[1] * dphidt[1] ) ) / ( 2.0 );
        alphidot1[1] = ( ( alpha1[0] * dphidt[0] ) + ( alpha1[2] * dphidt[2] ) ) / ( 2.0 );
        alphidot1[2] = ( ( alpha1[0] * dphidt[0] ) + ( alpha1[3] * dphidt[3] ) ) / ( 2.0 );
        alphidot1[3] = ( ( alpha1[0] * dphidt[0] ) + ( alpha1[4] * dphidt[4] ) ) / ( 2.0 );

        for (jj=0; jj<5; jj++) {
            modgradphi[jj] = sqrt( phix[jj]*phix[jj] + phiy[jj]*phiy[jj] );
        }

        jat1[0] = ( alphidot1[0] * phix[1] ) / ( modgradphi[1] );
        jat1[1] = ( alphidot1[1] * phix[2] ) / ( modgradphi[2] );
        jat1[2] = ( alphidot1[2] * phiy[3] ) / ( modgradphi[3] );
        jat1[3] = ( alphidot1[3] * phiy[4] ) / ( modgradphi[4] );

        for (jj=0; jj<4; jj++) {
            if ( modgradphi[jj+1] == 0.0 ) {
                jat1[jj] = 0.0;
            }
        }

        c1jatx = ( jat1[0] - jat1[1] ) / ( (pfmvar->deltax) );
        c1jaty = ( jat1[2] - jat1[3] ) / ( (pfmvar->deltay) );

        c1jat = c1jatx + c1jaty;
        

        for ( ii = 0; ii < 5; ii++ ) { 

          //printf("%d, %d, %d, %d, %le\n", tstep[0], i, j, ii, fhi[ii]);

          interface = 1;
          phaseliquid = 1; 
          phasesolid = 0; 
          
          if (fhi[ii] >= pfmdat->interfaceUplimit) { 
            bulkphase = 0;
            interface = 0;
          }
          else if (fhi[ii] <= pfmdat->interfaceDownlimit) { 
            bulkphase = 1;
            interface = 0;
          }

          if ( interface ) { 

            dc_dmu[ii][0][0][0] = propf3->cmu[0][0][0];
            dc_dmu[ii][1][0][0] = propf3->cmu[1][0][0];
            
            //dc_dmu[ii][0][0][0] = 1.0 / (2.0*propf3->A[0][0][0]);
            //dc_dmu[ii][1][0][0] = 1.0 / (2.0*propf3->A[1][0][0]);

            Da[ii] = pfmdat->D11l * (1.0-hphid[ii]) * dc_dmu[ii][1][0][0] + pfmdat->D11s * (hphid[ii]) * dc_dmu[ii][0][0][0];
            
          }
          else { 
            if ( bulkphase == 0 ) { 
              
              Da[ii] = pfmdat->D11s * propf3->cmu[0][0][0];
              
              //Da[ii] = pfmdat->D11s * ( 1.0 / (2.0*propf3->A[0][0][0]) );
              
            }
            else if ( bulkphase == 1 ) { 
              
              Da[ii] = pfmdat->D11l * propf3->cmu[1][0][0];
              
              //Da[ii] = pfmdat->D11l * ( 1.0 / (2.0*propf3->A[1][0][0]) );
              
            }
          }
        }

        Damidx1 = ( Da[2] + Da[4] ) / 2.0;
        Damidx2 = ( Da[4] + Da[1] ) / 2.0;
        Damidy1 = ( Da[3] + Da[4] ) / 2.0;
        Damidy2 = ( Da[4] + Da[0] ) / 2.0;
        
        gradmux1[0] = ( stgridO[5].mu[0] - stgridO[4].mu[0] ) / pfmvar->deltax;
        gradmux2[0] = ( stgridO[4].mu[0] - stgridO[3].mu[0] ) / pfmvar->deltax;
        gradmuy1[0] = ( stgridO[7].mu[0] - stgridO[4].mu[0] ) / pfmvar->deltay;
        gradmuy2[0] = ( stgridO[4].mu[0] - stgridO[1].mu[0] ) / pfmvar->deltay;

        divflux = ( Damidx1*gradmux1[0] - Damidx2*gradmux2[0] ) / pfmvar->deltax + ( Damidy1*gradmuy1[0] - Damidy2*gradmuy2[0] ) / pfmvar->deltay;
        
        interface = 1;
        phaseliquid = 1; 
        phasesolid = 0; 

        if (stgridO[4].phi >= pfmdat->interfaceUplimit) {
          bulkphase=0;
          interface = 0;
        }
        else if (stgridO[4].phi <= pfmdat->interfaceDownlimit) {
          bulkphase = 1;
          interface = 0; 
        }

        if ( interface == 0) { 

          gridNew[index].c1 = stgridO[4].c1 + pfmvar->deltat * (divflux + c1jat);
          //printf("%d, %d, %d, %d, %le %le, %le\n", tstep[0], i, j, interface, stgridO[4].c1, gridNew[index].c1, divflux);

          if ( bulkphase == 0 ) { 
            retmu[0] = 2.0*propf3->A[0][0][0]*gridNew[index].c1 + (propf3->Beq[0][0] + propf3->dBbdT[0][0]*(Ti-pfmdat->Teq));
          }
          else if ( bulkphase == 1 ) { 
            retmu[0] = 2.0*propf3->A[1][0][0]*gridNew[index].c1 + (propf3->Beq[1][0] + propf3->dBbdT[1][0]*(Ti-pfmdat->Teq));
          }

          deltamu = retmu[0] - stgridO[4].mu[0];
          gridNew[index].mu[0] = retmu[0];
          //printf("0# %d, %d, %d, %d, %le, %le, %le, %le\n", tstep[0], i, j, interface, stgridO[4].phi, stgridO[4].mu[0], gridNew[index].mu[0], deltamu);

        }
        else if ( interface ) {  
          
          if ( !pfmdat->ISOTHERMAL ) {
            
            DELTAT = pfmvar->deltat * ( -pfmdat->TGRADIENT * pfmdat->velocity ); 
            
            c_tdt[0] = propf3->cmu[0][0][0] * ( stgridO[4].mu[0] - ( propf3->Beq[0][0] + propf3->dBbdT[0][0] * (Ti+DELTAT-pfmdat->Teq) ) );
            
            dcbdT_phase[0][0] = c_tdt[0] - stcscl[4].c1s;
            
            c_tdt[0] = propf3->cmu[1][0][0] * ( stgridO[4].mu[0] - ( propf3->Beq[1][0] + propf3->dBbdT[1][0] * (Ti+DELTAT-pfmdat->Teq) ) );
            
            dcbdT_phase[1][0] = c_tdt[0] - stcscl[4].c1l;
          
          }
          
          suma = ( stcscl[4].c1s - stcscl[4].c1l ) * (1.0) * ( gridNew[index].phi - stgridO[4].phi );
          
          sum_dcbdT = fhi[4] * dcbdT_phase[0][0] + (1.0-fhi[4]) * dcbdT_phase[1][0];

          deltac = pfmvar->deltat * (divflux + c1jat); 

          gridNew[index].c1 = stgridO[4].c1 + deltac;
          //printf("%d, %d, %d, %d, %le %le, %le\n", tstep[0], i, j, interface, stgridO[4].c1, gridNew[index].c1, divflux);
          
          if ( pfmdat->ISOTHERMAL ) { 
            deltac = deltac - suma;
          }
          else { 
            deltac = deltac - suma  - sum_dcbdT;
          }
          
          //hphi[4] = 3.0*fhi[4]*fhi[4] - 2.0*fhi[4]*fhi[4]*fhi[4];
          hphi[4] = fhi[4];
          dcdmudenom = dc_dmu[4][0][0][0]*hphi[4] + dc_dmu[4][1][0][0]*(1.0-hphi[4]);
          
          gridNew[index].mu[0] = stgridO[4].mu[0] + deltac / dcdmudenom;
          //printf("1# %d, %d, %d, %d, %le, %le, %le, %le, %le, %le, %le, %le\n", tstep[0], i, j, interface, stgridO[4].phi, stgridO[4].mu[0], dcdmu[0][0][0], hphi[4], dcdmu[1][0][0], (1.0-hphi[4]), gridNew[index].mu[0], dcdmudenom);
          //printf("1# %d, %d, %d, %d, %le, %le, %le, %le, %le, %le, %le\n", tstep[0], i, j, interface, stgridO[4].phi, stgridO[4].mu[0],  gridNew[index].mu[0], deltac / dcdmudenom, suma, ( gridNew[index].phi - stgridO[4].phi ), (6.0*fhi[4]*(1.0-fhi[4])) );
        }
        
    }
}

__kernel void SolverCatr_4(__global struct grid *gridOld, __global struct grid *gridNew, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf4 *propf4, __constant struct propmatf4spline *propf4spline, __constant struct propmatf4spline *propf4spline1) { 

    int i;
    int j;
    int k;
    int jj, ii, i1, i2, j1, j2;
    int nx;
    int ny;
    int nz;
    int index;
    int indij[5];

    double hphi[5];
    double A1;
    double alpha1[5];
    double dphidt[5];
    double gradx_phi[3];
    double grady_phi[3];
    double phix[5];
    double phiy[5];
    double modgradphi[5];
    double alphidot1[4];
    double jat1[4];
    double c1jatx;
    double c1jaty;
    double c1jat;
    double term1lx;
    double term1ly;
    double term1l;
    double c1dot;
    
    double fhi[5], hphid[5];
    double gradmux1[1], gradmux2[1], gradmuy1[1], gradmuy2[1], mu[1]; 
    double Ti, ddgldx1ldx1l, ddgsdx1sdx1s, dcdmu[2][1][1], dc_dmu[5][2][1][1];
    double Da[5], Damidx1, Damidx2, Damidy1, Damidy2, divflux, deltamu; 
    double suma, deltac, dcdmudenom, call[2], retmu[1], Tij[5]; 
    int interface, phaseliquid, phasesolid, bulkphase;
    double retdmuphase[1], DELTAT, cg[2], dcbdT_phase[2][1], c_tdt[1], sum_dcbdT, muc1[1], cmu1[1]; 
    
    struct grid stgridO[9];
    struct grid stgridN[5];
    struct csle stcscl[5];

    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;

    if ( i != 0 && i != ny-1 && j != 0 && j != nx-1 ) {

      Ti = temp[j];

        index = (i*nx + j);

        stgridO[0] = gridOld[(i-1)*nx + (j-1)];
        stgridO[1] = gridOld[(i-1)*nx + ( j )];
        stgridO[2] = gridOld[(i-1)*nx + (j+1)];
        stgridO[3] = gridOld[( i )*nx + (j-1)];
        stgridO[4] = gridOld[( i )*nx + ( j )];
        stgridO[5] = gridOld[( i )*nx + (j+1)];
        stgridO[6] = gridOld[(i+1)*nx + (j-1)];
        stgridO[7] = gridOld[(i+1)*nx + ( j )];
        stgridO[8] = gridOld[(i+1)*nx + (j+1)];

        fhi[0] = stgridO[1].phi;
        fhi[1] = stgridO[3].phi;
        fhi[4] = stgridO[4].phi;
        fhi[2] = stgridO[5].phi;
        fhi[3] = stgridO[7].phi;

        hphi[0] = stgridO[1].phi;
        hphi[1] = stgridO[3].phi;
        hphi[4] = stgridO[4].phi;
        hphi[2] = stgridO[5].phi;
        hphi[3] = stgridO[7].phi;

        hphid[0] = fhi[0]; //stgridO[1].phi;
        hphid[1] = fhi[1]; //stgridO[3].phi;
        hphid[4] = fhi[4]; //stgridO[4].phi;
        hphid[2] = fhi[2]; //stgridO[5].phi;
        hphid[3] = fhi[3]; //stgridO[7].phi;

        stcscl[0] = cscl[(i-1)*nx + ( j )];
        stcscl[1] = cscl[( i )*nx + (j-1)];
        stcscl[4] = cscl[( i )*nx + ( j )];
        stcscl[2] = cscl[( i )*nx + (j+1)];
        stcscl[3] = cscl[(i+1)*nx + ( j )];

        Tij[0] = temp[(i-1)];
        Tij[1] = temp[( i )];
        Tij[4] = temp[( i )];
        Tij[2] = temp[( i )];
        Tij[3] = temp[(i+1)];
        
        indij[0] = ( j );
        indij[1] = (j-1);
        indij[4] = ( j );
        indij[2] = (j+1);
        indij[3] = ( j );

        stgridN[0].phi = gridNew[(i-1)*nx + ( j )].phi;
	    stgridN[1].phi = gridNew[( i )*nx + (j-1)].phi;
	    stgridN[4].phi = gridNew[( i )*nx + ( j )].phi;
	    stgridN[2].phi = gridNew[( i )*nx + (j+1)].phi;
	    stgridN[3].phi = gridNew[(i+1)*nx + ( j )].phi;

        dphidt[0] = ( stgridN[4].phi - stgridO[4].phi ) / ( (pfmvar->deltat) );
	    dphidt[1] = ( stgridN[2].phi - stgridO[5].phi ) / ( (pfmvar->deltat) );
	    dphidt[2] = ( stgridN[1].phi - stgridO[3].phi ) / ( (pfmvar->deltat) );
	    dphidt[3] = ( stgridN[3].phi - stgridO[7].phi ) / ( (pfmvar->deltat) );
	    dphidt[4] = ( stgridN[0].phi - stgridO[1].phi ) / ( (pfmvar->deltat) );

        A1 = (sqrt((pfmvar->ee)))/(sqrt(2.0*(pfmvar->w)));

        alpha1[0] = (A1)*( stcscl[4].c1l - stcscl[4].c1s );
        alpha1[1] = (A1)*( stcscl[2].c1l - stcscl[2].c1s );
        alpha1[2] = (A1)*( stcscl[1].c1l - stcscl[1].c1s );
        alpha1[3] = (A1)*( stcscl[3].c1l - stcscl[3].c1s );
        alpha1[4] = (A1)*( stcscl[0].c1l - stcscl[0].c1s );

        gradx_phi[0] = ( stgridO[5].phi - stgridO[3].phi ) / ( 2.0*(pfmvar->deltax) );
        gradx_phi[1] = ( stgridO[8].phi - stgridO[6].phi ) / ( 2.0*(pfmvar->deltax) );
        gradx_phi[2] = ( stgridO[2].phi - stgridO[0].phi ) / ( 2.0*(pfmvar->deltax) );

        grady_phi[0] = ( stgridO[7].phi - stgridO[1].phi ) / ( 2.0*(pfmvar->deltay) );
        grady_phi[1] = ( stgridO[8].phi - stgridO[2].phi ) / ( 2.0*(pfmvar->deltay) );
        grady_phi[2] = ( stgridO[6].phi - stgridO[0].phi ) / ( 2.0*(pfmvar->deltay) );

        phix[0] = gradx_phi[0];
        phix[1] = ( stgridO[5].phi - stgridO[4].phi ) / ( (pfmvar->deltax) );
        phix[2] = ( stgridO[4].phi - stgridO[3].phi ) / ( (pfmvar->deltax) );
        phix[3] = ( gradx_phi[0] + gradx_phi[1] ) / ( 2.0 );
        phix[4] = ( gradx_phi[0] + gradx_phi[2] ) / ( 2.0 );

        phiy[0] = grady_phi[0];
        phiy[1] = ( grady_phi[0] + grady_phi[1] ) / ( 2.0 );
        phiy[2] = ( grady_phi[0] + grady_phi[2] ) / ( 2.0 );
        phiy[3] = ( stgridO[7].phi - stgridO[4].phi ) / ( (pfmvar->deltay) );
        phiy[4] = ( stgridO[4].phi - stgridO[1].phi ) / ( (pfmvar->deltay) );

        alphidot1[0] = ( ( alpha1[0] * dphidt[0] ) + ( alpha1[1] * dphidt[1] ) ) / ( 2.0 );
        alphidot1[1] = ( ( alpha1[0] * dphidt[0] ) + ( alpha1[2] * dphidt[2] ) ) / ( 2.0 );
        alphidot1[2] = ( ( alpha1[0] * dphidt[0] ) + ( alpha1[3] * dphidt[3] ) ) / ( 2.0 );
        alphidot1[3] = ( ( alpha1[0] * dphidt[0] ) + ( alpha1[4] * dphidt[4] ) ) / ( 2.0 );

        for (jj=0; jj<5; jj++) {
            modgradphi[jj] = sqrt( phix[jj]*phix[jj] + phiy[jj]*phiy[jj] );
        }

        jat1[0] = ( alphidot1[0] * phix[1] ) / ( modgradphi[1] );
        jat1[1] = ( alphidot1[1] * phix[2] ) / ( modgradphi[2] );
        jat1[2] = ( alphidot1[2] * phiy[3] ) / ( modgradphi[3] );
        jat1[3] = ( alphidot1[3] * phiy[4] ) / ( modgradphi[4] );

        for (jj=0; jj<4; jj++) {
            if ( modgradphi[jj+1] == 0.0 ) {
                jat1[jj] = 0.0;
            }
        }

        c1jatx = ( jat1[0] - jat1[1] ) / ( (pfmvar->deltax) );
        c1jaty = ( jat1[2] - jat1[3] ) / ( (pfmvar->deltay) );

        c1jat = c1jatx + c1jaty;
        

        for ( ii = 0; ii < 5; ii++ ) { 

          //printf("%d, %d, %d, %d, %le\n", tstep[0], i, j, ii, fhi[ii]);

          interface = 1;
          phaseliquid = 1; 
          phasesolid = 0; 
          
          if (fhi[ii] >= pfmdat->interfaceUplimit) { 
            bulkphase = 0;
            interface = 0;
          }
          else if (fhi[ii] <= pfmdat->interfaceDownlimit) { 
            bulkphase = 1;
            interface = 0;
          }

          if ( interface ) { 

            if ( pfmdat->ISOTHERMAL ) { 
               dc_dmu[ii][0][0][0] = propf4->cmu[0][0][0];
               dc_dmu[ii][1][0][0] = propf4->cmu[1][0][0];
             }
             else { 
               dc_dmu[ii][0][0][0] = 1 / (2.0 * propf4spline[indij[ii]].A[0][0][0]);
               dc_dmu[ii][1][0][0] = 1 / (2.0 * propf4spline[indij[ii]].A[1][0][0]);
             }

            Da[ii] = pfmdat->D11l * (1.0-hphid[ii]) * dc_dmu[ii][1][0][0] + pfmdat->D11s * (hphid[ii]) * dc_dmu[ii][0][0][0];
            
          }
          else { 
           if ( pfmdat->ISOTHERMAL ) { 
            if ( bulkphase == 0 ) { 
              
              dc_dmu[ii][0][0][0] = propf4->cmu[bulkphase][0][0];
              
              
              Da[ii] = pfmdat->D11s * dc_dmu[ii][0][0][0];
              
            }
            else if ( bulkphase == 1 ) { 
              
              dc_dmu[ii][1][0][0] = propf4->cmu[bulkphase][0][0];
              
              Da[ii] = pfmdat->D11l * dc_dmu[ii][1][0][0];
              
            }
           }
           else {
              
            if ( bulkphase == 0 ) { 
              
              dc_dmu[ii][0][0][0] = 1.0 / (2.0 * propf4spline[indij[ii]].A[bulkphase][0][0]);
              
              
              Da[ii] = pfmdat->D11s * dc_dmu[ii][0][0][0];
              
            }
            else if ( bulkphase == 1 ) { 
              
              dc_dmu[ii][1][0][0] = 1.0 / (2.0 * propf4spline[indij[ii]].A[bulkphase][0][0]);
              
              Da[ii] = pfmdat->D11l * dc_dmu[ii][1][0][0];
              
            }
           }
           
          }
        }

        Damidx1 = ( Da[2] + Da[4] ) / 2.0;
        Damidx2 = ( Da[4] + Da[1] ) / 2.0;
        Damidy1 = ( Da[3] + Da[4] ) / 2.0;
        Damidy2 = ( Da[4] + Da[0] ) / 2.0;
        
        gradmux1[0] = ( stgridO[5].mu[0] - stgridO[4].mu[0] ) / pfmvar->deltax;
        gradmux2[0] = ( stgridO[4].mu[0] - stgridO[3].mu[0] ) / pfmvar->deltax;
        gradmuy1[0] = ( stgridO[7].mu[0] - stgridO[4].mu[0] ) / pfmvar->deltay;
        gradmuy2[0] = ( stgridO[4].mu[0] - stgridO[1].mu[0] ) / pfmvar->deltay;

        divflux = ( Damidx1*gradmux1[0] - Damidx2*gradmux2[0] ) / pfmvar->deltax + ( Damidy1*gradmuy1[0] - Damidy2*gradmuy2[0] ) / pfmvar->deltay;
        
        interface = 1;
        phaseliquid = 1; 
        phasesolid = 0; 

        if (stgridO[4].phi >= pfmdat->interfaceUplimit) {
          bulkphase=0;
          interface = 0;
        }
        else if (stgridO[4].phi <= pfmdat->interfaceDownlimit) {
          bulkphase = 1;
          interface = 0; 
        }

        if ( interface == 0) { 

          gridNew[index].c1 = stgridO[4].c1 + pfmvar->deltat * (divflux + c1jat);
          //printf("%d, %d, %d, %d, %le %le, %le\n", tstep[0], i, j, interface, stgridO[4].c1, gridNew[index].c1, divflux);

          if ( pfmdat->ISOTHERMAL ) { 
          
          if ( bulkphase == 0 ) { 
            retmu[0] = 2.0*propf4->A[bulkphase][0][0]*gridNew[index].c1 + propf4->B[bulkphase][0];
          }
          else if ( bulkphase == 1 ) { 
            retmu[0] = 2.0*propf4->A[bulkphase][0][0]*gridNew[index].c1 + propf4->B[bulkphase][0];
          }
          
         }
         else { 
           
           if ( bulkphase == 0 ) { 
            retmu[0] = 2.0*propf4spline[j].A[bulkphase][0][0]*gridNew[index].c1 + propf4spline[j].B[bulkphase][0];
          }
          else if ( bulkphase == 1 ) { 
            retmu[0] = 2.0*propf4spline[j].A[bulkphase][0][0]*gridNew[index].c1 + propf4spline[j].B[bulkphase][0];
          }
           
         }
          

          deltamu = retmu[0] - stgridO[4].mu[0];
          gridNew[index].mu[0] = retmu[0];
          //printf("0# %d, %d, %d, %d, %le, %le, %le, %le\n", tstep[0], i, j, interface, stgridO[4].phi, stgridO[4].mu[0], gridNew[index].mu[0], deltamu);

        }
        else if ( interface ) {  
          
          if ( pfmdat->ISOTHERMAL ) { 
               dcbdT_phase[0][0] = 0.0;
               dcbdT_phase[1][0] = 0.0;
          }
          else { 
            
            DELTAT = ( pfmvar->deltat ) * ( -pfmdat->TGRADIENT * pfmdat->velocity ); 
            
            muc1[0] = 2.0 * propf4spline1[j].A[0][0][0]; 
            
            cmu1[0] = 1.0 / muc1[0];
            
            c_tdt[0] = cmu1[0] * ( gridOld[index].mu[0] - propf4spline1[j].B[0][0] );
            
            dcbdT_phase[0][0] = c_tdt[0] - stcscl[4].c1s;
            
            muc1[0] = 2.0 * propf4spline1[j].A[1][0][0]; 
            
            cmu1[0] = 1.0 / muc1[0];
            
            c_tdt[0] = cmu1[0] * ( gridOld[index].mu[0] - propf4spline1[j].B[1][0] );
            
            dcbdT_phase[1][0] = c_tdt[0] - stcscl[4].c1l;
            
          }
          
          suma = ( stcscl[4].c1s - stcscl[4].c1l ) * (1.0) * ( gridNew[index].phi - stgridO[4].phi );
          
          sum_dcbdT = fhi[4] * dcbdT_phase[0][0] + (1.0-fhi[4]) * dcbdT_phase[1][0];

          deltac = pfmvar->deltat * (divflux + c1jat); 

          gridNew[index].c1 = stgridO[4].c1 + deltac;
          //printf("%d, %d, %d, %d, %le %le, %le\n", tstep[0], i, j, interface, stgridO[4].c1, gridNew[index].c1, divflux);
          
          if ( pfmdat->ISOTHERMAL ) { 
            deltac = deltac - suma;
          }
          else { 
            deltac = deltac - suma  - sum_dcbdT;
          }
          
          //hphi[4] = 3.0*fhi[4]*fhi[4] - 2.0*fhi[4]*fhi[4]*fhi[4];
          hphi[4] = fhi[4];
          dcdmudenom = dc_dmu[4][0][0][0]*hphi[4] + dc_dmu[4][1][0][0]*(1.0-hphi[4]);
          
          gridNew[index].mu[0] = stgridO[4].mu[0] + deltac / dcdmudenom;
          //printf("1# %d, %d, %d, %d, %le, %le, %le, %le, %le, %le, %le, %le\n", tstep[0], i, j, interface, stgridO[4].phi, stgridO[4].mu[0], dcdmu[0][0][0], hphi[4], dcdmu[1][0][0], (1.0-hphi[4]), gridNew[index].mu[0], dcdmudenom);
          //printf("1# %d, %d, %d, %d, %le, %le, %le, %le, %le, %le, %le\n", tstep[0], i, j, interface, stgridO[4].phi, stgridO[4].mu[0],  gridNew[index].mu[0], deltac / dcdmudenom, suma, ( gridNew[index].phi - stgridO[4].phi ), (6.0*fhi[4]*(1.0-fhi[4])) );
        }
        
    }
}

__kernel void apply_BC(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    
    int i;
    int j;
    int k;
    int nx;
    int ny;
    int nz;
    
    int index;
    int indexa;
    
    int me;
    
    int istart;
    int iend;
    int totny;
    int nxtotny;
    int totsliceny;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    istart = 0;
    iend = ny;
    
    me = pfmdat->myrank;
    
    index = (i*nx + j);
    
      if ( i == 0 ) {
	gridNew[(i*nx + j)] = gridNew[((i+1)*nx + j)];
      }  
    
      if ( i == (iend-1) ) {
	gridNew[(i*nx + j)] = gridNew[((i-1)*nx + j)];
      }  
    
    if ( j == 0 ) {
        gridNew[(i*nx + j)] = gridNew[(i*nx + (j+1))];
    }
    
    if ( j == (nx-1) ) {
        gridNew[(i*nx + j)] = gridNew[(i*nx + (j-1))];
    }
    
}

__kernel void apply_BC_phi(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    
    int i;
    int j;
    int k;
    int nx;
    int ny;
    int nz;
    
    int index;
    int indexa;
    
    int me;
    
    int istart;
    int iend;
    int totny;
    int nxtotny;
    int totsliceny;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    istart = 0;
    iend = ny;
    
    me = pfmdat->myrank;
    
    index = (i*nx + j);
    
    if ( i == 0 ) {
	gridNew[(i*nx + j)].phi = gridNew[((i+1)*nx + j)].phi;
    }
    
    if ( i == (iend-1) ) {
	gridNew[(i*nx + j)].phi = gridNew[((i-1)*nx + j)].phi;
    }
    
    if ( j == 0 ) {
        gridNew[(i*nx + j)].phi = gridNew[(i*nx + (j+1))].phi;
    }
    
    if ( j == (nx-1) ) {
        gridNew[(i*nx + j)].phi = gridNew[(i*nx + (j-1))].phi;
    }
    
}

__kernel void apply_BC_com(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    
    int i;
    int j;
    int k;
    int nx;
    int ny;
    int nz;
    
    int index;
    int indexa;
    
    int me;
    
    int istart;
    int iend;
    int totny;
    int nxtotny;
    int totsliceny;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    istart = 0;
    iend = ny;
    
    me = pfmdat->myrank;
    
    index = (i*nx + j);
    
      if ( i == 0 ) {
	gridNew[(i*nx + j)].c1 = gridNew[((i+1)*nx + j)].c1;
  gridNew[(i*nx + j)].mu[0] = gridNew[((i+1)*nx + j)].mu[0];
      }  
    
      if ( i == (iend-1) ) {
	gridNew[(i*nx + j)].c1 = gridNew[((i-1)*nx + j)].c1;
  gridNew[(i*nx + j)].mu[0] = gridNew[((i-1)*nx + j)].mu[0];
      }  
    
    if ( j == 0 ) {
        gridNew[(i*nx + j)].c1 = gridNew[(i*nx + (j+1))].c1;
        gridNew[(i*nx + j)].mu[0] = gridNew[(i*nx + (j+1))].mu[0];
    }
    
    if ( j == (nx-1) ) {
        gridNew[(i*nx + j)].c1 = gridNew[(i*nx + (j-1))].c1;
        gridNew[(i*nx + j)].mu[0] = gridNew[(i*nx + (j-1))].mu[0];
    }
    
}

__kernel void apply_BC_it_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;

    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == 0 ) {
      gridNew[(i*nx + j)] = gridNew[((i+1)*nx + j)];
    }
    
}

__kernel void apply_BC_ib_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == (ny-1) ) {
      gridNew[(i*nx + j)] = gridNew[((i-1)*nx + j)];  
    }
    
}

__kernel void apply_BC_it_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == 0 ) {
      gridNew[(i*nx + j)] = gridNew[((ny-2)*nx + j)];
    }
    
}

__kernel void apply_BC_ib_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == (ny-1) ) {
      gridNew[(i*nx + j)] = gridNew[((1)*nx + j)];  
    }
    
}

__kernel void apply_BC_jl_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == 0 ) {
        gridNew[(i*nx + j)] = gridNew[(i*nx + (j+1))];
    }
    
}

__kernel void apply_BC_jr_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == (nx-1) ) {
        gridNew[(i*nx + j)] = gridNew[(i*nx + (j-1))];
    }
    
}

__kernel void apply_BC_jl_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == 0 ) {
        gridNew[(i*nx + j)] = gridNew[(i*nx + (nx-2))];
    }
    
}

__kernel void apply_BC_jr_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == (nx-1) ) {
        gridNew[(i*nx + j)] = gridNew[(i*nx + (1))];
    }
    
}

__kernel void apply_BC_phi_it_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == 0 ) {
	gridNew[(i*nx + j)].phi = gridNew[((i+1)*nx + j)].phi;
    }
    
}

__kernel void apply_BC_phi_ib_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == (ny-1) ) {
	gridNew[(i*nx + j)].phi = gridNew[((i-1)*nx + j)].phi;
    }
    
}

__kernel void apply_BC_phi_it_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == 0 ) {
	gridNew[(i*nx + j)].phi = gridNew[((ny-2)*nx + j)].phi;
    }
    
}

__kernel void apply_BC_phi_ib_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == (ny-1) ) {
	gridNew[(i*nx + j)].phi = gridNew[((1)*nx + j)].phi;
    }
    
}

__kernel void apply_BC_phi_jl_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == 0 ) {
        gridNew[(i*nx + j)].phi = gridNew[(i*nx + (j+1))].phi;
    }
    
}

__kernel void apply_BC_phi_jr_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == (nx-1) ) {
        gridNew[(i*nx + j)].phi = gridNew[(i*nx + (j-1))].phi;
    }
    
}

__kernel void apply_BC_phi_jl_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == 0 ) {
        gridNew[(i*nx + j)].phi = gridNew[(i*nx + (nx-2))].phi;
    }
    
}

__kernel void apply_BC_phi_jr_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == (nx-1) ) {
        gridNew[(i*nx + j)].phi = gridNew[(i*nx + (1))].phi;
    }
    
}

__kernel void apply_BC_com_it_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == 0 ) {
      gridNew[(i*nx + j)].c1 = gridNew[((i+1)*nx + j)].c1;
      gridNew[(i*nx + j)].mu[0] = gridNew[((i+1)*nx + j)].mu[0];
    }
    
}

__kernel void apply_BC_com_ib_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == (ny-1) ) {
      gridNew[(i*nx + j)].c1 = gridNew[((i-1)*nx + j)].c1;
      gridNew[(i*nx + j)].mu[0] = gridNew[((i-1)*nx + j)].mu[0];
    }
    
}

__kernel void apply_BC_com_it_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == 0 ) {
      gridNew[(i*nx + j)].c1 = gridNew[((ny-2)*nx + j)].c1;
      gridNew[(i*nx + j)].mu[0] = gridNew[((ny-2)*nx + j)].mu[0];
    }
    
}

__kernel void apply_BC_com_ib_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == (ny-1) ) {
      gridNew[(i*nx + j)].c1 = gridNew[((1)*nx + j)].c1;
      gridNew[(i*nx + j)].mu[0] = gridNew[((1)*nx + j)].mu[0];
    }
    
}

__kernel void apply_BC_com_jl_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == 0 ) {
        gridNew[(i*nx + j)].c1 = gridNew[(i*nx + (j+1))].c1;
        gridNew[(i*nx + j)].mu[0] = gridNew[(i*nx + (j+1))].mu[0];
    }
    
}

__kernel void apply_BC_com_jr_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == (nx-1) ) {
        gridNew[(i*nx + j)].c1 = gridNew[(i*nx + (j-1))].c1;
        gridNew[(i*nx + j)].mu[0] = gridNew[(i*nx + (j-1))].mu[0];
    }
    
}

__kernel void apply_BC_com_jl_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == 0 ) {
        gridNew[(i*nx + j)].c1 = gridNew[(i*nx + (nx-2))].c1;
        gridNew[(i*nx + j)].mu[0] = gridNew[(i*nx + (nx-2))].mu[0];
    }
    
}

__kernel void apply_BC_com_jr_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == (nx-1) ) {
        gridNew[(i*nx + j)].c1 = gridNew[(i*nx + (1))].c1;
        gridNew[(i*nx + j)].mu[0] = gridNew[(i*nx + (1))].mu[0];
    }
    
}

__kernel void addNoise(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {

    int i;
    int j;
    int k;
    int nx;
    int ny;
    int nz;
    int index;

    j = get_global_id(0);
    i = get_global_id(1);

    nx = pfmdat->jNx;
    ny = pfmdat->iNy;

    index = (i*nx + j);

    double minr = -1.0;
    double maxr = 1.0;
    double rt, noise;

    if ( (gridNew[index].phi > 0.1) && (gridNew[index].phi < 0.5) ) {

        unsigned long x = get_local_id(0);
        unsigned long y = get_local_id(1);
        unsigned long z = get_global_id(0);
        unsigned long w = get_global_id(1);

        ulong t = x ^ (x << 11);
        x = y; 
        y = z; 
        z = w;
        w=index;
        w = w ^ (w >> 51) ^ (t ^ (t >> 8));
        rt = w%2222;
        noise = minr + rt/3333 * (maxr - minr + 1);

        gridNew[index].phi = gridNew[index].phi - ( pfmdat->NoiseFac * noise );

    }

    if (gridNew[index].phi < 0.0) {
        gridNew[index].phi = 0.0;
    }
    if (gridNew[index].phi > 1.0) {
        gridNew[index].phi = 1.0;
    }

}

__kernel void copy_New_To_Old(__global struct grid *gridOld, __global struct grid *gridNew, __constant struct pfmval *pfmdat) {

    int i;
    int j;
    int k;
    int nx;
    int ny;
    int nz;
    int index;

    j = get_global_id(0);
    i = get_global_id(1);

    nx = pfmdat->jNx;
    ny = pfmdat->iNy;

    index = (i*nx + j);

    gridOld[index] = gridNew[index];

}

__kernel void update_temp_UC(__global double *temp, __global int *tstep, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar) {

    int i;
    int j;
    int k;
    int nx;
    int ny;
    int nz;
    int index;

    j = get_global_id(0);
    i = get_global_id(1);

    nx = pfmdat->jNx;
    ny = pfmdat->iNy;

    temp[j] = pfmdat->T0;
    

}

__kernel void update_temp_DS(__global double *temp, __global int *tstep, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar) {

    int i;
    int j;
    int k;
    int nx;
    int ny;
    int nz;
    int index;
    int i0;
    int j0;

    j = get_global_id(0);
    i = get_global_id(1);

    nx = pfmdat->jNx;
    ny = pfmdat->iNy;

    //i0 = i + (pfmdat->myrank*(ny-2))-1; 
    //temp[i] = pfmdat->Toffset + pfmdat->TG*((i0-pfmdat->PosOffset)*pfmvar->dx - pfmdat->Vp*tstep[0]*pfmvar->dt);

    //i0 = i + (pfmdat->myrank*(ny-2));
    //temp[i] = pfmdat->Toffset + pfmdat->TGRADIENT*((i0-pfmdat->TPosOffset+pfmdat->shift_OFFSET)*pfmvar->deltax-(pfmdat->velocity*tstep[0]*pfmvar->deltat));

    j0 = j + (pfmdat->myrank*(nx-2));
    temp[j] = pfmdat->Toffset + pfmdat->TGRADIENT*((j0-pfmdat->TPosOffset+pfmdat->shift_OFFSET)*pfmvar->deltax-(pfmdat->velocity*tstep[0]*pfmvar->deltat));

}
  
  
__kernel void update_temp_CR(__global double *temp, __global int *tstep, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar) { 
    
    int i;
    int j;
    int k;
    int nx;
    int ny;
    int nz;
    int index;

    j = get_global_id(0);
    i = get_global_id(1);

    nx = pfmdat->jNx;
    ny = pfmdat->iNy;

    //	temp[i] = pfmdat->Tempstart - pfmdat->CR*tstep[0]*pfmvar->dt;

}
  
__kernel void apply_BC_temp_it_noflux(__global double *temp, __constant struct pfmval *pfmdat) {

    int i;
    int j;
    int k;
    int nx;
    int ny;
    int nz;
    int index;
    int i0;

    j = get_global_id(0);
    i = get_global_id(1);

    nx = pfmdat->jNx;
    ny = pfmdat->iNy;

    if ( i == 0 ) { 
        temp[i] = temp[i+1];
    }

}
  
__kernel void apply_BC_temp_ib_noflux(__global double *temp, __constant struct pfmval *pfmdat) {

    int i;
    int j;
    int k;
    int nx;
    int ny;
    int nz;
    int index;
    int i0;

    j = get_global_id(0);
    i = get_global_id(1);

    nx = pfmdat->jNx;
    ny = pfmdat->iNy;

    if ( i == ny-1 ) { 
        temp[i] = temp[i-1];
    }

}
