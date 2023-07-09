#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#include "solverloop/defines.h"
#include "solverloop/defines1.h"
#include "solverloop/CL_struct_var_kernels.h"
#include "solverloop/ThermoCL_2.h"
#include "solverloop/ThermoCL_2.c"
//#include "solverloop/GibbsEnergyData1.h"
#include "solverloop/GibbsEnergyData2.h"
#include "solverloop/MatrixInversion.h"
#include "solverloop/matrix_nsolCL.h"
#include "solverloop/FunctionF_2_c_mu_NR_GEData2.h"
#include "solverloop/functionH_CL_5th.h"
#include "solverloop/functionW_02_CL.h"
#include "solverloop/anisotropy_01_CL.h"
#include "solverloop/functionA_01_CL.h"
#include "solverloop/functionF_03_CL.h" 


__kernel void SolverCsClEq_F2(__global struct fields *gridinfoO, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf3 *propf3) {

    int x, y, z;
    int nx, ny, nz;
    int index;
    int is, ip, is1, is2;

    int interface, bulkphase;
    double cg[nsol], mu[nsol], cphas[nsol];
    double Ti;
    double tmp0;

    y = get_global_id(0);
    z = get_global_id(1);
    x = get_global_id(2);

    ny = get_global_size(0);
    nz = get_global_size(1);
    nx = get_global_size(2);

    index = y + ny*(z + x*nz);

    Ti = pfmdat->T0;

    interface = 1;
    bulkphase = 0;
    for ( ip = 0; ip < npha; ip++ ) {
      if ( gridinfoO[index].phi[ip] >= pfmdat->interfaceUplimit ) {
        bulkphase = ip;
        interface = 0;
        break;
      }
    }

    if ( interface ) {

      for ( ip = 0; ip < npha; ip++ ) {

        for ( is = 0; is < nsol; is++ ) {

          cg[is] = pfmdat->cguess[ip*npha+ip][is];
          mu[is] = gridinfoO[index].mu[is];
          //printf("%d, %d, %d, %d, %d, %le, %le, %le, %le, %d\n", tstep[0], i, j, ip, is, pfmdat->cguess[ip][ip][is], cg[is], gridinfoO[index].mu[is], mu[is],pfmdat->thermophase[ip]);

        }

        c_mu(mu, cphas, Ti, ip, cg, tstep[0], y, z, x, pfmdat->thermophase[ip]);

        for ( is = 0; is < nsol; is++ ) {
          cscl[index].comie[ip][is]= cphas[is];
          //printf("%d, %d, %d, %d, %d, %le\n", tstep[0], i, j, ip, is, cscl[index].comie[ip][is]);
        }

      }

    }

}

__kernel void SolverCsClEq_F3(__global struct fields *gridinfoO, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf3 *propf3) {
    
    int x, y, z;
    int nx, ny, nz;
    int index;
    int is, ip, is1, is2; ;
    
    int interface, bulkphase; 
    double cg[nsol], ci[nsol];
    double Ti;
    double tmp0;
    
    double F3_Beq[npha*nsol], F3_dBbdT[npha*nsol], F3_cmu[npha*nsol*nsol], Tequ, mu[nsol];
    
    // id0 = get_global_id(0);
    // id1 = get_global_id(1);
    // id2 = get_global_id(2);
    // 
    // sz0 = get_global_size(0);
    // sz1 = get_global_size(1);
    // sz2 = get_global_size(2);
    // 
    // index = id0 + sz0*(id1 + sz1*id2);
    
    y = get_global_id(0);
    z = get_global_id(1);
    x = get_global_id(2);
    
    ny = get_global_size(0);
    nz = get_global_size(1);
    nx = get_global_size(2);
    
    index = y + ny*(z + x*nz);
    
    Ti = pfmdat->T0;
    
    interface = 1; 
    bulkphase = 0;
    for ( ip = 0; ip < npha; ip++ ) {
      if ( gridinfoO[index].phi[ip] >= pfmdat->interfaceUplimit ) { 
        bulkphase = ip;
        interface = 0;
        break;
      }
    }
    
    for ( ip = 0; ip < npha; ip++ ) { 
      
      for ( is1 = 0; is1 < nsol; is1++ ) { 
        
        F3_Beq[ip*nsol+is1] = propf3->Beq[ip][is];
        F3_dBbdT[ip*nsol+is1] = propf3->dBbdT[ip][is1]; 
        
        cg[(ip*npha+ip)*nsol+is1] = pfmdat->cguess[ip*npha+ip][is1];
        
        for ( is2 = 0; is2 < nsol; is2++ ) { 
          F3_cmu[(ip*nsol+is1)*nsol+is2] = propf3->cmu[ip][is1][is2];
        }
        
      }
      
    } 
    
    Tequ = pfmdat->Teq;
    
    
    if ( interface ) {
      
      for ( is = 0; is < nsol; is++ ) { 
        mu[is] = gridinfoO[index].mu[is];
      }
      
      for ( ip = 0; ip < npha; ip++ ) {
        
        function_F_03_c_mu_CL(mu, ci, Ti, ip, cg, F3_cmu, F3_Beq, F3_dBbdT, Tequ);
        
        for ( is = 0; is < nsol; is++ ) { 
          cscl[index].comie[ip][is] = ci[is];
        }
        
      }
      
    }
    
 
}

__kernel void SolverCsClEq_F4(__global struct fields *gridinfoO, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf4 *propf4, __constant struct propmatf4spline *propf4spline) {

    int x, y, z;
    int nx, ny, nz;
    int index;
    int is, ip, is1, is2, js;

    int interface, bulkphase;
    double cg[nsol], ci[nsol];
    double cmu[nsol*nsol], muc[nsol*nsol];
    double Ti;
    double tmp0;

    y = get_global_id(0);
    z = get_global_id(1);
    x = get_global_id(2);

    ny = get_global_size(0);
    nz = get_global_size(1);
    nx = get_global_size(2);

    index = y + ny*(z + x*nz);

    Ti = pfmdat->T0;

    interface = 1;
    bulkphase = 0;
    for ( ip = 0; ip < npha; ip++ ) {
      if ( gridinfoO[index].phi[ip] >= pfmdat->interfaceUplimit ) {
        bulkphase = ip;
        interface = 0;
        break;
      }
    }
    //printf("! %d, %d, %d, %d, %le, %le, %le, %le, %le, %le\n", tstep[0], i, j, ip, cscl[index].comie[ip][is], propf3->cmu[ip][is1][is], gridinfoO[4].mu[is], propf3->Beq[ip][is], Ti, pfmdat->Teq);
    if ( interface ) {

      for ( ip = 0; ip < npha; ip++ ) {

        // if ( pfmdat->ISOTHERMAL ) {

          for ( is1 = 0; is1 < nsol; is1++ ) {
            tmp0 = 0.0;
            for ( is2 = 0; is2 < nsol; is2++ ) {
              tmp0 += propf4->cmu[ip][is1][is2] * ( gridinfoO[index].mu[is2] - propf4->B[ip][is2] );
            }
            ci[is1] = tmp0;
          }

        // }
        // else {
        //   for ( is = 0; is < nsol; is++ ) {
        //     for ( js = 0; js < nsol; js++ ) {
        //       if ( is == js ) {
        //         muc[is*nsol+js] = 2.0 * propf4spline[j].A[ip][is][js];
        //       }
        //       else {
        //         muc[is*nsol+js] = propf4spline[j].A[ip][is][js];
        //       }
        //     }
        //   }
        //
        //   matinvnew_nsol(muc, cmu);
        //
        //   for ( is1 = 0; is1 < nsol; is1++ ) {
        //     tmp0 = 0.0;
        //     for ( is2 = 0; is2 < nsol; is2++ ) {
        //       tmp0 += cmu[is1*nsol+is2] * ( gridinfoO[index].mu[is2] - propf4spline[j].B[ip][is2] );
        //     }
        //     ci[is1] = tmp0;
        //   }
        // }

        for ( is = 0; is < nsol; is++ ) {
          cscl[index].comie[ip][is] = ci[is];
          //printf("I %d, %d, %d, %d, %le, %le, %le\n", tstep[0], i, j, ip, cscl[index].comie[ip][is], gridinfoO[4].mu[is], gridinfoO[4].com[is]);
          //printf("I %d, %d, %d, %d, %le, %le, %le, %le, %le, %le\n", tstep[0], i, j, ip, cscl[index].comie[ip][is], propf3->cmu[ip][is1][is], gridinfoO[4].mu[is], propf3->Beq[ip][is], Ti, pfmdat->Teq);
        }
      }

    }
}

__kernel void SolverPhi_F2_smooth(__global struct fields *gridinfoO, __global struct fields *gridinfo, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf3 *propf3) {

  int x, y, z, ii, i1, j1, k1, i2;
  int nx, ny, nz;
  int index;

  int is, js, ks, il, il1, jl, jl1, ipx;
  int is1, is2;
  int ig, jg, ip, ip1, ip2, ip3;
  int interface,bulkphase;
  int activephases;

  double Ti;

  double tmp0, tmp1;

  double yc[nsol+1], retdmuphase[nsol*nsol], retGphase[1];
  double ge[npha], fhi[npha], psi;
  double dgedc[npha][nsol];
  double ddgedcdc[npha][nsol*nsol], InvD[npha][nsol*nsol], deltac[nsol], Ddc[nsol], ddgDdc[nsol];
  double cie[npha][nsol], min_tau, tau;
  double lap[npha], dff[npha], phidot, Tau[npha][npha];
  double zeta[npha][npha], TAU, df[npha], hphi[npha];
  double sumdf;
  double retG[npha], retmu[npha*nsol], retdmu[npha*nsol*nsol];
  double gphi[npha], dgphidphi[npha], dhphidphi[npha], mu[nsol], dhphi[npha][npha];
  double dwh[npha*npha], dwhabc[npha*npha*npha], laps[npha];
  double sumlambdaphi, lambaphi[npha];
  double dab[npha*npha], ee_ab[npha*npha];
  //double facecenter, edgecenter, corner;

  struct fields stgridO[27];

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  //index = y + ny*(z + x*nz);

  Ti = pfmdat->T0;

  if ( x != 0 && x != nx-1 && y != 0 && y != ny-1 && z != 0 && z != nz-1 ) {

    index = y + ny*(z + x*nz);

    stgridO[0]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[1]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[2]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[3]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[4]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[5]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[6]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[7]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[8]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[9]  = gridinfoO[ (y-1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[10] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[11] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[12] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[13] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[14] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[15] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[16] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[17] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[18] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[19] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[20] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[21] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[22] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[23] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[24] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[25] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[26] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z+1) ) ];

    for ( ip = 0; ip < npha; ip++ ) {
      // fhi[ip] = stgridO[13].phi[ip];
      fhi[ip] = gridinfoO[index].phi[ip];
    }

    // Calculate laplacian
    for ( ip = 0; ip < npha; ip++ ) {

      lap[ip]  = 0.0;

      // facecenter = stgridO[4].phi[ip] + stgridO[22].phi[ip] + stgridO[10].phi[ip] + stgridO[16].phi[ip] + stgridO[12].phi[ip] + stgridO[14].phi[ip];

      // edgecenter = stgridO[1].phi[ip] + stgridO[5].phi[ip] + stgridO[7].phi[ip] + stgridO[3].phi[ip] + stgridO[9].phi[ip] + stgridO[11].phi[ip] + stgridO[17].phi[ip] + stgridO[15].phi[ip] + stgridO[19].phi[ip] + stgridO[21].phi[ip] + stgridO[25].phi[ip] + stgridO[23].phi[ip];

      // corner = stgridO[0].phi[ip] + stgridO[2].phi[ip] + stgridO[8].phi[ip] + stgridO[6].phi[ip] + stgridO[18].phi[ip] + stgridO[20].phi[ip] + stgridO[26].phi[ip] + stgridO[24].phi[ip];

      // lap[ip] = 3.0*( facecenter + edgecenter/2.0 + corner/3.0 - 44.0*stgridO[13].phi[ip]/3.0 ) / (13.0 * (pfmvar->deltax) * (pfmvar->deltax));

      lap[ip]  = ( stgridO[14].phi[ip] - 2.0*stgridO[13].phi[ip] + stgridO[12].phi[ip] ) / (pfmvar->deltay*pfmvar->deltay);
      lap[ip] += ( stgridO[16].phi[ip] - 2.0*stgridO[13].phi[ip] + stgridO[10].phi[ip] ) / (pfmvar->deltaz*pfmvar->deltaz);
      lap[ip] += ( stgridO[22].phi[ip] - 2.0*stgridO[13].phi[ip] + stgridO[ 4].phi[ip] ) / (pfmvar->deltax*pfmvar->deltax);



      //printf("%le\n", lap[ip]);
    }

    // Get phase compositions
    for ( ip = 0; ip < npha; ip++ ) {
      for ( is = 0; is < nsol; is++ ) {
        cie[ip][is] = cscl[index].comie[ip][is];
      }
    }

    for ( ip = 0; ip < npha; ip++ ) {

      tmp0 = 0.0;
      for ( is = 0; is < nsol; is++ ) {
        yc[is] = cie[ip][is];
        tmp0 += yc[is];
      }
      yc[nsol] = 1.0 - tmp0;

      // Calculate free energy
      Ge(Ti, yc, retGphase, pfmdat->thermophase[ip]);

      ge[ip] = retGphase[0] / pfmdat->Vm;

      // Calculate second derivative of free energy .i.e ddgedcdc and get diffusion coefficient to local variable
      dMudc(Ti, yc, retdmuphase, pfmdat->thermophase[ip]);

      for ( is1 = 0; is1 < nsol; is1++ ) {
        for ( is2 = 0; is2 < nsol; is2++ ) {
          ddgedcdc[ip][is1*nsol+is2] = retdmuphase[is1*nsol+is2] / pfmdat->Vm;
          InvD[ip][is1*nsol+is2] = pfmdat->DInv[ip][is1*nsol+is2];
          //printf("%d, %d, %d, %le, %le, %le, %le, %le\n", tstep[0], i, j, ddgedcdc[ip][is1*nsol+is2], retdmuphase[is1*nsol+is2], InvD[ip][is1*nsol+is2], pfmdat->DInv[ip][is1][is2], y[is2]);
        }
      }

    }

    for ( is = 0; is < nsol; is++ ) {
      //mu[is] = stgridO[13].mu[is] / pfmdat->Vm;
      mu[is] = gridinfoO[index].mu[is] / pfmdat->Vm;
    }

    // Calculate TAU
    for ( ip1 = 0; ip1 < npha-1; ip1++ ) {
      for ( is = 0; is < nsol; is++ ) {
        deltac[is] = pfmdat->c_eq[(npha-1)*npha+(npha-1)][is] - pfmdat->c_eq[ip1*npha+ip1][is];
      }
      multiply_nsol(InvD[npha-1], deltac, Ddc);
      multiply_nsol(ddgedcdc[npha-1], Ddc, ddgDdc);

      tmp0 = 0.0;
      for ( is = 0; is < nsol; is++ ) {
        tmp0 += ddgDdc[is] * deltac[is];
      }

      zeta[ip1][npha-1] = tmp0;

      Tau[ip1][npha-1] = ( 6.0 * pfmvar->ee[ip1*npha+(npha-1)] * pfmdat->a2 * zeta[ip1][npha-1] ) / pfmvar->w[ip1*npha+(npha-1)];
      Tau[npha-1][ip1] = Tau[ip1][npha-1];

      if ( ip1 == 0 ) {
        min_tau = Tau[ip1][npha-1];
      }
      if ( Tau[ip1][npha-1] < min_tau ) {
        min_tau = Tau[ip1][npha-1];
      }
    }
    for ( ip1 = 0; ip1 < npha-1; ip1++ ) {
      for ( ip2 = 0; ip2 < npha-1; ip2++ ) {
        Tau[ip1][ip2] = min_tau;
      }
    }
    tau = min_tau;

    tmp0 = 0.0;
    tmp1 = 0.0;
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      for ( ip2 = 0; ip2 < npha; ip2++ ) {
        if ( ip1 < ip2 ) {
          tmp0 += Tau[ip1][ip2]*fhi[ip1]*fhi[ip2];
          tmp1 += fhi[ip1]*fhi[ip2];
        }
      }
    }
    if ( tmp1 ) {
      TAU = tmp0 / tmp1;
    } else {
      TAU = tau;
    }

    // get w_ab w_abc ee dab ti local linear indexed variable
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      for ( ip2 = 0; ip2 < npha; ip2++ ) {
        dwh[ip1*npha+ip2] = pfmvar->w[ip1*npha+ip2];
        dab[ip1*npha+ip2] = pfmdat->epsc[ip1*npha+ip2];
        ee_ab[ip1*npha+ip2] = pfmvar->ee[ip1*npha+ip2];
        for ( ip3 = 0; ip3 < npha; ip3++ ) {
          dwhabc[(ip1*npha+ip2)*npha+ip3] = pfmvar->w_abc[(ip1*npha+ip2)*npha+ip3];
        }
      }
    }

    // Calculate derivative of well potential
    for ( ip = 0; ip < npha; ip++ ) {
      dgphidphi[ip] = F_W_02_dwdphi(fhi, lap, dwh, dwhabc, ip);
    }

    // Get the laplacian for all betas
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      laps[ip1] = 0.0;
      for ( ip2 = 0; ip2 < npha; ip2++ ) {
        if ( ip1 != ip2 ) {
          laps[ip1] += pfmvar->ee[ip1*npha+ip2] * lap[ip2];
        }
      }
      // printf("    =%le\n", laps[ip]);
    }

    // Calculate -dF/dphi for all betas and check for active phases
    sumlambdaphi = 0.0;
    activephases = 0;
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      if ( fabs(lap[ip1]) > 0.0 ) {
         lambaphi[ip1] = -laps[ip1] - dgphidphi[ip1];
        sumlambdaphi += lambaphi[ip1];
        activephases++;
      }
    }

    // Divide the lambda with active phases
    if ( activephases ) {
      sumlambdaphi /= activephases;
    }

    // Update the phi by applying lagrangian multiplier
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      if ( fabs(lap[ip1]) > 0.0 ) {
        phidot = (2.0/(TAU)) * ( lambaphi[ip1] - sumlambdaphi );
        gridinfo[ index ].phi[ip1] = gridinfoO[index].phi[ip1] + (pfmvar->deltat)*phidot;

        //printf("t=%d, i=%d, j=%d, lap=%le,%le, mphi=%le, vard=%le, phi[0]=%le, ip=%d, phiold=%le\n", tstep[0], i, j, lap[0],lap[1],  1.0/TAU, 2.0*( lambaphi[ip1] - sumlambdaphi ), gridinfo[ index ].phi[ip1], ip1, stgridO[4].phi[ip1]);
      }
    }

  }

}

__kernel void SolverPhi_F2(__global struct fields *gridinfoO, __global struct fields *gridinfo, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf3 *propf3) {

  int x, y, z, ii, i1, j1, k1, i2;
  int nx, ny, nz;
  int index;

  int is, js, ks, il, il1, jl, jl1, ipx;
  int is1, is2;
  int ig, jg, ip, ip1, ip2, ip3;
  int ilmax;
  int interface,bulkphase;
  int activephases;

  double Ti;

  double tmp0, tmp1;

  double yc[nsol+1], retdmuphase[nsol*nsol], retGphase[1];
  double ge[npha], fhi[npha], psi;
  double dgedc[npha][nsol];
  double ddgedcdc[npha][nsol*nsol], InvD[npha][nsol*nsol], deltac[nsol], Ddc[nsol], ddgDdc[nsol];
  double cie[npha][nsol], min_tau, tau;
  double lap[npha], dff[npha], phidot, Tau[npha][npha];
  double zeta[npha][npha], TAU, df[npha], hphi[npha];
  double sumdf;
  double retG[npha], retmu[npha*nsol], retdmu[npha*nsol*nsol];
  double gphi[npha], dgphidphi[npha], dhphidphi[npha], mu[nsol], dhphi[npha][npha];
  double dwh[npha*npha], dwhabc[npha*npha*npha], laps[npha];
  double sumlambdaphi, lambaphi[npha];
  double stphi[npha*3*3*3], dab[npha*npha], ee_ab[npha*npha];
  double Rot_mat[npha*npha*3*3], Inv_Rot_mat[npha*npha*3*3];
  //double facecenter, edgecenter, corner;

  struct fields stgridO[27];

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  //index = y + ny*(z + x*nz);

  Ti = pfmdat->T0;

  if ( x != 0 && x != nx-1 && y != 0 && y != ny-1 && z != 0 && z != nz-1 ) {

    index = y + ny*(z + x*nz);

    stgridO[0]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[1]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[2]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[3]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[4]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[5]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[6]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[7]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[8]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[9]  = gridinfoO[ (y-1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[10] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[11] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[12] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[13] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[14] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[15] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[16] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[17] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[18] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[19] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[20] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[21] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[22] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[23] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[24] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[25] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[26] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z+1) ) ];

    for (ip = 0; ip < npha; ip++) {
      ii = 0;
        for (i1 = 0; i1 < 3; i1++) {
          for (j1 = 0; j1 < 3; j1++) {
            for (k1 = 0; k1 < 3; k1++) {
              //stphi[ (ip*3 + i1)*3 + j1 ] = stphi_1[ip][i1][j1];
              stphi[ ( (ip*3 + i1)*3 + j1 )*3 + k1 ] = stgridO[ii].phi[ip];
              ii++;
            }
          }
        }
    }


    for ( ip = 0; ip < npha; ip++ ) {
      // fhi[ip] = stgridO[13].phi[ip];
      fhi[ip] = gridinfoO[index].phi[ip];
    }

    // Calculate laplacian
    for ( ip = 0; ip < npha; ip++ ) {

      lap[ip]  = 0.0;

      // facecenter = stgridO[4].phi[ip] + stgridO[22].phi[ip] + stgridO[10].phi[ip] + stgridO[16].phi[ip] + stgridO[12].phi[ip] + stgridO[14].phi[ip];

      // edgecenter = stgridO[1].phi[ip] + stgridO[5].phi[ip] + stgridO[7].phi[ip] + stgridO[3].phi[ip] + stgridO[9].phi[ip] + stgridO[11].phi[ip] + stgridO[17].phi[ip] + stgridO[15].phi[ip] + stgridO[19].phi[ip] + stgridO[21].phi[ip] + stgridO[25].phi[ip] + stgridO[23].phi[ip];

      // corner = stgridO[0].phi[ip] + stgridO[2].phi[ip] + stgridO[8].phi[ip] + stgridO[6].phi[ip] + stgridO[18].phi[ip] + stgridO[20].phi[ip] + stgridO[26].phi[ip] + stgridO[24].phi[ip];

      // lap[ip] = 3.0*( facecenter + edgecenter/2.0 + corner/3.0 - 44.0*stgridO[13].phi[ip]/3.0 ) / (13.0 * (pfmvar->deltax) * (pfmvar->deltax));

      lap[ip]  = ( stgridO[14].phi[ip] - 2.0*stgridO[13].phi[ip] + stgridO[12].phi[ip] ) / (pfmvar->deltay*pfmvar->deltay);
      lap[ip] += ( stgridO[16].phi[ip] - 2.0*stgridO[13].phi[ip] + stgridO[10].phi[ip] ) / (pfmvar->deltaz*pfmvar->deltaz);
      lap[ip] += ( stgridO[22].phi[ip] - 2.0*stgridO[13].phi[ip] + stgridO[ 4].phi[ip] ) / (pfmvar->deltax*pfmvar->deltax);



      //printf("%le\n", lap[ip]);
    }

    // Get phase compositions
    for ( ip = 0; ip < npha; ip++ ) {
      for ( is = 0; is < nsol; is++ ) {
        cie[ip][is] = cscl[index].comie[ip][is];
      }
    }

    for ( ip = 0; ip < npha; ip++ ) {

      tmp0 = 0.0;
      for ( is = 0; is < nsol; is++ ) {
        yc[is] = cie[ip][is];
        tmp0 += yc[is];
      }
      yc[nsol] = 1.0 - tmp0;

      // Calculate free energy
      Ge(Ti, yc, retGphase, pfmdat->thermophase[ip]);

      ge[ip] = retGphase[0] / pfmdat->Vm;

      // Calculate second derivative of free energy .i.e ddgedcdc and get diffusion coefficient to local variable
      dMudc(Ti, yc, retdmuphase, pfmdat->thermophase[ip]);

      for ( is1 = 0; is1 < nsol; is1++ ) {
        for ( is2 = 0; is2 < nsol; is2++ ) {
          ddgedcdc[ip][is1*nsol+is2] = retdmuphase[is1*nsol+is2] / pfmdat->Vm;
          InvD[ip][is1*nsol+is2] = pfmdat->DInv[ip][is1*nsol+is2];
          //printf("%d, %d, %d, %le, %le, %le, %le, %le\n", tstep[0], i, j, ddgedcdc[ip][is1*nsol+is2], retdmuphase[is1*nsol+is2], InvD[ip][is1*nsol+is2], pfmdat->DInv[ip][is1][is2], y[is2]);
        }
      }

    }

    for ( is = 0; is < nsol; is++ ) {
      //mu[is] = stgridO[13].mu[is] / pfmdat->Vm;
      mu[is] = gridinfoO[index].mu[is] / pfmdat->Vm;
    }

    // Calculate TAU
    for ( ip1 = 0; ip1 < npha-1; ip1++ ) {
      for ( is = 0; is < nsol; is++ ) {
        deltac[is] = pfmdat->c_eq[(npha-1)*npha+(npha-1)][is] - pfmdat->c_eq[ip1*npha+ip1][is];
      }
      multiply_nsol(InvD[npha-1], deltac, Ddc);
      multiply_nsol(ddgedcdc[npha-1], Ddc, ddgDdc);

      tmp0 = 0.0;
      for ( is = 0; is < nsol; is++ ) {
        tmp0 += ddgDdc[is] * deltac[is];
      }

      zeta[ip1][npha-1] = tmp0;

      Tau[ip1][npha-1] = ( 6.0 * pfmvar->ee[ip1*npha+(npha-1)] * pfmdat->a2 * zeta[ip1][npha-1] ) / pfmvar->w[ip1*npha+(npha-1)];
      Tau[npha-1][ip1] = Tau[ip1][npha-1];

      if ( ip1 == 0 ) {
        min_tau = Tau[ip1][npha-1];
      }
      if ( Tau[ip1][npha-1] < min_tau ) {
        min_tau = Tau[ip1][npha-1];
      }
    }
    for ( ip1 = 0; ip1 < npha-1; ip1++ ) {
      for ( ip2 = 0; ip2 < npha-1; ip2++ ) {
        Tau[ip1][ip2] = min_tau;
      }
    }
    tau = min_tau;

    tmp0 = 0.0;
    tmp1 = 0.0;
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      for ( ip2 = 0; ip2 < npha; ip2++ ) {
        if ( ip1 < ip2 ) {
          tmp0 += Tau[ip1][ip2]*fhi[ip1]*fhi[ip2];
          tmp1 += fhi[ip1]*fhi[ip2];
        }
      }
    }
    if ( tmp1 ) {
      TAU = tmp0 / tmp1;
    } else {
      TAU = tau;
    }


    // Calculate df_TD/dphi
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      dff[ip1] = 0.0;
      for ( ip2 = 0; ip2 < npha; ip2++ ) {

          psi = 0.0;
          psi += ge[ip2];
          for ( is = 0; is < nsol; is++ ) {
            psi -= mu[is]*cie[ip2][is];
          }
          psi *= dhfhi(fhi, ip2, ip1);

          dff[ip1] +=  psi;
      }
    }

    // get w_ab w_abc ee dab ti local linear indexed variable
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      for ( ip2 = 0; ip2 < npha; ip2++ ) {
        dwh[ip1*npha+ip2] = pfmvar->w[ip1*npha+ip2];
        dab[ip1*npha+ip2] = pfmdat->epsc[ip1*npha+ip2];
        ee_ab[ip1*npha+ip2] = pfmvar->ee[ip1*npha+ip2];
        for ( ip3 = 0; ip3 < npha; ip3++ ) {
          dwhabc[(ip1*npha+ip2)*npha+ip3] = pfmvar->w_abc[(ip1*npha+ip2)*npha+ip3];
        }
      }
    }

    // Calculate derivative of well potential
    for ( ip = 0; ip < npha; ip++ ) {
      dgphidphi[ip] = F_W_02_dwdphi(fhi, lap, dwh, dwhabc, ip);
    }

    // Get the laplacian for all betas
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      laps[ip1] = 0.0;
      for ( ip2 = 0; ip2 < npha; ip2++ ) {
        if ( ip1 != ip2 ) {
          laps[ip1] += pfmvar->ee[ip1*npha+ip2] * lap[ip2];
        }
      }
      // printf("    =%le\n", laps[ip]);
    }

    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      for ( ip2 = 0; ip2 < npha; ip2++ ) {
        for ( i1 = 0; i1 < 3; i1++ ) {
          for ( i2 = 0; i2 < 3; i2++ ) {
                Rot_mat[((ip1*npha+ip2)*npha+i1)*3+i2] = pfmdat->Rotation_matrix[ip1][ip2][i1][i2];
            Inv_Rot_mat[((ip1*npha+ip2)*npha+i1)*3+i2] = pfmdat->Inv_Rotation_matrix[ip1][ip2][i1][i2];
          }
        }
      }
    }

    for ( ip = 0; ip < npha; ip++ ) {
      //printf("%le, %le\n", laps[ip], calcAnisotropy_01(stphi, dab, ee_ab, ip, pfmvar->deltax, pfmvar->deltay, dz));
      laps[ip] = calcAnisotropy_01(stphi, dab, ee_ab, ip, pfmvar->deltax, pfmvar->deltay, pfmvar->deltaz, Rot_mat, Inv_Rot_mat, y, z, x);
    }

    // Calculate -dF/dphi for all betas and check for active phases
    sumlambdaphi = 0.0;
    activephases = 0;
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      if ( fabs(lap[ip1]) > 0.0 ) {
         lambaphi[ip1] = -laps[ip1] - dgphidphi[ip1] - dff[ip1];
        sumlambdaphi += lambaphi[ip1];
        activephases++;
      }
    }

    // Divide the lambda with active phases
    if ( activephases ) {
      sumlambdaphi /= activephases;
    }

    // Update the phi by applying lagrangian multiplier
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      if ( fabs(lap[ip1]) > 0.0 ) {
        phidot = (2.0/(TAU)) * ( lambaphi[ip1] - sumlambdaphi );
        gridinfo[ index ].phi[ip1] = gridinfoO[index].phi[ip1] + (pfmvar->deltat)*phidot;

        //printf("t=%d, i=%d, j=%d, lap=%le,%le, mphi=%le, vard=%le, phi[0]=%le, ip=%d, phiold=%le\n", tstep[0], i, j, lap[0],lap[1],  1.0/TAU, 2.0*( lambaphi[ip1] - sumlambdaphi ), gridinfo[ index ].phi[ip1], ip1, stgridO[4].phi[ip1]);
      }
    }

  }

}

__kernel void SolverPhi_F3_smooth(__global struct fields *gridinfoO, __global struct fields *gridinfo, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf3 *propf3) {

  int x, y, z, ii, i1, j1, k1, i2;
  int nx, ny, nz;
  int index;

  int is, js, ks, il, il1, jl, jl1, ipx;
  int is1, is2;
  int ig, jg, ip, ip1, ip2, ip3;
  int ilmax;
  int interface,bulkphase;
  int activephases;

  double Ti;

  double tmp0, tmp1, tmp2, tmp3, tmp4;

  double ge[npha], fhi[npha], psi, F3_A[npha*nsol*nsol], F3_B[npha*nsol], F3_C[npha], dcdmu[npha][nsol*nsol], F3_cmu[npha*nsol*nsol];
  double ddgedcdc[npha][nsol*nsol], InvD[npha][nsol*nsol], deltac[nsol], Ddc[nsol], ddgDdc[nsol];
  double cie[npha][nsol], min_tau, tau;
  double cdiff[npha][nsol], lap[npha], dff[npha], phidot, Tau[npha][npha];
  double zeta[npha][npha], TAU, df[npha], hphi[npha];
  double sumdf;
  double retG[npha], retmu[npha*nsol], retdmu[npha*nsol*nsol];
  double gphi[npha], dgphidphi[npha], dhphidphi[npha], mu[nsol], dhphi[npha][npha];
  double dwh[npha*npha], dwhabc[npha*npha*npha], laps[npha];
  double sumlambdaphi, lambaphi[npha];
  double dab[npha*npha], ee_ab[npha*npha];
  //double facecenter, edgecenter, corner;

  struct fields stgridO[27];

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  //index = y + ny*(z + x*nz);

  if ( x != 0 && x != nx-1 && y != 0 && y != ny-1 && z != 0 && z != nz-1 ) {

    index = y + ny*(z + x*nz);

    Ti = pfmdat->T0;

    stgridO[0]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[1]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[2]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[3]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[4]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[5]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[6]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[7]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[8]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[9]  = gridinfoO[ (y-1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[10] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[11] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[12] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[13] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[14] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[15] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[16] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[17] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[18] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[19] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[20] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[21] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[22] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[23] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[24] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[25] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[26] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z+1) ) ];

    for ( ip = 0; ip < npha; ip++ ) {
      // fhi[ip] = stgridO[13].phi[ip];
      fhi[ip] = gridinfoO[index].phi[ip];
    }

    // Calculate laplacian
    for ( ip = 0; ip < npha; ip++ ) {

      lap[ip]  = 0.0;

      // facecenter = stgridO[4].phi[ip] + stgridO[22].phi[ip] + stgridO[10].phi[ip] + stgridO[16].phi[ip] + stgridO[12].phi[ip] + stgridO[14].phi[ip];

      // edgecenter = stgridO[1].phi[ip] + stgridO[5].phi[ip] + stgridO[7].phi[ip] + stgridO[3].phi[ip] + stgridO[9].phi[ip] + stgridO[11].phi[ip] + stgridO[17].phi[ip] + stgridO[15].phi[ip] + stgridO[19].phi[ip] + stgridO[21].phi[ip] + stgridO[25].phi[ip] + stgridO[23].phi[ip];

      // corner = stgridO[0].phi[ip] + stgridO[2].phi[ip] + stgridO[8].phi[ip] + stgridO[6].phi[ip] + stgridO[18].phi[ip] + stgridO[20].phi[ip] + stgridO[26].phi[ip] + stgridO[24].phi[ip];

      // lap[ip] = 3.0*( facecenter + edgecenter/2.0 + corner/3.0 - 44.0*stgridO[13].phi[ip]/3.0 ) / (13.0 * (pfmvar->deltax) * (pfmvar->deltax));

      lap[ip]  = ( stgridO[14].phi[ip] - 2.0*stgridO[13].phi[ip] + stgridO[12].phi[ip] ) / (pfmvar->deltay*pfmvar->deltay);
      lap[ip] += ( stgridO[16].phi[ip] - 2.0*stgridO[13].phi[ip] + stgridO[10].phi[ip] ) / (pfmvar->deltaz*pfmvar->deltaz);
      lap[ip] += ( stgridO[22].phi[ip] - 2.0*stgridO[13].phi[ip] + stgridO[ 4].phi[ip] ) / (pfmvar->deltax*pfmvar->deltax);



      //printf("%le\n", lap[ip]);
    }



    // Get phase compositions
    for ( ip = 0; ip < npha; ip++ ) {
      for ( is = 0; is < nsol; is++ ) {
        cie[ip][is] = cscl[index].comie[ip][is];
      }
    }

    for ( ip = 0; ip < npha; ip++ ) {
      F3_C[ip] = propf3->C[ip];
      for ( is1 = 0; is1 < nsol; is1++ ) {
        F3_B[ip*nsol+is1] = propf3->B[ip][is1];
        for ( is2 = 0; is2 < nsol; is2++ ) {
          F3_A[(ip*nsol+is1)*nsol+is2] = propf3->A[ip][is1][is2];
          F3_cmu[(ip*nsol+is1)*nsol+is2] = propf3->cmu[ip][is1][is2];
        }
      }
    }

    // Get mu
    for ( is = 0; is < nsol; is++ ) {
      //mu[is] = stgridO[13].mu[is] / pfmdat->Vm;
      mu[is] = gridinfoO[index].mu[is] / pfmdat->Vm;
    }

    for ( ip = 0; ip < npha; ip++ ) {
      function_F_03_dc_dmu_CL(mu, cie[ip], Ti, ip, F3_cmu, dcdmu[ip]);
    }

    for ( ip = 0; ip < npha; ip++ ) {
      matinvnew_nsol(dcdmu[ip], ddgedcdc[ip]);
    }


    // Calculate second derivative of free energy .i.e ddgedcdc and get diffusion coefficient to local variable
    for ( ip = 0; ip < npha; ip++ ) {
      for ( is1 = 0; is1 < nsol; is1++ ) {
        for ( is2 = 0; is2 < nsol; is2++ ) {

          ddgedcdc[ip][is1*nsol+is2] /= pfmdat->Vm;

          InvD[ip][is1*nsol+is2] = pfmdat->DInv[ip][is1*nsol+is2];
        }
      }
    }

    // Calculate TAU
    for ( ip1 = 0; ip1 < npha-1; ip1++ ) {
      for ( is = 0; is < nsol; is++ ) {
        deltac[is] = pfmdat->c_eq[(npha-1)*npha+(npha-1)][is] - pfmdat->c_eq[ip1*npha+ip1][is];
      }
      multiply_nsol(InvD[npha-1], deltac, Ddc);
      multiply_nsol(ddgedcdc[npha-1], Ddc, ddgDdc);

      tmp0 = 0.0;
      for ( is = 0; is < nsol; is++ ) {
        tmp0 += ddgDdc[is] * deltac[is];
      }

      zeta[ip1][npha-1] = tmp0;

      Tau[ip1][npha-1] = ( 6.0 * pfmvar->ee[ip1*npha+(npha-1)] * pfmdat->a2 * zeta[ip1][npha-1] ) / pfmvar->w[ip1*npha+(npha-1)];
      Tau[npha-1][ip1] = Tau[ip1][npha-1];

      if ( ip1 == 0 ) {
        min_tau = Tau[ip1][npha-1];
      }
      if ( Tau[ip1][npha-1] < min_tau ) {
        min_tau = Tau[ip1][npha-1];
      }
    }
    for ( ip1 = 0; ip1 < npha-1; ip1++ ) {
      for ( ip2 = 0; ip2 < npha-1; ip2++ ) {
        Tau[ip1][ip2] = min_tau;
      }
    }
    tau = min_tau;

    tmp0 = 0.0;
    tmp1 = 0.0;
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      for ( ip2 = 0; ip2 < npha; ip2++ ) {
        if ( ip1 < ip2 ) {
          tmp0 += Tau[ip1][ip2]*fhi[ip1]*fhi[ip2];
          tmp1 += fhi[ip1]*fhi[ip2];
        }
      }
    }
    if ( tmp1 ) {
      TAU = tmp0 / tmp1;
    } else {
      TAU = tau;
    }

    // get w_ab w_abc ee dab ti local linear indexed variable
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      for ( ip2 = 0; ip2 < npha; ip2++ ) {
        dwh[ip1*npha+ip2] = pfmvar->w[ip1*npha+ip2];
        dab[ip1*npha+ip2] = pfmdat->epsc[ip1*npha+ip2];
        ee_ab[ip1*npha+ip2] = pfmvar->ee[ip1*npha+ip2];
        for ( ip3 = 0; ip3 < npha; ip3++ ) {
          dwhabc[(ip1*npha+ip2)*npha+ip3] = pfmvar->w_abc[(ip1*npha+ip2)*npha+ip3];
        }
      }
    }

    // Calculate derivative of well potential
    for ( ip = 0; ip < npha; ip++ ) {
      dgphidphi[ip] = F_W_02_dwdphi(fhi, lap, dwh, dwhabc, ip);
    }

    // Get the laplacian for all betas
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      laps[ip1] = 0.0;
      for ( ip2 = 0; ip2 < npha; ip2++ ) {
        if ( ip1 != ip2 ) {
          laps[ip1] += pfmvar->ee[ip1*npha+ip2] * lap[ip2];
        }
      }
      // printf("    =%le\n", laps[ip]);
    }

    // Calculate -dF/dphi for all betas and check for active phases
    sumlambdaphi = 0.0;
    activephases = 0;
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      if ( fabs(lap[ip1]) > 0.0 ) {
         lambaphi[ip1] = -laps[ip1] - dgphidphi[ip1];
        sumlambdaphi += lambaphi[ip1];
        activephases++;
      }
    }

    // Divide the lambda with active phases
    if ( activephases ) {
      sumlambdaphi /= activephases;
    }

    // Update the phi by applying lagrangian multiplier
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      if ( fabs(lap[ip1]) > 0.0 ) {
        phidot = (2.0/(TAU)) * ( lambaphi[ip1] - sumlambdaphi );
        gridinfo[ index ].phi[ip1] = gridinfoO[index].phi[ip1] + (pfmvar->deltat)*phidot;

        //printf("t=%d, i=%d, j=%d, lap=%le,%le, mphi=%le, vard=%le, phi[0]=%le, ip=%d, phiold=%le\n", tstep[0], i, j, lap[0],lap[1],  1.0/TAU, 2.0*( lambaphi[ip1] - sumlambdaphi ), gridinfo[ index ].phi[ip1], ip1, stgridO[4].phi[ip1]);
      }
    }

  }

}

__kernel void SolverPhi_F3(__global struct fields *gridinfoO, __global struct fields *gridinfo, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf3 *propf3) {
  
  int x, y, z, ii, i1, j1, k1, i2;
  int nx, ny, nz;
  int index;

  int is, js, ks, il, il1, jl, jl1, ipx;
  int is1, is2;
  int ig, jg, ip, ip1, ip2, ip3;
  int ilmax;
  int interface,bulkphase;
  int activephases;

  double Ti;

  double tmp0, tmp1;
  
  double ge[npha], fhi[npha], psi, F3_A[npha*nsol*nsol], F3_B[npha*nsol], F3_C[npha], dcdmu[npha][nsol*nsol], F3_cmu[npha*nsol*nsol];
  double ddgedcdc[npha][nsol*nsol], InvD[npha][nsol*nsol], deltac[nsol], Ddc[nsol], ddgDdc[nsol];
  double cie[npha][nsol], min_tau, tau;
  double cdiff[npha][nsol], lap[npha], dff[npha], phidot, Tau[npha][npha];
  double zeta[npha][npha], TAU, df[npha], hphi[npha];
  double sumdf;
  double retG[npha], retmu[npha*nsol], retdmu[npha*nsol*nsol];
  double gphi[npha], dgphidphi[npha], dhphidphi[npha], mu[nsol], dhphi[npha][npha];
  double dwh[npha*npha], dwhabc[npha*npha*npha], laps[npha]; 
  double sumlambdaphi, lambaphi[npha];
  double stphi[npha*3*3*3], dab[npha*npha], ee_ab[npha*npha];
  double Rot_mat[npha*npha*3*3], Inv_Rot_mat[npha*npha*3*3];
  //double facecenter, edgecenter, corner;
  
  struct fields stgridO[27];
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  //index = y + ny*(z + x*nz);
  
  if ( x != 0 && x != nx-1 && y != 0 && y != ny-1 && z != 0 && z != nz-1 ) {
    
    index = y + ny*(z + x*nz);
    
    Ti = pfmdat->T0;
    
    stgridO[0]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[1]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[2]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[3]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[4]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[5]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[6]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[7]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[8]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[9]  = gridinfoO[ (y-1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[10] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[11] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[12] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[13] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[14] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[15] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[16] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[17] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[18] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[19] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[20] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[21] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[22] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[23] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[24] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[25] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[26] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z+1) ) ];
/*
    if (pfmdat->ELASTICITY) {
      stgO[0] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];

      st_it_gO[0] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
      st_it_gO[1] = it_gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
      st_it_gO[2] = it_gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
      st_it_gO[3] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
      st_it_gO[4] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
      st_it_gO[5] = it_gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
      st_it_gO[6] = it_gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];

      eigen_strain[0].xx = 0.0;
      eigen_strain[0].yy = 0.0;
      eigen_strain[0].zz = 0.0;
      eigen_strain[0].yz = 0.0;
      eigen_strain[0].xz = 0.0;
      eigen_strain[0].xy = 0.0;
      for (ip = 0; ip < npha; ip++) {
        eigen_strain[X].xx += eigen_strain_phase[ip].xx*stgO[0].phi[ip];
        eigen_strain[X].yy += eigen_strain_phase[ip].yy*stgO[0].phi[ip];
        eigen_strain[X].zz += eigen_strain_phase[ip].zz*stgO[0].phi[ip];
        eigen_strain[X].yz += eigen_strain_phase[ip].yz*stgO[0].phi[ip];
        eigen_strain[X].xz += eigen_strain_phase[ip].xz*stgO[0].phi[ip];
        eigen_strain[X].xy += eigen_strain_phase[ip].xy*stgO[0].phi[ip];
      }
      eigen_strain[Y] = eigen_strain[X];
      eigen_strain[Z] = eigen_strain[X];

      strain[X].xx = 0.5*(st_it_gO[6].disp[X][2] - st_it_gO[5].disp[X][2]) - eigen_strain[X].xx;
      strain[X].yy = 0.5*(st_it_gO[2].disp[Y][2] - st_it_gO[1].disp[Y][2]) - eigen_strain[X].yy;
      strain[X].zz = 0.5*(st_it_gO[4].disp[Z][2] - st_it_gO[3].disp[Z][2]) - eigen_strain[X].zz;

      strain[Y].xx = strain[X].xx;
      strain[Y].yy = strain[X].yy;
      strain[Y].zz = strain[X].zz;
      strain[Z].yy = strain[X].yy;
      strain[Z].xx = strain[X].xx;
      strain[Z].zz = strain[X].zz;

      strain[X].xy = 0.25*((st_it_gO[2].disp[X][2] - st_it_gO[1].disp[X][2]) + (st_it_gO[6].disp[Y][2] - st_it_gO[5].disp[Y][2]));
      strain[X].xz = 0.25*((st_it_gO[4].disp[X][2] - st_it_gO[3].disp[X][2]) + (st_it_gO[6].disp[Z][2] - st_it_gO[5].disp[Z][2]));
      strain[X].yz = 0.25*((st_it_gO[2].disp[Z][2] - st_it_gO[1].disp[Z][2]) + (st_it_gO[4].disp[Y][2] - st_it_gO[3].disp[Y][2]));

      strain[Y].xy = strain[X].xy;
      strain[Y].xz = strain[X].xz;
      strain[Y].yz = strain[X].yz;
      strain[Z].xy = strain[X].xy;
      strain[Z].xz = strain[X].xz;
      strain[Z].yz = strain[X].yz;

      stiffness_c.C11 = 0.0;
      stiffness_c.C12 = 0.0;
      stiffness_c.C44 = 0.0;
      for (ip = 0; ip < npha; ip++) {
        stiffness_c[X].C11 += (stiffness_phase_n[ip].C11)*stgO[0].phi[ip];
        stiffness_c[X].C12 += (stiffness_phase_n[ip].C12)*stgO[0].phi[ip];
        stiffness_c[X].C44 += (stiffness_phase_n[ip].C44)*stgO[0].phi[ip];
      }
      stiffness_c[Y] = stiffness_c[X];
      stiffness_c[Z] = stiffness_c[X];


      sigma.xx  = stiffness_c[X].C11*(strain[X].xx) + stiffness_c[X].C12*(strain[X].yy + strain[X].zz);

      sigma.yy  = stiffness_c[X].C12*(strain[X].xx  + strain[X].zz) + stiffness_c[X].C11*(strain[X].yy);

      sigma.zz  = stiffness_c[X].C12*(strain[X].xx  + strain[X].yy) + stiffness_c[X].C11*(strain[X].zz);

      sigma.xy  = 2.0*stiffness_c[X].C44*(strain[X].xy);

      sigma.xz  = 2.0*stiffness_c[X].C44*(strain[X].xz);

      sigma.yz  = 2.0*stiffness_c[X].C44*(strain[X].yz);
    }
*/
    for (ip = 0; ip < npha; ip++) {
      ii = 0; 
        for (i1 = 0; i1 < 3; i1++) { 
          for (j1 = 0; j1 < 3; j1++) { 
            for (k1 = 0; k1 < 3; k1++) {
              //stphi[ (ip*3 + i1)*3 + j1 ] = stphi_1[ip][i1][j1];
              stphi[ ( (ip*3 + i1)*3 + j1 )*3 + k1 ] = stgridO[ii].phi[ip];
              //if (stphi[ ( (ip*3 + i1)*3 + j1 )*3 + k1 ] > 0.5) {
              //      printf("%d, %d, %d, %d, %d, %d, %d, %d, %lf, %lf\n", y, z, x, ip, i1, j1, k1, ( (ip*3 + i1)*3 + j1 )*3 + k1, stphi[ ( (ip*3 + i1)*3 + j1 )*3 + k1 ], stgridO[ii].phi[ip]);
              //}
              ii++;
            }
          }
        }
    }
    
    
    for ( ip = 0; ip < npha; ip++ ) { 
      // fhi[ip] = stgridO[13].phi[ip];
      fhi[ip] = gridinfoO[index].phi[ip];
    }
    
    // Calculate laplacian 
    for ( ip = 0; ip < npha; ip++ ) { 

      lap[ip]  = 0.0;

      // facecenter = stgridO[4].phi[ip] + stgridO[22].phi[ip] + stgridO[10].phi[ip] + stgridO[16].phi[ip] + stgridO[12].phi[ip] + stgridO[14].phi[ip];

      // edgecenter = stgridO[1].phi[ip] + stgridO[5].phi[ip] + stgridO[7].phi[ip] + stgridO[3].phi[ip] + stgridO[9].phi[ip] + stgridO[11].phi[ip] + stgridO[17].phi[ip] + stgridO[15].phi[ip] + stgridO[19].phi[ip] + stgridO[21].phi[ip] + stgridO[25].phi[ip] + stgridO[23].phi[ip];

      // corner = stgridO[0].phi[ip] + stgridO[2].phi[ip] + stgridO[8].phi[ip] + stgridO[6].phi[ip] + stgridO[18].phi[ip] + stgridO[20].phi[ip] + stgridO[26].phi[ip] + stgridO[24].phi[ip];

      // lap[ip] = 3.0*( facecenter + edgecenter/2.0 + corner/3.0 - 44.0*stgridO[13].phi[ip]/3.0 ) / (13.0 * (pfmvar->deltax) * (pfmvar->deltax));
      
      lap[ip]  = ( stgridO[14].phi[ip] - 2.0*stgridO[13].phi[ip] + stgridO[12].phi[ip] ) / (pfmvar->deltay*pfmvar->deltay);
      lap[ip] += ( stgridO[16].phi[ip] - 2.0*stgridO[13].phi[ip] + stgridO[10].phi[ip] ) / (pfmvar->deltaz*pfmvar->deltaz);
      lap[ip] += ( stgridO[22].phi[ip] - 2.0*stgridO[13].phi[ip] + stgridO[ 4].phi[ip] ) / (pfmvar->deltax*pfmvar->deltax);


      
      //printf("%le\n", lap[ip]);
    }
    
    
    
    // Get phase compositions
    for ( ip = 0; ip < npha; ip++ ) { 
      for ( is = 0; is < nsol; is++ ) { 
        cie[ip][is] = cscl[index].comie[ip][is];
      }
    }
    
    for ( ip = 0; ip < npha; ip++ ) { 
      F3_C[ip] = propf3->C[ip];
      for ( is1 = 0; is1 < nsol; is1++ ) { 
        F3_B[ip*nsol+is1] = propf3->B[ip][is1];
        for ( is2 = 0; is2 < nsol; is2++ ) { 
          F3_A[(ip*nsol+is1)*nsol+is2] = propf3->A[ip][is1][is2]; 
          F3_cmu[(ip*nsol+is1)*nsol+is2] = propf3->cmu[ip][is1][is2]; 
        }
      }
    }
    
    // Calculate free energy 
    for ( ip = 0; ip < npha; ip ++ ) { 
      ge[ip] = function_F_03_free_energy_CL(cie[ip], Ti, ip, F3_A, F3_B, F3_C) / pfmdat->Vm;
    }
    
    // Get mu 
    for ( is = 0; is < nsol; is++ ) { 
      //mu[is] = stgridO[13].mu[is] / pfmdat->Vm;
      mu[is] = gridinfoO[index].mu[is] / pfmdat->Vm;
    }
    
    for ( ip = 0; ip < npha; ip++ ) { 
      function_F_03_dc_dmu_CL(mu, cie[ip], Ti, ip, F3_cmu, dcdmu[ip]);
    }
    
    for ( ip = 0; ip < npha; ip++ ) { 
      matinvnew_nsol(dcdmu[ip], ddgedcdc[ip]); 
    }
    
    
    // Calculate second derivative of free energy .i.e ddgedcdc and get diffusion coefficient to local variable 
    for ( ip = 0; ip < npha; ip++ ) { 
      for ( is1 = 0; is1 < nsol; is1++ ) { 
        for ( is2 = 0; is2 < nsol; is2++ ) { 
            
          ddgedcdc[ip][is1*nsol+is2] /= pfmdat->Vm;
          
          InvD[ip][is1*nsol+is2] = pfmdat->DInv[ip][is1*nsol+is2];
        }
      }
    }
    
    // Calculate TAU
    for ( ip1 = 0; ip1 < npha-1; ip1++ ) { 
      for ( is = 0; is < nsol; is++ ) { 
        deltac[is] = pfmdat->c_eq[(npha-1)*npha+(npha-1)][is] - pfmdat->c_eq[ip1*npha+ip1][is];
      }
      multiply_nsol(InvD[npha-1], deltac, Ddc);
      multiply_nsol(ddgedcdc[npha-1], Ddc, ddgDdc);
      
      tmp0 = 0.0; 
      for ( is = 0; is < nsol; is++ ) { 
        tmp0 += ddgDdc[is] * deltac[is];
      }

      zeta[ip1][npha-1] = tmp0;

      Tau[ip1][npha-1] = ( 6.0 * pfmvar->ee[ip1*npha+(npha-1)] * pfmdat->a2 * zeta[ip1][npha-1] ) / pfmvar->w[ip1*npha+(npha-1)];
      Tau[npha-1][ip1] = Tau[ip1][npha-1];

      if ( ip1 == 0 ) { 
        min_tau = Tau[ip1][npha-1];
      }
      if ( Tau[ip1][npha-1] < min_tau ) { 
        min_tau = Tau[ip1][npha-1];
      }
    }
    for ( ip1 = 0; ip1 < npha-1; ip1++ ) { 
      for ( ip2 = 0; ip2 < npha-1; ip2++ ) { 
        Tau[ip1][ip2] = min_tau;
      }
    }
    tau = min_tau;
    
    tmp0 = 0.0;
    tmp1 = 0.0;
    for ( ip1 = 0; ip1 < npha; ip1++ ) { 
      for ( ip2 = 0; ip2 < npha; ip2++ ) { 
        if ( ip1 < ip2 ) {
          tmp0 += Tau[ip1][ip2]*fhi[ip1]*fhi[ip2];
          tmp1 += fhi[ip1]*fhi[ip2];
        }
      }
    }
    if ( tmp1 ) {
      TAU = tmp0 / tmp1;
    } else {
      TAU = tau;
    }
    
    
    // Calculate df_TD/dphi
    for ( ip1 = 0; ip1 < npha; ip1++ ) { 
      dff[ip1] = 0.0;
      for ( ip2 = 0; ip2 < npha; ip2++ ) { 
          
          psi = 0.0; 
          psi += ge[ip2]; 
          for ( is = 0; is < nsol; is++ ) { 
            psi -= mu[is]*cie[ip2][is];
          }
          psi *= dhfhi(fhi, ip2, ip1); 
          
          dff[ip1] +=  psi; 
      }
    }
    
    // get w_ab w_abc ee dab ti local linear indexed variable
    for ( ip1 = 0; ip1 < npha; ip1++ ) { 
      for ( ip2 = 0; ip2 < npha; ip2++ ) { 
        dwh[ip1*npha+ip2] = pfmvar->w[ip1*npha+ip2];
        dab[ip1*npha+ip2] = pfmdat->epsc[ip1*npha+ip2];
        ee_ab[ip1*npha+ip2] = pfmvar->ee[ip1*npha+ip2];
        for ( ip3 = 0; ip3 < npha; ip3++ ) { 
          dwhabc[(ip1*npha+ip2)*npha+ip3] = pfmvar->w_abc[(ip1*npha+ip2)*npha+ip3];
        }
      }
    }
    
    // Calculate derivative of well potential 
    for ( ip = 0; ip < npha; ip++ ) { 
      dgphidphi[ip] = F_W_02_dwdphi(fhi, lap, dwh, dwhabc, ip);
    }
    
    // Get the laplacian for all betas
    for ( ip1 = 0; ip1 < npha; ip1++ ) { 
      laps[ip1] = 0.0;
      for ( ip2 = 0; ip2 < npha; ip2++ ) { 
        if ( ip1 != ip2 ) { 
          laps[ip1] += pfmvar->ee[ip1*npha+ip2] * lap[ip2];
        }
      }
      // printf("    =%le\n", laps[ip]);
    }

    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      for ( ip2 = 0; ip2 < npha; ip2++ ) {
        for ( i1 = 0; i1 < 3; i1++ ) {
          for ( i2 = 0; i2 < 3; i2++ ) {
                Rot_mat[((ip1*npha+ip2)*npha+i1)*3+i2] = pfmdat->Rotation_matrix[ip1][ip2][i1][i2];
            Inv_Rot_mat[((ip1*npha+ip2)*npha+i1)*3+i2] = pfmdat->Inv_Rotation_matrix[ip1][ip2][i1][i2];
          }
        }
      }
    }

    for ( ip = 0; ip < npha; ip++ ) {
      //printf("%le, %le\n", laps[ip], calcAnisotropy_01(stphi, dab, ee_ab, ip, pfmvar->deltax, pfmvar->deltay, dz));
      //printf("%d, %d, %d, %d, %le, %le, %le\n", y, z, x, ip, pfmvar->deltax, pfmvar->deltay, pfmvar->deltaz);
      laps[ip] = calcAnisotropy_01(stphi, dab, ee_ab, ip, pfmvar->deltax, pfmvar->deltay, pfmvar->deltaz, Rot_mat, Inv_Rot_mat, y, z, x);
    }
    //printf("%d, %d, %d, %d, %le, %le, %le, %le, %le\n", y, z, x, ip, pfmvar->deltax, pfmvar->deltay, pfmvar->deltaz, laps[0], laps[1]);


    
    // Calculate -dF/dphi for all betas and check for active phases 
    sumlambdaphi = 0.0;
    activephases = 0;
    for ( ip1 = 0; ip1 < npha; ip1++ ) { 
      if ( fabs(lap[ip1]) > 0.0 ) {
        lambaphi[ip1] = -laps[ip1] - dgphidphi[ip1] - dff[ip1];

        /*
        if (pfmdat->ELASTICITY) {
          //lambda_phi[ip1] -= df_elast(grad, sigma, gridinfo_w[center].phi, a);

          delast = -(sigma.xx*eigen_strain_phase[ip1].xx + sigma.yy*eigen_strain_phase[ip1].yy + sigma.zz*eigen_strain_phase[ip1].zz + 2.0*sigma.xy*eigen_strain_phase[ip1].xy + 2.0*sigma.xz*eigen_strain_phase[ip1].xz + 2.0*sigma.yz*eigen_strain_phase[ip1].yz);

          sigma_phase.xx = stiffness_phase[ip1].C11*strain[X].xx  + stiffness_phase[ip1].C12*(strain[X].yy + strain[X].zz);
          sigma_phase.yy = stiffness_phase[ip1].C12*(strain[X].xx + strain[X].zz) + stiffness_phase[ip1].C11*strain[X].yy;
          sigma_phase.zz = stiffness_phase[ip1].C12*(strain[X].xx + strain[X].yy) + stiffness_phase[ip1].C11*strain[X].zz;

          sigma_phase.xy = 2.0*stiffness_phase[ip1].C44*strain[X].xy;
          sigma_phase.xz = 2.0*stiffness_phase[ip1].C44*strain[X].xz;
          sigma_phase.yz = 2.0*stiffness_phase[ip1].C44*strain[X].yz;

          delast += 0.5*(sigma_phase.xx*strain[X].xx + sigma_phase.yy*strain[X].yy + sigma_phase.zz*strain[X].zz + 2.0*sigma_phase.xy*strain[X].xy + 2.0*sigma_phase.xz*strain[X].xz + 2.0*sigma_phase.yz*strain[X].yz);


          lambda_phi[ip1] -= delast;

        }
        */
        sumlambdaphi += lambaphi[ip1];
        activephases++;
      }
    }
    
    // Divide the lambda with active phases
    if ( activephases ) { 
      sumlambdaphi /= activephases; 
    }
    
    // Update the phi by applying lagrangian multiplier
    for ( ip1 = 0; ip1 < npha; ip1++ ) { 
      if ( fabs(lap[ip1]) > 0.0 ) { 
        phidot = (2.0/(TAU)) * ( lambaphi[ip1] - sumlambdaphi );
        gridinfo[ index ].phi[ip1] = gridinfoO[index].phi[ip1] + (pfmvar->deltat)*phidot;
        
        //printf("t=%d, i=%d, j=%d, lap=%le,%le, mphi=%le, vard=%le, phi[0]=%le, ip=%d, phiold=%le\n", tstep[0], i, j, lap[0],lap[1],  1.0/TAU, 2.0*( lambaphi[ip1] - sumlambdaphi ), gridinfo[ index ].phi[ip1], ip1, stgridO[4].phi[ip1]);
      }
    }

  }

}

__kernel void SolverPhi_F4_smooth(__global struct fields *gridinfoO, __global struct fields *gridinfo, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf4 *propf4, __constant struct propmatf4spline *propf4spline) {


  int x, y, z, ii, i1, j1, k1, i2;
  int nx, ny, nz;
  int index;

  int is, js, ks, il, il1, jl, jl1, ipx;
  int is1, is2;
  int ig, jg, ip, ip1, ip2, ip3;
  int ilmax;
  int interface,bulkphase;
  int activephases;

  double Ti;

  double tmp0, tmp1;
  double facecenter, edgecenter, corner;

  double ge[npha], fhi[npha], psi;
  double dgedc[npha][nsol];
  double ddgedcdc[npha][nsol*nsol], InvD[npha][nsol*nsol], deltac[nsol], Ddc[nsol], ddgDdc[nsol];
  double cie[npha][nsol], min_tau, tau;
  double lap[npha], dff[npha], phidot, Tau[npha][npha];
  double zeta[npha][npha], TAU, df[npha], hphi[npha];
  double sumdf;
  double retG[npha], retmu[npha*nsol], retdmu[npha*nsol*nsol];
  double gphi[npha], dgphidphi[npha], dhphidphi[npha], mu[nsol], dhphi[npha][npha];
  double dwh[npha*npha], dwhabc[npha*npha*npha], laps[npha];
  double sumlambdaphi, lambaphi[npha];
  double dab[npha*npha], ee_ab[npha*npha];

  struct fields stgridO[27];

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  //index = y + ny*(z + x*nz);

  Ti = pfmdat->T0;

  if ( x != 0 && x != nx-1 && y != 0 && y != ny-1 && z != 0 && z != nz-1 ) {

    index = y + ny*(z + x*nz);

    stgridO[0]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[1]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[2]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[3]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[4]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[5]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[6]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[7]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[8]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[9]  = gridinfoO[ (y-1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[10] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[11] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[12] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[13] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[14] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[15] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[16] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[17] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[18] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[19] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[20] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[21] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[22] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[23] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[24] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[25] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[26] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z+1) ) ];

    for ( ip = 0; ip < npha; ip++ ) {
      // fhi[ip] = stgridO[13].phi[ip];
      fhi[ip] = gridinfoO[index].phi[ip];
    }

    // Calculate laplacian
    for ( ip = 0; ip < npha; ip++ ) {

      lap[ip]  = 0.0;

      facecenter = stgridO[4].phi[ip] + stgridO[22].phi[ip] + stgridO[10].phi[ip] + stgridO[16].phi[ip] + stgridO[12].phi[ip] + stgridO[14].phi[ip];

      edgecenter = stgridO[1].phi[ip] + stgridO[5].phi[ip] + stgridO[7].phi[ip] + stgridO[3].phi[ip] + stgridO[9].phi[ip] + stgridO[11].phi[ip] + stgridO[17].phi[ip] + stgridO[15].phi[ip] + stgridO[19].phi[ip] + stgridO[21].phi[ip] + stgridO[25].phi[ip] + stgridO[23].phi[ip];

      corner = stgridO[0].phi[ip] + stgridO[2].phi[ip] + stgridO[8].phi[ip] + stgridO[6].phi[ip] + stgridO[18].phi[ip] + stgridO[20].phi[ip] + stgridO[26].phi[ip] + stgridO[24].phi[ip];

      lap[ip] = 3.0*( facecenter + edgecenter/2.0 + corner/3.0 - 44.0*stgridO[13].phi[ip]/3.0 ) / (13.0 * (pfmvar->deltax) * (pfmvar->deltax));

//       lap[ip]  = ( stgridO[14].phi[ip] - 2.0*stgridO[13].phi[ip] + stgridO[12].phi[ip] ) / (pfmvar->deltay*pfmvar->deltay);
//       lap[ip] += ( stgridO[16].phi[ip] - 2.0*stgridO[13].phi[ip] + stgridO[10].phi[ip] ) / (pfmvar->deltaz*pfmvar->deltaz);
//       lap[ip] += ( stgridO[22].phi[ip] - 2.0*stgridO[13].phi[ip] + stgridO[ 4].phi[ip] ) / (pfmvar->deltax*pfmvar->deltax);



      //printf("%le\n", lap[ip]);
    }


    for ( ip = 0; ip < npha; ip++ ) {
      for ( is = 0; is < nsol; is++ ) {
        cie[ip][is] = cscl[index].comie[ip][is];
      }
    }

    // if ( pfmdat->ISOTHERMAL ) {
      for ( ip = 0; ip < npha; ip ++ ) {
        tmp0 = 0.0;
        tmp1 = 0.0;
        for ( is1 = 0; is1 < nsol; is1++ ) {
          for ( is2 = 0; is2 < nsol; is2++ ) {
            if ( is1 <= is2 ) {
              tmp0 += propf4->A[ip][is1][is2]*cie[ip][is1]*cie[ip][is2];
            }
          }
          tmp1 += propf4->B[ip][is1]*cie[ip][is1];
        }
        ge[ip] = ( tmp0 + tmp1 + propf4->C[ip] ) / pfmdat->Vm;
      }
    // }
    // else {
    //   for ( ip = 0; ip < npha; ip ++ ) {
    //     tmp0 = 0.0;
    //     tmp1 = 0.0;
    //     for ( is1 = 0; is1 < nsol; is1++ ) {
    //       for ( is2 = 0; is2 < nsol; is2++ ) {
    //         if ( is1 <= is2 ) {
    //           tmp0 += propf4spline[i].A[ip][is1][is2]*cie[ip][is1]*cie[ip][is2];
    //         }
    //       }
    //       tmp1 += propf4spline[i].B[ip][is1]*cie[ip][is1];
    //     }
    //     ge[ip] = ( tmp0 + tmp1 + propf4spline[i].C[ip] ) / pfmdat->Vm;
    //   }
    // }


    for ( is = 0; is < nsol; is++ ) {
      //mu[is] = stgridO[13].mu[is] / pfmdat->Vm;
      mu[is] = gridinfoO[index].mu[is] / pfmdat->Vm;
    }

    // if ( pfmdat->ISOTHERMAL ) {
      for ( ip = 0; ip < npha; ip++ ) {
        for ( is1 = 0; is1 < nsol; is1++ ) {
          for ( is2 = 0; is2 < nsol; is2++ ) {
            if ( is1 == is2 ) {
              ddgedcdc[ip][is1*nsol+is2] = 2.0*propf4->A[ip][is1][is2] / pfmdat->Vm;
            }
            else {
              ddgedcdc[ip][is1*nsol+is2] = propf4->A[ip][is1][is2] / pfmdat->Vm;
            }
            InvD[ip][is1*nsol+is2] = pfmdat->DInv[ip][is1*nsol+is2];
          }
        }
      }
    // }
    // else {
    //   for ( ip = 0; ip < npha; ip++ ) {
    //     for ( is1 = 0; is1 < nsol; is1++ ) {
    //       for ( is2 = 0; is2 < nsol; is2++ ) {
    //         if ( is1 == is2 ) {
    //           ddgedcdc[ip][is1*nsol+is2] = 2.0*propf4spline[j].A[ip][is1][is2] / pfmdat->Vm;
    //         }
    //         else {
    //           ddgedcdc[ip][is1*nsol+is2] = propf4spline[j].A[ip][is1][is2] / pfmdat->Vm;
    //         }
    //         InvD[ip][is1*nsol+is2] = pfmdat->DInv[ip][is1*nsol+is2];
    //       }
    //     }
    //   }
    // }

    // Calculate TAU
    for ( ip1 = 0; ip1 < npha-1; ip1++ ) {
      for ( is = 0; is < nsol; is++ ) {
        deltac[is] = pfmdat->c_eq[(npha-1)*npha+(npha-1)][is] - pfmdat->c_eq[ip1*npha+ip1][is];
      }
      multiply_nsol(InvD[npha-1], deltac, Ddc);
      multiply_nsol(ddgedcdc[npha-1], Ddc, ddgDdc);

      tmp0 = 0.0;
      for ( is = 0; is < nsol; is++ ) {
        tmp0 += ddgDdc[is] * deltac[is];
      }

      zeta[ip1][npha-1] = tmp0;

      Tau[ip1][npha-1] = ( 6.0 * pfmvar->ee[ip1*npha+(npha-1)] * pfmdat->a2 * zeta[ip1][npha-1] ) / pfmvar->w[ip1*npha+(npha-1)];
      Tau[npha-1][ip1] = Tau[ip1][npha-1];

      if ( ip1 == 0 ) {
        min_tau = Tau[ip1][npha-1];
      }
      if ( Tau[ip1][npha-1] < min_tau ) {
        min_tau = Tau[ip1][npha-1];
      }
    }
    for ( ip1 = 0; ip1 < npha-1; ip1++ ) {
      for ( ip2 = 0; ip2 < npha-1; ip2++ ) {
        Tau[ip1][ip2] = min_tau;
      }
    }
    tau = min_tau;

    tmp0 = 0.0;
    tmp1 = 0.0;
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      for ( ip2 = 0; ip2 < npha; ip2++ ) {
        if ( ip1 < ip2 ) {
          tmp0 += Tau[ip1][ip2]*fhi[ip1]*fhi[ip2];
          tmp1 += fhi[ip1]*fhi[ip2];
        }
      }
    }
    if ( tmp1 ) {
      TAU = tmp0 / tmp1;
    } else {
      TAU = tau;
    }

    // get w_ab w_abc ee dab ti local linear indexed variable
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      for ( ip2 = 0; ip2 < npha; ip2++ ) {
        dwh[ip1*npha+ip2] = pfmvar->w[ip1*npha+ip2];
        dab[ip1*npha+ip2] = pfmdat->epsc[ip1*npha+ip2];
        ee_ab[ip1*npha+ip2] = pfmvar->ee[ip1*npha+ip2];
        for ( ip3 = 0; ip3 < npha; ip3++ ) {
          dwhabc[(ip1*npha+ip2)*npha+ip3] = pfmvar->w_abc[(ip1*npha+ip2)*npha+ip3];
        }
      }
    }

    // Calculate derivative of well potential
    for ( ip = 0; ip < npha; ip++ ) {
      dgphidphi[ip] = F_W_02_dwdphi(fhi, lap, dwh, dwhabc, ip);
    }

    // Get the laplacian for all betas
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      laps[ip1] = 0.0;
      for ( ip2 = 0; ip2 < npha; ip2++ ) {
        if ( ip1 != ip2 ) {
          laps[ip1] += pfmvar->ee[ip1*npha+ip2] * lap[ip2];
        }
      }
      // printf("    =%le\n", laps[ip]);
    }

    // Calculate -dF/dphi for all betas and check for active phases
    sumlambdaphi = 0.0;
    activephases = 0;
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      if ( fabs(lap[ip1]) > 0.0 ) {
         lambaphi[ip1] = -laps[ip1] - dgphidphi[ip1];
        sumlambdaphi += lambaphi[ip1];
        activephases++;
      }
    }

    // Divide the lambda with active phases
    if ( activephases ) {
      sumlambdaphi /= activephases;
    }

    // Update the phi by applying lagrangian multiplier
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      if ( fabs(lap[ip1]) > 0.0 ) {
        phidot = (2.0/(TAU)) * ( lambaphi[ip1] - sumlambdaphi );
        gridinfo[ index ].phi[ip1] = gridinfoO[index].phi[ip1] + (pfmvar->deltat)*phidot;

        //printf("t=%d, i=%d, j=%d, lap=%le,%le, mphi=%le, vard=%le, phi[0]=%le, ip=%d, phiold=%le\n", tstep[0], i, j, lap[0],lap[1],  1.0/TAU, 2.0*( lambaphi[ip1] - sumlambdaphi ), gridinfo[ index ].phi[ip1], ip1, stgridO[4].phi[ip1]);
      }
    }

  }
}

__kernel void SolverPhi_F4(__global struct fields *gridinfoO, __global struct fields *gridinfo, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf4 *propf4, __constant struct propmatf4spline *propf4spline, __global struct iter_variables *it_gridinfoO, __global struct symmetric_tensor *eigen_strain_phase, __global struct Stiffness_cubic *stiffness_phase, __global struct Stiffness_cubic *stiffness_phase_n) {

  int x, y, z, ii, i1, j1, k1, i2;
  int nx, ny, nz;
  int index;
  int X, Y, Z;

  int is, js, ks, il, il1, jl, jl1, ipx;
  int is1, is2;
  int ig, jg, ip, ip1, ip2, ip3;
  int ilmax;
  int interface,bulkphase;
  int activephases;

  double Ti;

  double tmp0, tmp1;
  double facecenter, edgecenter, corner;

  double ge[npha], fhi[npha], psi;
  double dgedc[npha][nsol];
  double ddgedcdc[npha][nsol*nsol], InvD[npha][nsol*nsol], deltac[nsol], Ddc[nsol], ddgDdc[nsol];
  double cie[npha][nsol], min_tau, tau;
  double lap[npha], dff[npha], phidot, Tau[npha][npha];
  double zeta[npha][npha], TAU, df[npha], hphi[npha];
  double sumdf;
  double retG[npha], retmu[npha*nsol], retdmu[npha*nsol*nsol];
  double gphi[npha], dgphidphi[npha], dhphidphi[npha], mu[nsol], dhphi[npha][npha];
  double dwh[npha*npha], dwhabc[npha*npha*npha], laps[npha];
  double sumlambdaphi, lambaphi[npha];
  double stphi[npha*3*3*3], dab[npha*npha], ee_ab[npha*npha];
  double Rot_mat[npha*npha*3*3], Inv_Rot_mat[npha*npha*3*3];

  double delast;

  struct fields stgridO[27], stgO[1];
  struct symmetric_tensor eigen_strain[3];
  struct symmetric_tensor strain[3];
  struct Stiffness_cubic stiffness_c[3];
  struct symmetric_tensor sigma_phase;
  struct symmetric_tensor sigma;

  struct iter_variables st_it_gO[7];

  X = 0;
  Y = 1;
  Z = 2;

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  //index = y + ny*(z + x*nz);

  Ti = pfmdat->T0;

  if ( x != 0 && x != nx-1 && y != 0 && y != ny-1 && z != 0 && z != nz-1 ) {

    index = y + ny*(z + x*nz);

    stgridO[0]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[1]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[2]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[3]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[4]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[5]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[6]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[7]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[8]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[9]  = gridinfoO[ (y-1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[10] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[11] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[12] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[13] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[14] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[15] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[16] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[17] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[18] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[19] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[20] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[21] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[22] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[23] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[24] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[25] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[26] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z+1) ) ];


    if (pfmdat->ELASTICITY) {
      stgO[0] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];

      st_it_gO[0] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
      st_it_gO[1] = it_gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
      st_it_gO[2] = it_gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
      st_it_gO[3] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
      st_it_gO[4] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
      st_it_gO[5] = it_gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
      st_it_gO[6] = it_gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];

      eigen_strain[0].xx = 0.0;
      eigen_strain[0].yy = 0.0;
      eigen_strain[0].zz = 0.0;
      eigen_strain[0].yz = 0.0;
      eigen_strain[0].xz = 0.0;
      eigen_strain[0].xy = 0.0;
      for (ip = 0; ip < npha; ip++) {
        eigen_strain[X].xx += eigen_strain_phase[ip].xx*stgO[0].phi[ip];
        eigen_strain[X].yy += eigen_strain_phase[ip].yy*stgO[0].phi[ip];
        eigen_strain[X].zz += eigen_strain_phase[ip].zz*stgO[0].phi[ip];
        eigen_strain[X].yz += eigen_strain_phase[ip].yz*stgO[0].phi[ip];
        eigen_strain[X].xz += eigen_strain_phase[ip].xz*stgO[0].phi[ip];
        eigen_strain[X].xy += eigen_strain_phase[ip].xy*stgO[0].phi[ip];
      }
      eigen_strain[Y] = eigen_strain[X];
      eigen_strain[Z] = eigen_strain[X];

      strain[X].xx = 0.5*(st_it_gO[6].disp[X][2] - st_it_gO[5].disp[X][2]) - eigen_strain[X].xx;
      strain[X].yy = 0.5*(st_it_gO[2].disp[Y][2] - st_it_gO[1].disp[Y][2]) - eigen_strain[X].yy;
      strain[X].zz = 0.5*(st_it_gO[4].disp[Z][2] - st_it_gO[3].disp[Z][2]) - eigen_strain[X].zz;

      strain[Y].xx = strain[X].xx;
      strain[Y].yy = strain[X].yy;
      strain[Y].zz = strain[X].zz;
      strain[Z].yy = strain[X].yy;
      strain[Z].xx = strain[X].xx;
      strain[Z].zz = strain[X].zz;

      strain[X].xy = 0.25*((st_it_gO[2].disp[X][2] - st_it_gO[1].disp[X][2]) + (st_it_gO[6].disp[Y][2] - st_it_gO[5].disp[Y][2]));
      strain[X].xz = 0.25*((st_it_gO[4].disp[X][2] - st_it_gO[3].disp[X][2]) + (st_it_gO[6].disp[Z][2] - st_it_gO[5].disp[Z][2]));
      strain[X].yz = 0.25*((st_it_gO[2].disp[Z][2] - st_it_gO[1].disp[Z][2]) + (st_it_gO[4].disp[Y][2] - st_it_gO[3].disp[Y][2]));

      strain[Y].xy = strain[X].xy;
      strain[Y].xz = strain[X].xz;
      strain[Y].yz = strain[X].yz;
      strain[Z].xy = strain[X].xy;
      strain[Z].xz = strain[X].xz;
      strain[Z].yz = strain[X].yz;

      stiffness_c[0].C11 = 0.0;
      stiffness_c[0].C12 = 0.0;
      stiffness_c[0].C44 = 0.0;
      for (ip = 0; ip < npha; ip++) {
        stiffness_c[X].C11 += (stiffness_phase[ip].C11)*stgO[0].phi[ip];
        stiffness_c[X].C12 += (stiffness_phase[ip].C12)*stgO[0].phi[ip];
        stiffness_c[X].C44 += (stiffness_phase[ip].C44)*stgO[0].phi[ip];
      }
      stiffness_c[Y] = stiffness_c[X];
      stiffness_c[Z] = stiffness_c[X];


      sigma.xx  = stiffness_c[X].C11*(strain[X].xx) + stiffness_c[X].C12*(strain[X].yy + strain[X].zz);

      sigma.yy  = stiffness_c[X].C12*(strain[X].xx  + strain[X].zz) + stiffness_c[X].C11*(strain[X].yy);

      sigma.zz  = stiffness_c[X].C12*(strain[X].xx  + strain[X].yy) + stiffness_c[X].C11*(strain[X].zz);

      sigma.xy  = 2.0*stiffness_c[X].C44*(strain[X].xy);

      sigma.xz  = 2.0*stiffness_c[X].C44*(strain[X].xz);

      sigma.yz  = 2.0*stiffness_c[X].C44*(strain[X].yz);

      //printf("%d, %d, %d, %d, %d, %d, %le, %le, %le\n", tstep[0], x, y, z, index, ip1, stiffness_c[X].C11, stiffness_c[X].C12, stiffness_c[X].C44);

      //printf("%d, %d, %d, %d, %d, %d, %le, %le, %le\n", tstep[0], x, y, z, index, ip1, stiffness_phase[1].C11, stiffness_phase[1].C12, stiffness_phase[1].C44);

    }

    for (ip = 0; ip < npha; ip++) {
      ii = 0;
        for (i1 = 0; i1 < 3; i1++) {
          for (j1 = 0; j1 < 3; j1++) {
            for (k1 = 0; k1 < 3; k1++) {
              //stphi[ (ip*3 + i1)*3 + j1 ] = stphi_1[ip][i1][j1];
              stphi[ ( (ip*3 + i1)*3 + j1 )*3 + k1 ] = stgridO[ii].phi[ip];
              ii++;
            }
          }
        }
    }

    for ( ip = 0; ip < npha; ip++ ) {
      // fhi[ip] = stgridO[13].phi[ip];
      fhi[ip] = gridinfoO[index].phi[ip];
    }

    // Calculate laplacian
    for ( ip = 0; ip < npha; ip++ ) {

      lap[ip]  = 0.0;

      facecenter = stgridO[4].phi[ip] + stgridO[22].phi[ip] + stgridO[10].phi[ip] + stgridO[16].phi[ip] + stgridO[12].phi[ip] + stgridO[14].phi[ip];

      edgecenter = stgridO[1].phi[ip] + stgridO[5].phi[ip] + stgridO[7].phi[ip] + stgridO[3].phi[ip] + stgridO[9].phi[ip] + stgridO[11].phi[ip] + stgridO[17].phi[ip] + stgridO[15].phi[ip] + stgridO[19].phi[ip] + stgridO[21].phi[ip] + stgridO[25].phi[ip] + stgridO[23].phi[ip];

      corner = stgridO[0].phi[ip] + stgridO[2].phi[ip] + stgridO[8].phi[ip] + stgridO[6].phi[ip] + stgridO[18].phi[ip] + stgridO[20].phi[ip] + stgridO[26].phi[ip] + stgridO[24].phi[ip];

      lap[ip] = 3.0*( facecenter + edgecenter/2.0 + corner/3.0 - 44.0*stgridO[13].phi[ip]/3.0 ) / (13.0 * (pfmvar->deltax) * (pfmvar->deltax));

//       lap[ip]  = ( stgridO[14].phi[ip] - 2.0*stgridO[13].phi[ip] + stgridO[12].phi[ip] ) / (pfmvar->deltay*pfmvar->deltay);
//       lap[ip] += ( stgridO[16].phi[ip] - 2.0*stgridO[13].phi[ip] + stgridO[10].phi[ip] ) / (pfmvar->deltaz*pfmvar->deltaz);
//       lap[ip] += ( stgridO[22].phi[ip] - 2.0*stgridO[13].phi[ip] + stgridO[ 4].phi[ip] ) / (pfmvar->deltax*pfmvar->deltax);



      //printf("%le\n", lap[ip]);
    }


    for ( ip = 0; ip < npha; ip++ ) {
      for ( is = 0; is < nsol; is++ ) {
        cie[ip][is] = cscl[index].comie[ip][is];
      }
    }

    // if ( pfmdat->ISOTHERMAL ) {
      for ( ip = 0; ip < npha; ip ++ ) {
        tmp0 = 0.0;
        tmp1 = 0.0;
        for ( is1 = 0; is1 < nsol; is1++ ) {
          for ( is2 = 0; is2 < nsol; is2++ ) {
            if ( is1 <= is2 ) {
              tmp0 += propf4->A[ip][is1][is2]*cie[ip][is1]*cie[ip][is2];
            }
          }
          tmp1 += propf4->B[ip][is1]*cie[ip][is1];
        }
        ge[ip] = ( tmp0 + tmp1 + propf4->C[ip] ) / pfmdat->Vm;
      }
    // }
    // else {
    //   for ( ip = 0; ip < npha; ip ++ ) {
    //     tmp0 = 0.0;
    //     tmp1 = 0.0;
    //     for ( is1 = 0; is1 < nsol; is1++ ) {
    //       for ( is2 = 0; is2 < nsol; is2++ ) {
    //         if ( is1 <= is2 ) {
    //           tmp0 += propf4spline[i].A[ip][is1][is2]*cie[ip][is1]*cie[ip][is2];
    //         }
    //       }
    //       tmp1 += propf4spline[i].B[ip][is1]*cie[ip][is1];
    //     }
    //     ge[ip] = ( tmp0 + tmp1 + propf4spline[i].C[ip] ) / pfmdat->Vm;
    //   }
    // }


    for ( is = 0; is < nsol; is++ ) {
      //mu[is] = stgridO[13].mu[is] / pfmdat->Vm;
      mu[is] = gridinfoO[index].mu[is] / pfmdat->Vm;
    }

    // if ( pfmdat->ISOTHERMAL ) {
      for ( ip = 0; ip < npha; ip++ ) {
        for ( is1 = 0; is1 < nsol; is1++ ) {
          for ( is2 = 0; is2 < nsol; is2++ ) {
            if ( is1 == is2 ) {
              ddgedcdc[ip][is1*nsol+is2] = 2.0*propf4->A[ip][is1][is2] / pfmdat->Vm;
            }
            else {
              ddgedcdc[ip][is1*nsol+is2] = propf4->A[ip][is1][is2] / pfmdat->Vm;
            }
            InvD[ip][is1*nsol+is2] = pfmdat->DInv[ip][is1*nsol+is2];
          }
        }
      }
    // }
    // else {
    //   for ( ip = 0; ip < npha; ip++ ) {
    //     for ( is1 = 0; is1 < nsol; is1++ ) {
    //       for ( is2 = 0; is2 < nsol; is2++ ) {
    //         if ( is1 == is2 ) {
    //           ddgedcdc[ip][is1*nsol+is2] = 2.0*propf4spline[j].A[ip][is1][is2] / pfmdat->Vm;
    //         }
    //         else {
    //           ddgedcdc[ip][is1*nsol+is2] = propf4spline[j].A[ip][is1][is2] / pfmdat->Vm;
    //         }
    //         InvD[ip][is1*nsol+is2] = pfmdat->DInv[ip][is1*nsol+is2];
    //       }
    //     }
    //   }
    // }

    // Calculate TAU
    for ( ip1 = 0; ip1 < npha-1; ip1++ ) {
      for ( is = 0; is < nsol; is++ ) {
        deltac[is] = pfmdat->c_eq[(npha-1)*npha+(npha-1)][is] - pfmdat->c_eq[ip1*npha+ip1][is];
      }
      multiply_nsol(InvD[npha-1], deltac, Ddc);
      multiply_nsol(ddgedcdc[npha-1], Ddc, ddgDdc);

      tmp0 = 0.0;
      for ( is = 0; is < nsol; is++ ) {
        tmp0 += ddgDdc[is] * deltac[is];
      }

      zeta[ip1][npha-1] = tmp0;

      Tau[ip1][npha-1] = ( 6.0 * pfmvar->ee[ip1*npha+(npha-1)] * pfmdat->a2 * zeta[ip1][npha-1] ) / pfmvar->w[ip1*npha+(npha-1)];
      Tau[npha-1][ip1] = Tau[ip1][npha-1];

      if ( ip1 == 0 ) {
        min_tau = Tau[ip1][npha-1];
      }
      if ( Tau[ip1][npha-1] < min_tau ) {
        min_tau = Tau[ip1][npha-1];
      }
    }
    for ( ip1 = 0; ip1 < npha-1; ip1++ ) {
      for ( ip2 = 0; ip2 < npha-1; ip2++ ) {
        Tau[ip1][ip2] = min_tau;
      }
    }
    tau = min_tau;

    tmp0 = 0.0;
    tmp1 = 0.0;
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      for ( ip2 = 0; ip2 < npha; ip2++ ) {
        if ( ip1 < ip2 ) {
          tmp0 += Tau[ip1][ip2]*fhi[ip1]*fhi[ip2];
          tmp1 += fhi[ip1]*fhi[ip2];
        }
      }
    }
    if ( tmp1 ) {
      TAU = tmp0 / tmp1;
    } else {
      TAU = tau;
    }


    // Calculate df_TD/dphi
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      dff[ip1] = 0.0;
      for ( ip2 = 0; ip2 < npha; ip2++ ) {

          psi = 0.0;
          psi += ge[ip2];
          for ( is = 0; is < nsol; is++ ) {
            psi -= mu[is]*cie[ip2][is];
          }
          psi *= dhfhi(fhi, ip2, ip1);

          dff[ip1] +=  psi;
      }
    }
      //printf("%d, %d, %d, %d, %d, %le, %le\n", tstep[0], x, y, z, index, dff[0], dff[1]);

    // get w_ab w_abc ee dab ti local linear indexed variable
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      for ( ip2 = 0; ip2 < npha; ip2++ ) {
        dwh[ip1*npha+ip2] = pfmvar->w[ip1*npha+ip2];
        dab[ip1*npha+ip2] = pfmdat->epsc[ip1*npha+ip2];
        ee_ab[ip1*npha+ip2] = pfmvar->ee[ip1*npha+ip2];
        for ( ip3 = 0; ip3 < npha; ip3++ ) {
          dwhabc[(ip1*npha+ip2)*npha+ip3] = pfmvar->w_abc[(ip1*npha+ip2)*npha+ip3];
        }
      }
    }

    // Calculate derivative of well potential
    for ( ip = 0; ip < npha; ip++ ) {
      dgphidphi[ip] = F_W_02_dwdphi(fhi, lap, dwh, dwhabc, ip);
      //printf("%d, %d, %d, %d, %d, %d, %le\n", tstep[0], x, y, z, index, ip, dgphidphi[ip]);
    }

    // Get the laplacian for all betas
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      laps[ip1] = 0.0;
      for ( ip2 = 0; ip2 < npha; ip2++ ) {
        if ( ip1 != ip2 ) {
          laps[ip1] += pfmvar->ee[ip1*npha+ip2] * lap[ip2];
        }
      }
      // printf("    =%le\n", laps[ip]);
    }

    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      for ( ip2 = 0; ip2 < npha; ip2++ ) {
        for ( i1 = 0; i1 < 3; i1++ ) {
          for ( i2 = 0; i2 < 3; i2++ ) {
                Rot_mat[((ip1*npha+ip2)*3+i1)*3+i2] = pfmdat->Rotation_matrix[ip1][ip2][i1][i2];
            Inv_Rot_mat[((ip1*npha+ip2)*3+i1)*3+i2] = pfmdat->Inv_Rotation_matrix[ip1][ip2][i1][i2];
          }
        }
      }
    }

//     for ( ip = 0; ip < npha; ip++ ) {
//       //printf("%le, %le\n", laps[ip], calcAnisotropy_01(stphi, dab, ee_ab, ip, pfmvar->deltax, pfmvar->deltay, pfmvar->deltaz, Rot_mat, Inv_Rot_mat, y, z, x));
//       laps[ip] = calcAnisotropy_01(stphi, dab, ee_ab, ip, pfmvar->deltax, pfmvar->deltay, pfmvar->deltaz, Rot_mat, Inv_Rot_mat, y, z, x);
//     }
    //printf("%d, %d, %d, %d, %d, %le, %le, %le, %le, %le, %le, %le\n", tstep[0], y, z, x, index, pfmvar->deltax, pfmvar->deltay, pfmvar->deltaz, Rot_mat[1], Inv_Rot_mat[1], laps[0], laps[1]);

    // Calculate -dF/dphi for all betas and check for active phases
    sumlambdaphi = 0.0;
    activephases = 0;
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      if ( fabs(lap[ip1]) > 0.0 ) {
         lambaphi[ip1] = -laps[ip1] - dgphidphi[ip1] - dff[ip1];


        if (pfmdat->ELASTICITY) {
          //lambda_phi[ip1] -= df_elast(grad, sigma, gridinfo_w[center].phi, a);

          delast = -(sigma.xx*eigen_strain_phase[ip1].xx + sigma.yy*eigen_strain_phase[ip1].yy + sigma.zz*eigen_strain_phase[ip1].zz + 2.0*sigma.xy*eigen_strain_phase[ip1].xy + 2.0*sigma.xz*eigen_strain_phase[ip1].xz + 2.0*sigma.yz*eigen_strain_phase[ip1].yz);

          //printf("%d, %d, %d, %d, %d, %d, %le\n", tstep[0], x, y, z, index, ip1, delast);
          //printf("%d, %d, %d, %d, %d, %d, %le, %le, %le, %le, %le, %le\n", tstep[0], x, y, z, index, ip1, eigen_strain_phase[ip1].xx, eigen_strain_phase[ip1].yy, eigen_strain_phase[ip1].zz, eigen_strain_phase[ip1].xy, eigen_strain_phase[ip1].xz, eigen_strain_phase[ip1].yz);
          //printf("%d, %d, %d, %d, %d, %d, %le, %le, %le, %le, %le, %le\n", tstep[0], x, y, z, index, ip1, sigma.xx, sigma.yy, sigma.zz, sigma.xy, sigma.xz, sigma.yz);

          sigma_phase.xx = stiffness_phase[ip1].C11*strain[X].xx  + stiffness_phase[ip1].C12*(strain[X].yy + strain[X].zz);
          sigma_phase.yy = stiffness_phase[ip1].C12*(strain[X].xx + strain[X].zz) + stiffness_phase[ip1].C11*strain[X].yy;
          sigma_phase.zz = stiffness_phase[ip1].C12*(strain[X].xx + strain[X].yy) + stiffness_phase[ip1].C11*strain[X].zz;

          sigma_phase.xy = 2.0*stiffness_phase[ip1].C44*strain[X].xy;
          sigma_phase.xz = 2.0*stiffness_phase[ip1].C44*strain[X].xz;
          sigma_phase.yz = 2.0*stiffness_phase[ip1].C44*strain[X].yz;

          delast += 0.5*(sigma_phase.xx*strain[X].xx + sigma_phase.yy*strain[X].yy + sigma_phase.zz*strain[X].zz + 2.0*sigma_phase.xy*strain[X].xy + 2.0*sigma_phase.xz*strain[X].xz + 2.0*sigma_phase.yz*strain[X].yz);


          lambaphi[ip1] -= delast;

          //printf("%d, %d, %d, %d, %d, %d, %le\n", tstep[0], x, y, z, index, ip1, delast);

          //printf("%d, %d, %d, %d, %d, %d, %le, %le, %le, %le, %le, %le\n", tstep[0], x, y, z, index, ip1, sigma_phase.xx, sigma_phase.yy, sigma_phase.zz, sigma_phase.xy, sigma_phase.xz, sigma_phase.yz);

          //printf("%d, %d, %d, %d, %d, %d, %le, %le, %le, %le, %le, %le\n", tstep[0], x, y, z, index, ip1, strain[X].xx, strain[X].yy, strain[X].zz, strain[X].xy, strain[X].xz, strain[X].yz);

        }

        //printf("%d, %d, %d, %d, %d, %d, %le, %d, %d\n", tstep[0], x, y, z, index, ip1, delast, pfmdat->ELASTICITY, pfmdat->ISOTHERMAL);


        sumlambdaphi += lambaphi[ip1];
        activephases++;
        //printf("%d, %d, %d, %d, %d, %d, %le, %le\n", tstep[0], x, y, z, index, ip1, dff[ip1], -laps[ip1]);
        //printf("t=%d, y=%d, z=%d, x=%d, lap=%le, %le, mphi=%le, vard=%le, phi[0]=%le, ip=%d, phiold=%le\n", tstep[0], y, z, x, laps[0], laps[1], 1.0/TAU,  lambaphi[ip1], gridinfo[ index ].phi[ip1], ip1, stgridO[13].phi[ip1]);
      }
    }

    // Divide the lambda with active phases
    if ( activephases ) {
      sumlambdaphi /= activephases;
    }

    // Update the phi by applying lagrangian multiplier
    for ( ip1 = 0; ip1 < npha; ip1++ ) {
      if ( fabs(lap[ip1]) > 0.0 ) {
        phidot = (2.0/(TAU)) * ( lambaphi[ip1] - sumlambdaphi );
        gridinfo[ index ].phi[ip1] = gridinfoO[index].phi[ip1] + (pfmvar->deltat)*phidot;

        // printf("t=%d, y=%d, z=%d, x=%d, lap=%le, %le, mphi=%le, vard=%le, phi[0]=%le, ip=%d, phiold=%le\n", tstep[0], y, z, x, lap[0], lap[1],  1.0/TAU, 2.0*( lambaphi[ip1] - sumlambdaphi ), gridinfo[ index ].phi[ip1], ip1, stgridO[13].phi[ip1]);
      }
    }

  }
}

__kernel void SolverCatr_F2_smooth(__global struct fields *gridinfoO, __global struct fields *gridinfo, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf3 *propf3) {

  int x, y, z, jj, ii;
  int nx, ny, nz;
  int index;

  int ig, is, js, ip, i1, j1, is1, is2, ip1, ip2, ks;
  int interface, bulkphase;

  double A1[npha-1];
  double dphidt[npha-1][7], gradx_phi[npha][5], grady_phi[npha][5], gradz_phi[npha][5], phix[npha][7], phiy[npha][7], phiz[npha][7], modgradphi[npha][7], scalprodct[npha][7];
  double cjatx, cjaty, cjatz;
  double jatc[npha-1][nsol], jatr[nsol];

  //double hphid[5][npha], hphi[npha], dhphi[npha][npha], hphicii[5][npha];
  double tmp0, tmp1;
  double DELTAT, sum_dhphi, Ti;
  double sum[nsol], sum_dcbdT[nsol], dcdmu[nsol*nsol], dc_dmu[7][npha][nsol*nsol], Da[7][nsol][nsol];
  double Damidx[2][nsol][nsol], Damidy[2][nsol][nsol], Damidz[2][nsol][nsol], gradmux[2][nsol], gradmuy[2][nsol], gradmuz[2][nsol];
  double divflux[nsol], mu[nsol], deltamu[nsol], c_tdt[nsol], dcbdT_phase[npha][nsol];
  double deltac[nsol], inv_dcdmu[nsol*nsol], cg[nsol], cphas[nsol], cie[npha][nsol];
  double cas[npha*(nsol+1)], retdmu[npha*nsol*nsol], ddgedcdc[5][npha][nsol*nsol], retmu[npha*nsol], dgedc[npha][nsol];
  double yc[nsol+1], retdmuphase[nsol*nsol];
  double alpha[npha-1][nsol][7], alphidot[npha-1][nsol][6], jat[npha-1][nsol][6];

  struct fields stgridO[27];
  struct fields stgN[7];
  struct csle stcscl[7];
  struct fields stgO[7];

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  //index = y + ny*(z + x*nz);

  Ti = pfmdat->T0;

  if ( x != 0 && x != nx-1 && y != 0 && y != ny-1 && z != 0 && z != nz-1 ) {

    index = y + ny*(z + x*nz);

    stgridO[0]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[1]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[2]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[3]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[4]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[5]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[6]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[7]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[8]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[9]  = gridinfoO[ (y-1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[10] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[11] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[12] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[13] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[14] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[15] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[16] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[17] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[18] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[19] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[20] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[21] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[22] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[23] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[24] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[25] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[26] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z+1) ) ];

    stcscl[0] = cscl[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stcscl[1] = cscl[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stcscl[2] = cscl[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stcscl[3] = cscl[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stcscl[4] = cscl[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stcscl[5] = cscl[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stcscl[6] = cscl[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];

    stgO[0] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stgO[1] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stgO[2] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stgO[3] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stgO[4] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stgO[5] = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stgO[6] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];


      for ( is1 = 0; is1 < nsol; is1++ ) {
        sum[is1] = 0.0;
        sum_dcbdT[is1] = 0.0;
        for ( is2 = 0; is2 < nsol; is2++ ) {
          dcdmu[is1*nsol+is2] = 0.0;
        }
      }

      for ( ip = 0; ip < npha; ip++ ) {
        for ( is = 0; is < nsol; is++ ) {
          dcbdT_phase[ip][is] = 0.0;
        }
      }

      // Calculations are carried out at stencil points of 5 point stencil
      for ( ii = 0; ii < 7; ii++ ) {

        interface = 1;
        bulkphase = 0;
        for ( ip = 0; ip < npha; ip++ ) {
          if ( stgO[ii].phi[ip] >= pfmdat->interfaceUplimit ) {
            bulkphase = ip;
            interface = 0;
            break;
          }
        }

        if ( interface ) {

           for ( ip = 0; ip < npha; ip++ ) {

            // Get dc_dmu
            tmp0 = 0.0;
            for ( is = 0; is < nsol; is++ ) {
              yc[is] = stcscl[ii].comie[ip][is];
              tmp0 +=  yc[is];
            }
            yc[nsol] = 1.0 - tmp0;

            //printf("t = %d, i = %d, j = %d, ii = %d, ip = %d, tdbip = %d", tstep[0], i, j, ii, ip, pfmdat->thermophase[ip]);
            dMudc(Ti, yc, retdmuphase, pfmdat->thermophase[ip]);

            matinvnew_nsol(retdmuphase, dc_dmu[ii][ip]);

          }

          // Calculate mobility at stencil points
          for ( is = 0; is < nsol; is++ ) {
            for ( js =0; js < nsol; js++ ) {
              Da[ii][is][js] = 0.0;
              for ( ip = 0; ip < npha; ip++ ) {
                for ( ks = 0; ks < nsol; ks++ ) {
                  Da[ii][is][js] += pfmdat->D[ip][is*nsol+ks] * dc_dmu[ii][ip][ks*nsol+js] * stgO[ii].phi[ip];
                }
              }
            }
          }

        }
        else {

          // Get dc_dmu for bulk phase
          tmp0 = 0.0;
          for ( is = 0; is < nsol; is++ ) {
            yc[is] = stgO[ii].com[is];
            tmp0 +=  yc[is];
          }
          yc[nsol] = 1.0 - tmp0;

          dMudc(Ti, yc, retdmuphase, pfmdat->thermophase[bulkphase]);

          matinvnew_nsol(retdmuphase, dc_dmu[ii][bulkphase]);

          //printf("dBC %d, %d, %d, %d, %le\n", tstep[0], i, j, ii, dc_dmu[ii][ip][0]);
          for ( is = 0; is < nsol; is++ ) {
            for ( js =0; js < nsol; js++ ) {
              Da[ii][is][js] = 0.0;
              //for ( ip = 0; ip < npha; ip++ ) {
                for ( ks = 0; ks < nsol; ks++ ) {
                  Da[ii][is][js] += pfmdat->D[bulkphase][is*nsol+ks] * dc_dmu[ii][bulkphase][ks*nsol+js];
                }
              //}
            }
          }

        }

      }

      // Calculating mobility at half grid positions
      for ( is = 0; is < nsol; is++ ) {
        for ( js =0; js < nsol; js++ ) {

          Damidy[0][is][js] = ( Da[0][is][js] + Da[1][is][js] ) / 2.0;
          Damidy[1][is][js] = ( Da[2][is][js] + Da[0][is][js] ) / 2.0;

          Damidz[0][is][js] = ( Da[0][is][js] + Da[3][is][js] ) / 2.0;
          Damidz[1][is][js] = ( Da[4][is][js] + Da[0][is][js] ) / 2.0;

          Damidx[0][is][js] = ( Da[0][is][js] + Da[5][is][js] ) / 2.0;
          Damidx[1][is][js] = ( Da[6][is][js] + Da[0][is][js] ) / 2.0;
        }
      }

      // Calculating grad of mu
      for ( is = 0; is < nsol; is++ ) {
        gradmuy[0][is] = ( stgO[0].mu[is] - stgO[1].mu[is] ) / pfmvar->deltay;
        gradmuy[1][is] = ( stgO[2].mu[is] - stgO[0].mu[is] ) / pfmvar->deltay;

        gradmuz[0][is] = ( stgO[0].mu[is] - stgO[3].mu[is] ) / pfmvar->deltaz;
        gradmuz[1][is] = ( stgO[4].mu[is] - stgO[0].mu[is] ) / pfmvar->deltaz;

        gradmux[0][is] = ( stgO[0].mu[is] - stgO[5].mu[is] ) / pfmvar->deltax;
        gradmux[1][is] = ( stgO[6].mu[is] - stgO[0].mu[is] ) / pfmvar->deltax;
      }

      // Calculating divergence of flux
      for ( is = 0; is < nsol; is++ ) {
        divflux[is] = 0.0;
        for ( js =0; js < nsol; js++ ) {

          divflux[is] =  ( Damidy[1][is][js] * gradmuy[1][js] - Damidy[0][is][js] * gradmuy[0][js] ) / pfmvar->deltay;

          divflux[is] += ( Damidz[1][is][js] * gradmuz[1][js] - Damidz[0][is][js] * gradmuz[0][js] ) / pfmvar->deltaz;

          divflux[is] += ( Damidx[1][is][js] * gradmux[1][js] - Damidx[0][is][js] * gradmux[0][js] ) / pfmvar->deltax;

        }
      }


      // Claculation are done at i,j for updating composition and mu

      // Check for interface
      interface = 1;
      bulkphase = 0;
      for ( ip = 0; ip < npha; ip++ ) {
        if ( stgO[0].phi[ip] >= pfmdat->interfaceUplimit ) {
          bulkphase = ip;
          interface = 0;
          break;
        }
      }

      if ( !interface ) {

        for ( is = 0; is < nsol; is++ ) {
          gridinfo[index].com[is] = gridinfoO[index].com[is] + pfmvar->deltat * ( divflux[is] );
        }

        tmp0 = 0.0;
        for ( is = 0; is < nsol; is++ ) {
          yc[is] = gridinfo[index].com[is];
          tmp0 +=  yc[is];
        }
        yc[nsol] = 1.0 - tmp0;

        Mu(Ti, yc, mu, pfmdat->thermophase[bulkphase]);


        for ( is = 0; is < nsol; is++ ) {
          //deltamu[is] = mu[is] - gridinfoO[index].mu[is];
          gridinfo[index].mu[is] = mu[is];
        }

      }
      else {

        // Calculation are carried out in iterface region

        /*
        if ( !pfmdat->ISOTHERMAL ) {
          DELTAT = pfmvar->deltat * ( -pfmdat->TGRADIENT * pfmdat->velocity );
          for ( ip = 0; ip < npha; ip++ ) {

            for ( is = 0; is < nsol; is++ ) {
              cg[is] = pfmdat->cguess[ip*npha+ip][is];
              mu[is] = stgO[4].mu[is];
            }

            c_mu(mu, c_tdt, Ti+DELTAT, ip, cg, tstep[0], i, j, pfmdat->thermophase[ip]);

            for ( is = 0; is < nsol; is++ ) {
              dcbdT_phase[ip][is] = c_tdt[is] - stcscl[4].comie[ip][is];
            }
          }
        }
        */


        for ( is1 = 0; is1 < nsol; is1++ ) {

          // Calculate deltac
          deltac[is1] = pfmvar->deltat * ( divflux[is1] );

          // // Update composition using deltac
          // gridinfo[index].com[is1] = gridinfoO[index].com[is1] + deltac[is1];

          // // // Update deltac with c_alpha * dh_alpha/dt
          // // if ( pfmdat->ISOTHERMAL ) {
          //   deltac[is1] += -sum[is1];
          // // }
          // // else {
          // //  deltac[is1] += -sum[is] - sum_dcbdT[is1];
          // // }
          for ( is2 = 0; is2 < nsol; is2++ ) {
            for ( ip = 0; ip < npha; ip++ ) {
              dcdmu[is1*nsol+is2] += dc_dmu[0][ip][is1*nsol+is2] * hfhi(stgO[0].phi, ip);
            }
          }
        }

        matinvnew_nsol(dcdmu, inv_dcdmu);

        multiply_nsol(inv_dcdmu, deltac, deltamu);
        // for ( is1 = 0; is1 < nsol; is1++ ) {
        //   tmp0 = 0.0;
        //   for ( is2 = 0; is2 < nsol; is2++ ) {
        //     tmp0 += inv_dcdmu[is1*nsol+is2] * deltac[is2];
        //   }
        //   deltamu[is1] = tmp0;
        // }

        //vectorsum(deltamu, stgO[index].mu, gridinfo[index].mu, nsol);
        for ( is = 0; is < nsol; is++ ) {
          gridinfo[index].mu[is] = gridinfoO[index].mu[is] + deltamu[is];
        }
        //printf("IC %d, %d, %d, %le, %le\n", tstep[0], i, j, gridinfo[index].mu[0], gridinfo[index].com[0]);
        //printf("IC %d, %d, %d, %le, %le, %le, %le\n", tstep[0], i, j, gridinfoO[index].mu[0], gridinfo[index].mu[0], gridinfoO[index].com[0], gridinfo[index].com[0]);
      }

  }
}

__kernel void SolverCatr_F2(__global struct fields *gridinfoO, __global struct fields *gridinfo, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf3 *propf3) {

  int x, y, z, jj, ii;
  int nx, ny, nz;
  int index;

  int ig, is, js, ip, i1, j1, is1, is2, ip1, ip2, ks;
  int interface, bulkphase;

  double A1[npha-1];
  double dphidt[npha-1][7], gradx_phi[npha][5], grady_phi[npha][5], gradz_phi[npha][5], phix[npha][7], phiy[npha][7], phiz[npha][7], modgradphi[npha][7], scalprodct[npha][7];
  double cjatx, cjaty, cjatz;
  double jatc[npha-1][nsol], jatr[nsol];

  //double hphid[5][npha], hphi[npha], dhphi[npha][npha], hphicii[5][npha];
  double tmp0, tmp1;
  double DELTAT, sum_dhphi, Ti;
  double sum[nsol], sum_dcbdT[nsol], dcdmu[nsol*nsol], dc_dmu[7][npha][nsol*nsol], Da[7][nsol][nsol];
  double Damidx[2][nsol][nsol], Damidy[2][nsol][nsol], Damidz[2][nsol][nsol], gradmux[2][nsol], gradmuy[2][nsol], gradmuz[2][nsol];
  double divflux[nsol], mu[nsol], deltamu[nsol], c_tdt[nsol], dcbdT_phase[npha][nsol];
  double deltac[nsol], inv_dcdmu[nsol*nsol], cg[nsol], cphas[nsol], cie[npha][nsol];
  double cas[npha*(nsol+1)], retdmu[npha*nsol*nsol], ddgedcdc[5][npha][nsol*nsol], retmu[npha*nsol], dgedc[npha][nsol];
  double yc[nsol+1], retdmuphase[nsol*nsol];
  double alpha[npha-1][nsol][7], alphidot[npha-1][nsol][6], jat[npha-1][nsol][6];

  struct fields stgridO[27];
  struct fields stgN[7];
  struct csle stcscl[7];
  struct fields stgO[7];

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  //index = y + ny*(z + x*nz);

  Ti = pfmdat->T0;

  if ( x != 0 && x != nx-1 && y != 0 && y != ny-1 && z != 0 && z != nz-1 ) {

    index = y + ny*(z + x*nz);

    stgridO[0]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[1]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[2]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[3]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[4]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[5]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[6]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[7]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[8]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[9]  = gridinfoO[ (y-1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[10] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[11] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[12] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[13] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[14] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[15] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[16] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[17] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[18] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[19] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[20] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[21] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[22] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[23] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[24] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[25] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[26] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z+1) ) ];

    stcscl[0] = cscl[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stcscl[1] = cscl[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stcscl[2] = cscl[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stcscl[3] = cscl[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stcscl[4] = cscl[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stcscl[5] = cscl[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stcscl[6] = cscl[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];

    stgO[0] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stgO[1] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stgO[2] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stgO[3] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stgO[4] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stgO[5] = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stgO[6] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];

    for ( ip = 0; ip < npha-1; ip++ ) {
        stgN[0].phi[ip] = gridinfo[ (y  ) + ny*( (x  )*nz + (z  ) ) ].phi[ip];
        stgN[1].phi[ip] = gridinfo[ (y-1) + ny*( (x  )*nz + (z  ) ) ].phi[ip];
        stgN[2].phi[ip] = gridinfo[ (y+1) + ny*( (x  )*nz + (z  ) ) ].phi[ip];
        stgN[3].phi[ip] = gridinfo[ (y  ) + ny*( (x  )*nz + (z-1) ) ].phi[ip];
        stgN[4].phi[ip] = gridinfo[ (y  ) + ny*( (x  )*nz + (z+1) ) ].phi[ip];
        stgN[5].phi[ip] = gridinfo[ (y  ) + ny*( (x-1)*nz + (z  ) ) ].phi[ip];
        stgN[6].phi[ip] = gridinfo[ (y  ) + ny*( (x+1)*nz + (z  ) ) ].phi[ip];
    }

    for ( ip = 0; ip < npha-1; ip++ ) {
        dphidt[ip][0] = ( stgN[0].phi[ip] - stgO[0].phi[ip] ) / ( pfmvar->deltat );
        dphidt[ip][1] = ( stgN[1].phi[ip] - stgO[1].phi[ip] ) / ( pfmvar->deltat );
        dphidt[ip][2] = ( stgN[2].phi[ip] - stgO[2].phi[ip] ) / ( pfmvar->deltat );
        dphidt[ip][3] = ( stgN[3].phi[ip] - stgO[3].phi[ip] ) / ( pfmvar->deltat );
        dphidt[ip][4] = ( stgN[4].phi[ip] - stgO[4].phi[ip] ) / ( pfmvar->deltat );
        dphidt[ip][5] = ( stgN[5].phi[ip] - stgO[5].phi[ip] ) / ( pfmvar->deltat );
        dphidt[ip][6] = ( stgN[6].phi[ip] - stgO[6].phi[ip] ) / ( pfmvar->deltat );

      A1[ip] = sqrt ( pfmvar->ee[ip*npha+(npha-1)] ) / sqrt( pfmvar->w[ip*npha+(npha-1)] );

    }


    for ( ip = 0; ip < npha-1; ip++ ) {
        for ( is = 0; is < nsol; is++ ) {
          alpha[ip][is][0] = (A1[ip])*( stcscl[0].comie[npha-1][is] - stcscl[0].comie[ip][is] );
          alpha[ip][is][1] = (A1[ip])*( stcscl[1].comie[npha-1][is] - stcscl[1].comie[ip][is] );
          alpha[ip][is][2] = (A1[ip])*( stcscl[2].comie[npha-1][is] - stcscl[2].comie[ip][is] );
          alpha[ip][is][3] = (A1[ip])*( stcscl[3].comie[npha-1][is] - stcscl[3].comie[ip][is] );
          alpha[ip][is][4] = (A1[ip])*( stcscl[4].comie[npha-1][is] - stcscl[4].comie[ip][is] );
          alpha[ip][is][5] = (A1[ip])*( stcscl[5].comie[npha-1][is] - stcscl[5].comie[ip][is] );
          alpha[ip][is][6] = (A1[ip])*( stcscl[6].comie[npha-1][is] - stcscl[6].comie[ip][is] );
      }
    }

    for ( ip = 0; ip < npha; ip++ ) {
        grady_phi[ip][0] = ( stgridO[14].phi[ip] - stgridO[12].phi[ip] ) / ( 2.0*pfmvar->deltay );
        grady_phi[ip][1] = ( stgridO[11].phi[ip] - stgridO[9 ].phi[ip] ) / ( 2.0*pfmvar->deltay );
        grady_phi[ip][2] = ( stgridO[17].phi[ip] - stgridO[15].phi[ip] ) / ( 2.0*pfmvar->deltay );
        grady_phi[ip][3] = ( stgridO[5 ].phi[ip] - stgridO[3 ].phi[ip] ) / ( 2.0*pfmvar->deltay );
        grady_phi[ip][4] = ( stgridO[23].phi[ip] - stgridO[21].phi[ip] ) / ( 2.0*pfmvar->deltay );

        gradz_phi[ip][0] = ( stgridO[16].phi[ip] - stgridO[10].phi[ip] ) / ( 2.0*pfmvar->deltaz );
        gradz_phi[ip][1] = ( stgridO[15].phi[ip] - stgridO[9 ].phi[ip] ) / ( 2.0*pfmvar->deltaz );
        gradz_phi[ip][2] = ( stgridO[17].phi[ip] - stgridO[11].phi[ip] ) / ( 2.0*pfmvar->deltaz );
        gradz_phi[ip][3] = ( stgridO[7 ].phi[ip] - stgridO[1 ].phi[ip] ) / ( 2.0*pfmvar->deltaz );
        gradz_phi[ip][4] = ( stgridO[25].phi[ip] - stgridO[19].phi[ip] ) / ( 2.0*pfmvar->deltaz );

        gradx_phi[ip][0] = ( stgridO[22].phi[ip] - stgridO[4 ].phi[ip] ) / ( 2.0*pfmvar->deltax );
        gradx_phi[ip][1] = ( stgridO[21].phi[ip] - stgridO[3 ].phi[ip] ) / ( 2.0*pfmvar->deltax );
        gradx_phi[ip][2] = ( stgridO[23].phi[ip] - stgridO[5 ].phi[ip] ) / ( 2.0*pfmvar->deltax );
        gradx_phi[ip][3] = ( stgridO[19].phi[ip] - stgridO[1 ].phi[ip] ) / ( 2.0*pfmvar->deltax );
        gradx_phi[ip][4] = ( stgridO[25].phi[ip] - stgridO[7 ].phi[ip] ) / ( 2.0*pfmvar->deltax );
    }

    for ( ip = 0; ip < npha; ip++ ) {

        phiy[ip][0] = grady_phi[ip][0];
        phiy[ip][1] = ( stgridO[13].phi[ip] - stgridO[12].phi[ip] ) / ( pfmvar->deltay );
        phiy[ip][2] = ( stgridO[14].phi[ip] - stgridO[13].phi[ip] ) / ( pfmvar->deltay );
        phiy[ip][3] = ( grady_phi[ip][0] + grady_phi[ip][1] ) / ( 2.0 );
        phiy[ip][4] = ( grady_phi[ip][0] + grady_phi[ip][2] ) / ( 2.0 );
        phiy[ip][5] = ( grady_phi[ip][0] + grady_phi[ip][3] ) / ( 2.0 );
        phiy[ip][6] = ( grady_phi[ip][0] + grady_phi[ip][4] ) / ( 2.0 );

        phiz[ip][0] = gradz_phi[ip][0];
        phiz[ip][1] = ( gradz_phi[ip][0] + gradz_phi[ip][1] ) / ( 2.0 );
        phiz[ip][2] = ( gradz_phi[ip][0] + gradz_phi[ip][2] ) / ( 2.0 );
        phiz[ip][3] = ( stgridO[13].phi[ip] - stgridO[10].phi[ip] ) / ( pfmvar->deltaz );
        phiz[ip][4] = ( stgridO[16].phi[ip] - stgridO[13].phi[ip] ) / ( pfmvar->deltaz );
        phiz[ip][5] = ( gradz_phi[ip][0] + gradz_phi[ip][3] ) / ( 2.0 );
        phiz[ip][6] = ( gradz_phi[ip][0] + gradz_phi[ip][4] ) / ( 2.0 );

        phix[ip][0] = gradx_phi[ip][0];
        phix[ip][1] = ( gradx_phi[ip][0] + gradx_phi[ip][1] ) / ( 2.0 );
        phix[ip][2] = ( gradx_phi[ip][0] + gradx_phi[ip][2] ) / ( 2.0 );
        phix[ip][3] = ( gradx_phi[ip][0] + gradx_phi[ip][3] ) / ( 2.0 );
        phix[ip][4] = ( gradx_phi[ip][0] + gradx_phi[ip][4] ) / ( 2.0 );
        phix[ip][5] = ( stgridO[13].phi[ip] - stgridO[4 ].phi[ip] ) / ( pfmvar->deltax );
        phix[ip][6] = ( stgridO[22].phi[ip] - stgridO[13].phi[ip] ) / ( pfmvar->deltax );
    }

    for ( ip = 0; ip < npha-1; ip++ ) {
	      for ( is = 0; is < nsol; is++ ) {
          alphidot[ip][is][0] = ( ( alpha[ip][is][0] * dphidt[ip][0] ) + ( alpha[ip][is][1] * dphidt[ip][1] ) ) / ( 2.0 );
          alphidot[ip][is][1] = ( ( alpha[ip][is][0] * dphidt[ip][0] ) + ( alpha[ip][is][2] * dphidt[ip][2] ) ) / ( 2.0 );
          alphidot[ip][is][2] = ( ( alpha[ip][is][0] * dphidt[ip][0] ) + ( alpha[ip][is][3] * dphidt[ip][3] ) ) / ( 2.0 );
          alphidot[ip][is][3] = ( ( alpha[ip][is][0] * dphidt[ip][0] ) + ( alpha[ip][is][4] * dphidt[ip][4] ) ) / ( 2.0 );
          alphidot[ip][is][4] = ( ( alpha[ip][is][0] * dphidt[ip][0] ) + ( alpha[ip][is][5] * dphidt[ip][5] ) ) / ( 2.0 );
          alphidot[ip][is][5] = ( ( alpha[ip][is][0] * dphidt[ip][0] ) + ( alpha[ip][is][6] * dphidt[ip][6] ) ) / ( 2.0 );
	      }
    }

    for ( ip = 0; ip < npha; ip++ ) {
      for ( jj = 0; jj < 7; jj++ ) {
          modgradphi[ip][jj] = sqrt( phiy[ip][jj]*phiy[ip][jj] + phiz[ip][jj]*phiz[ip][jj] + phix[ip][jj]*phix[ip][jj]);
      }
    }

    for ( ip = 0; ip < npha-1; ip++ ) {
      for ( jj = 0; jj < 7; jj++ ) {
          scalprodct[ip][jj] = -1.0*( phiy[ip][jj]*phiy[npha-1][jj] + phiz[ip][jj]*phiz[npha-1][jj] + phix[ip][jj]*phix[npha-1][jj] );
          if ( modgradphi[npha-1][jj] > 0.0 ) {
            scalprodct[ip][jj] /= ( modgradphi[ip][jj] * modgradphi[npha-1][jj] );
        }
      }
    }



    for ( ip = 0; ip < npha-1; ip++ ) {
	    for ( is = 0; is < nsol; is++ ) {
          jat[ip][is][0] = ( ( alphidot[ip][is][0] * phiy[ip][1] ) / ( modgradphi[ip][1] ) )   * ( ( 1.0 - pfmdat->D[ip][is*nsol+is] /   pfmdat->D[npha-1][is*nsol+is] ) * fabs(scalprodct[ip][1]) );
          jat[ip][is][1] = ( ( alphidot[ip][is][1] * phiy[ip][2] ) / ( modgradphi[ip][2] ) )   * ( ( 1.0 - pfmdat->D[ip][is*nsol+is] /   pfmdat->D[npha-1][is*nsol+is] ) * fabs(scalprodct[ip][2]) );
          jat[ip][is][2] = ( ( alphidot[ip][is][2] * phiz[ip][3] ) / ( modgradphi[ip][3] ) )   * ( ( 1.0 - pfmdat->D[ip][is*nsol+is] /   pfmdat->D[npha-1][is*nsol+is] ) * fabs(scalprodct[ip][3]) );
          jat[ip][is][3] = ( ( alphidot[ip][is][3] * phiz[ip][4] ) / ( modgradphi[ip][4] ) )   * ( ( 1.0 - pfmdat->D[ip][is*nsol+is] /   pfmdat->D[npha-1][is*nsol+is] ) * fabs(scalprodct[ip][4]) );
          jat[ip][is][4] = ( ( alphidot[ip][is][4] * phix[ip][5] ) / ( modgradphi[ip][5] ) )   * ( ( 1.0 - pfmdat->D[ip][is*nsol+is] /   pfmdat->D[npha-1][is*nsol+is] ) * fabs(scalprodct[ip][5]) );
          jat[ip][is][5] = ( ( alphidot[ip][is][5] * phix[ip][6] ) / ( modgradphi[ip][6] ) )   * ( ( 1.0 - pfmdat->D[ip][is*nsol+is] /   pfmdat->D[npha-1][is*nsol+is] ) * fabs(scalprodct[ip][6]) );
      }
    }

    for ( ip = 0; ip < npha-1; ip++ ) {
      for ( is = 0; is < nsol; is++ ) {
          for (jj=0; jj< 6; jj++) {
            if ( modgradphi[ip][jj+1] == 0.0 ) {
              jat[ip][is][jj] = 0.0;
          }
        }
      }
    }

    for ( ip = 0; ip < npha-1; ip++ ) {
      for ( is = 0; is < nsol; is++ ) {
          cjaty = ( jat[ip][is][1] - jat[ip][is][0] ) / ( pfmvar->deltay );
          cjatz = ( jat[ip][is][3] - jat[ip][is][2] ) / ( pfmvar->deltaz );
          cjatx = ( jat[ip][is][5] - jat[ip][is][4] ) / ( pfmvar->deltax );
          jatc[ip][is] = cjatx + cjaty + cjatz;
	    }
    }

    for ( is = 0; is < nsol; is++ ) {
      jatr[is] = 0.0;
        for ( ip = 0; ip < npha-1; ip++ ) {
          jatr[is] += jatc[ip][is];
      }
    }



      for ( is1 = 0; is1 < nsol; is1++ ) {
        sum[is1] = 0.0;
        sum_dcbdT[is1] = 0.0;
        for ( is2 = 0; is2 < nsol; is2++ ) {
          dcdmu[is1*nsol+is2] = 0.0;
        }
      }

      for ( ip = 0; ip < npha; ip++ ) {
        for ( is = 0; is < nsol; is++ ) {
          dcbdT_phase[ip][is] = 0.0;
        }
      }

      // Calculations are carried out at stencil points of 5 point stencil
      for ( ii = 0; ii < 7; ii++ ) {

        interface = 1;
        bulkphase = 0;
        for ( ip = 0; ip < npha; ip++ ) {
          if ( stgO[ii].phi[ip] >= pfmdat->interfaceUplimit ) {
            bulkphase = ip;
            interface = 0;
            break;
          }
        }

        if ( interface ) {

           for ( ip = 0; ip < npha; ip++ ) {

            // Get dc_dmu
            tmp0 = 0.0;
            for ( is = 0; is < nsol; is++ ) {
              yc[is] = stcscl[ii].comie[ip][is];
              tmp0 +=  yc[is];
            }
            yc[nsol] = 1.0 - tmp0;

            //printf("t = %d, i = %d, j = %d, ii = %d, ip = %d, tdbip = %d", tstep[0], i, j, ii, ip, pfmdat->thermophase[ip]);
            dMudc(Ti, yc, retdmuphase, pfmdat->thermophase[ip]);

            matinvnew_nsol(retdmuphase, dc_dmu[ii][ip]);

          }

          // Calculate mobility at stencil points
          for ( is = 0; is < nsol; is++ ) {
            for ( js =0; js < nsol; js++ ) {
              Da[ii][is][js] = 0.0;
              for ( ip = 0; ip < npha; ip++ ) {
                for ( ks = 0; ks < nsol; ks++ ) {
                  Da[ii][is][js] += pfmdat->D[ip][is*nsol+ks] * dc_dmu[ii][ip][ks*nsol+js] * stgO[ii].phi[ip];
                }
              }
            }
          }

        }
        else {

          // Get dc_dmu for bulk phase
          tmp0 = 0.0;
          for ( is = 0; is < nsol; is++ ) {
            yc[is] = stgO[ii].com[is];
            tmp0 +=  yc[is];
          }
          yc[nsol] = 1.0 - tmp0;

          dMudc(Ti, yc, retdmuphase, pfmdat->thermophase[bulkphase]);

          matinvnew_nsol(retdmuphase, dc_dmu[ii][bulkphase]);

          //printf("dBC %d, %d, %d, %d, %le\n", tstep[0], i, j, ii, dc_dmu[ii][ip][0]);
          for ( is = 0; is < nsol; is++ ) {
            for ( js =0; js < nsol; js++ ) {
              Da[ii][is][js] = 0.0;
              //for ( ip = 0; ip < npha; ip++ ) {
                for ( ks = 0; ks < nsol; ks++ ) {
                  Da[ii][is][js] += pfmdat->D[bulkphase][is*nsol+ks] * dc_dmu[ii][bulkphase][ks*nsol+js];
                }
              //}
            }
          }

        }

      }

      // Calculating mobility at half grid positions
      for ( is = 0; is < nsol; is++ ) {
        for ( js =0; js < nsol; js++ ) {

          Damidy[0][is][js] = ( Da[0][is][js] + Da[1][is][js] ) / 2.0;
          Damidy[1][is][js] = ( Da[2][is][js] + Da[0][is][js] ) / 2.0;

          Damidz[0][is][js] = ( Da[0][is][js] + Da[3][is][js] ) / 2.0;
          Damidz[1][is][js] = ( Da[4][is][js] + Da[0][is][js] ) / 2.0;

          Damidx[0][is][js] = ( Da[0][is][js] + Da[5][is][js] ) / 2.0;
          Damidx[1][is][js] = ( Da[6][is][js] + Da[0][is][js] ) / 2.0;
        }
      }

      // Calculating grad of mu
      for ( is = 0; is < nsol; is++ ) {
        gradmuy[0][is] = ( stgO[0].mu[is] - stgO[1].mu[is] ) / pfmvar->deltay;
        gradmuy[1][is] = ( stgO[2].mu[is] - stgO[0].mu[is] ) / pfmvar->deltay;

        gradmuz[0][is] = ( stgO[0].mu[is] - stgO[3].mu[is] ) / pfmvar->deltaz;
        gradmuz[1][is] = ( stgO[4].mu[is] - stgO[0].mu[is] ) / pfmvar->deltaz;

        gradmux[0][is] = ( stgO[0].mu[is] - stgO[5].mu[is] ) / pfmvar->deltax;
        gradmux[1][is] = ( stgO[6].mu[is] - stgO[0].mu[is] ) / pfmvar->deltax;
      }

      // Calculating divergence of flux
      for ( is = 0; is < nsol; is++ ) {
        divflux[is] = 0.0;
        for ( js =0; js < nsol; js++ ) {

          divflux[is] =  ( Damidy[1][is][js] * gradmuy[1][js] - Damidy[0][is][js] * gradmuy[0][js] ) / pfmvar->deltay;

          divflux[is] += ( Damidz[1][is][js] * gradmuz[1][js] - Damidz[0][is][js] * gradmuz[0][js] ) / pfmvar->deltaz;

          divflux[is] += ( Damidx[1][is][js] * gradmux[1][js] - Damidx[0][is][js] * gradmux[0][js] ) / pfmvar->deltax;

        }
      }


      // Claculation are done at i,j for updating composition and mu

      // Check for interface
      interface = 1;
      bulkphase = 0;
      for ( ip = 0; ip < npha; ip++ ) {
        if ( stgO[0].phi[ip] >= pfmdat->interfaceUplimit ) {
          bulkphase = ip;
          interface = 0;
          break;
        }
      }

      if ( !interface ) {

        for ( is = 0; is < nsol; is++ ) {
          gridinfo[index].com[is] = gridinfoO[index].com[is] + pfmvar->deltat * ( divflux[is] + jatr[is] );
        }

        tmp0 = 0.0;
        for ( is = 0; is < nsol; is++ ) {
          yc[is] = gridinfo[index].com[is];
          tmp0 +=  yc[is];
        }
        yc[nsol] = 1.0 - tmp0;

        Mu(Ti, yc, mu, pfmdat->thermophase[bulkphase]);


        for ( is = 0; is < nsol; is++ ) {
          //deltamu[is] = mu[is] - gridinfoO[index].mu[is];
          gridinfo[index].mu[is] = mu[is];
        }

      }
      else {

        // Calculation are carried out in iterface region

        /*
        if ( !pfmdat->ISOTHERMAL ) {
          DELTAT = pfmvar->deltat * ( -pfmdat->TGRADIENT * pfmdat->velocity );
          for ( ip = 0; ip < npha; ip++ ) {

            for ( is = 0; is < nsol; is++ ) {
              cg[is] = pfmdat->cguess[ip*npha+ip][is];
              mu[is] = stgO[4].mu[is];
            }

            c_mu(mu, c_tdt, Ti+DELTAT, ip, cg, tstep[0], i, j, pfmdat->thermophase[ip]);

            for ( is = 0; is < nsol; is++ ) {
              dcbdT_phase[ip][is] = c_tdt[is] - stcscl[4].comie[ip][is];
            }
          }
        }
        */

        // Calculate Sum ( c_alpha * dh_alpha/ dt )
        for ( ip1 = 0; ip1 < npha; ip1++ ) {
          sum_dhphi = 0.0;
          for ( ip2 = 0; ip2 < npha; ip2++ ) {
            sum_dhphi += dhfhi(stgO[0].phi, ip1, ip2) * ( gridinfo[index].phi[ip2] - gridinfoO[index].phi[ip2] );
          }

          for ( is = 0; is < nsol; is++ ) {
            sum[is]       += stcscl[0].comie[ip1][is] * sum_dhphi;
            // sum_dcbdT[is] += dcbdT_phase[ip1][is] * hfhi(stgO[0].phi, ip1);
          }
        }


        for ( is1 = 0; is1 < nsol; is1++ ) {

          // Calculate deltac
          deltac[is1] = pfmvar->deltat * ( divflux[is1] + jatr[is1] );

          // Update composition using deltac
          gridinfo[index].com[is1] = gridinfoO[index].com[is1] + deltac[is1];

          // // Update deltac with c_alpha * dh_alpha/dt
          // if ( pfmdat->ISOTHERMAL ) {
            deltac[is1] += -sum[is1];
          // }
          // else {
          //  deltac[is1] += -sum[is] - sum_dcbdT[is1];
          // }
          for ( is2 = 0; is2 < nsol; is2++ ) {
            for ( ip = 0; ip < npha; ip++ ) {
              dcdmu[is1*nsol+is2] += dc_dmu[0][ip][is1*nsol+is2] * hfhi(stgO[0].phi, ip);
            }
          }
        }

        matinvnew_nsol(dcdmu, inv_dcdmu);

        multiply_nsol(inv_dcdmu, deltac, deltamu);
        // for ( is1 = 0; is1 < nsol; is1++ ) {
        //   tmp0 = 0.0;
        //   for ( is2 = 0; is2 < nsol; is2++ ) {
        //     tmp0 += inv_dcdmu[is1*nsol+is2] * deltac[is2];
        //   }
        //   deltamu[is1] = tmp0;
        // }

        //vectorsum(deltamu, stgO[index].mu, gridinfo[index].mu, nsol);
        for ( is = 0; is < nsol; is++ ) {
          gridinfo[index].mu[is] = gridinfoO[index].mu[is] + deltamu[is];
        }
        //printf("IC %d, %d, %d, %le, %le\n", tstep[0], i, j, gridinfo[index].mu[0], gridinfo[index].com[0]);
        //printf("IC %d, %d, %d, %le, %le, %le, %le\n", tstep[0], i, j, gridinfoO[index].mu[0], gridinfo[index].mu[0], gridinfoO[index].com[0], gridinfo[index].com[0]);
      }

  }
}

__kernel void SolverCatr_F3_smooth(__global struct fields *gridinfoO, __global struct fields *gridinfo, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf3 *propf3) {

  int x, y, z, ii, jj;
  int nx, ny, nz;
  int index;

  int is, js, ks, ip, ip1, ip2, is1, is2;
  int interface, bulkphase;
  double tmp0, tmp1;
  double DELTAT, sum_dhphi, Ti;
  double sum[nsol], sum_dcbdT[nsol], dcdmu[nsol*nsol], dc_dmu[7][npha][nsol][nsol], Da[7][nsol][nsol];
  double Damidx[2][nsol][nsol], Damidy[2][nsol][nsol], Damidz[2][nsol][nsol], gradmux[2][nsol], gradmuy[2][nsol], gradmuz[2][nsol];
  double divflux[nsol], mu[nsol], deltamu[nsol], c_tdt[nsol], dcbdT_phase[npha][nsol];
  double deltac[nsol], inv_dcdmu[nsol*nsol], cg[nsol];
  double yc[nsol+1], retdmuphase[nsol*nsol];
  double F3_A[npha*nsol*nsol], F3_Beq[npha*nsol], F3_dBbdT[npha*nsol], Tequ, comp[nsol];

  double alpha[npha-1][nsol][7], alphidot[npha-1][nsol][6], jat[npha-1][nsol][6];
  double A1[npha-1];
  double dphidt[npha-1][7], gradx_phi[npha][5], grady_phi[npha][5], gradz_phi[npha][5], phix[npha][7], phiy[npha][7], phiz[npha][7], modgradphi[npha][7], scalprodct[npha][7];
  double cjatx, cjaty, cjatz;
  double jatc[npha-1][nsol], jatr[nsol];

  struct fields stgridO[27];
  struct fields stgN[7];
  struct csle stcscl[7];
  struct fields stgO[7];

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  //index = y + ny*(z + x*nz);

  if ( x != 0 && x != nx-1 && y != 0 && y != ny-1 && z != 0 && z != nz-1 ) {

    index = y + ny*(z + x*nz);

    Ti = pfmdat->T0;

    stgridO[0]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[1]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[2]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[3]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[4]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[5]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[6]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[7]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[8]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[9]  = gridinfoO[ (y-1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[10] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[11] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[12] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[13] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[14] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[15] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[16] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[17] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[18] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[19] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[20] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[21] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[22] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[23] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[24] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[25] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[26] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z+1) ) ];

    stcscl[0] = cscl[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stcscl[1] = cscl[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stcscl[2] = cscl[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stcscl[3] = cscl[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stcscl[4] = cscl[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stcscl[5] = cscl[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stcscl[6] = cscl[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];

    stgO[0] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stgO[1] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stgO[2] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stgO[3] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stgO[4] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stgO[5] = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stgO[6] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];


      for ( is1 = 0; is1 < nsol; is1++ ) {
        sum[is1] = 0.0;
        sum_dcbdT[is1] = 0.0;
        for ( is2 = 0; is2 < nsol; is2++ ) {
          dcdmu[is1*nsol+is2] = 0.0;
        }
      }

      for ( ip = 0; ip < npha; ip++ ) {
        for ( is = 0; is < nsol; is++ ) {
          dcbdT_phase[ip][is] = 0;
        }
      }

      // Calculations are carried out at stencil points of 7 point stencil
      for ( ii = 0; ii < 7; ii++ ) {

        interface = 1;
        bulkphase = 0;
        for ( ip = 0; ip < npha; ip++ ) {
          if ( stgO[ii].phi[ip] >= pfmdat->interfaceUplimit ) {
            bulkphase = ip;
            interface = 0;
            break;
          }
        }

        if ( interface ) {

          // Get dc_dmu
          for ( ip = 0; ip < npha; ip++ ) {
            for ( is1 = 0; is1 < nsol; is1++ ) {
              for ( is2 =0; is2 < nsol; is2++ ) {
                dc_dmu[ii][ip][is1][is2] = propf3->cmu[ip][is1][is2];
              }
            }
          }

          // Calculate mobility at stencil points
          for ( is = 0; is < nsol; is++ ) {
            for ( js =0; js < nsol; js++ ) {
              Da[ii][is][js] = 0.0;
              for ( ip = 0; ip < npha; ip++ ) {
                for ( ks = 0; ks < nsol; ks++ ) {
                  //Da[ii][is][js] += pfmdat->D[ip][is*nsol+ks] * dc_dmu[ii][ip][ks][js] * hfhi(stgO[ii].phi, ip);
                  Da[ii][is][js] += pfmdat->D[ip][is*nsol+ks] * dc_dmu[ii][ip][ks][js] * stgO[ii].phi[ip];
                }
              }
            }
          }

        }
        else {

          // Get dc_dmu for bulk phase
          for ( is1 = 0; is1 < nsol; is1++ ) {
            for ( is2 = 0; is2 < nsol; is2++ ) {
              dc_dmu[ii][bulkphase][is1][is2] = propf3->cmu[bulkphase][is1][is2];
            }
          }

          // Calculate mobility for bulk phase at stencil points
          for ( is = 0; is < nsol; is++ ) {
            for ( js =0; js < nsol; js++ ) {
              Da[ii][is][js] = 0.0;

                for ( ks = 0; ks < nsol; ks++ ) {
                  Da[ii][is][js] += pfmdat->D[bulkphase][is*nsol+ks] * dc_dmu[ii][bulkphase][ks][js];
                }

            }
          }

        }

      }

      // Calculating mobility at half grid positions
      for ( is = 0; is < nsol; is++ ) {
        for ( js =0; js < nsol; js++ ) {

          Damidy[0][is][js] = ( Da[0][is][js] + Da[1][is][js] ) / 2.0;
          Damidy[1][is][js] = ( Da[2][is][js] + Da[0][is][js] ) / 2.0;

          Damidz[0][is][js] = ( Da[0][is][js] + Da[3][is][js] ) / 2.0;
          Damidz[1][is][js] = ( Da[4][is][js] + Da[0][is][js] ) / 2.0;

          Damidx[0][is][js] = ( Da[0][is][js] + Da[5][is][js] ) / 2.0;
          Damidx[1][is][js] = ( Da[6][is][js] + Da[0][is][js] ) / 2.0;

          //Damidx[0][is][js] = ( Da[2][is][js] + Da[4][is][js] ) / 2.0;
          //Damidx[1][is][js] = ( Da[4][is][js] + Da[1][is][js] ) / 2.0;

          //Damidy[0][is][js] = ( Da[3][is][js] + Da[4][is][js] ) / 2.0;
          //Damidy[1][is][js] = ( Da[4][is][js] + Da[0][is][js] ) / 2.0;
        }
      }

      // Calculating grad of mu
      for ( is = 0; is < nsol; is++ ) {
        gradmuy[0][is] = ( stgO[0].mu[is] - stgO[1].mu[is] ) / pfmvar->deltay;
        gradmuy[1][is] = ( stgO[2].mu[is] - stgO[0].mu[is] ) / pfmvar->deltay;

        gradmuz[0][is] = ( stgO[0].mu[is] - stgO[3].mu[is] ) / pfmvar->deltaz;
        gradmuz[1][is] = ( stgO[4].mu[is] - stgO[0].mu[is] ) / pfmvar->deltaz;

        gradmux[0][is] = ( stgO[0].mu[is] - stgO[5].mu[is] ) / pfmvar->deltax;
        gradmux[1][is] = ( stgO[6].mu[is] - stgO[0].mu[is] ) / pfmvar->deltax;

        //gradmux[0][is] = ( stgO[2].mu[is] - stgO[4].mu[is] ) / pfmvar->deltax;
        //gradmux[1][is] = ( stgO[4].mu[is] - stgO[1].mu[is] ) / pfmvar->deltax;
        //gradmuy[0][is] = ( stgO[3].mu[is] - stgO[4].mu[is] ) / pfmvar->deltay;
        //gradmuy[1][is] = ( stgO[4].mu[is] - stgO[0].mu[is] ) / pfmvar->deltay;
      }

      // Calculating divergence of flux
      for ( is = 0; is < nsol; is++ ) {
        divflux[is] = 0.0;
        for ( js =0; js < nsol; js++ ) {

          divflux[is] =  ( Damidy[1][is][js] * gradmuy[1][js] - Damidy[0][is][js] * gradmuy[0][js] ) / pfmvar->deltay;

          divflux[is] += ( Damidz[1][is][js] * gradmuz[1][js] - Damidz[0][is][js] * gradmuz[0][js] ) / pfmvar->deltaz;

          divflux[is] += ( Damidx[1][is][js] * gradmux[1][js] - Damidx[0][is][js] * gradmux[0][js] ) / pfmvar->deltax;

          //divflux[is] += ( Damidx[0][is][js] * gradmux[0][js] - Damidx[1][is][js] * gradmux[1][js] ) / pfmvar->deltax;

          //divflux[is] += ( Damidy[0][is][js] * gradmuy[0][js] - Damidy[1][is][js] * gradmuy[1][js] ) / pfmvar->deltay;

        }
      }


      // Claculation are done at x,y,z for updating composition and mu

      // Check for interface
      interface = 1;
      bulkphase = 0;
      for ( ip = 0; ip < npha; ip++ ) {
        if ( stgO[0].phi[ip] >= pfmdat->interfaceUplimit ) {
          bulkphase = ip;
          interface = 0;
          break;
        }
      }



          for ( ip = 0; ip < npha; ip++ ) {

            for ( is1 = 0; is1 < nsol; is1++ ) {

              F3_Beq[ip*nsol+is1] = propf3->Beq[ip][is1];

              F3_dBbdT[ip*nsol+is1] = propf3->dBbdT[ip][is1];

              for ( is2 = 0; is2 < nsol; is2++ ) {
                F3_A[(ip*nsol+is1)*nsol+is2] = propf3->A[ip][is1][is2];
              }

            }
          }
          Tequ = pfmdat->Teq;

      // Update composition and mu in bulk phase
      if ( !interface ) {

        // update composition
        for ( is = 0; is < nsol; is++ ) {
          // gridinfo[index].com[is] = gridinfoO[index].com[is] + pfmvar->deltat * ( divflux[is] + jatr[is] );
          gridinfo[index].com[is] = gridinfoO[index].com[is] + pfmvar->deltat * ( divflux[is] );
          comp[is] = gridinfo[index].com[is];
        }


          function_F_03_Mu_CL(comp, Ti, bulkphase, mu, F3_A, F3_Beq, F3_dBbdT, Tequ);

        // Updating mu
        for ( is = 0; is < nsol; is++ ) {
          gridinfo[index].mu[is] = mu[is];
        }

      }
      else {

        // Calculation are carried out in iterface region

        for ( is1 = 0; is1 < nsol; is1++ ) {

          // Calculate deltac
          deltac[is1] = pfmvar->deltat * ( divflux[is1] );

          //// Update composition using deltac
          //gridinfo[index].com[is1] = stgO[0].com[is1] + deltac[is1];

          //// Update deltac with c_alpha * dh_alpha/dt
          //deltac[is1] += -sum[is1];


          // Calculating h_alpha * dc_dmu
          for ( is2 = 0; is2 < nsol; is2++ ) {
            for ( ip = 0; ip < npha; ip++ ) {
              dcdmu[is1*nsol+is2] += dc_dmu[0][ip][is1][is2] * hfhi(stgO[0].phi, ip);
            }
          }

        }

        // printf("t=%d, %d=i, %d=j, sum=%le\n", tstep[0], i, j, sum[0]);

        //printf("t=%d, %d=i, %d=j, phi0=%le, phi1=%le, h0=%le, h1=%le, dh00=%le, dh01=%le, dh10=%le, dh11=%le\n", tstep[0], i, j, stgO[4].phi[0], stgO[4].phi[1], hfhi(stgO[4].phi, 0), hfhi(stgO[4].phi, 1), dhfhi(stgO[4].phi, 0, 0), dhfhi(stgO[4].phi, 0, 1), dhfhi(stgO[4].phi, 1, 0), dhfhi(stgO[4].phi, 1, 1));

        // Inverse of dcdmu
        matinvnew_nsol(dcdmu, inv_dcdmu);

        multiply_nsol(inv_dcdmu, deltac, deltamu);

        //vectorsum(deltamu, stgO[index].mu, gridinfo[index].mu, nsol);
        for ( is = 0; is < nsol; is++ ) {
          gridinfo[index].mu[is] = gridinfoO[index].mu[is] + deltamu[is];
        }
        //printf("IC %d, %d, %d, %le, %le, %le, %le\n", tstep[0], i, j, gridinfoO[index].mu[0], gridinfo[index].mu[0], gridinfoO[index].com[0], gridinfo[index].com[0]);
      }

    }
}

__kernel void SolverCatr_F3(__global struct fields *gridinfoO, __global struct fields *gridinfo, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf3 *propf3) { 
  
  int x, y, z, ii, jj;
  int nx, ny, nz;
  int index;

  int is, js, ks, ip, ip1, ip2, is1, is2;
  int interface, bulkphase;
  double tmp0, tmp1;
  double DELTAT, sum_dhphi, Ti;
  double sum[nsol], sum_dcbdT[nsol], dcdmu[nsol*nsol], dc_dmu[7][npha][nsol][nsol], Da[7][nsol][nsol];
  double Damidx[2][nsol][nsol], Damidy[2][nsol][nsol], Damidz[2][nsol][nsol], gradmux[2][nsol], gradmuy[2][nsol], gradmuz[2][nsol];
  double divflux[nsol], mu[nsol], deltamu[nsol], c_tdt[nsol], dcbdT_phase[npha][nsol];
  double deltac[nsol], inv_dcdmu[nsol*nsol], cg[nsol];
  double yc[nsol+1], retdmuphase[nsol*nsol];
  double F3_A[npha*nsol*nsol], F3_Beq[npha*nsol], F3_dBbdT[npha*nsol], Tequ, comp[nsol];

  double alpha[npha-1][nsol][7], alphidot[npha-1][nsol][6], jat[npha-1][nsol][6];
  double A1[npha-1];
  double dphidt[npha-1][7], gradx_phi[npha][5], grady_phi[npha][5], gradz_phi[npha][5], phix[npha][7], phiy[npha][7], phiz[npha][7], modgradphi[npha][7], scalprodct[npha][7];
  double cjatx, cjaty, cjatz;
  double jatc[npha-1][nsol], jatr[nsol];

  struct fields stgridO[27];
  struct fields stgN[7];
  struct csle stcscl[7];
  struct fields stgO[7];
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  //index = y + ny*(z + x*nz);

  if ( x != 0 && x != nx-1 && y != 0 && y != ny-1 && z != 0 && z != nz-1 ) {

    index = y + ny*(z + x*nz);

    Ti = pfmdat->T0;
    
    stgridO[0]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[1]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[2]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[3]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[4]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[5]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[6]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[7]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[8]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[9]  = gridinfoO[ (y-1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[10] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[11] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[12] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[13] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[14] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[15] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[16] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[17] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[18] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[19] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[20] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[21] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[22] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[23] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[24] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[25] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[26] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z+1) ) ];
    
    // stgridO[0] = gridinfoO[(i-1)*ny + (j-1)];
    // stgridO[1] = gridinfoO[(i-1)*ny + ( j )];
    // stgridO[2] = gridinfoO[(i-1)*ny + (j+1)];
    // stgridO[3] = gridinfoO[( i )*ny + (j-1)];
    // stgridO[4] = gridinfoO[( i )*ny + ( j )];
    // stgridO[5] = gridinfoO[( i )*ny + (j+1)];
    // stgridO[6] = gridinfoO[(i+1)*ny + (j-1)];
    // stgridO[7] = gridinfoO[(i+1)*ny + ( j )];
    // stgridO[8] = gridinfoO[(i+1)*ny + (j+1)];
    
    stcscl[0] = cscl[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stcscl[1] = cscl[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stcscl[2] = cscl[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stcscl[3] = cscl[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stcscl[4] = cscl[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stcscl[5] = cscl[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stcscl[6] = cscl[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];
    
    // stcscl[0] = cscl[(i-1)*ny + ( j )];
    // stcscl[1] = cscl[( i )*ny + (j-1)];
    // stcscl[4] = cscl[( i )*ny + ( j )];
    // stcscl[2] = cscl[( i )*ny + (j+1)];
    // stcscl[3] = cscl[(i+1)*ny + ( j )];

    stgO[0] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stgO[1] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stgO[2] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stgO[3] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stgO[4] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stgO[5] = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stgO[6] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];

    // stgO[0] = gridinfoO[(i-1)*ny + ( j )];
    // stgO[1] = gridinfoO[( i )*ny + (j-1)];
    // stgO[4] = gridinfoO[( i )*ny + ( j )];
    // stgO[2] = gridinfoO[( i )*ny + (j+1)];
    // stgO[3] = gridinfoO[(i+1)*ny + ( j )];

    //if (pfmdat->atr) {
      
      for ( ip = 0; ip < npha-1; ip++ ) {
        stgN[0].phi[ip] = gridinfo[ (y  ) + ny*( (x  )*nz + (z  ) ) ].phi[ip];
        stgN[1].phi[ip] = gridinfo[ (y-1) + ny*( (x  )*nz + (z  ) ) ].phi[ip];
        stgN[2].phi[ip] = gridinfo[ (y+1) + ny*( (x  )*nz + (z  ) ) ].phi[ip];
        stgN[3].phi[ip] = gridinfo[ (y  ) + ny*( (x  )*nz + (z-1) ) ].phi[ip];
        stgN[4].phi[ip] = gridinfo[ (y  ) + ny*( (x  )*nz + (z+1) ) ].phi[ip];
        stgN[5].phi[ip] = gridinfo[ (y  ) + ny*( (x-1)*nz + (z  ) ) ].phi[ip];
        stgN[6].phi[ip] = gridinfo[ (y  ) + ny*( (x+1)*nz + (z  ) ) ].phi[ip];
        //stgridN[0].phi[ip] = gridinfo[(i-1)*ny + ( j )].phi[ip];
        //stgridN[1].phi[ip] = gridinfo[( i )*ny + (j-1)].phi[ip];
        //stgridN[4].phi[ip] = gridinfo[( i )*ny + ( j )].phi[ip];
        //stgridN[2].phi[ip] = gridinfo[( i )*ny + (j+1)].phi[ip];
        //stgridN[3].phi[ip] = gridinfo[(i+1)*ny + ( j )].phi[ip];
      }
      
      for ( ip = 0; ip < npha-1; ip++ ) {
        dphidt[ip][0] = ( stgN[0].phi[ip] - stgO[0].phi[ip] ) / ( pfmvar->deltat );
        dphidt[ip][1] = ( stgN[1].phi[ip] - stgO[1].phi[ip] ) / ( pfmvar->deltat );
        dphidt[ip][2] = ( stgN[2].phi[ip] - stgO[2].phi[ip] ) / ( pfmvar->deltat );
        dphidt[ip][3] = ( stgN[3].phi[ip] - stgO[3].phi[ip] ) / ( pfmvar->deltat );
        dphidt[ip][4] = ( stgN[4].phi[ip] - stgO[4].phi[ip] ) / ( pfmvar->deltat );
        dphidt[ip][5] = ( stgN[5].phi[ip] - stgO[5].phi[ip] ) / ( pfmvar->deltat );
        dphidt[ip][6] = ( stgN[6].phi[ip] - stgO[6].phi[ip] ) / ( pfmvar->deltat );
        
        // dphidt[ip][0] = ( stgridN[4].phi[ip] - stgridO[4].phi[ip] ) / ( pfmvar->deltat );
        // dphidt[ip][1] = ( stgridN[2].phi[ip] - stgridO[5].phi[ip] ) / ( pfmvar->deltat );
        // dphidt[ip][2] = ( stgridN[1].phi[ip] - stgridO[3].phi[ip] ) / ( pfmvar->deltat );
        // dphidt[ip][3] = ( stgridN[3].phi[ip] - stgridO[7].phi[ip] ) / ( pfmvar->deltat );
        // dphidt[ip][4] = ( stgridN[0].phi[ip] - stgridO[1].phi[ip] ) / ( pfmvar->deltat );
  
        A1[ip] = sqrt ( pfmvar->ee[ip*npha+(npha-1)] ) / sqrt( pfmvar->w[ip*npha+(npha-1)] );
      }
    


    
      for ( ip = 0; ip < npha-1; ip++ ) {
        for ( is = 0; is < nsol; is++ ) {
          alpha[ip][is][0] = (A1[ip])*( stcscl[0].comie[npha-1][is] - stcscl[0].comie[ip][is] );
          alpha[ip][is][1] = (A1[ip])*( stcscl[1].comie[npha-1][is] - stcscl[1].comie[ip][is] );
          alpha[ip][is][2] = (A1[ip])*( stcscl[2].comie[npha-1][is] - stcscl[2].comie[ip][is] );
          alpha[ip][is][3] = (A1[ip])*( stcscl[3].comie[npha-1][is] - stcscl[3].comie[ip][is] );
          alpha[ip][is][4] = (A1[ip])*( stcscl[4].comie[npha-1][is] - stcscl[4].comie[ip][is] );
          alpha[ip][is][5] = (A1[ip])*( stcscl[5].comie[npha-1][is] - stcscl[5].comie[ip][is] );
          alpha[ip][is][6] = (A1[ip])*( stcscl[6].comie[npha-1][is] - stcscl[6].comie[ip][is] );
          
          //alpha[ip][is][0] = (A1[ip])*( stcscl[4].comie[npha-1][is] - stcscl[4].comie[ip][is] );
          //alpha[ip][is][1] = (A1[ip])*( stcscl[2].comie[npha-1][is] - stcscl[2].comie[ip][is] );
          //alpha[ip][is][2] = (A1[ip])*( stcscl[1].comie[npha-1][is] - stcscl[1].comie[ip][is] );
          //alpha[ip][is][3] = (A1[ip])*( stcscl[3].comie[npha-1][is] - stcscl[3].comie[ip][is] );
          //alpha[ip][is][4] = (A1[ip])*( stcscl[0].comie[npha-1][is] - stcscl[0].comie[ip][is] );
        }
      }
  
      for ( ip = 0; ip < npha; ip++ ) {
        grady_phi[ip][0] = ( stgridO[14].phi[ip] - stgridO[12].phi[ip] ) / ( 2.0*pfmvar->deltay );
        grady_phi[ip][1] = ( stgridO[11].phi[ip] - stgridO[9 ].phi[ip] ) / ( 2.0*pfmvar->deltay );
        grady_phi[ip][2] = ( stgridO[17].phi[ip] - stgridO[15].phi[ip] ) / ( 2.0*pfmvar->deltay );
        grady_phi[ip][3] = ( stgridO[5 ].phi[ip] - stgridO[3 ].phi[ip] ) / ( 2.0*pfmvar->deltay );
        grady_phi[ip][4] = ( stgridO[23].phi[ip] - stgridO[21].phi[ip] ) / ( 2.0*pfmvar->deltay );
        
        gradz_phi[ip][0] = ( stgridO[16].phi[ip] - stgridO[10].phi[ip] ) / ( 2.0*pfmvar->deltaz );
        gradz_phi[ip][1] = ( stgridO[15].phi[ip] - stgridO[9 ].phi[ip] ) / ( 2.0*pfmvar->deltaz );
        gradz_phi[ip][2] = ( stgridO[17].phi[ip] - stgridO[11].phi[ip] ) / ( 2.0*pfmvar->deltaz );
        gradz_phi[ip][3] = ( stgridO[7 ].phi[ip] - stgridO[1 ].phi[ip] ) / ( 2.0*pfmvar->deltaz );
        gradz_phi[ip][4] = ( stgridO[25].phi[ip] - stgridO[19].phi[ip] ) / ( 2.0*pfmvar->deltaz );
        
        gradx_phi[ip][0] = ( stgridO[22].phi[ip] - stgridO[4 ].phi[ip] ) / ( 2.0*pfmvar->deltax );
        gradx_phi[ip][1] = ( stgridO[21].phi[ip] - stgridO[3 ].phi[ip] ) / ( 2.0*pfmvar->deltax );
        gradx_phi[ip][2] = ( stgridO[23].phi[ip] - stgridO[5 ].phi[ip] ) / ( 2.0*pfmvar->deltax );
        gradx_phi[ip][3] = ( stgridO[19].phi[ip] - stgridO[1 ].phi[ip] ) / ( 2.0*pfmvar->deltax );
        gradx_phi[ip][4] = ( stgridO[25].phi[ip] - stgridO[7 ].phi[ip] ) / ( 2.0*pfmvar->deltax );
      }
  
      for ( ip = 0; ip < npha; ip++ ) {
  
        phiy[ip][0] = grady_phi[ip][0];
        phiy[ip][1] = ( stgridO[13].phi[ip] - stgridO[12].phi[ip] ) / ( pfmvar->deltay );
        phiy[ip][2] = ( stgridO[14].phi[ip] - stgridO[13].phi[ip] ) / ( pfmvar->deltay );
        phiy[ip][3] = ( grady_phi[ip][0] + grady_phi[ip][1] ) / ( 2.0 );
        phiy[ip][4] = ( grady_phi[ip][0] + grady_phi[ip][2] ) / ( 2.0 );
        phiy[ip][5] = ( grady_phi[ip][0] + grady_phi[ip][3] ) / ( 2.0 );
        phiy[ip][6] = ( grady_phi[ip][0] + grady_phi[ip][4] ) / ( 2.0 );
        
        phiz[ip][0] = gradz_phi[ip][0];
        phiz[ip][1] = ( gradz_phi[ip][0] + gradz_phi[ip][1] ) / ( 2.0 );
        phiz[ip][2] = ( gradz_phi[ip][0] + gradz_phi[ip][2] ) / ( 2.0 );
        phiz[ip][3] = ( stgridO[13].phi[ip] - stgridO[10].phi[ip] ) / ( pfmvar->deltaz );
        phiz[ip][4] = ( stgridO[16].phi[ip] - stgridO[13].phi[ip] ) / ( pfmvar->deltaz );
        phiz[ip][5] = ( gradz_phi[ip][0] + gradz_phi[ip][3] ) / ( 2.0 );
        phiz[ip][6] = ( gradz_phi[ip][0] + gradz_phi[ip][4] ) / ( 2.0 );
        
        phix[ip][0] = gradx_phi[ip][0];
        phix[ip][1] = ( gradx_phi[ip][0] + gradx_phi[ip][1] ) / ( 2.0 );
        phix[ip][2] = ( gradx_phi[ip][0] + gradx_phi[ip][2] ) / ( 2.0 );
        phix[ip][3] = ( gradx_phi[ip][0] + gradx_phi[ip][3] ) / ( 2.0 );
        phix[ip][4] = ( gradx_phi[ip][0] + gradx_phi[ip][4] ) / ( 2.0 );
        phix[ip][5] = ( stgridO[13].phi[ip] - stgridO[4 ].phi[ip] ) / ( pfmvar->deltax );
        phix[ip][6] = ( stgridO[22].phi[ip] - stgridO[13].phi[ip] ) / ( pfmvar->deltax );
        
        
        // phix[ip][0] = gradx_phi[ip][0];
        // phix[ip][1] = ( stgridO[5].phi[ip] - stgridO[4].phi[ip] ) / ( pfmvar->deltax );
        // phix[ip][2] = ( stgridO[4].phi[ip] - stgridO[3].phi[ip] ) / ( pfmvar->deltax );
        // phix[ip][3] = ( gradx_phi[ip][0] + gradx_phi[ip][1] ) / ( 2.0 );
        // phix[ip][4] = ( gradx_phi[ip][0] + gradx_phi[ip][2] ) / ( 2.0 );
        // 
        // phiy[ip][0] = grady_phi[ip][0];
        // phiy[ip][1] = ( grady_phi[ip][0] + grady_phi[ip][1] ) / ( 2.0 );
        // phiy[ip][2] = ( grady_phi[ip][0] + grady_phi[ip][2] ) / ( 2.0 );
        // phiy[ip][3] = ( stgridO[7].phi[ip] - stgridO[4].phi[ip] ) / ( pfmvar->deltay );
        // phiy[ip][4] = ( stgridO[4].phi[ip] - stgridO[1].phi[ip] ) / ( pfmvar->deltay );
      }
  
      for ( ip = 0; ip < npha-1; ip++ ) {
	      for ( is = 0; is < nsol; is++ ) {
          alphidot[ip][is][0] = ( ( alpha[ip][is][0] * dphidt[ip][0] ) + ( alpha[ip][is][1] * dphidt[ip][1] ) ) / ( 2.0 );
          alphidot[ip][is][1] = ( ( alpha[ip][is][0] * dphidt[ip][0] ) + ( alpha[ip][is][2] * dphidt[ip][2] ) ) / ( 2.0 );
          alphidot[ip][is][2] = ( ( alpha[ip][is][0] * dphidt[ip][0] ) + ( alpha[ip][is][3] * dphidt[ip][3] ) ) / ( 2.0 );
          alphidot[ip][is][3] = ( ( alpha[ip][is][0] * dphidt[ip][0] ) + ( alpha[ip][is][4] * dphidt[ip][4] ) ) / ( 2.0 );
          alphidot[ip][is][4] = ( ( alpha[ip][is][0] * dphidt[ip][0] ) + ( alpha[ip][is][5] * dphidt[ip][5] ) ) / ( 2.0 );
          alphidot[ip][is][5] = ( ( alpha[ip][is][0] * dphidt[ip][0] ) + ( alpha[ip][is][6] * dphidt[ip][6] ) ) / ( 2.0 );
	      }
      }
  
      for ( ip = 0; ip < npha; ip++ ) {
        for ( jj = 0; jj < 7; jj++ ) {
          modgradphi[ip][jj] = sqrt( phiy[ip][jj]*phiy[ip][jj] + phiz[ip][jj]*phiz[ip][jj] + phix[ip][jj]*phix[ip][jj]);
        }
      }
  
      for ( ip = 0; ip < npha-1; ip++ ) {
        for ( jj = 0; jj < 7; jj++ ) {
          scalprodct[ip][jj] = -1.0*( phiy[ip][jj]*phiy[npha-1][jj] + phiz[ip][jj]*phiz[npha-1][jj] + phix[ip][jj]*phix[npha-1][jj] );
          if ( modgradphi[npha-1][jj] > 0.0 ) {
            scalprodct[ip][jj] /= ( modgradphi[ip][jj] * modgradphi[npha-1][jj] );
          }
        }
      }
  
  
  
      for ( ip = 0; ip < npha-1; ip++ ) {
	      for ( is = 0; is < nsol; is++ ) {
          jat[ip][is][0] = ( ( alphidot[ip][is][0] * phiy[ip][1] ) / ( modgradphi[ip][1] ) )   * ( ( 1.0 - pfmdat->D[ip][is*nsol+is] /   pfmdat->D[npha-1][is*nsol+is] ) * fabs(scalprodct[ip][1]) );
          jat[ip][is][1] = ( ( alphidot[ip][is][1] * phiy[ip][2] ) / ( modgradphi[ip][2] ) )   * ( ( 1.0 - pfmdat->D[ip][is*nsol+is] /   pfmdat->D[npha-1][is*nsol+is] ) * fabs(scalprodct[ip][2]) );
          jat[ip][is][2] = ( ( alphidot[ip][is][2] * phiz[ip][3] ) / ( modgradphi[ip][3] ) )   * ( ( 1.0 - pfmdat->D[ip][is*nsol+is] /   pfmdat->D[npha-1][is*nsol+is] ) * fabs(scalprodct[ip][3]) );
          jat[ip][is][3] = ( ( alphidot[ip][is][3] * phiz[ip][4] ) / ( modgradphi[ip][4] ) )   * ( ( 1.0 - pfmdat->D[ip][is*nsol+is] /   pfmdat->D[npha-1][is*nsol+is] ) * fabs(scalprodct[ip][4]) );
          jat[ip][is][4] = ( ( alphidot[ip][is][4] * phix[ip][5] ) / ( modgradphi[ip][5] ) )   * ( ( 1.0 - pfmdat->D[ip][is*nsol+is] /   pfmdat->D[npha-1][is*nsol+is] ) * fabs(scalprodct[ip][5]) );
          jat[ip][is][5] = ( ( alphidot[ip][is][5] * phix[ip][6] ) / ( modgradphi[ip][6] ) )   * ( ( 1.0 - pfmdat->D[ip][is*nsol+is] /   pfmdat->D[npha-1][is*nsol+is] ) * fabs(scalprodct[ip][6]) );
        
          //jat[ip][is][0] = ( ( alphidot[ip][is][0] * phix[ip][1] ) / ( modgradphi[ip][1] ) )   * ( ( 1.0 - pfmdat->D[ip][is*nsol+is] /   pfmdat->D[npha-1][is*nsol+is] ) * fabs(scalprodct[ip][1]) );
          //jat[ip][is][1] = ( ( alphidot[ip][is][1] * phix[ip][2] ) / ( modgradphi[ip][2] ) )   * ( ( 1.0 - pfmdat->D[ip][is*nsol+is] /   pfmdat->D[npha-1][is*nsol+is] ) * fabs(scalprodct[ip][2]) );
          //jat[ip][is][2] = ( ( alphidot[ip][is][2] * phiy[ip][3] ) / ( modgradphi[ip][3] ) )   * ( ( 1.0 - pfmdat->D[ip][is*nsol+is] /   pfmdat->D[npha-1][is*nsol+is] ) * fabs(scalprodct[ip][3]) );
          //jat[ip][is][3] = ( ( alphidot[ip][is][3] * phiy[ip][4] ) / ( modgradphi[ip][4] ) )   * ( ( 1.0 - pfmdat->D[ip][is*nsol+is] /   pfmdat->D[npha-1][is*nsol+is] ) * fabs(scalprodct[ip][4]) );
        }
      }
  
      for ( ip = 0; ip < npha-1; ip++ ) {
        for ( is = 0; is < nsol; is++ ) {
          for (jj=0; jj< 6; jj++) {
            if ( modgradphi[ip][jj+1] == 0.0 ) {
              jat[ip][is][jj] = 0.0;
            }
          }
        }
      }
  
      for ( ip = 0; ip < npha-1; ip++ ) {
        for ( is = 0; is < nsol; is++ ) {
          cjaty = ( jat[ip][is][1] - jat[ip][is][0] ) / ( pfmvar->deltay );
          cjatz = ( jat[ip][is][3] - jat[ip][is][2] ) / ( pfmvar->deltaz );
          cjatx = ( jat[ip][is][5] - jat[ip][is][4] ) / ( pfmvar->deltax );
          jatc[ip][is] = cjatx + cjaty + cjatz;
          
          //cjatx = ( jat[ip][is][0] - jat[ip][is][1] ) / ( pfmvar->deltax );
          //cjaty = ( jat[ip][is][2] - jat[ip][is][3] ) / ( pfmvar->deltay );
          //jatc[ip][is] = cjatx + cjaty;
          
        }
      }
  
      for ( is = 0; is < nsol; is++ ) {
        jatr[is] = 0.0;
        for ( ip = 0; ip < npha-1; ip++ ) {
          jatr[is] += jatc[ip][is];
        }
      }
    //}
  

      for ( is1 = 0; is1 < nsol; is1++ ) {
        sum[is1] = 0.0;
        sum_dcbdT[is1] = 0.0;
        for ( is2 = 0; is2 < nsol; is2++ ) {
          dcdmu[is1*nsol+is2] = 0.0;
        }
      }

      for ( ip = 0; ip < npha; ip++ ) {
        for ( is = 0; is < nsol; is++ ) {
          dcbdT_phase[ip][is] = 0;
        }
      }

      // Calculations are carried out at stencil points of 7 point stencil
      for ( ii = 0; ii < 7; ii++ ) {

        interface = 1;
        bulkphase = 0;
        for ( ip = 0; ip < npha; ip++ ) {
          if ( stgO[ii].phi[ip] >= pfmdat->interfaceUplimit ) {
            bulkphase = ip;
            interface = 0;
            break;
          }
        }

        if ( interface ) {

          // Get dc_dmu
          for ( ip = 0; ip < npha; ip++ ) {
            for ( is1 = 0; is1 < nsol; is1++ ) {
              for ( is2 =0; is2 < nsol; is2++ ) {
                dc_dmu[ii][ip][is1][is2] = propf3->cmu[ip][is1][is2];
              }
            }
          }

          // Calculate mobility at stencil points
          for ( is = 0; is < nsol; is++ ) {
            for ( js =0; js < nsol; js++ ) {
              Da[ii][is][js] = 0.0;
              for ( ip = 0; ip < npha; ip++ ) {
                for ( ks = 0; ks < nsol; ks++ ) {
                  //Da[ii][is][js] += pfmdat->D[ip][is*nsol+ks] * dc_dmu[ii][ip][ks][js] * hfhi(stgO[ii].phi, ip);
                  Da[ii][is][js] += pfmdat->D[ip][is*nsol+ks] * dc_dmu[ii][ip][ks][js] * stgO[ii].phi[ip];
                }
              }
            }
          }

        }
        else {

          // Get dc_dmu for bulk phase
          for ( is1 = 0; is1 < nsol; is1++ ) {
            for ( is2 = 0; is2 < nsol; is2++ ) {
              dc_dmu[ii][bulkphase][is1][is2] = propf3->cmu[bulkphase][is1][is2];
            }
          }

          // Calculate mobility for bulk phase at stencil points
          for ( is = 0; is < nsol; is++ ) {
            for ( js =0; js < nsol; js++ ) {
              Da[ii][is][js] = 0.0;

                for ( ks = 0; ks < nsol; ks++ ) {
                  //Da[ii][is][js] += pfmdat->D[bulkphase][is*nsol+ks] * dc_dmu[ii][bulkphase][ks][js];
                  Da[ii][is][js] += pfmdat->D[bulkphase][is*nsol+ks] * dc_dmu[ii][bulkphase][ks][js];
                }

            }
          }

        }

      }

      // Calculating mobility at half grid positions
      for ( is = 0; is < nsol; is++ ) {
        for ( js =0; js < nsol; js++ ) {

          Damidy[0][is][js] = ( Da[0][is][js] + Da[1][is][js] ) / 2.0;
          Damidy[1][is][js] = ( Da[2][is][js] + Da[0][is][js] ) / 2.0;
          
          Damidz[0][is][js] = ( Da[0][is][js] + Da[3][is][js] ) / 2.0;
          Damidz[1][is][js] = ( Da[4][is][js] + Da[0][is][js] ) / 2.0;
          
          Damidx[0][is][js] = ( Da[0][is][js] + Da[5][is][js] ) / 2.0;
          Damidx[1][is][js] = ( Da[6][is][js] + Da[0][is][js] ) / 2.0;
          
          //Damidx[0][is][js] = ( Da[2][is][js] + Da[4][is][js] ) / 2.0;
          //Damidx[1][is][js] = ( Da[4][is][js] + Da[1][is][js] ) / 2.0;

          //Damidy[0][is][js] = ( Da[3][is][js] + Da[4][is][js] ) / 2.0;
          //Damidy[1][is][js] = ( Da[4][is][js] + Da[0][is][js] ) / 2.0;
        }
      }

      // Calculating grad of mu
      for ( is = 0; is < nsol; is++ ) {
        gradmuy[0][is] = ( stgO[0].mu[is] - stgO[1].mu[is] ) / pfmvar->deltay;
        gradmuy[1][is] = ( stgO[2].mu[is] - stgO[0].mu[is] ) / pfmvar->deltay;

        gradmuz[0][is] = ( stgO[0].mu[is] - stgO[3].mu[is] ) / pfmvar->deltaz;
        gradmuz[1][is] = ( stgO[4].mu[is] - stgO[0].mu[is] ) / pfmvar->deltaz;

        gradmux[0][is] = ( stgO[0].mu[is] - stgO[5].mu[is] ) / pfmvar->deltax;
        gradmux[1][is] = ( stgO[6].mu[is] - stgO[0].mu[is] ) / pfmvar->deltax;
        
        //gradmux[0][is] = ( stgO[2].mu[is] - stgO[4].mu[is] ) / pfmvar->deltax;
        //gradmux[1][is] = ( stgO[4].mu[is] - stgO[1].mu[is] ) / pfmvar->deltax;
        //gradmuy[0][is] = ( stgO[3].mu[is] - stgO[4].mu[is] ) / pfmvar->deltay;
        //gradmuy[1][is] = ( stgO[4].mu[is] - stgO[0].mu[is] ) / pfmvar->deltay;
      }

      // Calculating divergence of flux
      for ( is = 0; is < nsol; is++ ) {
        divflux[is] = 0.0;
        for ( js =0; js < nsol; js++ ) {

          divflux[is] =  ( Damidy[1][is][js] * gradmuy[1][js] - Damidy[0][is][js] * gradmuy[0][js] ) / pfmvar->deltay;

          divflux[is] += ( Damidz[1][is][js] * gradmuz[1][js] - Damidz[0][is][js] * gradmuz[0][js] ) / pfmvar->deltaz;

          divflux[is] += ( Damidx[1][is][js] * gradmux[1][js] - Damidx[0][is][js] * gradmux[0][js] ) / pfmvar->deltax;

          //divflux[is] += ( Damidx[0][is][js] * gradmux[0][js] - Damidx[1][is][js] * gradmux[1][js] ) / pfmvar->deltax;

          //divflux[is] += ( Damidy[0][is][js] * gradmuy[0][js] - Damidy[1][is][js] * gradmuy[1][js] ) / pfmvar->deltay;

        }
      }


      // Claculation are done at x,y,z for updating composition and mu

      // Check for interface
      interface = 1;
      bulkphase = 0;
      for ( ip = 0; ip < npha; ip++ ) {
        if ( stgO[0].phi[ip] >= pfmdat->interfaceUplimit ) {
          bulkphase = ip;
          interface = 0;
          break;
        }
      }



          for ( ip = 0; ip < npha; ip++ ) {

            for ( is1 = 0; is1 < nsol; is1++ ) {

              F3_Beq[ip*nsol+is1] = propf3->Beq[ip][is1];

              F3_dBbdT[ip*nsol+is1] = propf3->dBbdT[ip][is1];

              for ( is2 = 0; is2 < nsol; is2++ ) {
                F3_A[(ip*nsol+is1)*nsol+is2] = propf3->A[ip][is1][is2];
              }

            }
          }
          Tequ = pfmdat->Teq;

      // Update composition and mu in bulk phase
      if ( !interface ) {

        // update composition
        for ( is = 0; is < nsol; is++ ) {
          gridinfo[index].com[is] = gridinfoO[index].com[is] + pfmvar->deltat * ( divflux[is] + jatr[is] );
          comp[is] = gridinfo[index].com[is];
        }


          function_F_03_Mu_CL(comp, Ti, bulkphase, mu, F3_A, F3_Beq, F3_dBbdT, Tequ);

        // Updating mu
        for ( is = 0; is < nsol; is++ ) {
          gridinfo[index].mu[is] = mu[is];
        }

      }
      else {

        // Calculation are carried out in iterface region

        // Calculate Sum ( c_alpha * dh_alpha/ dt )
        for ( ip1 = 0; ip1 < npha; ip1++ ) {
          sum_dhphi = 0.0;
          for ( ip2 = 0; ip2 < npha; ip2++ ) {
            sum_dhphi += dhfhi(stgO[0].phi, ip1, ip2) * ( gridinfo[index].phi[ip2] - gridinfoO[index].phi[ip2] );
          }

          for ( is = 0; is < nsol; is++ ) {
            sum[is]       += stcscl[0].comie[ip1][is] * sum_dhphi;
          }
        }

        for ( is1 = 0; is1 < nsol; is1++ ) {

          // Calculate deltac
          deltac[is1] = pfmvar->deltat * ( divflux[is1] + jatr[is1] );

          //if (tstep[0] > 10) {
          // Update composition using deltac
          gridinfo[index].com[is1] = stgO[0].com[is1] + deltac[is1];

          // Update deltac with c_alpha * dh_alpha/dt
          deltac[is1] += -sum[is1];
          //}


          // Calculating h_alpha * dc_dmu
          for ( is2 = 0; is2 < nsol; is2++ ) {
            for ( ip = 0; ip < npha; ip++ ) {
              dcdmu[is1*nsol+is2] += dc_dmu[0][ip][is1][is2] * hfhi(stgO[0].phi, ip);
            }
          }

        }

        // printf("t=%d, %d=i, %d=j, sum=%le\n", tstep[0], i, j, sum[0]);

        //printf("t=%d, %d=i, %d=j, phi0=%le, phi1=%le, h0=%le, h1=%le, dh00=%le, dh01=%le, dh10=%le, dh11=%le\n", tstep[0], i, j, stgO[4].phi[0], stgO[4].phi[1], hfhi(stgO[4].phi, 0), hfhi(stgO[4].phi, 1), dhfhi(stgO[4].phi, 0, 0), dhfhi(stgO[4].phi, 0, 1), dhfhi(stgO[4].phi, 1, 0), dhfhi(stgO[4].phi, 1, 1));

        // Inverse of dcdmu
        matinvnew_nsol(dcdmu, inv_dcdmu);

        multiply_nsol(inv_dcdmu, deltac, deltamu);
        // for ( is1 = 0; is1 < nsol; is1++ ) {
        //   tmp0 = 0.0;
        //   for ( is2 = 0; is2 < nsol; is2++ ) {
        //     tmp0 += inv_dcdmu[is1*nsol+is2] * deltac[is2];
        //   }
        //   deltamu[is1] = tmp0;
        // }

        //vectorsum(deltamu, stgO[index].mu, gridinfo[index].mu, nsol);
        for ( is = 0; is < nsol; is++ ) {
          gridinfo[index].mu[is] = gridinfoO[index].mu[is] + deltamu[is];
        }
        //printf("IC %d, %d, %d, %le, %le, %le, %le\n", tstep[0], i, j, gridinfoO[index].mu[0], gridinfo[index].mu[0], gridinfoO[index].com[0], gridinfo[index].com[0]);
      }

    }
}

__kernel void SolverCatr_F4_smooth(__global struct fields *gridinfoO, __global struct fields *gridinfo, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf4 *propf4, __constant struct propmatf4spline *propf4spline, __constant struct propmatf4spline *propf4spline1) {

  int x, y, z;
  int jj, ii;
  int nx, ny, nz;
  int index;

  int ig, is, js, ip, i1, j1, is1, is2, ip1, ip2, ks;
  int interface, bulkphase, indij[5];

  double A1[npha-1];
  double dphidt[npha-1][7], gradx_phi[npha][5], grady_phi[npha][5], gradz_phi[npha][5], phix[npha][7], phiy[npha][7], phiz[npha][7], modgradphi[npha][7], scalprodct[npha][7];
  double cjatx, cjaty, cjatz;
  double jatc[npha-1][nsol], jatr[nsol];

  //double hphid[5][npha], hphi[npha], dhphi[npha][npha], hphicii[5][npha];
  double tmp0, tmp1;
  double DELTAT, sum_dhphi, Ti;
  double sum[nsol], sum_dcbdT[nsol], dcdmu[nsol*nsol], dc_dmu[7][npha][nsol*nsol], Da[7][nsol][nsol];
  double Damidx[2][nsol][nsol], Damidy[2][nsol][nsol], Damidz[2][nsol][nsol], gradmux[2][nsol], gradmuy[2][nsol], gradmuz[2][nsol];
  double divflux[nsol], mu[nsol], deltamu[nsol], c_tdt[nsol], dcbdT_phase[npha][nsol];
  double deltac[nsol], inv_dcdmu[nsol*nsol], cg[nsol];

  double alpha[npha-1][nsol][7], alphidot[npha-1][nsol][6], jat[npha-1][nsol][6], muc[npha][nsol*nsol], muc1[nsol*nsol], cmu1[nsol*nsol];

  struct fields stgridO[27];
  struct fields stgN[7];
  struct csle stcscl[7];
  struct fields stgO[7];

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  //index = y + ny*(z + x*nz);

  if ( x != 0 && x != nx-1 && y != 0 && y != ny-1 && z != 0 && z != nz-1 ) {

    index = y + ny*(z + x*nz);

    Ti = pfmdat->T0;

    stgridO[0]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[1]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[2]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[3]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[4]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[5]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[6]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[7]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[8]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[9]  = gridinfoO[ (y-1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[10] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[11] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[12] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[13] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[14] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[15] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[16] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[17] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[18] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[19] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[20] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[21] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[22] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[23] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[24] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[25] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[26] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z+1) ) ];

    stcscl[0] = cscl[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stcscl[1] = cscl[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stcscl[2] = cscl[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stcscl[3] = cscl[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stcscl[4] = cscl[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stcscl[5] = cscl[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stcscl[6] = cscl[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];

    stgO[0] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stgO[1] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stgO[2] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stgO[3] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stgO[4] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stgO[5] = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stgO[6] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];

    for ( is1 = 0; is1 < nsol; is1++ ) {
        sum[is1] = 0.0;
        sum_dcbdT[is1] = 0.0;
        for ( is2 = 0; is2 < nsol; is2++ ) {
          dcdmu[is1*nsol+is2] = 0.0;
        }
      }

      for ( ip = 0; ip < npha; ip++ ) {
        for ( is = 0; is < nsol; is++ ) {
          dcbdT_phase[ip][is] = 0;
        }
      }

      for ( ii = 0; ii < 7; ii++ ) {

        interface = 1;
        bulkphase = 0;
        for ( ip = 0; ip < npha; ip++ ) {
          if ( stgO[ii].phi[ip] >= pfmdat->interfaceUplimit ) {
            bulkphase = ip;
            interface = 0;
            break;
          }
        }

        if ( interface ) {

          //stcscl[2].comie[ip][is]

          //if ( pfmdat->ISOTHERMAL ) {
            for ( ip = 0; ip < npha; ip++ ) {
              for ( is1 = 0; is1 < nsol; is1++ ) {
                for ( is2 =0; is2 < nsol; is2++ ) {
                  dc_dmu[ii][ip][is1*nsol+is2] = propf4->cmu[ip][is1][is2];
                }
              }
            }
          //}
          /*else {
            for ( ip = 0; ip < npha; ip++ ) {
              for ( is1 = 0; is1 < nsol; is1++ ) {
                for ( is2 =0; is2 < nsol; is2++ ) {
                  if ( is1 == is2 ) {
                    muc[ip][is1*nsol+is2] = 2.0 * propf4spline[indij[ii]].A[ip][is1][is2];
                  }
                  else {
                    muc[ip][is1*nsol+is2] = propf4spline[indij[ii]].A[ip][is1][is2];
                  }
                }
              }
              matinvnew_nsol(muc[ip], dc_dmu[ii][ip]);
            }
          }*/

          //printf("%d, %d, %d\n", tstep[0], i, j);
          for ( is = 0; is < nsol; is++ ) {
            for ( js =0; js < nsol; js++ ) {
              Da[ii][is][js] = 0.0;
              for ( ip = 0; ip < npha; ip++ ) {
                for ( ks = 0; ks < nsol; ks++ ) {
                  Da[ii][is][js] += pfmdat->D[ip][is*nsol+ks] * dc_dmu[ii][ip][ks*nsol+js] * stgO[ii].phi[ip];
                  //Da[ii][is][js] += pfmdat->D[ip][is*nsol+ks] * dc_dmu[ii][ip][ks+nsol*js] * stgO[ii].phi[ip];
                }
              }
            }
          }

        }
        else {

          // bulk compositon

          //if ( pfmdat->ISOTHERMAL ) {
            for ( is1 = 0; is1 < nsol; is1++ ) {
              for ( is2 = 0; is2 < nsol; is2++ ) {
                dc_dmu[ii][bulkphase][is1*nsol+is2] = propf4->cmu[bulkphase][is1][is2];
              }
            }
          //}
          /*else {
            for ( is1 = 0; is1 < nsol; is1++ ) {
              for ( is2 =0; is2 < nsol; is2++ ) {
                if ( is1 == is2 ) {
                  muc[bulkphase][is1*nsol+is2] = 2.0 * propf4spline[indij[ii]].A[bulkphase][is1][is2];
                }
                else {
                  muc[bulkphase][is1*nsol+is2] = propf4spline[indij[ii]].A[bulkphase][is1][is2];
                }
              }
              matinvnew_nsol(muc[bulkphase], dc_dmu[ii][bulkphase]);
            }
          }*/


          for ( is = 0; is < nsol; is++ ) {
            for ( js =0; js < nsol; js++ ) {
              Da[ii][is][js] = 0.0;
              //for ( ip = 0; ip < npha; ip++ ) {
                for ( ks = 0; ks < nsol; ks++ ) {
                  Da[ii][is][js] += pfmdat->D[bulkphase][is*nsol+ks] * dc_dmu[ii][bulkphase][ks*nsol+js];
                  //Da[ii][is][js] += pfmdat->D[bulkphase][is*nsol+ks] * dc_dmu[ii][bulkphase][ks+nsol*js];
                }
              //}
            }
          }

        }

      }

      for ( is = 0; is < nsol; is++ ) {
        for ( js =0; js < nsol; js++ ) {

          Damidy[0][is][js] = ( Da[0][is][js] + Da[1][is][js] ) / 2.0;
          Damidy[1][is][js] = ( Da[2][is][js] + Da[0][is][js] ) / 2.0;

          Damidz[0][is][js] = ( Da[0][is][js] + Da[3][is][js] ) / 2.0;
          Damidz[1][is][js] = ( Da[4][is][js] + Da[0][is][js] ) / 2.0;

          Damidx[0][is][js] = ( Da[0][is][js] + Da[5][is][js] ) / 2.0;
          Damidx[1][is][js] = ( Da[6][is][js] + Da[0][is][js] ) / 2.0;

        }
      }

      for ( is = 0; is < nsol; is++ ) {
        gradmuy[0][is] = ( stgO[0].mu[is] - stgO[1].mu[is] ) / pfmvar->deltay;
        gradmuy[1][is] = ( stgO[2].mu[is] - stgO[0].mu[is] ) / pfmvar->deltay;

        gradmuz[0][is] = ( stgO[0].mu[is] - stgO[3].mu[is] ) / pfmvar->deltaz;
        gradmuz[1][is] = ( stgO[4].mu[is] - stgO[0].mu[is] ) / pfmvar->deltaz;

        gradmux[0][is] = ( stgO[0].mu[is] - stgO[5].mu[is] ) / pfmvar->deltax;
        gradmux[1][is] = ( stgO[6].mu[is] - stgO[0].mu[is] ) / pfmvar->deltax;
      }

      for ( is = 0; is < nsol; is++ ) {
        divflux[is] = 0.0;
        for ( js =0; js < nsol; js++ ) {

          divflux[is] +=  ( Damidy[1][is][js] * gradmuy[1][js] - Damidy[0][is][js] * gradmuy[0][js] ) / pfmvar->deltay;

          divflux[is] += ( Damidz[1][is][js] * gradmuz[1][js] - Damidz[0][is][js] * gradmuz[0][js] ) / pfmvar->deltaz;

          divflux[is] += ( Damidx[1][is][js] * gradmux[1][js] - Damidx[0][is][js] * gradmux[0][js] ) / pfmvar->deltax;

        }
      }

      interface = 1;
      bulkphase = 0;
      for ( ip = 0; ip < npha; ip++ ) {
        if ( stgO[0].phi[ip] >= pfmdat->interfaceUplimit ) {
          bulkphase = ip;
          interface = 0;
          break;
        }
      }

      if ( !interface ) {

        for ( is = 0; is < nsol; is++ ) {
          gridinfo[index].com[is] = gridinfoO[index].com[is] + pfmvar->deltat * ( divflux[is] );
        }

        ip = bulkphase;
        //if ( pfmdat->ISOTHERMAL ) {
          tmp0 = 0.0;
          for ( is1 = 0; is1 < nsol; is1++ ) {
            tmp0  = 2.0*propf4->A[bulkphase][is1][is1]*gridinfo[index].com[is1] + propf4->B[bulkphase][is1];
            for ( is2 = 0; is2 < nsol; is2++ ) {
              if ( is1 != is2 ) {
                tmp0 += propf4->A[bulkphase][is1][is2]*gridinfo[index].com[is2];
              }
            }
            mu[is1] = tmp0;
          }
        //}
        /*else {
          tmp0 = 0.0;
          for ( is1 = 0; is1 < nsol; is1++ ) {
            tmp0  = 2.0*propf4spline[j].A[bulkphase][is1][is1]*gridinfo[index].com[is1] + propf4spline[j].B[bulkphase][is1];
            for ( is2 = 0; is2 < nsol; is2++ ) {
              if ( is1 != is2 ) {
                tmp0 += propf4spline[j].A[bulkphase][is1][is2]*gridinfo[index].com[is2];
              }
            }
            mu[is1] = tmp0;
          }
        }*/

        for ( is = 0; is < nsol; is++ ) {
          //deltamu[is] = mu[is] - stgO[0].mu[is];
          gridinfo[index].mu[is] = mu[is];
        }
        //printf("BC %d, %d, %d, %le, %le, %le, %le\n", tstep[0], i, j, gridOld[index].mu[0], gridinfo[index].mu[0], gridOld[index].com[0], gridinfo[index].com[0]);
      }
      else {

        /*if ( pfmdat->ISOTHERMAL ) {
          for ( ip = 0; ip < npha; ip++ ) {
            for ( is = 0; is < nsol; is++ ) {
              dcbdT_phase[ip][is] = 0.0;
            }
          }
        }*/
        /*else {

          DELTAT = ( pfmvar->deltat ) * ( -pfmdat->TGRADIENT * pfmdat->velocity );

          for ( ip = 0; ip < npha; ip++ ) {
            for ( is1 = 0; is1 < nsol; is1++ ) {
              for ( is2 = 0; is2 < nsol; is2++ ) {
                if ( is1 == is2 ) {
                  muc1[is1*nsol+is2] = 2.0 * propf4spline1[j].A[ip][is1][is2];
                }
                else {
                  muc1[is1*nsol+is2] = propf4spline1[j].A[ip][is1][is2];
                }
              }
            }
            matinvnew_nsol(muc1, cmu1);

            for ( is1 = 0; is1 < nsol; is1++ ) {
              tmp0 = 0.0;
              for ( is2 = 0; is2 < nsol; is2++ ) {
                tmp0 += cmu1[is1*nsol+is2] * ( gridOld[index].mu[is2] - propf4spline1[j].B[ip][is2] );
              }
              c_tdt[is1] = tmp0;
            }

            for ( is = 0; is < nsol; is++ ) {
              dcbdT_phase[ip][is] = c_tdt[is] - stcscl[4].comie[ip][is];
            }

          }
        }*/


        for ( is1 = 0; is1 < nsol; is1++ ) {

          deltac[is1] = pfmvar->deltat * ( divflux[is1] + jatr[is1] );

          // gridinfo[index].com[is1] = stgO[0].com[is1] + deltac[is1];
          //
          // //if ( pfmdat->ISOTHERMAL ) {
          //   deltac[is1] += -sum[is1];
          // //}
          // //else {
          // //  deltac[is1] += -sum[is1] - sum_dcbdT[is1];
          // //}
          for ( is2 = 0; is2 < nsol; is2++ ) {
            for ( ip = 0; ip < npha; ip++ ) {
              dcdmu[is1*nsol+is2] += dc_dmu[0][ip][is1*nsol+is2] * hfhi(stgO[0].phi, ip);
            }
          }
        }

        matinvnew_nsol(dcdmu, inv_dcdmu);

        multiply_nsol(inv_dcdmu, deltac, deltamu);
        // for ( is1 = 0; is1 < nsol; is1++ ) {
        //   tmp0 = 0.0;
        //   for ( is2 = 0; is2 < nsol; is2++ ) {
        //     tmp0 += inv_dcdmu[is1*nsol+is2] * deltac[is2];
        //   }
        //   deltamu[is1] = tmp0;
        // }

        //vectorsum(deltamu, stgO[index].mu, gridinfo[index].mu, nsol);
        for ( is = 0; is < nsol; is++ ) {
          gridinfo[index].mu[is] = gridinfoO[index].mu[is] + deltamu[is];
        }
        //printf("IC %d, %d, %d, %le, %le\n", tstep[0], i, j, gridinfo[index].mu[0], gridinfo[index].com[0]);
        //printf("IC %d, %d, %d, %le, %le, %le, %le\n", tstep[0], i, j, gridOld[index].mu[0], gridinfo[index].mu[0], gridOld[index].com[0], gridinfo[index].com[0]);
      }
      //printf("%d, %d, %d, %le, %le\n", tstep[0], i, j, gridinfo[index].mu[0], gridinfo[index].com[0]);

    }


}

__kernel void SolverCatr_F4(__global struct fields *gridinfoO, __global struct fields *gridinfo, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp, __global int *tstep, __constant struct propmatf4 *propf4, __constant struct propmatf4spline *propf4spline, __constant struct propmatf4spline *propf4spline1) { 

  int x, y, z;
  int jj, ii;
  int nx, ny, nz;
  int index;

  int ig, is, js, ip, i1, j1, is1, is2, ip1, ip2, ks;
  int interface, bulkphase, indij[5];

  double A1[npha-1];
  double dphidt[npha-1][7], gradx_phi[npha][5], grady_phi[npha][5], gradz_phi[npha][5], phix[npha][7], phiy[npha][7], phiz[npha][7], modgradphi[npha][7], scalprodct[npha][7];
  double cjatx, cjaty, cjatz;
  double jatc[npha-1][nsol], jatr[nsol];

  //double hphid[5][npha], hphi[npha], dhphi[npha][npha], hphicii[5][npha];
  double tmp0, tmp1;
  double DELTAT, sum_dhphi, Ti;
  double sum[nsol], sum_dcbdT[nsol], dcdmu[nsol*nsol], dc_dmu[7][npha][nsol*nsol], Da[7][nsol][nsol];
  double Damidx[2][nsol][nsol], Damidy[2][nsol][nsol], Damidz[2][nsol][nsol], gradmux[2][nsol], gradmuy[2][nsol], gradmuz[2][nsol];
  double divflux[nsol], mu[nsol], deltamu[nsol], c_tdt[nsol], dcbdT_phase[npha][nsol];
  double deltac[nsol], inv_dcdmu[nsol*nsol], cg[nsol];

  double alpha[npha-1][nsol][7], alphidot[npha-1][nsol][6], jat[npha-1][nsol][6], muc[npha][nsol*nsol], muc1[nsol*nsol], cmu1[nsol*nsol];

  struct fields stgridO[27];
  struct fields stgN[7];
  struct csle stcscl[7];
  struct fields stgO[7];

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  //index = y + ny*(z + x*nz);

  if ( x != 0 && x != nx-1 && y != 0 && y != ny-1 && z != 0 && z != nz-1 ) {

    index = y + ny*(z + x*nz);

    Ti = pfmdat->T0;

    stgridO[0]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[1]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[2]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z-1) ) ];
    stgridO[3]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[4]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[5]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z  ) ) ];
    stgridO[6]  = gridinfoO[ (y-1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[7]  = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[8]  = gridinfoO[ (y+1) + ny*( (x-1)*nz + (z+1) ) ];
    stgridO[9]  = gridinfoO[ (y-1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[10] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[11] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z-1) ) ];
    stgridO[12] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[13] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[14] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stgridO[15] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[16] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[17] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z+1) ) ];
    stgridO[18] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[19] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[20] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z-1) ) ];
    stgridO[21] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[22] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[23] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z  ) ) ];
    stgridO[24] = gridinfoO[ (y-1) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[25] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z+1) ) ];
    stgridO[26] = gridinfoO[ (y+1) + ny*( (x+1)*nz + (z+1) ) ];

    stcscl[0] = cscl[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stcscl[1] = cscl[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stcscl[2] = cscl[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stcscl[3] = cscl[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stcscl[4] = cscl[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stcscl[5] = cscl[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stcscl[6] = cscl[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];

    stgO[0] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stgO[1] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stgO[2] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stgO[3] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stgO[4] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stgO[5] = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stgO[6] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];

    for ( ip = 0; ip < npha-1; ip++ ) {
        stgN[0].phi[ip] = gridinfo[ (y  ) + ny*( (x  )*nz + (z  ) ) ].phi[ip];
        stgN[1].phi[ip] = gridinfo[ (y-1) + ny*( (x  )*nz + (z  ) ) ].phi[ip];
        stgN[2].phi[ip] = gridinfo[ (y+1) + ny*( (x  )*nz + (z  ) ) ].phi[ip];
        stgN[3].phi[ip] = gridinfo[ (y  ) + ny*( (x  )*nz + (z-1) ) ].phi[ip];
        stgN[4].phi[ip] = gridinfo[ (y  ) + ny*( (x  )*nz + (z+1) ) ].phi[ip];
        stgN[5].phi[ip] = gridinfo[ (y  ) + ny*( (x-1)*nz + (z  ) ) ].phi[ip];
        stgN[6].phi[ip] = gridinfo[ (y  ) + ny*( (x+1)*nz + (z  ) ) ].phi[ip];
    }

    for ( ip = 0; ip < npha-1; ip++ ) {
        dphidt[ip][0] = ( stgN[0].phi[ip] - stgO[0].phi[ip] ) / ( pfmvar->deltat );
        dphidt[ip][1] = ( stgN[1].phi[ip] - stgO[1].phi[ip] ) / ( pfmvar->deltat );
        dphidt[ip][2] = ( stgN[2].phi[ip] - stgO[2].phi[ip] ) / ( pfmvar->deltat );
        dphidt[ip][3] = ( stgN[3].phi[ip] - stgO[3].phi[ip] ) / ( pfmvar->deltat );
        dphidt[ip][4] = ( stgN[4].phi[ip] - stgO[4].phi[ip] ) / ( pfmvar->deltat );
        dphidt[ip][5] = ( stgN[5].phi[ip] - stgO[5].phi[ip] ) / ( pfmvar->deltat );
        dphidt[ip][6] = ( stgN[6].phi[ip] - stgO[6].phi[ip] ) / ( pfmvar->deltat );

      A1[ip] = sqrt ( pfmvar->ee[ip*npha+(npha-1)] ) / sqrt( pfmvar->w[ip*npha+(npha-1)] );

    }


    for ( ip = 0; ip < npha-1; ip++ ) {
      for ( is = 0; is < nsol; is++ ) {
          alpha[ip][is][0] = (A1[ip])*( stcscl[0].comie[npha-1][is] - stcscl[0].comie[ip][is] );
          alpha[ip][is][1] = (A1[ip])*( stcscl[1].comie[npha-1][is] - stcscl[1].comie[ip][is] );
          alpha[ip][is][2] = (A1[ip])*( stcscl[2].comie[npha-1][is] - stcscl[2].comie[ip][is] );
          alpha[ip][is][3] = (A1[ip])*( stcscl[3].comie[npha-1][is] - stcscl[3].comie[ip][is] );
          alpha[ip][is][4] = (A1[ip])*( stcscl[4].comie[npha-1][is] - stcscl[4].comie[ip][is] );
          alpha[ip][is][5] = (A1[ip])*( stcscl[5].comie[npha-1][is] - stcscl[5].comie[ip][is] );
          alpha[ip][is][6] = (A1[ip])*( stcscl[6].comie[npha-1][is] - stcscl[6].comie[ip][is] );
      }
    }

    for ( ip = 0; ip < npha; ip++ ) {
        grady_phi[ip][0] = ( stgridO[14].phi[ip] - stgridO[12].phi[ip] ) / ( 2.0*pfmvar->deltay );
        grady_phi[ip][1] = ( stgridO[11].phi[ip] - stgridO[9 ].phi[ip] ) / ( 2.0*pfmvar->deltay );
        grady_phi[ip][2] = ( stgridO[17].phi[ip] - stgridO[15].phi[ip] ) / ( 2.0*pfmvar->deltay );
        grady_phi[ip][3] = ( stgridO[5 ].phi[ip] - stgridO[3 ].phi[ip] ) / ( 2.0*pfmvar->deltay );
        grady_phi[ip][4] = ( stgridO[23].phi[ip] - stgridO[21].phi[ip] ) / ( 2.0*pfmvar->deltay );

        gradz_phi[ip][0] = ( stgridO[16].phi[ip] - stgridO[10].phi[ip] ) / ( 2.0*pfmvar->deltaz );
        gradz_phi[ip][1] = ( stgridO[15].phi[ip] - stgridO[9 ].phi[ip] ) / ( 2.0*pfmvar->deltaz );
        gradz_phi[ip][2] = ( stgridO[17].phi[ip] - stgridO[11].phi[ip] ) / ( 2.0*pfmvar->deltaz );
        gradz_phi[ip][3] = ( stgridO[7 ].phi[ip] - stgridO[1 ].phi[ip] ) / ( 2.0*pfmvar->deltaz );
        gradz_phi[ip][4] = ( stgridO[25].phi[ip] - stgridO[19].phi[ip] ) / ( 2.0*pfmvar->deltaz );

        gradx_phi[ip][0] = ( stgridO[22].phi[ip] - stgridO[4 ].phi[ip] ) / ( 2.0*pfmvar->deltax );
        gradx_phi[ip][1] = ( stgridO[21].phi[ip] - stgridO[3 ].phi[ip] ) / ( 2.0*pfmvar->deltax );
        gradx_phi[ip][2] = ( stgridO[23].phi[ip] - stgridO[5 ].phi[ip] ) / ( 2.0*pfmvar->deltax );
        gradx_phi[ip][3] = ( stgridO[19].phi[ip] - stgridO[1 ].phi[ip] ) / ( 2.0*pfmvar->deltax );
        gradx_phi[ip][4] = ( stgridO[25].phi[ip] - stgridO[7 ].phi[ip] ) / ( 2.0*pfmvar->deltax );
    }

    for ( ip = 0; ip < npha; ip++ ) {

        phiy[ip][0] = grady_phi[ip][0];
        phiy[ip][1] = ( stgridO[13].phi[ip] - stgridO[12].phi[ip] ) / ( pfmvar->deltay );
        phiy[ip][2] = ( stgridO[14].phi[ip] - stgridO[13].phi[ip] ) / ( pfmvar->deltay );
        phiy[ip][3] = ( grady_phi[ip][0] + grady_phi[ip][1] ) / ( 2.0 );
        phiy[ip][4] = ( grady_phi[ip][0] + grady_phi[ip][2] ) / ( 2.0 );
        phiy[ip][5] = ( grady_phi[ip][0] + grady_phi[ip][3] ) / ( 2.0 );
        phiy[ip][6] = ( grady_phi[ip][0] + grady_phi[ip][4] ) / ( 2.0 );

        phiz[ip][0] = gradz_phi[ip][0];
        phiz[ip][1] = ( gradz_phi[ip][0] + gradz_phi[ip][1] ) / ( 2.0 );
        phiz[ip][2] = ( gradz_phi[ip][0] + gradz_phi[ip][2] ) / ( 2.0 );
        phiz[ip][3] = ( stgridO[13].phi[ip] - stgridO[10].phi[ip] ) / ( pfmvar->deltaz );
        phiz[ip][4] = ( stgridO[16].phi[ip] - stgridO[13].phi[ip] ) / ( pfmvar->deltaz );
        phiz[ip][5] = ( gradz_phi[ip][0] + gradz_phi[ip][3] ) / ( 2.0 );
        phiz[ip][6] = ( gradz_phi[ip][0] + gradz_phi[ip][4] ) / ( 2.0 );

        phix[ip][0] = gradx_phi[ip][0];
        phix[ip][1] = ( gradx_phi[ip][0] + gradx_phi[ip][1] ) / ( 2.0 );
        phix[ip][2] = ( gradx_phi[ip][0] + gradx_phi[ip][2] ) / ( 2.0 );
        phix[ip][3] = ( gradx_phi[ip][0] + gradx_phi[ip][3] ) / ( 2.0 );
        phix[ip][4] = ( gradx_phi[ip][0] + gradx_phi[ip][4] ) / ( 2.0 );
        phix[ip][5] = ( stgridO[13].phi[ip] - stgridO[4 ].phi[ip] ) / ( pfmvar->deltax );
        phix[ip][6] = ( stgridO[22].phi[ip] - stgridO[13].phi[ip] ) / ( pfmvar->deltax );
    }

    for ( ip = 0; ip < npha-1; ip++ ) {
	      for ( is = 0; is < nsol; is++ ) {
          alphidot[ip][is][0] = ( ( alpha[ip][is][0] * dphidt[ip][0] ) + ( alpha[ip][is][1] * dphidt[ip][1] ) ) / ( 2.0 );
          alphidot[ip][is][1] = ( ( alpha[ip][is][0] * dphidt[ip][0] ) + ( alpha[ip][is][2] * dphidt[ip][2] ) ) / ( 2.0 );
          alphidot[ip][is][2] = ( ( alpha[ip][is][0] * dphidt[ip][0] ) + ( alpha[ip][is][3] * dphidt[ip][3] ) ) / ( 2.0 );
          alphidot[ip][is][3] = ( ( alpha[ip][is][0] * dphidt[ip][0] ) + ( alpha[ip][is][4] * dphidt[ip][4] ) ) / ( 2.0 );
          alphidot[ip][is][4] = ( ( alpha[ip][is][0] * dphidt[ip][0] ) + ( alpha[ip][is][5] * dphidt[ip][5] ) ) / ( 2.0 );
          alphidot[ip][is][5] = ( ( alpha[ip][is][0] * dphidt[ip][0] ) + ( alpha[ip][is][6] * dphidt[ip][6] ) ) / ( 2.0 );
	      }
    }

    for ( ip = 0; ip < npha; ip++ ) {
      for ( jj = 0; jj < 7; jj++ ) {
          modgradphi[ip][jj] = sqrt( phiy[ip][jj]*phiy[ip][jj] + phiz[ip][jj]*phiz[ip][jj] + phix[ip][jj]*phix[ip][jj]);
      }
    }

    for ( ip = 0; ip < npha-1; ip++ ) {
      for ( jj = 0; jj < 7; jj++ ) {
          scalprodct[ip][jj] = -1.0*( phiy[ip][jj]*phiy[npha-1][jj] + phiz[ip][jj]*phiz[npha-1][jj] + phix[ip][jj]*phix[npha-1][jj] );
          if ( modgradphi[npha-1][jj] > 0.0 ) {
            scalprodct[ip][jj] /= ( modgradphi[ip][jj] * modgradphi[npha-1][jj] );
        }
      }
    }



    for ( ip = 0; ip < npha-1; ip++ ) {
	    for ( is = 0; is < nsol; is++ ) {
          jat[ip][is][0] = ( ( alphidot[ip][is][0] * phiy[ip][1] ) / ( modgradphi[ip][1] ) )   * ( ( 1.0 - pfmdat->D[ip][is*nsol+is] /   pfmdat->D[npha-1][is*nsol+is] ) * fabs(scalprodct[ip][1]) );
          jat[ip][is][1] = ( ( alphidot[ip][is][1] * phiy[ip][2] ) / ( modgradphi[ip][2] ) )   * ( ( 1.0 - pfmdat->D[ip][is*nsol+is] /   pfmdat->D[npha-1][is*nsol+is] ) * fabs(scalprodct[ip][2]) );
          jat[ip][is][2] = ( ( alphidot[ip][is][2] * phiz[ip][3] ) / ( modgradphi[ip][3] ) )   * ( ( 1.0 - pfmdat->D[ip][is*nsol+is] /   pfmdat->D[npha-1][is*nsol+is] ) * fabs(scalprodct[ip][3]) );
          jat[ip][is][3] = ( ( alphidot[ip][is][3] * phiz[ip][4] ) / ( modgradphi[ip][4] ) )   * ( ( 1.0 - pfmdat->D[ip][is*nsol+is] /   pfmdat->D[npha-1][is*nsol+is] ) * fabs(scalprodct[ip][4]) );
          jat[ip][is][4] = ( ( alphidot[ip][is][4] * phix[ip][5] ) / ( modgradphi[ip][5] ) )   * ( ( 1.0 - pfmdat->D[ip][is*nsol+is] /   pfmdat->D[npha-1][is*nsol+is] ) * fabs(scalprodct[ip][5]) );
          jat[ip][is][5] = ( ( alphidot[ip][is][5] * phix[ip][6] ) / ( modgradphi[ip][6] ) )   * ( ( 1.0 - pfmdat->D[ip][is*nsol+is] /   pfmdat->D[npha-1][is*nsol+is] ) * fabs(scalprodct[ip][6]) );
      }
    }

    for ( ip = 0; ip < npha-1; ip++ ) {
      for ( is = 0; is < nsol; is++ ) {
          for (jj=0; jj< 6; jj++) {
            if ( modgradphi[ip][jj+1] == 0.0 ) {
              jat[ip][is][jj] = 0.0;
          }
        }
      }
    }

    for ( ip = 0; ip < npha-1; ip++ ) {
      for ( is = 0; is < nsol; is++ ) {
          cjaty = ( jat[ip][is][1] - jat[ip][is][0] ) / ( pfmvar->deltay );
          cjatz = ( jat[ip][is][3] - jat[ip][is][2] ) / ( pfmvar->deltaz );
          cjatx = ( jat[ip][is][5] - jat[ip][is][4] ) / ( pfmvar->deltax );
          jatc[ip][is] = cjatx + cjaty + cjatz;
	    }
    }

    for ( is = 0; is < nsol; is++ ) {
      jatr[is] = 0.0;
      for ( ip = 0; ip < npha-1; ip++ ) {
        jatr[is] += jatc[ip][is];
      }
    }

    for ( is1 = 0; is1 < nsol; is1++ ) {
        sum[is1] = 0.0;
        sum_dcbdT[is1] = 0.0;
        for ( is2 = 0; is2 < nsol; is2++ ) {
          dcdmu[is1*nsol+is2] = 0.0;
        }
      }

      for ( ip = 0; ip < npha; ip++ ) {
        for ( is = 0; is < nsol; is++ ) {
          dcbdT_phase[ip][is] = 0;
        }
      }

      for ( ii = 0; ii < 7; ii++ ) {

        interface = 1;
        bulkphase = 0;
        for ( ip = 0; ip < npha; ip++ ) {
          if ( stgO[ii].phi[ip] >= pfmdat->interfaceUplimit ) {
            bulkphase = ip;
            interface = 0;
            break;
          }
        }

        if ( interface ) {

          //stcscl[2].comie[ip][is]

          //if ( pfmdat->ISOTHERMAL ) {
            for ( ip = 0; ip < npha; ip++ ) {
              for ( is1 = 0; is1 < nsol; is1++ ) {
                for ( is2 =0; is2 < nsol; is2++ ) {
                  dc_dmu[ii][ip][is1*nsol+is2] = propf4->cmu[ip][is1][is2];
                  //printf("IC: %d, %d, %d, %d, %d, %d, %d, %d, %le, %le\n", tstep[0], x, y, z, index, ip, is1, is2, dc_dmu[ii][ip][is1*nsol+is2], propf4->cmu[ip][is1][is2]);
                }
              }
            }
          //}
          /*else {
            for ( ip = 0; ip < npha; ip++ ) {
              for ( is1 = 0; is1 < nsol; is1++ ) {
                for ( is2 =0; is2 < nsol; is2++ ) {
                  if ( is1 == is2 ) {
                    muc[ip][is1*nsol+is2] = 2.0 * propf4spline[indij[ii]].A[ip][is1][is2];
                  }
                  else {
                    muc[ip][is1*nsol+is2] = propf4spline[indij[ii]].A[ip][is1][is2];
                  }
                }
              }
              matinvnew_nsol(muc[ip], dc_dmu[ii][ip]);
            }
          }*/

          //printf("%d, %d, %d\n", tstep[0], i, j);
          for ( is = 0; is < nsol; is++ ) {
            for ( js =0; js < nsol; js++ ) {
              Da[ii][is][js] = 0.0;
              for ( ip = 0; ip < npha; ip++ ) {
                for ( ks = 0; ks < nsol; ks++ ) {
                  Da[ii][is][js] += pfmdat->D[ip][is*nsol+ks] * dc_dmu[ii][ip][ks*nsol+js] * stgO[ii].phi[ip];
                  //Da[ii][is][js] += pfmdat->D[ip][is*nsol+ks] * dc_dmu[ii][ip][ks+nsol*js] * stgO[ii].phi[ip];
                  //printf("BC: %d, %d, %d, %d, %d, %d, %d, %d, %d, %le\n", tstep[0], x, y, z, index, is, js, ip, ks, Da[ii][is][js]);
                }
              }
            }
          }

        }
        else {

          // bulk compositon

          //if ( pfmdat->ISOTHERMAL ) {
            for ( is1 = 0; is1 < nsol; is1++ ) {
              for ( is2 = 0; is2 < nsol; is2++ ) {
                dc_dmu[ii][bulkphase][is1*nsol+is2] = propf4->cmu[bulkphase][is1][is2];
                //printf("BC: %d, %d, %d, %d, %d, %d, %d, %d, %le, %le\n", tstep[0], x, y, z, index, ip, is1, is2, dc_dmu[ii][ip][is1*nsol+is2], propf4->cmu[ip][is1][is2]);
              }
            }
          //}
          /*else {
            for ( is1 = 0; is1 < nsol; is1++ ) {
              for ( is2 =0; is2 < nsol; is2++ ) {
                if ( is1 == is2 ) {
                  muc[bulkphase][is1*nsol+is2] = 2.0 * propf4spline[indij[ii]].A[bulkphase][is1][is2];
                }
                else {
                  muc[bulkphase][is1*nsol+is2] = propf4spline[indij[ii]].A[bulkphase][is1][is2];
                }
              }
              matinvnew_nsol(muc[bulkphase], dc_dmu[ii][bulkphase]);
            }
          }*/


          for ( is = 0; is < nsol; is++ ) {
            for ( js =0; js < nsol; js++ ) {
              Da[ii][is][js] = 0.0;
              //for ( ip = 0; ip < npha; ip++ ) {
                for ( ks = 0; ks < nsol; ks++ ) {
                  Da[ii][is][js] += pfmdat->D[bulkphase][is*nsol+ks] * dc_dmu[ii][bulkphase][ks*nsol+js];
                  //Da[ii][is][js] += pfmdat->D[bulkphase][is*nsol+ks] * dc_dmu[ii][bulkphase][ks+nsol*js];
                  //printf("BC: %d, %d, %d, %d, %d, %d, %d, %d, %d, %le\n", tstep[0], x, y, z, index, is, js, bulkphase, ks, Da[ii][is][js]);
                }
              //}
            }
          }

        }

        //for ( is = 0; is < nsol; is++ ) {
        //  for ( js =0; js < nsol; js++ ) {
        //    printf("%d, %d, %d, %d, %d, %d, %d, %d, %le\n", tstep[0], x, y, z, index, is, js, interface, Da[ii][is][js]);
        //  }
        //}


      }

      for ( is = 0; is < nsol; is++ ) {
        for ( js =0; js < nsol; js++ ) {

          Damidy[0][is][js] = ( Da[0][is][js] + Da[1][is][js] ) / 2.0;
          Damidy[1][is][js] = ( Da[2][is][js] + Da[0][is][js] ) / 2.0;

          Damidz[0][is][js] = ( Da[0][is][js] + Da[3][is][js] ) / 2.0;
          Damidz[1][is][js] = ( Da[4][is][js] + Da[0][is][js] ) / 2.0;

          Damidx[0][is][js] = ( Da[0][is][js] + Da[5][is][js] ) / 2.0;
          Damidx[1][is][js] = ( Da[6][is][js] + Da[0][is][js] ) / 2.0;

        }
      }

        // for ( is = 0; is < nsol; is++ ) {
        //  for ( js =0; js < nsol; js++ ) {
        //    printf("%d, %d, %d, %d, %d, %d, %d, %le, %le, %le, %le, %le, %le\n", tstep[0], x, y, z, index, is, js, Damidy[0][is][js], Damidy[1][is][js], Damidz[0][is][js], Damidz[1][is][js], Damidx[0][is][js], Damidx[1][is][js]);
        //  }
        // }

      for ( is = 0; is < nsol; is++ ) {
        gradmuy[0][is] = ( stgO[0].mu[is] - stgO[1].mu[is] ) / pfmvar->deltay;
        gradmuy[1][is] = ( stgO[2].mu[is] - stgO[0].mu[is] ) / pfmvar->deltay;

        gradmuz[0][is] = ( stgO[0].mu[is] - stgO[3].mu[is] ) / pfmvar->deltaz;
        gradmuz[1][is] = ( stgO[4].mu[is] - stgO[0].mu[is] ) / pfmvar->deltaz;

        gradmux[0][is] = ( stgO[0].mu[is] - stgO[5].mu[is] ) / pfmvar->deltax;
        gradmux[1][is] = ( stgO[6].mu[is] - stgO[0].mu[is] ) / pfmvar->deltax;
      }

        //for ( is = 0; is < nsol; is++ ) {
        //   printf("%d, %d, %d, %d, %d, %d, %le, %le, %le, %le, %le, %le\n", tstep[0], x, y, z, index, is, gradmuy[0][is], gradmuy[1][is], gradmuz[0][is], gradmuz[1][is], gradmux[0][is], gradmux[1][is]);
        //}

      for ( is = 0; is < nsol; is++ ) {
        divflux[is] = 0.0;
        for ( js = 0; js < nsol; js++ ) {
          //printf("x:%d, %d, %d, %d, %d, %d, %d, %le, %le, %le, %le, %le, %le\n", tstep[0], x, y, z, index, is, js, Damidy[1][is][js], Damidy[0][is][js], Damidz[1][is][js], Damidz[0][is][js], Damidx[1][is][js], Damidx[0][is][js]);
          //printf("x:%d, %d, %d, %d, %d, %d, %d, %le, %le, %le, %le, %le, %le\n", tstep[0], x, y, z, index, is, js, gradmuy[1][js], gradmuy[0][js], gradmuz[1][js], gradmuz[0][js], gradmux[1][js], gradmux[0][js]);

          divflux[is] +=  ( Damidy[1][is][js] * gradmuy[1][js] - Damidy[0][is][js] * gradmuy[0][js] ) / pfmvar->deltay;

          divflux[is] += ( Damidz[1][is][js] * gradmuz[1][js] - Damidz[0][is][js] * gradmuz[0][js] ) / pfmvar->deltaz;

          divflux[is] += ( Damidx[1][is][js] * gradmux[1][js] - Damidx[0][is][js] * gradmux[0][js] ) / pfmvar->deltax;

        }
        //printf("x:%d, %d, %d, %d, %d, %le, %le, %le, %le, %le\n", tstep[0], x, y, z, index, ( Damidy[1][is][js] * gradmuy[1][js] - Damidy[0][is][js] * gradmuy[0][js] ), ( Damidz[1][is][js] * gradmuz[1][js] - Damidz[0][is][js] * gradmuz[0][js] ), ( Damidx[1][is][js] * gradmux[1][js] - Damidx[0][is][js] * gradmux[0][js] ), jatr[1], pfmvar->deltat);
      }

      //printf("x:%d, %d, %d, %d, %d, %le, %le, %le, %le, %le\n", tstep[0], x, y, z, index, divflux[0], divflux[1], jatr[0], jatr[1], pfmvar->deltat);
      //printf("%d, %d, %d, %d, %d, %le, %le\n", tstep[0], x, y, z, index, divflux[0], divflux[1]);

      interface = 1;
      bulkphase = 0;
      for ( ip = 0; ip < npha; ip++ ) {
        if ( stgO[0].phi[ip] >= pfmdat->interfaceUplimit ) {
          bulkphase = ip;
          interface = 0;
          break;
        }
      }

      if ( !interface ) {

        for ( is = 0; is < nsol; is++ ) {
          gridinfo[index].com[is] = gridinfoO[index].com[is] + pfmvar->deltat * ( divflux[is] + jatr[is] );
        }
        //printf("BC:%d, %d, %d, %d, %d, %le, %le, %le, %le, %le\n", tstep[0], x, y, z, index, divflux[0], divflux[1], jatr[0], jatr[1], pfmvar->deltat);

        ip = bulkphase;
        //if ( pfmdat->ISOTHERMAL ) {
          tmp0 = 0.0;
          for ( is1 = 0; is1 < nsol; is1++ ) {
            tmp0  = 2.0*propf4->A[bulkphase][is1][is1]*gridinfo[index].com[is1] + propf4->B[bulkphase][is1];
            for ( is2 = 0; is2 < nsol; is2++ ) {
              if ( is1 != is2 ) {
                tmp0 += propf4->A[bulkphase][is1][is2]*gridinfo[index].com[is2];
              }
            }
            mu[is1] = tmp0;
          }
        //}
        /*else {
          tmp0 = 0.0;
          for ( is1 = 0; is1 < nsol; is1++ ) {
            tmp0  = 2.0*propf4spline[j].A[bulkphase][is1][is1]*gridinfo[index].com[is1] + propf4spline[j].B[bulkphase][is1];
            for ( is2 = 0; is2 < nsol; is2++ ) {
              if ( is1 != is2 ) {
                tmp0 += propf4spline[j].A[bulkphase][is1][is2]*gridinfo[index].com[is2];
              }
            }
            mu[is1] = tmp0;
          }
        }*/

        for ( is = 0; is < nsol; is++ ) {
          //deltamu[is] = mu[is] - stgO[0].mu[is];
          gridinfo[index].mu[is] = mu[is];
        }
        //printf("BC %d, %d, %d, %le, %le, %le, %le\n", tstep[0], i, j, gridOld[index].mu[0], gridinfo[index].mu[0], gridOld[index].com[0], gridinfo[index].com[0]);
      }
      else {

        /*if ( pfmdat->ISOTHERMAL ) {
          for ( ip = 0; ip < npha; ip++ ) {
            for ( is = 0; is < nsol; is++ ) {
              dcbdT_phase[ip][is] = 0.0;
            }
          }
        }*/
        /*else {

          DELTAT = ( pfmvar->deltat ) * ( -pfmdat->TGRADIENT * pfmdat->velocity );

          for ( ip = 0; ip < npha; ip++ ) {
            for ( is1 = 0; is1 < nsol; is1++ ) {
              for ( is2 = 0; is2 < nsol; is2++ ) {
                if ( is1 == is2 ) {
                  muc1[is1*nsol+is2] = 2.0 * propf4spline1[j].A[ip][is1][is2];
                }
                else {
                  muc1[is1*nsol+is2] = propf4spline1[j].A[ip][is1][is2];
                }
              }
            }
            matinvnew_nsol(muc1, cmu1);

            for ( is1 = 0; is1 < nsol; is1++ ) {
              tmp0 = 0.0;
              for ( is2 = 0; is2 < nsol; is2++ ) {
                tmp0 += cmu1[is1*nsol+is2] * ( gridOld[index].mu[is2] - propf4spline1[j].B[ip][is2] );
              }
              c_tdt[is1] = tmp0;
            }

            for ( is = 0; is < nsol; is++ ) {
              dcbdT_phase[ip][is] = c_tdt[is] - stcscl[4].comie[ip][is];
            }

          }
        }*/

        for ( ip1 = 0; ip1 < npha; ip1++ ) {
          sum_dhphi = 0.0;
          for ( ip2 = 0; ip2 < npha; ip2++ ) {
            sum_dhphi += dhfhi(stgO[0].phi, ip1, ip2) * ( gridinfo[index].phi[ip2] - gridinfoO[index].phi[ip2] );
          }

          for ( is = 0; is < nsol; is++ ) {
            sum[is]       += stcscl[0].comie[ip1][is] * sum_dhphi;
            // sum_dcbdT[is] += dcbdT_phase[ip1][is] * hfhi(stgO[0].phi, ip1);
          }
        }

        for ( is1 = 0; is1 < nsol; is1++ ) {

          deltac[is1] = pfmvar->deltat * ( divflux[is1] + jatr[is1] );



          gridinfo[index].com[is1] = stgO[0].com[is1] + deltac[is1];

          //if ( pfmdat->ISOTHERMAL ) {
            deltac[is1] += -sum[is1];
          //}
          //else {
          //  deltac[is1] += -sum[is1] - sum_dcbdT[is1];
          //}
          for ( is2 = 0; is2 < nsol; is2++ ) {
            for ( ip = 0; ip < npha; ip++ ) {
              dcdmu[is1*nsol+is2] += dc_dmu[0][ip][is1*nsol+is2] * hfhi(stgO[0].phi, ip);
            }
          }
        }
        //printf("IC:%d, %d, %d, %d, %d, %le, %le, %le, %le, %le\n", tstep[0], x, y, z, index, divflux[0], divflux[1], jatr[0], jatr[1], pfmvar->deltat);

        matinvnew_nsol(dcdmu, inv_dcdmu);

        multiply_nsol(inv_dcdmu, deltac, deltamu);
        // for ( is1 = 0; is1 < nsol; is1++ ) {
        //   tmp0 = 0.0;
        //   for ( is2 = 0; is2 < nsol; is2++ ) {
        //     tmp0 += inv_dcdmu[is1*nsol+is2] * deltac[is2];
        //   }
        //   deltamu[is1] = tmp0;
        // }

        //vectorsum(deltamu, stgO[index].mu, gridinfo[index].mu, nsol);
        for ( is = 0; is < nsol; is++ ) {
          gridinfo[index].mu[is] = gridinfoO[index].mu[is] + deltamu[is];
        }
        //printf("IC %d, %d, %d, %le, %le\n", tstep[0], i, j, gridinfo[index].mu[0], gridinfo[index].com[0]);
        //printf("IC %d, %d, %d, %le, %le, %le, %le\n", tstep[0], i, j, gridOld[index].mu[0], gridinfo[index].mu[0], gridOld[index].com[0], gridinfo[index].com[0]);
      }
      //printf("%d, %d, %d, %d, %d, %le, %le, %le, %le\n", tstep[0], x, y, z, index, divflux[0], divflux[1], jatr[0], jatr[1]);
      //printf("%d, %d, %d, %d, %d, %le, %le, %le, %le\n", tstep[0], x, y, z, index, gridinfo[index].com[0], gridinfoO[index].com[0], gridinfo[index].com[1], gridinfoO[index].com[1]);
      //printf("%d, %d, %d, %d, %d, %le, %le, %le, %le\n", tstep[0], x, y, z, index, stgO[0].com[0], gridinfoO[index].com[0], stgO[0].com[1], gridinfoO[index].com[1]);

    }


}



__kernel void SolverStress_iterative(__global struct fields *gridinfoO, __global struct iter_variables *it_gridinfoO, __global struct symmetric_tensor *eigen_strain_phase, __global struct Stiffness_cubic *stiffness_phase, __global struct Stiffness_cubic *stiffness_phase_n, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global int *tstep) {

  int x, y, z, ii, i1, j1, k1, i2;
  int nx, ny, nz;
  int index;
  int X, Y, Z;
  int center;

  int is, js, ks, il, il1, jl, jl1, ipx;
  int is1, is2;
  int ig, jg, ip, ip1, ip2, ip3;

  double deltat_e, damping_factor,  rho;
 double trace_strain_right;
 double trace_strain_left;
 double trace_strain_back;
 double trace_strain_front;
 double trace_strain_top;
 double trace_strain_bottom;

  double lambda_front   ;
  double lambda_back    ;
  double mu_front       ;
  double mu_back        ;
  double mu_prime_front ;
  double mu_prime_back  ;
  double lambda_right   ;
  double lambda_left    ;
  double mu_right       ;
  double mu_left        ;
  double mu_prime_right ;
  double mu_prime_left  ;
  double lambda_top     ;
  double lambda_bottom  ;
  double mu_top         ;
  double mu_bottom      ;
  double mu_prime_top   ;
  double mu_prime_bottom;

  double sigma_xx_front;
  double sigma_xx_back;
  double sigma_yx_right;
  double sigma_yx_left;
  double sigma_xy_front;
  double sigma_xy_back;
  double sigma_yy_right;
  double sigma_yy_left;

 double forceX;
 double forceY;
 double forceZ;

 double div_phi_front;
 double div_phi_back;
 double div_phi_right;
 double div_phi_left;
 double div_phi_top;
 double div_phi_bottom;


  struct symmetric_tensor eigen_strain[3], eigen_strain_lf[3], eigen_strain_rt[3], eigen_strain_bt[3], eigen_strain_tp[3], eigen_strain_bk[3], eigen_strain_ft[3];

  //struct symmetric_tensor eigen_strain_phase[npha];

  struct symmetric_tensor strain[3], strain_lf[3], strain_rt[3], strain_bt[3], strain_tp[3], strain_bk[3], strain_ft[3];

  struct Stiffness_cubic stiffness_c[3], stiffness_c_lf[3], stiffness_c_rt[3], stiffness_c_bt[3], stiffness_c_tp[3], stiffness_c_bk[3], stiffness_c_ft[3];

  //struct Stiffness_cubic stiffness_phase_n[npha];

 struct symmetric_tensor sigma_front;
 struct symmetric_tensor sigma_back;
 struct symmetric_tensor sigma_right;
 struct symmetric_tensor sigma_left;
 struct symmetric_tensor sigma_top;
 struct symmetric_tensor sigma_bottom;

  struct fields stgO[7], stgO_lf[7], stgO_rt[7], stgO_bt[7], stgO_tp[7], stgO_bk[7], stgO_ft[7];
  struct iter_variables st_it_gO[7], st_it_gO_lf[7], st_it_gO_rt[7], st_it_gO_bt[7], st_it_gO_tp[7], st_it_gO_bk[7], st_it_gO_ft[7];

  X = 0;
  Y = 1;
  Z = 2;

  deltat_e = pfmdat->deltat_e;
  rho = pfmdat->rho;
  damping_factor = pfmdat->damping_factor;

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  //index = y + ny*(z + x*nz);

  if ( nz > 4 ) {
  if ( ( x > 1 ) && ( x < (nx-2) ) && ( y > 1 ) && ( y < (ny-2) ) && ( z > 1 ) && ( z < (nz-2) ) ) {
  //if ( x > 1 && x < nx-2 && y > 1 && y < ny-2 && z > 1 && z < nz-2 ) {

    index = y + ny*(z + x*nz);
    center = index;

    printf("%d, %d, %d, %d, %le, %le, %le\n", x, y, z , index, it_gridinfoO[center].disp[X][2], it_gridinfoO[center].disp[Y][2], it_gridinfoO[center].disp[Z][2]);

    stgO[0] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stgO[1] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stgO[2] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stgO[3] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stgO[4] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stgO[5] = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stgO[6] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];

    st_it_gO[0] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO[1] = it_gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO[2] = it_gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO[3] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    st_it_gO[4] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    st_it_gO[5] = it_gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    st_it_gO[6] = it_gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];

    eigen_strain[0].xx = 0.0;
    eigen_strain[0].yy = 0.0;
    eigen_strain[0].zz = 0.0;
    eigen_strain[0].yz = 0.0;
    eigen_strain[0].xz = 0.0;
    eigen_strain[0].xy = 0.0;

    for (ip = 0; ip < npha; ip++) {

      eigen_strain[X].xx += eigen_strain_phase[ip].xx*stgO[0].phi[ip];
      eigen_strain[X].yy += eigen_strain_phase[ip].yy*stgO[0].phi[ip];
      eigen_strain[X].zz += eigen_strain_phase[ip].zz*stgO[0].phi[ip];
      eigen_strain[X].yz += eigen_strain_phase[ip].yz*stgO[0].phi[ip];
      eigen_strain[X].xz += eigen_strain_phase[ip].xz*stgO[0].phi[ip];
      eigen_strain[X].xy += eigen_strain_phase[ip].xy*stgO[0].phi[ip];

    }

    for (i1 = 1; i1 < 3; i1++) {

      eigen_strain[i1].xx = eigen_strain[X].xx;
      eigen_strain[i1].yy = eigen_strain[X].yy;
      eigen_strain[i1].zz = eigen_strain[X].zz;
      eigen_strain[i1].yz = eigen_strain[X].yz;
      eigen_strain[i1].xz = eigen_strain[X].xz;
      eigen_strain[i1].xy = eigen_strain[X].xy;

    }
    // or
    eigen_strain[Y] = eigen_strain[X];
    eigen_strain[Z] = eigen_strain[X];

    strain[X].xx = 0.5*(st_it_gO[6].disp[X][2] - st_it_gO[5].disp[X][2]) - eigen_strain[X].xx;
    strain[X].yy = 0.5*(st_it_gO[2].disp[Y][2] - st_it_gO[1].disp[Y][2]) - eigen_strain[X].yy;
    strain[X].zz = 0.5*(st_it_gO[4].disp[Z][2] - st_it_gO[3].disp[Z][2]) - eigen_strain[X].zz;

    strain[Y].xx = strain[X].xx;
    strain[Y].yy = strain[X].yy;
    strain[Y].zz = strain[X].zz;
    strain[Z].yy = strain[X].yy;
    strain[Z].xx = strain[X].xx;
    strain[Z].zz = strain[X].zz;

    strain[X].xy = 0.25*((st_it_gO[2].disp[X][2] - st_it_gO[1].disp[X][2]) + (st_it_gO[6].disp[Y][2] - st_it_gO[5].disp[Y][2]));
    strain[X].xz = 0.25*((st_it_gO[4].disp[X][2] - st_it_gO[3].disp[X][2]) + (st_it_gO[6].disp[Z][2] - st_it_gO[5].disp[Z][2]));
    strain[X].yz = 0.25*((st_it_gO[2].disp[Z][2] - st_it_gO[1].disp[Z][2]) + (st_it_gO[4].disp[Y][2] - st_it_gO[3].disp[Y][2]));

    strain[Y].xy = strain[X].xy;
    strain[Y].xz = strain[X].xz;
    strain[Y].yz = strain[X].yz;
    strain[Z].xy = strain[X].xy;
    strain[Z].xz = strain[X].xz;
    strain[Z].yz = strain[X].yz;

    stiffness_c[0].C11 = 0.0;
    stiffness_c[0].C12 = 0.0;
    stiffness_c[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c[X].C11 += (stiffness_phase_n[ip].C11)*stgO[0].phi[ip];
      stiffness_c[X].C12 += (stiffness_phase_n[ip].C12)*stgO[0].phi[ip];
      stiffness_c[X].C44 += (stiffness_phase_n[ip].C44)*stgO[0].phi[ip];
    }
    stiffness_c[Y] = stiffness_c[X];
    stiffness_c[Z] = stiffness_c[X];


    stgO_lf[0] = gridinfoO[ (y  -1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_lf[1] = gridinfoO[ (y-1-1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_lf[2] = gridinfoO[ (y+1-1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_lf[3] = gridinfoO[ (y  -1) + ny*( (x  )*nz + (z-1) ) ];
    stgO_lf[4] = gridinfoO[ (y  -1) + ny*( (x  )*nz + (z+1) ) ];
    stgO_lf[5] = gridinfoO[ (y  -1) + ny*( (x-1)*nz + (z  ) ) ];
    stgO_lf[6] = gridinfoO[ (y  -1) + ny*( (x+1)*nz + (z  ) ) ];

    st_it_gO_lf[0] = it_gridinfoO[ (y  -1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_lf[1] = it_gridinfoO[ (y-1-1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_lf[2] = it_gridinfoO[ (y+1-1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_lf[3] = it_gridinfoO[ (y  -1) + ny*( (x  )*nz + (z-1) ) ];
    st_it_gO_lf[4] = it_gridinfoO[ (y  -1) + ny*( (x  )*nz + (z+1) ) ];
    st_it_gO_lf[5] = it_gridinfoO[ (y  -1) + ny*( (x-1)*nz + (z  ) ) ];
    st_it_gO_lf[6] = it_gridinfoO[ (y  -1) + ny*( (x+1)*nz + (z  ) ) ];


    stgO_rt[0] = gridinfoO[ (y  +1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_rt[1] = gridinfoO[ (y-1+1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_rt[2] = gridinfoO[ (y+1+1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_rt[3] = gridinfoO[ (y  +1) + ny*( (x  )*nz + (z-1) ) ];
    stgO_rt[4] = gridinfoO[ (y  +1) + ny*( (x  )*nz + (z+1) ) ];
    stgO_rt[5] = gridinfoO[ (y  +1) + ny*( (x-1)*nz + (z  ) ) ];
    stgO_rt[6] = gridinfoO[ (y  +1) + ny*( (x+1)*nz + (z  ) ) ];

    st_it_gO_rt[0] = it_gridinfoO[ (y  +1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_rt[1] = it_gridinfoO[ (y-1+1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_rt[2] = it_gridinfoO[ (y+1+1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_rt[3] = it_gridinfoO[ (y  +1) + ny*( (x  )*nz + (z-1) ) ];
    st_it_gO_rt[4] = it_gridinfoO[ (y  +1) + ny*( (x  )*nz + (z+1) ) ];
    st_it_gO_rt[5] = it_gridinfoO[ (y  +1) + ny*( (x-1)*nz + (z  ) ) ];
    st_it_gO_rt[6] = it_gridinfoO[ (y  +1) + ny*( (x+1)*nz + (z  ) ) ];


    stgO_bt[0] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  -1) ) ];
    stgO_bt[1] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  -1) ) ];
    stgO_bt[2] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  -1) ) ];
    stgO_bt[3] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1-1) ) ];
    stgO_bt[4] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1-1) ) ];
    stgO_bt[5] = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  -1) ) ];
    stgO_bt[6] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  -1) ) ];

    st_it_gO_bt[0] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  -1) ) ];
    st_it_gO_bt[1] = it_gridinfoO[ (y-1) + ny*( (x  )*nz + (z  -1) ) ];
    st_it_gO_bt[2] = it_gridinfoO[ (y+1) + ny*( (x  )*nz + (z  -1) ) ];
    st_it_gO_bt[3] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1-1) ) ];
    st_it_gO_bt[4] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1-1) ) ];
    st_it_gO_bt[5] = it_gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  -1) ) ];
    st_it_gO_bt[6] = it_gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  -1) ) ];


    stgO_tp[0] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  +1) ) ];
    stgO_tp[1] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  +1) ) ];
    stgO_tp[2] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  +1) ) ];
    stgO_tp[3] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1+1) ) ];
    stgO_tp[4] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1+1) ) ];
    stgO_tp[5] = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  +1) ) ];
    stgO_tp[6] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  +1) ) ];

    st_it_gO_tp[0] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  +1) ) ];
    st_it_gO_tp[1] = it_gridinfoO[ (y-1) + ny*( (x  )*nz + (z  +1) ) ];
    st_it_gO_tp[2] = it_gridinfoO[ (y+1) + ny*( (x  )*nz + (z  +1) ) ];
    st_it_gO_tp[3] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1+1) ) ];
    st_it_gO_tp[4] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1+1) ) ];
    st_it_gO_tp[5] = it_gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  +1) ) ];
    st_it_gO_tp[6] = it_gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  +1) ) ];


    stgO_bk[0] = gridinfoO[ (y  ) + ny*( (x  -1)*nz + (z  ) ) ];
    stgO_bk[1] = gridinfoO[ (y-1) + ny*( (x  -1)*nz + (z  ) ) ];
    stgO_bk[2] = gridinfoO[ (y+1) + ny*( (x  -1)*nz + (z  ) ) ];
    stgO_bk[3] = gridinfoO[ (y  ) + ny*( (x  -1)*nz + (z-1) ) ];
    stgO_bk[4] = gridinfoO[ (y  ) + ny*( (x  -1)*nz + (z+1) ) ];
    stgO_bk[5] = gridinfoO[ (y  ) + ny*( (x-1-1)*nz + (z  ) ) ];
    stgO_bk[6] = gridinfoO[ (y  ) + ny*( (x+1-1)*nz + (z  ) ) ];

    st_it_gO_bk[0] = it_gridinfoO[ (y  ) + ny*( (x  -1)*nz + (z  ) ) ];
    st_it_gO_bk[1] = it_gridinfoO[ (y-1) + ny*( (x  -1)*nz + (z  ) ) ];
    st_it_gO_bk[2] = it_gridinfoO[ (y+1) + ny*( (x  -1)*nz + (z  ) ) ];
    st_it_gO_bk[3] = it_gridinfoO[ (y  ) + ny*( (x  -1)*nz + (z-1) ) ];
    st_it_gO_bk[4] = it_gridinfoO[ (y  ) + ny*( (x  -1)*nz + (z+1) ) ];
    st_it_gO_bk[5] = it_gridinfoO[ (y  ) + ny*( (x-1-1)*nz + (z  ) ) ];
    st_it_gO_bk[6] = it_gridinfoO[ (y  ) + ny*( (x+1-1)*nz + (z  ) ) ];


    stgO_ft[0] = gridinfoO[ (y  ) + ny*( (x  +1)*nz + (z  ) ) ];
    stgO_ft[1] = gridinfoO[ (y-1) + ny*( (x  +1)*nz + (z  ) ) ];
    stgO_ft[2] = gridinfoO[ (y+1) + ny*( (x  +1)*nz + (z  ) ) ];
    stgO_ft[3] = gridinfoO[ (y  ) + ny*( (x  +1)*nz + (z-1) ) ];
    stgO_ft[4] = gridinfoO[ (y  ) + ny*( (x  +1)*nz + (z+1) ) ];
    stgO_ft[5] = gridinfoO[ (y  ) + ny*( (x-1+1)*nz + (z  ) ) ];
    stgO_ft[6] = gridinfoO[ (y  ) + ny*( (x+1+1)*nz + (z  ) ) ];

    st_it_gO_ft[0] = it_gridinfoO[ (y  ) + ny*( (x  +1)*nz + (z  ) ) ];
    st_it_gO_ft[1] = it_gridinfoO[ (y-1) + ny*( (x  +1)*nz + (z  ) ) ];
    st_it_gO_ft[2] = it_gridinfoO[ (y+1) + ny*( (x  +1)*nz + (z  ) ) ];
    st_it_gO_ft[3] = it_gridinfoO[ (y  ) + ny*( (x  +1)*nz + (z-1) ) ];
    st_it_gO_ft[4] = it_gridinfoO[ (y  ) + ny*( (x  +1)*nz + (z+1) ) ];
    st_it_gO_ft[5] = it_gridinfoO[ (y  ) + ny*( (x-1+1)*nz + (z  ) ) ];
    st_it_gO_ft[6] = it_gridinfoO[ (y  ) + ny*( (x+1+1)*nz + (z  ) ) ];

    eigen_strain_lf[0].xx = 0.0;
    eigen_strain_lf[0].yy = 0.0;
    eigen_strain_lf[0].zz = 0.0;
    eigen_strain_lf[0].yz = 0.0;
    eigen_strain_lf[0].xz = 0.0;
    eigen_strain_lf[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_lf[X].xx += eigen_strain_phase[ip].xx*stgO_lf[0].phi[ip];
      eigen_strain_lf[X].yy += eigen_strain_phase[ip].yy*stgO_lf[0].phi[ip];
      eigen_strain_lf[X].zz += eigen_strain_phase[ip].zz*stgO_lf[0].phi[ip];
      eigen_strain_lf[X].yz += eigen_strain_phase[ip].yz*stgO_lf[0].phi[ip];
      eigen_strain_lf[X].xz += eigen_strain_phase[ip].xz*stgO_lf[0].phi[ip];
      eigen_strain_lf[X].xy += eigen_strain_phase[ip].xy*stgO_lf[0].phi[ip];
    }
    eigen_strain_lf[Y] = eigen_strain_lf[X];
    eigen_strain_lf[Z] = eigen_strain_lf[X];

    eigen_strain_rt[0].xx = 0.0;
    eigen_strain_rt[0].yy = 0.0;
    eigen_strain_rt[0].zz = 0.0;
    eigen_strain_rt[0].yz = 0.0;
    eigen_strain_rt[0].xz = 0.0;
    eigen_strain_rt[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_rt[X].xx += eigen_strain_phase[ip].xx*stgO_rt[0].phi[ip];
      eigen_strain_rt[X].yy += eigen_strain_phase[ip].yy*stgO_rt[0].phi[ip];
      eigen_strain_rt[X].zz += eigen_strain_phase[ip].zz*stgO_rt[0].phi[ip];
      eigen_strain_rt[X].yz += eigen_strain_phase[ip].yz*stgO_rt[0].phi[ip];
      eigen_strain_rt[X].xz += eigen_strain_phase[ip].xz*stgO_rt[0].phi[ip];
      eigen_strain_rt[X].xy += eigen_strain_phase[ip].xy*stgO_rt[0].phi[ip];
    }
    eigen_strain_rt[Y] = eigen_strain_rt[X];
    eigen_strain_rt[Z] = eigen_strain_rt[X];

    eigen_strain_bt[0].xx = 0.0;
    eigen_strain_bt[0].yy = 0.0;
    eigen_strain_bt[0].zz = 0.0;
    eigen_strain_bt[0].yz = 0.0;
    eigen_strain_bt[0].xz = 0.0;
    eigen_strain_bt[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_bt[X].xx += eigen_strain_phase[ip].xx*stgO_bt[0].phi[ip];
      eigen_strain_bt[X].yy += eigen_strain_phase[ip].yy*stgO_bt[0].phi[ip];
      eigen_strain_bt[X].zz += eigen_strain_phase[ip].zz*stgO_bt[0].phi[ip];
      eigen_strain_bt[X].yz += eigen_strain_phase[ip].yz*stgO_bt[0].phi[ip];
      eigen_strain_bt[X].xz += eigen_strain_phase[ip].xz*stgO_bt[0].phi[ip];
      eigen_strain_bt[X].xy += eigen_strain_phase[ip].xy*stgO_bt[0].phi[ip];
    }
    eigen_strain_bt[Y] = eigen_strain_bt[X];
    eigen_strain_bt[Z] = eigen_strain_bt[X];

    eigen_strain_tp[0].xx = 0.0;
    eigen_strain_tp[0].yy = 0.0;
    eigen_strain_tp[0].zz = 0.0;
    eigen_strain_tp[0].yz = 0.0;
    eigen_strain_tp[0].xz = 0.0;
    eigen_strain_tp[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_tp[X].xx += eigen_strain_phase[ip].xx*stgO_tp[0].phi[ip];
      eigen_strain_tp[X].yy += eigen_strain_phase[ip].yy*stgO_tp[0].phi[ip];
      eigen_strain_tp[X].zz += eigen_strain_phase[ip].zz*stgO_tp[0].phi[ip];
      eigen_strain_tp[X].yz += eigen_strain_phase[ip].yz*stgO_tp[0].phi[ip];
      eigen_strain_tp[X].xz += eigen_strain_phase[ip].xz*stgO_tp[0].phi[ip];
      eigen_strain_tp[X].xy += eigen_strain_phase[ip].xy*stgO_tp[0].phi[ip];
    }
    eigen_strain_tp[Y] = eigen_strain_tp[X];
    eigen_strain_tp[Z] = eigen_strain_tp[X];

    eigen_strain_bk[0].xx = 0.0;
    eigen_strain_bk[0].yy = 0.0;
    eigen_strain_bk[0].zz = 0.0;
    eigen_strain_bk[0].yz = 0.0;
    eigen_strain_bk[0].xz = 0.0;
    eigen_strain_bk[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_bk[X].xx += eigen_strain_phase[ip].xx*stgO_bk[0].phi[ip];
      eigen_strain_bk[X].yy += eigen_strain_phase[ip].yy*stgO_bk[0].phi[ip];
      eigen_strain_bk[X].zz += eigen_strain_phase[ip].zz*stgO_bk[0].phi[ip];
      eigen_strain_bk[X].yz += eigen_strain_phase[ip].yz*stgO_bk[0].phi[ip];
      eigen_strain_bk[X].xz += eigen_strain_phase[ip].xz*stgO_bk[0].phi[ip];
      eigen_strain_bk[X].xy += eigen_strain_phase[ip].xy*stgO_bk[0].phi[ip];
    }
    eigen_strain_bk[Y] = eigen_strain_bk[X];
    eigen_strain_bk[Z] = eigen_strain_bk[X];

    eigen_strain_ft[0].xx = 0.0;
    eigen_strain_ft[0].yy = 0.0;
    eigen_strain_ft[0].zz = 0.0;
    eigen_strain_ft[0].yz = 0.0;
    eigen_strain_ft[0].xz = 0.0;
    eigen_strain_ft[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_ft[X].xx += eigen_strain_phase[ip].xx*stgO_ft[0].phi[ip];
      eigen_strain_ft[X].yy += eigen_strain_phase[ip].yy*stgO_ft[0].phi[ip];
      eigen_strain_ft[X].zz += eigen_strain_phase[ip].zz*stgO_ft[0].phi[ip];
      eigen_strain_ft[X].yz += eigen_strain_phase[ip].yz*stgO_ft[0].phi[ip];
      eigen_strain_ft[X].xz += eigen_strain_phase[ip].xz*stgO_ft[0].phi[ip];
      eigen_strain_ft[X].xy += eigen_strain_phase[ip].xy*stgO_ft[0].phi[ip];
    }
    eigen_strain_ft[Y] = eigen_strain_ft[X];
    eigen_strain_ft[Z] = eigen_strain_ft[X];

    strain_lf[X].xx = 0.5*(st_it_gO_lf[6].disp[X][2] - st_it_gO_lf[5].disp[X][2]) - eigen_strain_lf[X].xx;
    strain_lf[X].yy = 0.5*(st_it_gO_lf[2].disp[Y][2] - st_it_gO_lf[1].disp[Y][2]) - eigen_strain_lf[X].yy;
    strain_lf[X].zz = 0.5*(st_it_gO_lf[4].disp[Z][2] - st_it_gO_lf[3].disp[Z][2]) - eigen_strain_lf[X].zz;
    strain_lf[Y].xx = strain_lf[X].xx;
    strain_lf[Y].yy = strain_lf[X].yy;
    strain_lf[Y].zz = strain_lf[X].zz;
    strain_lf[Z].yy = strain_lf[X].yy;
    strain_lf[Z].xx = strain_lf[X].xx;
    strain_lf[Z].zz = strain_lf[X].zz;

    strain_rt[X].xx = 0.5*(st_it_gO_rt[6].disp[X][2] - st_it_gO_rt[5].disp[X][2]) - eigen_strain_rt[X].xx;
    strain_rt[X].yy = 0.5*(st_it_gO_rt[2].disp[Y][2] - st_it_gO_rt[1].disp[Y][2]) - eigen_strain_rt[X].yy;
    strain_rt[X].zz = 0.5*(st_it_gO_rt[4].disp[Z][2] - st_it_gO_rt[3].disp[Z][2]) - eigen_strain_rt[X].zz;
    strain_rt[Y].xx = strain_rt[X].xx;
    strain_rt[Y].yy = strain_rt[X].yy;
    strain_rt[Y].zz = strain_rt[X].zz;
    strain_rt[Z].yy = strain_rt[X].yy;
    strain_rt[Z].xx = strain_rt[X].xx;
    strain_rt[Z].zz = strain_rt[X].zz;

    strain_bt[X].xx = 0.5*(st_it_gO_bt[6].disp[X][2] - st_it_gO_bt[5].disp[X][2]) - eigen_strain_bt[X].xx;
    strain_bt[X].yy = 0.5*(st_it_gO_bt[2].disp[Y][2] - st_it_gO_bt[1].disp[Y][2]) - eigen_strain_bt[X].yy;
    strain_bt[X].zz = 0.5*(st_it_gO_bt[4].disp[Z][2] - st_it_gO_bt[3].disp[Z][2]) - eigen_strain_bt[X].zz;
    strain_bt[Y].xx = strain_bt[X].xx;
    strain_bt[Y].yy = strain_bt[X].yy;
    strain_bt[Y].zz = strain_bt[X].zz;
    strain_bt[Z].yy = strain_bt[X].yy;
    strain_bt[Z].xx = strain_bt[X].xx;
    strain_bt[Z].zz = strain_bt[X].zz;

    strain_tp[X].xx = 0.5*(st_it_gO_tp[6].disp[X][2] - st_it_gO_tp[5].disp[X][2]) - eigen_strain_tp[X].xx;
    strain_tp[X].yy = 0.5*(st_it_gO_tp[2].disp[Y][2] - st_it_gO_tp[1].disp[Y][2]) - eigen_strain_tp[X].yy;
    strain_tp[X].zz = 0.5*(st_it_gO_tp[4].disp[Z][2] - st_it_gO_tp[3].disp[Z][2]) - eigen_strain_tp[X].zz;
    strain_tp[Y].xx = strain_tp[X].xx;
    strain_tp[Y].yy = strain_tp[X].yy;
    strain_tp[Y].zz = strain_tp[X].zz;
    strain_tp[Z].yy = strain_tp[X].yy;
    strain_tp[Z].xx = strain_tp[X].xx;
    strain_tp[Z].zz = strain_tp[X].zz;

    strain_bk[X].xx = 0.5*(st_it_gO_bk[6].disp[X][2] - st_it_gO_bk[5].disp[X][2]) - eigen_strain_bk[X].xx;
    strain_bk[X].yy = 0.5*(st_it_gO_bk[2].disp[Y][2] - st_it_gO_bk[1].disp[Y][2]) - eigen_strain_bk[X].yy;
    strain_bk[X].zz = 0.5*(st_it_gO_bk[4].disp[Z][2] - st_it_gO_bk[3].disp[Z][2]) - eigen_strain_bk[X].zz;
    strain_bk[Y].xx = strain_bk[X].xx;
    strain_bk[Y].yy = strain_bk[X].yy;
    strain_bk[Y].zz = strain_bk[X].zz;
    strain_bk[Z].yy = strain_bk[X].yy;
    strain_bk[Z].xx = strain_bk[X].xx;
    strain_bk[Z].zz = strain_bk[X].zz;

    strain_ft[X].xx = 0.5*(st_it_gO_ft[6].disp[X][2] - st_it_gO_ft[5].disp[X][2]) - eigen_strain_ft[X].xx;
    strain_ft[X].yy = 0.5*(st_it_gO_ft[2].disp[Y][2] - st_it_gO_ft[1].disp[Y][2]) - eigen_strain_ft[X].yy;
    strain_ft[X].zz = 0.5*(st_it_gO_ft[4].disp[Z][2] - st_it_gO_ft[3].disp[Z][2]) - eigen_strain_ft[X].zz;
    strain_ft[Y].xx = strain_ft[X].xx;
    strain_ft[Y].yy = strain_ft[X].yy;
    strain_ft[Y].zz = strain_ft[X].zz;
    strain_ft[Z].yy = strain_ft[X].yy;
    strain_ft[Z].xx = strain_ft[X].xx;
    strain_ft[Z].zz = strain_ft[X].zz;

    strain_lf[X].xy = 0.25*((st_it_gO_lf[2].disp[X][2] - st_it_gO_lf[1].disp[X][2]) + (st_it_gO_lf[6].disp[Y][2] - st_it_gO_lf[5].disp[Y][2]));
    strain_lf[X].xz = 0.25*((st_it_gO_lf[4].disp[X][2] - st_it_gO_lf[3].disp[X][2]) + (st_it_gO_lf[6].disp[Z][2] - st_it_gO_lf[5].disp[Z][2]));
    strain_lf[X].yz = 0.25*((st_it_gO_lf[2].disp[Z][2] - st_it_gO_lf[1].disp[Z][2]) + (st_it_gO_lf[4].disp[Y][2] - st_it_gO_lf[3].disp[Y][2]));
    strain_lf[Y].xy = strain_lf[X].xy;
    strain_lf[Y].xz = strain_lf[X].xz;
    strain_lf[Y].yz = strain_lf[X].yz;
    strain_lf[Z].xy = strain_lf[X].xy;
    strain_lf[Z].xz = strain_lf[X].xz;
    strain_lf[Z].yz = strain_lf[X].yz;

    strain_rt[X].xy = 0.25*((st_it_gO_rt[2].disp[X][2] - st_it_gO_rt[1].disp[X][2]) + (st_it_gO_rt[6].disp[Y][2] - st_it_gO_rt[5].disp[Y][2]));
    strain_rt[X].xz = 0.25*((st_it_gO_rt[4].disp[X][2] - st_it_gO_rt[3].disp[X][2]) + (st_it_gO_rt[6].disp[Z][2] - st_it_gO_rt[5].disp[Z][2]));
    strain_rt[X].yz = 0.25*((st_it_gO_rt[2].disp[Z][2] - st_it_gO_rt[1].disp[Z][2]) + (st_it_gO_rt[4].disp[Y][2] - st_it_gO_rt[3].disp[Y][2]));
    strain_rt[Y].xy = strain_rt[X].xy;
    strain_rt[Y].xz = strain_rt[X].xz;
    strain_rt[Y].yz = strain_rt[X].yz;
    strain_rt[Z].xy = strain_rt[X].xy;
    strain_rt[Z].xz = strain_rt[X].xz;
    strain_rt[Z].yz = strain_rt[X].yz;

    strain_bt[X].xy = 0.25*((st_it_gO_bt[2].disp[X][2] - st_it_gO_bt[1].disp[X][2]) + (st_it_gO_bt[6].disp[Y][2] - st_it_gO_bt[5].disp[Y][2]));
    strain_bt[X].xz = 0.25*((st_it_gO_bt[4].disp[X][2] - st_it_gO_bt[3].disp[X][2]) + (st_it_gO_bt[6].disp[Z][2] - st_it_gO_bt[5].disp[Z][2]));
    strain_bt[X].yz = 0.25*((st_it_gO_bt[2].disp[Z][2] - st_it_gO_bt[1].disp[Z][2]) + (st_it_gO_bt[4].disp[Y][2] - st_it_gO_bt[3].disp[Y][2]));
    strain_bt[Y].xy = strain_bt[X].xy;
    strain_bt[Y].xz = strain_bt[X].xz;
    strain_bt[Y].yz = strain_bt[X].yz;
    strain_bt[Z].xy = strain_bt[X].xy;
    strain_bt[Z].xz = strain_bt[X].xz;
    strain_bt[Z].yz = strain_bt[X].yz;

    strain_tp[X].xy = 0.25*((st_it_gO_tp[2].disp[X][2] - st_it_gO_tp[1].disp[X][2]) + (st_it_gO_tp[6].disp[Y][2] - st_it_gO_tp[5].disp[Y][2]));
    strain_tp[X].xz = 0.25*((st_it_gO_tp[4].disp[X][2] - st_it_gO_tp[3].disp[X][2]) + (st_it_gO_tp[6].disp[Z][2] - st_it_gO_tp[5].disp[Z][2]));
    strain_tp[X].yz = 0.25*((st_it_gO_tp[2].disp[Z][2] - st_it_gO_tp[1].disp[Z][2]) + (st_it_gO_tp[4].disp[Y][2] - st_it_gO_tp[3].disp[Y][2]));
    strain_tp[Y].xy = strain_tp[X].xy;
    strain_tp[Y].xz = strain_tp[X].xz;
    strain_tp[Y].yz = strain_tp[X].yz;
    strain_tp[Z].xy = strain_tp[X].xy;
    strain_tp[Z].xz = strain_tp[X].xz;
    strain_tp[Z].yz = strain_tp[X].yz;

    strain_bk[X].xy = 0.25*((st_it_gO_bk[2].disp[X][2] - st_it_gO_bk[1].disp[X][2]) + (st_it_gO_bk[6].disp[Y][2] - st_it_gO_bk[5].disp[Y][2]));
    strain_bk[X].xz = 0.25*((st_it_gO_bk[4].disp[X][2] - st_it_gO_bk[3].disp[X][2]) + (st_it_gO_bk[6].disp[Z][2] - st_it_gO_bk[5].disp[Z][2]));
    strain_bk[X].yz = 0.25*((st_it_gO_bk[2].disp[Z][2] - st_it_gO_bk[1].disp[Z][2]) + (st_it_gO_bk[4].disp[Y][2] - st_it_gO_bk[3].disp[Y][2]));
    strain_bk[Y].xy = strain_bk[X].xy;
    strain_bk[Y].xz = strain_bk[X].xz;
    strain_bk[Y].yz = strain_bk[X].yz;
    strain_bk[Z].xy = strain_bk[X].xy;
    strain_bk[Z].xz = strain_bk[X].xz;
    strain_bk[Z].yz = strain_bk[X].yz;

    strain_ft[X].xy = 0.25*((st_it_gO_ft[2].disp[X][2] - st_it_gO_ft[1].disp[X][2]) + (st_it_gO_ft[6].disp[Y][2] - st_it_gO_ft[5].disp[Y][2]));
    strain_ft[X].xz = 0.25*((st_it_gO_ft[4].disp[X][2] - st_it_gO_ft[3].disp[X][2]) + (st_it_gO_ft[6].disp[Z][2] - st_it_gO_ft[5].disp[Z][2]));
    strain_ft[X].yz = 0.25*((st_it_gO_ft[2].disp[Z][2] - st_it_gO_ft[1].disp[Z][2]) + (st_it_gO_ft[4].disp[Y][2] - st_it_gO_ft[3].disp[Y][2]));
    strain_ft[Y].xy = strain_ft[X].xy;
    strain_ft[Y].xz = strain_ft[X].xz;
    strain_ft[Y].yz = strain_ft[X].yz;
    strain_ft[Z].xy = strain_ft[X].xy;
    strain_ft[Z].xz = strain_ft[X].xz;
    strain_ft[Z].yz = strain_ft[X].yz;

    stiffness_c_lf[0].C11 = 0.0;
    stiffness_c_lf[0].C12 = 0.0;
    stiffness_c_lf[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_lf[X].C11 += (stiffness_phase_n[ip].C11)*stgO_lf[0].phi[ip];
      stiffness_c_lf[X].C12 += (stiffness_phase_n[ip].C12)*stgO_lf[0].phi[ip];
      stiffness_c_lf[X].C44 += (stiffness_phase_n[ip].C44)*stgO_lf[0].phi[ip];
    }
    stiffness_c_lf[Y] = stiffness_c_lf[X];
    stiffness_c_lf[Z] = stiffness_c_lf[X];

    stiffness_c_rt[0].C11 = 0.0;
    stiffness_c_rt[0].C12 = 0.0;
    stiffness_c_rt[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_rt[X].C11 += (stiffness_phase_n[ip].C11)*stgO_rt[0].phi[ip];
      stiffness_c_rt[X].C12 += (stiffness_phase_n[ip].C12)*stgO_rt[0].phi[ip];
      stiffness_c_rt[X].C44 += (stiffness_phase_n[ip].C44)*stgO_rt[0].phi[ip];
    }
    stiffness_c_rt[Y] = stiffness_c_rt[X];
    stiffness_c_rt[Z] = stiffness_c_rt[X];

    stiffness_c_bt[0].C11 = 0.0;
    stiffness_c_bt[0].C12 = 0.0;
    stiffness_c_bt[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_bt[X].C11 += (stiffness_phase_n[ip].C11)*stgO_bt[0].phi[ip];
      stiffness_c_bt[X].C12 += (stiffness_phase_n[ip].C12)*stgO_bt[0].phi[ip];
      stiffness_c_bt[X].C44 += (stiffness_phase_n[ip].C44)*stgO_bt[0].phi[ip];
    }
    stiffness_c_bt[Y] = stiffness_c_bt[X];
    stiffness_c_bt[Z] = stiffness_c_bt[X];

    stiffness_c_tp[0].C11 = 0.0;
    stiffness_c_tp[0].C12 = 0.0;
    stiffness_c_tp[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_tp[X].C11 += (stiffness_phase_n[ip].C11)*stgO_tp[0].phi[ip];
      stiffness_c_tp[X].C12 += (stiffness_phase_n[ip].C12)*stgO_tp[0].phi[ip];
      stiffness_c_tp[X].C44 += (stiffness_phase_n[ip].C44)*stgO_tp[0].phi[ip];
    }
    stiffness_c_tp[Y] = stiffness_c_tp[X];
    stiffness_c_tp[Z] = stiffness_c_tp[X];

    stiffness_c_bk[0].C11 = 0.0;
    stiffness_c_bk[0].C12 = 0.0;
    stiffness_c_bk[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_bk[X].C11 += (stiffness_phase_n[ip].C11)*stgO_bk[0].phi[ip];
      stiffness_c_bk[X].C12 += (stiffness_phase_n[ip].C12)*stgO_bk[0].phi[ip];
      stiffness_c_bk[X].C44 += (stiffness_phase_n[ip].C44)*stgO_bk[0].phi[ip];
    }
    stiffness_c_bk[Y] = stiffness_c_bk[X];
    stiffness_c_bk[Z] = stiffness_c_bk[X];

    stiffness_c_ft[0].C11 = 0.0;
    stiffness_c_ft[0].C12 = 0.0;
    stiffness_c_ft[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_ft[X].C11 += (stiffness_phase_n[ip].C11)*stgO_ft[0].phi[ip];
      stiffness_c_ft[X].C12 += (stiffness_phase_n[ip].C12)*stgO_ft[0].phi[ip];
      stiffness_c_ft[X].C44 += (stiffness_phase_n[ip].C44)*stgO_ft[0].phi[ip];
    }
    stiffness_c_ft[Y] = stiffness_c_ft[X];
    stiffness_c_ft[Z] = stiffness_c_ft[X];

    trace_strain_front     = strain_ft[X].xx + strain_ft[X].yy + strain_ft[X].zz;
    trace_strain_back      = strain_bk[X].xx + strain_bk[X].yy + strain_bk[X].zz;
    trace_strain_right     = strain_rt[Y].xx + strain_rt[Y].yy + strain_rt[Y].zz;
    trace_strain_left      = strain_lf[Y].xx + strain_lf[Y].yy + strain_lf[Y].zz;
    trace_strain_top       = strain_tp[Z].xx + strain_tp[Z].yy + strain_tp[Z].zz;
    trace_strain_bottom    = strain_bt[Z].xx + strain_bt[Z].yy + strain_bt[Z].zz;


    lambda_front    = stiffness_c_ft[X].C12;  //C12
    lambda_back     = stiffness_c_bk[X].C12;   //C12

    mu_front        = stiffness_c_ft[X].C44;  //C44
    mu_back         = stiffness_c_bk[X].C44;   //C44

    mu_prime_front  = stiffness_c_ft[X].C11 - stiffness_c_ft[X].C12 - 2.0*stiffness_c_ft[X].C44;
    mu_prime_back   = stiffness_c_bk[X].C11 - stiffness_c_bk[X].C12 - 2.0*stiffness_c_bk[X].C44;

    lambda_right    = stiffness_c_rt[Y].C12;  //C12
    lambda_left     = stiffness_c_lf[Y].C12;   //C12

    mu_right        = stiffness_c_rt[Y].C44;  //C44
    mu_left         = stiffness_c_lf[Y].C44;   //C44

    mu_prime_right  = stiffness_c_rt[Y].C11 - stiffness_c_rt[Y].C12 - 2.0*stiffness_c_rt[Y].C44;
    mu_prime_left   = stiffness_c_lf[Y].C11 - stiffness_c_lf[Y].C12 - 2.0*stiffness_c_lf[Y].C44;

    lambda_top      = stiffness_c_tp[Z].C12;  //C12
    lambda_bottom   = stiffness_c_bt[Z].C12;   //C12

    mu_top          = stiffness_c_tp[Z].C44;  //C44
    mu_bottom       = stiffness_c_bt[Z].C44;   //C44

    mu_prime_top    = stiffness_c_tp[Z].C11 - stiffness_c_tp[Z].C12 - 2.0*stiffness_c_tp[Z].C44;
    mu_prime_bottom = stiffness_c_bt[Z].C11 - stiffness_c_bt[Z].C12 - 2.0*stiffness_c_bt[Z].C44;

    sigma_front.xx  = lambda_front*trace_strain_front  + 2.0*mu_front*strain_ft[X].xx		//Cubic
                      + mu_prime_front*strain_ft[X].xx;
    sigma_back.xx   = lambda_back*trace_strain_back  + 2.0*mu_back*strain_bk[X].xx
                      + mu_prime_back*strain_bk[X].xx;

    sigma_right.xy  = 2.0*mu_right*strain_rt[X].xy;
    sigma_left.xy   = 2.0*mu_left*strain_lf[X].xy;

    sigma_right.yz  = 2.0*mu_right*strain_rt[Z].yz;
    sigma_left.yz   = 2.0*mu_left*strain_lf[Z].yz;

    sigma_top.xz    = 2.0*mu_top*strain_tp[X].xz;
    sigma_bottom.xz = 2.0*mu_bottom*strain_bt[X].xz;

    sigma_top.yz    = 2.0*mu_top*strain_tp[Y].yz;
    sigma_bottom.yz = 2.0*mu_bottom*strain_bt[Y].yz;

    sigma_front.xy  = 2.0*mu_front*strain_ft[Y].xy;
    sigma_back.xy   = 2.0*mu_back*strain_bk[Y].xy;

    sigma_front.xz  = 2.0*mu_front*strain_ft[Z].xz;
    sigma_back.xz   = 2.0*mu_back*strain_bk[Z].xz;

    sigma_right.yy  = lambda_right*trace_strain_right   + 2.0*mu_right*strain_rt[Y].yy + mu_prime_right*strain_rt[Y].yy;

    sigma_left.yy   = lambda_left*trace_strain_left     + 2.0*mu_left*strain_lf[Y].yy + mu_prime_left*strain_lf[Y].yy;

    sigma_top.zz    = lambda_top*trace_strain_top       + 2.0*mu_top*strain_tp[Z].zz + mu_prime_top*strain_tp[Z].zz;

    sigma_bottom.zz = lambda_bottom*trace_strain_bottom + 2.0*mu_bottom*strain_bt[Z].zz + mu_prime_bottom*strain_bt[Z].zz;

    forceX          = (sigma_front.xx - sigma_back.xx)  + (sigma_right.xy - sigma_left.xy) + (sigma_top.xz   - sigma_bottom.xz);
    forceX         *= 0.5;

    forceY          = (sigma_front.xy - sigma_back.xy)  + (sigma_top.yz - sigma_bottom.yz) + (sigma_right.yy - sigma_left.yy);
    forceY         *= 0.5;

    forceZ          = (sigma_front.xz - sigma_back.xz)  + (sigma_right.yz - sigma_left.yz) + (sigma_top.zz   - sigma_bottom.zz);
    forceZ         *= 0.5;

    //printf("%d, %d, %d, %d, %le, %le, %le, %le, %le, %le, %le\n", x, y, z , index, sigma_front.xx, sigma_back.xx, sigma_right.xy, sigma_left.xy, sigma_top.xz, sigma_bottom.xz, forceX);

    it_gridinfoO[center].disp[X][0] = it_gridinfoO[center].disp[X][1];
    it_gridinfoO[center].disp[Y][0] = it_gridinfoO[center].disp[Y][1];
    it_gridinfoO[center].disp[Z][0] = it_gridinfoO[center].disp[Z][1];
    it_gridinfoO[center].disp[X][1] = it_gridinfoO[center].disp[X][2];
    it_gridinfoO[center].disp[Y][1] = it_gridinfoO[center].disp[Y][2];
    it_gridinfoO[center].disp[Z][1] = it_gridinfoO[center].disp[Z][2];

    it_gridinfoO[center].disp[X][2] = (((deltat_e*deltat_e)/rho)*forceX - (1 - damping_factor*deltat_e)*it_gridinfoO[center].disp[X][0] + 2*it_gridinfoO[center].disp[X][1])/(1.0 + damping_factor*deltat_e);
    it_gridinfoO[center].disp[Y][2] = (((deltat_e*deltat_e)/rho)*forceY - (1 - damping_factor*deltat_e)*it_gridinfoO[center].disp[Y][0] + 2*it_gridinfoO[center].disp[Y][1])/(1.0 + damping_factor*deltat_e);
    it_gridinfoO[center].disp[Z][2] = (((deltat_e*deltat_e)/rho)*forceZ - (1 - damping_factor*deltat_e)*it_gridinfoO[center].disp[Z][0] + 2*it_gridinfoO[center].disp[Z][1])/(1.0 + damping_factor*deltat_e);

    //printf("%d, %d, %d, %d, %le, %le, %le\n", x, y, z , index, damping_factor, rho, forceZ);

    //printf("%d, %d, %d, %d, %d, %le, %le, %le\n", tstep[0], x, y, z, index, it_gridinfoO[center].disp[X][2], it_gridinfoO[center].disp[Y][2], it_gridinfoO[center].disp[Z][2]);




  }
  }
  else {

  if ( ( x > 1 ) && ( x < (nx-2) ) && ( y > 1 ) && ( y < (ny-2) ) ) {
  //if ( x > 1 && x < nx-2 && y > 1 && y < ny-2 && z > 1 && z < nz-2 ) {

    index = y + ny*(z + x*nz);
    center = index;

    //printf("#%d, %d, %d, %d, %le, %le, %le\n", x, y, z , index, it_gridinfoO[center].disp[X][2], it_gridinfoO[center].disp[Y][2], it_gridinfoO[center].disp[Z][2]);

    stgO[0] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stgO[1] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stgO[2] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stgO[3] = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stgO[4] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];

    st_it_gO[0] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO[1] = it_gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO[2] = it_gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO[3] = it_gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    st_it_gO[4] = it_gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];

    eigen_strain[0].xx = 0.0;
    eigen_strain[0].yy = 0.0;
    eigen_strain[0].xy = 0.0;

    for (ip = 0; ip < npha; ip++) {

      eigen_strain[X].xx += eigen_strain_phase[ip].xx*stgO[0].phi[ip];
      eigen_strain[X].yy += eigen_strain_phase[ip].yy*stgO[0].phi[ip];
      eigen_strain[X].xy += eigen_strain_phase[ip].xy*stgO[0].phi[ip];

    }

    for (i1 = 1; i1 < 3; i1++) {

      eigen_strain[i1].xx = eigen_strain[X].xx;
      eigen_strain[i1].yy = eigen_strain[X].yy;
      eigen_strain[i1].xy = eigen_strain[X].xy;

    }
    // or
    eigen_strain[Y] = eigen_strain[X];

    strain[X].xx = 0.5*(st_it_gO[4].disp[X][2] - st_it_gO[3].disp[X][2]) - eigen_strain[X].xx;
    strain[X].yy = 0.5*(st_it_gO[2].disp[Y][2] - st_it_gO[1].disp[Y][2]) - eigen_strain[X].yy;

    strain[Y].xx = strain[X].xx;
    strain[Y].yy = strain[X].yy;

    strain[X].xy = 0.25*((st_it_gO[2].disp[X][2] - st_it_gO[1].disp[X][2]) + (st_it_gO[4].disp[Y][2] - st_it_gO[3].disp[Y][2]));

    strain[Y].xy = strain[X].xy;

    stiffness_c[0].C11 = 0.0;
    stiffness_c[0].C12 = 0.0;
    stiffness_c[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c[X].C11 += (stiffness_phase_n[ip].C11)*stgO[0].phi[ip];
      stiffness_c[X].C12 += (stiffness_phase_n[ip].C12)*stgO[0].phi[ip];
      stiffness_c[X].C44 += (stiffness_phase_n[ip].C44)*stgO[0].phi[ip];
    }
    stiffness_c[Y] = stiffness_c[X];


    stgO_lf[0] = gridinfoO[ (y  -1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_lf[1] = gridinfoO[ (y-1-1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_lf[2] = gridinfoO[ (y+1-1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_lf[3] = gridinfoO[ (y  -1) + ny*( (x-1)*nz + (z  ) ) ];
    stgO_lf[4] = gridinfoO[ (y  -1) + ny*( (x+1)*nz + (z  ) ) ];

    st_it_gO_lf[0] = it_gridinfoO[ (y  -1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_lf[1] = it_gridinfoO[ (y-1-1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_lf[2] = it_gridinfoO[ (y+1-1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_lf[3] = it_gridinfoO[ (y  -1) + ny*( (x-1)*nz + (z  ) ) ];
    st_it_gO_lf[4] = it_gridinfoO[ (y  -1) + ny*( (x+1)*nz + (z  ) ) ];




    stgO_rt[0] = gridinfoO[ (y  +1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_rt[1] = gridinfoO[ (y-1+1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_rt[2] = gridinfoO[ (y+1+1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_rt[3] = gridinfoO[ (y  +1) + ny*( (x-1)*nz + (z  ) ) ];
    stgO_rt[4] = gridinfoO[ (y  +1) + ny*( (x+1)*nz + (z  ) ) ];



    st_it_gO_rt[0] = it_gridinfoO[ (y  +1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_rt[1] = it_gridinfoO[ (y-1+1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_rt[2] = it_gridinfoO[ (y+1+1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_rt[3] = it_gridinfoO[ (y  +1) + ny*( (x-1)*nz + (z  ) ) ];
    st_it_gO_rt[4] = it_gridinfoO[ (y  +1) + ny*( (x+1)*nz + (z  ) ) ];




    stgO_bt[0] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  -1) ) ];
    stgO_bt[1] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  -1) ) ];
    stgO_bt[2] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  -1) ) ];
    stgO_bt[3] = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  -1) ) ];
    stgO_bt[4] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  -1) ) ];



    st_it_gO_bt[0] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  -1) ) ];
    st_it_gO_bt[1] = it_gridinfoO[ (y-1) + ny*( (x  )*nz + (z  -1) ) ];
    st_it_gO_bt[2] = it_gridinfoO[ (y+1) + ny*( (x  )*nz + (z  -1) ) ];
    st_it_gO_bt[3] = it_gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  -1) ) ];
    st_it_gO_bt[4] = it_gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  -1) ) ];




    stgO_tp[0] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  +1) ) ];
    stgO_tp[1] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  +1) ) ];
    stgO_tp[2] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  +1) ) ];
    stgO_tp[3] = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  +1) ) ];
    stgO_tp[4] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  +1) ) ];



    st_it_gO_tp[0] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  +1) ) ];
    st_it_gO_tp[1] = it_gridinfoO[ (y-1) + ny*( (x  )*nz + (z  +1) ) ];
    st_it_gO_tp[2] = it_gridinfoO[ (y+1) + ny*( (x  )*nz + (z  +1) ) ];
    st_it_gO_tp[3] = it_gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  +1) ) ];
    st_it_gO_tp[4] = it_gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  +1) ) ];




    stgO_bk[0] = gridinfoO[ (y  ) + ny*( (x  -1)*nz + (z  ) ) ];
    stgO_bk[1] = gridinfoO[ (y-1) + ny*( (x  -1)*nz + (z  ) ) ];
    stgO_bk[2] = gridinfoO[ (y+1) + ny*( (x  -1)*nz + (z  ) ) ];
    stgO_bk[3] = gridinfoO[ (y  ) + ny*( (x-1-1)*nz + (z  ) ) ];
    stgO_bk[4] = gridinfoO[ (y  ) + ny*( (x+1-1)*nz + (z  ) ) ];



    st_it_gO_bk[0] = it_gridinfoO[ (y  ) + ny*( (x  -1)*nz + (z  ) ) ];
    st_it_gO_bk[1] = it_gridinfoO[ (y-1) + ny*( (x  -1)*nz + (z  ) ) ];
    st_it_gO_bk[2] = it_gridinfoO[ (y+1) + ny*( (x  -1)*nz + (z  ) ) ];
    st_it_gO_bk[3] = it_gridinfoO[ (y  ) + ny*( (x-1-1)*nz + (z  ) ) ];
    st_it_gO_bk[4] = it_gridinfoO[ (y  ) + ny*( (x+1-1)*nz + (z  ) ) ];




    stgO_ft[0] = gridinfoO[ (y  ) + ny*( (x  +1)*nz + (z  ) ) ];
    stgO_ft[1] = gridinfoO[ (y-1) + ny*( (x  +1)*nz + (z  ) ) ];
    stgO_ft[2] = gridinfoO[ (y+1) + ny*( (x  +1)*nz + (z  ) ) ];
    stgO_ft[3] = gridinfoO[ (y  ) + ny*( (x-1+1)*nz + (z  ) ) ];
    stgO_ft[4] = gridinfoO[ (y  ) + ny*( (x+1+1)*nz + (z  ) ) ];



    st_it_gO_ft[0] = it_gridinfoO[ (y  ) + ny*( (x  +1)*nz + (z  ) ) ];
    st_it_gO_ft[1] = it_gridinfoO[ (y-1) + ny*( (x  +1)*nz + (z  ) ) ];
    st_it_gO_ft[2] = it_gridinfoO[ (y+1) + ny*( (x  +1)*nz + (z  ) ) ];
    st_it_gO_ft[3] = it_gridinfoO[ (y  ) + ny*( (x-1+1)*nz + (z  ) ) ];
    st_it_gO_ft[4] = it_gridinfoO[ (y  ) + ny*( (x+1+1)*nz + (z  ) ) ];



    eigen_strain_lf[0].xx = 0.0;
    eigen_strain_lf[0].yy = 0.0;
    eigen_strain_lf[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_lf[X].xx += eigen_strain_phase[ip].xx*stgO_lf[0].phi[ip];
      eigen_strain_lf[X].yy += eigen_strain_phase[ip].yy*stgO_lf[0].phi[ip];
      eigen_strain_lf[X].xy += eigen_strain_phase[ip].xy*stgO_lf[0].phi[ip];
    }
    eigen_strain_lf[Y] = eigen_strain_lf[X];

    eigen_strain_rt[0].xx = 0.0;
    eigen_strain_rt[0].yy = 0.0;
    eigen_strain_rt[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_rt[X].xx += eigen_strain_phase[ip].xx*stgO_rt[0].phi[ip];
      eigen_strain_rt[X].yy += eigen_strain_phase[ip].yy*stgO_rt[0].phi[ip];
      eigen_strain_rt[X].xy += eigen_strain_phase[ip].xy*stgO_rt[0].phi[ip];
    }
    eigen_strain_rt[Y] = eigen_strain_rt[X];

    eigen_strain_bt[0].xx = 0.0;
    eigen_strain_bt[0].yy = 0.0;
    eigen_strain_bt[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_bt[X].xx += eigen_strain_phase[ip].xx*stgO_bt[0].phi[ip];
      eigen_strain_bt[X].yy += eigen_strain_phase[ip].yy*stgO_bt[0].phi[ip];
      eigen_strain_bt[X].xy += eigen_strain_phase[ip].xy*stgO_bt[0].phi[ip];
    }
    eigen_strain_bt[Y] = eigen_strain_bt[X];

    eigen_strain_tp[0].xx = 0.0;
    eigen_strain_tp[0].yy = 0.0;
    eigen_strain_tp[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_tp[X].xx += eigen_strain_phase[ip].xx*stgO_tp[0].phi[ip];
      eigen_strain_tp[X].yy += eigen_strain_phase[ip].yy*stgO_tp[0].phi[ip];
      eigen_strain_tp[X].xy += eigen_strain_phase[ip].xy*stgO_tp[0].phi[ip];
    }
    eigen_strain_tp[Y] = eigen_strain_tp[X];

    eigen_strain_bk[0].xx = 0.0;
    eigen_strain_bk[0].yy = 0.0;
    eigen_strain_bk[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_bk[X].xx += eigen_strain_phase[ip].xx*stgO_bk[0].phi[ip];
      eigen_strain_bk[X].yy += eigen_strain_phase[ip].yy*stgO_bk[0].phi[ip];
      eigen_strain_bk[X].xy += eigen_strain_phase[ip].xy*stgO_bk[0].phi[ip];
    }
    eigen_strain_bk[Y] = eigen_strain_bk[X];

    eigen_strain_ft[0].xx = 0.0;
    eigen_strain_ft[0].yy = 0.0;
    eigen_strain_ft[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_ft[X].xx += eigen_strain_phase[ip].xx*stgO_ft[0].phi[ip];
      eigen_strain_ft[X].yy += eigen_strain_phase[ip].yy*stgO_ft[0].phi[ip];
      eigen_strain_ft[X].xy += eigen_strain_phase[ip].xy*stgO_ft[0].phi[ip];
    }
    eigen_strain_ft[Y] = eigen_strain_ft[X];

    strain_lf[X].xx = 0.5*(st_it_gO_lf[4].disp[X][2] - st_it_gO_lf[3].disp[X][2]) - eigen_strain_lf[X].xx;
    strain_lf[X].yy = 0.5*(st_it_gO_lf[2].disp[Y][2] - st_it_gO_lf[1].disp[Y][2]) - eigen_strain_lf[X].yy;
    strain_lf[Y].xx = strain_lf[X].xx;
    strain_lf[Y].yy = strain_lf[X].yy;

    strain_rt[X].xx = 0.5*(st_it_gO_rt[4].disp[X][2] - st_it_gO_rt[3].disp[X][2]) - eigen_strain_rt[X].xx;
    strain_rt[X].yy = 0.5*(st_it_gO_rt[2].disp[Y][2] - st_it_gO_rt[1].disp[Y][2]) - eigen_strain_rt[X].yy;
    strain_rt[Y].xx = strain_rt[X].xx;
    strain_rt[Y].yy = strain_rt[X].yy;

    strain_bt[X].xx = 0.5*(st_it_gO_bt[4].disp[X][2] - st_it_gO_bt[3].disp[X][2]) - eigen_strain_bt[X].xx;
    strain_bt[X].yy = 0.5*(st_it_gO_bt[2].disp[Y][2] - st_it_gO_bt[1].disp[Y][2]) - eigen_strain_bt[X].yy;
    strain_bt[Y].xx = strain_bt[X].xx;
    strain_bt[Y].yy = strain_bt[X].yy;

    strain_tp[X].xx = 0.5*(st_it_gO_tp[4].disp[X][2] - st_it_gO_tp[3].disp[X][2]) - eigen_strain_tp[X].xx;
    strain_tp[X].yy = 0.5*(st_it_gO_tp[2].disp[Y][2] - st_it_gO_tp[1].disp[Y][2]) - eigen_strain_tp[X].yy;
    strain_tp[Y].xx = strain_tp[X].xx;
    strain_tp[Y].yy = strain_tp[X].yy;

    strain_bk[X].xx = 0.5*(st_it_gO_bk[4].disp[X][2] - st_it_gO_bk[3].disp[X][2]) - eigen_strain_bk[X].xx;
    strain_bk[X].yy = 0.5*(st_it_gO_bk[2].disp[Y][2] - st_it_gO_bk[1].disp[Y][2]) - eigen_strain_bk[X].yy;
    strain_bk[Y].xx = strain_bk[X].xx;
    strain_bk[Y].yy = strain_bk[X].yy;

    strain_ft[X].xx = 0.5*(st_it_gO_ft[4].disp[X][2] - st_it_gO_ft[3].disp[X][2]) - eigen_strain_ft[X].xx;
    strain_ft[X].yy = 0.5*(st_it_gO_ft[2].disp[Y][2] - st_it_gO_ft[1].disp[Y][2]) - eigen_strain_ft[X].yy;
    strain_ft[Y].xx = strain_ft[X].xx;
    strain_ft[Y].yy = strain_ft[X].yy;

    strain_lf[X].xy = 0.25*((st_it_gO_lf[2].disp[X][2] - st_it_gO_lf[1].disp[X][2]) + (st_it_gO_lf[4].disp[Y][2] - st_it_gO_lf[3].disp[Y][2]));
    strain_lf[Y].xy = strain_lf[X].xy;

    strain_rt[X].xy = 0.25*((st_it_gO_rt[2].disp[X][2] - st_it_gO_rt[1].disp[X][2]) + (st_it_gO_rt[4].disp[Y][2] - st_it_gO_rt[3].disp[Y][2]));
    strain_rt[Y].xy = strain_rt[X].xy;

    strain_bt[X].xy = 0.25*((st_it_gO_bt[2].disp[X][2] - st_it_gO_bt[1].disp[X][2]) + (st_it_gO_bt[4].disp[Y][2] - st_it_gO_bt[3].disp[Y][2]));
    strain_bt[Y].xy = strain_bt[X].xy;

    strain_tp[X].xy = 0.25*((st_it_gO_tp[2].disp[X][2] - st_it_gO_tp[1].disp[X][2]) + (st_it_gO_tp[4].disp[Y][2] - st_it_gO_tp[3].disp[Y][2]));
    strain_tp[Y].xy = strain_tp[X].xy;

    strain_bk[X].xy = 0.25*((st_it_gO_bk[2].disp[X][2] - st_it_gO_bk[1].disp[X][2]) + (st_it_gO_bk[4].disp[Y][2] - st_it_gO_bk[3].disp[Y][2]));
    strain_bk[Y].xy = strain_bk[X].xy;

    strain_ft[X].xy = 0.25*((st_it_gO_ft[2].disp[X][2] - st_it_gO_ft[1].disp[X][2]) + (st_it_gO_ft[4].disp[Y][2] - st_it_gO_ft[3].disp[Y][2]));
    strain_ft[Y].xy = strain_ft[X].xy;
    strain_ft[Z].xy = strain_ft[X].xy;

    stiffness_c_lf[0].C11 = 0.0;
    stiffness_c_lf[0].C12 = 0.0;
    stiffness_c_lf[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_lf[X].C11 += (stiffness_phase_n[ip].C11)*stgO_lf[0].phi[ip];
      stiffness_c_lf[X].C12 += (stiffness_phase_n[ip].C12)*stgO_lf[0].phi[ip];
      stiffness_c_lf[X].C44 += (stiffness_phase_n[ip].C44)*stgO_lf[0].phi[ip];
    }
    stiffness_c_lf[Y] = stiffness_c_lf[X];

    stiffness_c_rt[0].C11 = 0.0;
    stiffness_c_rt[0].C12 = 0.0;
    stiffness_c_rt[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_rt[X].C11 += (stiffness_phase_n[ip].C11)*stgO_rt[0].phi[ip];
      stiffness_c_rt[X].C12 += (stiffness_phase_n[ip].C12)*stgO_rt[0].phi[ip];
      stiffness_c_rt[X].C44 += (stiffness_phase_n[ip].C44)*stgO_rt[0].phi[ip];
    }
    stiffness_c_rt[Y] = stiffness_c_rt[X];

    stiffness_c_bt[0].C11 = 0.0;
    stiffness_c_bt[0].C12 = 0.0;
    stiffness_c_bt[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_bt[X].C11 += (stiffness_phase_n[ip].C11)*stgO_bt[0].phi[ip];
      stiffness_c_bt[X].C12 += (stiffness_phase_n[ip].C12)*stgO_bt[0].phi[ip];
      stiffness_c_bt[X].C44 += (stiffness_phase_n[ip].C44)*stgO_bt[0].phi[ip];
    }
    stiffness_c_bt[Y] = stiffness_c_bt[X];

    stiffness_c_tp[0].C11 = 0.0;
    stiffness_c_tp[0].C12 = 0.0;
    stiffness_c_tp[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_tp[X].C11 += (stiffness_phase_n[ip].C11)*stgO_tp[0].phi[ip];
      stiffness_c_tp[X].C12 += (stiffness_phase_n[ip].C12)*stgO_tp[0].phi[ip];
      stiffness_c_tp[X].C44 += (stiffness_phase_n[ip].C44)*stgO_tp[0].phi[ip];
    }
    stiffness_c_tp[Y] = stiffness_c_tp[X];

    stiffness_c_bk[0].C11 = 0.0;
    stiffness_c_bk[0].C12 = 0.0;
    stiffness_c_bk[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_bk[X].C11 += (stiffness_phase_n[ip].C11)*stgO_bk[0].phi[ip];
      stiffness_c_bk[X].C12 += (stiffness_phase_n[ip].C12)*stgO_bk[0].phi[ip];
      stiffness_c_bk[X].C44 += (stiffness_phase_n[ip].C44)*stgO_bk[0].phi[ip];
    }
    stiffness_c_bk[Y] = stiffness_c_bk[X];

    stiffness_c_ft[0].C11 = 0.0;
    stiffness_c_ft[0].C12 = 0.0;
    stiffness_c_ft[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_ft[X].C11 += (stiffness_phase_n[ip].C11)*stgO_ft[0].phi[ip];
      stiffness_c_ft[X].C12 += (stiffness_phase_n[ip].C12)*stgO_ft[0].phi[ip];
      stiffness_c_ft[X].C44 += (stiffness_phase_n[ip].C44)*stgO_ft[0].phi[ip];
    }
    stiffness_c_ft[Y] = stiffness_c_ft[X];

    trace_strain_front     = strain_ft[X].xx + strain_ft[X].yy;// + strain_ft[X].zz;
    trace_strain_back      = strain_bk[X].xx + strain_bk[X].yy;// + strain_bk[X].zz;
    trace_strain_right     = strain_rt[Y].xx + strain_rt[Y].yy;// + strain_rt[Y].zz;
    trace_strain_left      = strain_lf[Y].xx + strain_lf[Y].yy;// + strain_lf[Y].zz;
    trace_strain_top       = 0.0;//strain_tp[Z].xx + strain_tp[Z].yy;// + strain_tp[Z].zz;
    trace_strain_bottom    = 0.0;//strain_bt[Z].xx + strain_bt[Z].yy;// + strain_bt[Z].zz;


    lambda_front    = stiffness_c_ft[X].C12;  //C12
    lambda_back     = stiffness_c_bk[X].C12;   //C12

    mu_front        = stiffness_c_ft[X].C44;  //C44
    mu_back         = stiffness_c_bk[X].C44;   //C44

    mu_prime_front  = stiffness_c_ft[X].C11 - stiffness_c_ft[X].C12 - 2.0*stiffness_c_ft[X].C44;
    mu_prime_back   = stiffness_c_bk[X].C11 - stiffness_c_bk[X].C12 - 2.0*stiffness_c_bk[X].C44;

    lambda_right    = stiffness_c_rt[Y].C12;  //C12
    lambda_left     = stiffness_c_lf[Y].C12;   //C12

    mu_right        = stiffness_c_rt[Y].C44;  //C44
    mu_left         = stiffness_c_lf[Y].C44;   //C44

    mu_prime_right  = stiffness_c_rt[Y].C11 - stiffness_c_rt[Y].C12 - 2.0*stiffness_c_rt[Y].C44;
    mu_prime_left   = stiffness_c_lf[Y].C11 - stiffness_c_lf[Y].C12 - 2.0*stiffness_c_lf[Y].C44;

    lambda_top      = 0.0;//stiffness_c_tp[Z].C12;  //C12
    lambda_bottom   = 0.0;//stiffness_c_bt[Z].C12;   //C12

    mu_top          = 0.0;//stiffness_c_tp[Z].C44;  //C44
    mu_bottom       = 0.0;//stiffness_c_bt[Z].C44;   //C44

    mu_prime_top    = 0.0;//stiffness_c_tp[Z].C11 - stiffness_c_tp[Z].C12 - 2.0*stiffness_c_tp[Z].C44;
    mu_prime_bottom = 0.0;//stiffness_c_bt[Z].C11 - stiffness_c_bt[Z].C12 - 2.0*stiffness_c_bt[Z].C44;

    sigma_front.xx  = lambda_front*trace_strain_front  + 2.0*mu_front*strain_ft[X].xx		//Cubic
                      + mu_prime_front*strain_ft[X].xx;
    sigma_back.xx   = lambda_back*trace_strain_back  + 2.0*mu_back*strain_bk[X].xx
                      + mu_prime_back*strain_bk[X].xx;

    sigma_right.xy  = 2.0*mu_right*strain_rt[X].xy;
    sigma_left.xy   = 2.0*mu_left*strain_lf[X].xy;

    sigma_right.yz  = 0.0;//2.0*mu_right*strain_rt[Z].yz;
    sigma_left.yz   = 0.0;//2.0*mu_left*strain_lf[Z].yz;

    sigma_top.xz    = 0.0;//2.0*mu_top*strain_tp[X].xz;
    sigma_bottom.xz = 0.0;//2.0*mu_bottom*strain_bt[X].xz;

    sigma_top.yz    = 0.0;//2.0*mu_top*strain_tp[Y].yz;
    sigma_bottom.yz = 0.0;//2.0*mu_bottom*strain_bt[Y].yz;

    sigma_front.xy  = 2.0*mu_front*strain_ft[Y].xy;
    sigma_back.xy   = 2.0*mu_back*strain_bk[Y].xy;

    sigma_front.xz  = 0.0;//2.0*mu_front*strain_ft[Z].xz;
    sigma_back.xz   = 0.0;//2.0*mu_back*strain_bk[Z].xz;

    sigma_right.yy  = lambda_right*trace_strain_right   + 2.0*mu_right*strain_rt[Y].yy + mu_prime_right*strain_rt[Y].yy;

    sigma_left.yy   = lambda_left*trace_strain_left     + 2.0*mu_left*strain_lf[Y].yy + mu_prime_left*strain_lf[Y].yy;

    sigma_top.zz    = 0.0;//lambda_top*trace_strain_top       + 2.0*mu_top*strain_tp[Z].zz + mu_prime_top*strain_tp[Z].zz;

    sigma_bottom.zz = 0.0;//lambda_bottom*trace_strain_bottom + 2.0*mu_bottom*strain_bt[Z].zz + mu_prime_bottom*strain_bt[Z].zz;

    forceX          = (sigma_front.xx - sigma_back.xx)  + (sigma_right.xy - sigma_left.xy) + (sigma_top.xz   - sigma_bottom.xz);
    forceX         *= 0.5;

    forceY          = (sigma_front.xy - sigma_back.xy)  + (sigma_top.yz - sigma_bottom.yz) + (sigma_right.yy - sigma_left.yy);
    forceY         *= 0.5;

    forceZ          = 0.0;//(sigma_front.xz - sigma_back.xz)  + (sigma_right.yz - sigma_left.yz) + (sigma_top.zz   - sigma_bottom.zz);
    forceZ         *= 0.0;//0.5;

    //printf("%d, %d, %d, %d, %le, %le, %le, %le, %le, %le, %le\n", x, y, z , index, sigma_front.xx, sigma_back.xx, sigma_right.xy, sigma_left.xy, sigma_top.xz, sigma_bottom.xz, forceX);

    it_gridinfoO[center].disp[X][0] = it_gridinfoO[center].disp[X][1];
    it_gridinfoO[center].disp[Y][0] = it_gridinfoO[center].disp[Y][1];
    it_gridinfoO[center].disp[Z][0] = 0.0;//it_gridinfoO[center].disp[Z][1];
    it_gridinfoO[center].disp[X][1] = it_gridinfoO[center].disp[X][2];
    it_gridinfoO[center].disp[Y][1] = it_gridinfoO[center].disp[Y][2];
    it_gridinfoO[center].disp[Z][1] = 0.0;//it_gridinfoO[center].disp[Z][2];

    it_gridinfoO[center].disp[X][2] = (((deltat_e*deltat_e)/rho)*forceX - (1 - damping_factor*deltat_e)*it_gridinfoO[center].disp[X][0] + 2*it_gridinfoO[center].disp[X][1])/(1.0 + damping_factor*deltat_e);
    it_gridinfoO[center].disp[Y][2] = (((deltat_e*deltat_e)/rho)*forceY - (1 - damping_factor*deltat_e)*it_gridinfoO[center].disp[Y][0] + 2*it_gridinfoO[center].disp[Y][1])/(1.0 + damping_factor*deltat_e);
    it_gridinfoO[center].disp[Z][2] = 0.0;//(((deltat_e*deltat_e)/rho)*forceZ - (1 - damping_factor*deltat_e)*it_gridinfoO[center].disp[Z][0] + 2*it_gridinfoO[center].disp[Z][1])/(1.0 + damping_factor*deltat_e);

    //printf("%d, %d, %d, %d, %le, %le, %le\n", x, y, z , index, damping_factor, rho, forceZ);

    //printf("%d, %d, %d, %d, %d, %le, %le, %le\n", tstep[0], x, y, z, index, it_gridinfoO[center].disp[X][2], it_gridinfoO[center].disp[Y][2], it_gridinfoO[center].disp[Z][2]);




  }
  }


}



__kernel void SolverStress_iterative_2D(__global struct fields *gridinfoO, __global struct iter_variables *it_gridinfoO, __global struct symmetric_tensor *eigen_strain_phase, __global struct Stiffness_cubic *stiffness_phase, __global struct Stiffness_cubic *stiffness_phase_n, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global int *tstep) {

  int x, y, z, ii, i1, j1, k1, i2;
  int nx, ny, nz;
  int index;
  int X, Y;
  int center;

  int is, js, ks, il, il1, jl, jl1, ipx;
  int is1, is2;
  int ig, jg, ip, ip1, ip2, ip3;

  double deltat_e, damping_factor,  rho;
 double trace_strain_right;
 double trace_strain_left;
 double trace_strain_back;
 double trace_strain_front;

  double lambda_front   ;
  double lambda_back    ;
  double mu_front       ;
  double mu_back        ;
  double mu_prime_front ;
  double mu_prime_back  ;
  double lambda_right   ;
  double lambda_left    ;
  double mu_right       ;
  double mu_left        ;
  double mu_prime_right ;
  double mu_prime_left

  double sigma_xx_front;
  double sigma_xx_back;
  double sigma_yx_right;
  double sigma_yx_left;
  double sigma_xy_front;
  double sigma_xy_back;
  double sigma_yy_right;
  double sigma_yy_left;

 double forceX;
 double forceY;

 double div_phi_front;
 double div_phi_back;
 double div_phi_right;
 double div_phi_left;


  struct symmetric_tensor eigen_strain[3], eigen_strain_lf[3], eigen_strain_rt[3], eigen_strain_bt[3], eigen_strain_tp[3], eigen_strain_bk[3], eigen_strain_ft[3];

  //struct symmetric_tensor eigen_strain_phase[npha];

  struct symmetric_tensor strain[3], strain_lf[3], strain_rt[3], strain_bt[3], strain_tp[3], strain_bk[3], strain_ft[3];

  struct Stiffness_cubic stiffness_c[3], stiffness_c_lf[3], stiffness_c_rt[3], stiffness_c_bt[3], stiffness_c_tp[3], stiffness_c_bk[3], stiffness_c_ft[3];

  //struct Stiffness_cubic stiffness_phase_n[npha];

 struct symmetric_tensor sigma_front;
 struct symmetric_tensor sigma_back;
 struct symmetric_tensor sigma_right;
 struct symmetric_tensor sigma_left;
 struct symmetric_tensor sigma_top;
 struct symmetric_tensor sigma_bottom;

  struct fields stgO[5], stgO_lf[5], stgO_rt[5], stgO_bt[5], stgO_tp[5], stgO_bk[5], stgO_ft[5];
  struct iter_variables st_it_gO[5], st_it_gO_lf[5], st_it_gO_rt[5], st_it_gO_bt[5], st_it_gO_tp[5], st_it_gO_bk[5], st_it_gO_ft[5];

  X = 0;
  Y = 1;

  deltat_e = pfmdat->deltat_e;
  rho = pfmdat->rho;
  damping_factor = pfmdat->damping_factor;

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  //index = y + ny*(z + x*nz);

  if ( nz > 4 ) {
  if ( ( x > 1 ) && ( x < (nx-2) ) && ( y > 1 ) && ( y < (ny-2) ) && ( z > 1 ) && ( z < (nz-2) ) ) {
  //if ( x > 1 && x < nx-2 && y > 1 && y < ny-2 && z > 1 && z < nz-2 ) {

    index = y + ny*(z + x*nz);
    center = index;

    printf("%d, %d, %d, %d, %le, %le, %le\n", x, y, z , index, it_gridinfoO[center].disp[X][2], it_gridinfoO[center].disp[Y][2], it_gridinfoO[center].disp[Z][2]);

    stgO[0] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stgO[1] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stgO[2] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stgO[3] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    stgO[4] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    stgO[5] = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stgO[6] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];

    st_it_gO[0] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO[1] = it_gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO[2] = it_gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO[3] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1) ) ];
    st_it_gO[4] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1) ) ];
    st_it_gO[5] = it_gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    st_it_gO[6] = it_gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];

    eigen_strain[0].xx = 0.0;
    eigen_strain[0].yy = 0.0;
    eigen_strain[0].zz = 0.0;
    eigen_strain[0].yz = 0.0;
    eigen_strain[0].xz = 0.0;
    eigen_strain[0].xy = 0.0;

    for (ip = 0; ip < npha; ip++) {

      eigen_strain[X].xx += eigen_strain_phase[ip].xx*stgO[0].phi[ip];
      eigen_strain[X].yy += eigen_strain_phase[ip].yy*stgO[0].phi[ip];
      eigen_strain[X].zz += eigen_strain_phase[ip].zz*stgO[0].phi[ip];
      eigen_strain[X].yz += eigen_strain_phase[ip].yz*stgO[0].phi[ip];
      eigen_strain[X].xz += eigen_strain_phase[ip].xz*stgO[0].phi[ip];
      eigen_strain[X].xy += eigen_strain_phase[ip].xy*stgO[0].phi[ip];

    }

    for (i1 = 1; i1 < 3; i1++) {

      eigen_strain[i1].xx = eigen_strain[X].xx;
      eigen_strain[i1].yy = eigen_strain[X].yy;
      eigen_strain[i1].zz = eigen_strain[X].zz;
      eigen_strain[i1].yz = eigen_strain[X].yz;
      eigen_strain[i1].xz = eigen_strain[X].xz;
      eigen_strain[i1].xy = eigen_strain[X].xy;

    }
    // or
    eigen_strain[Y] = eigen_strain[X];
    eigen_strain[Z] = eigen_strain[X];

    strain[X].xx = 0.5*(st_it_gO[6].disp[X][2] - st_it_gO[5].disp[X][2]) - eigen_strain[X].xx;
    strain[X].yy = 0.5*(st_it_gO[2].disp[Y][2] - st_it_gO[1].disp[Y][2]) - eigen_strain[X].yy;
    strain[X].zz = 0.5*(st_it_gO[4].disp[Z][2] - st_it_gO[3].disp[Z][2]) - eigen_strain[X].zz;

    strain[Y].xx = strain[X].xx;
    strain[Y].yy = strain[X].yy;
    strain[Y].zz = strain[X].zz;
    strain[Z].yy = strain[X].yy;
    strain[Z].xx = strain[X].xx;
    strain[Z].zz = strain[X].zz;

    strain[X].xy = 0.25*((st_it_gO[2].disp[X][2] - st_it_gO[1].disp[X][2]) + (st_it_gO[6].disp[Y][2] - st_it_gO[5].disp[Y][2]));
    strain[X].xz = 0.25*((st_it_gO[4].disp[X][2] - st_it_gO[3].disp[X][2]) + (st_it_gO[6].disp[Z][2] - st_it_gO[5].disp[Z][2]));
    strain[X].yz = 0.25*((st_it_gO[2].disp[Z][2] - st_it_gO[1].disp[Z][2]) + (st_it_gO[4].disp[Y][2] - st_it_gO[3].disp[Y][2]));

    strain[Y].xy = strain[X].xy;
    strain[Y].xz = strain[X].xz;
    strain[Y].yz = strain[X].yz;
    strain[Z].xy = strain[X].xy;
    strain[Z].xz = strain[X].xz;
    strain[Z].yz = strain[X].yz;

    stiffness_c[0].C11 = 0.0;
    stiffness_c[0].C12 = 0.0;
    stiffness_c[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c[X].C11 += (stiffness_phase_n[ip].C11)*stgO[0].phi[ip];
      stiffness_c[X].C12 += (stiffness_phase_n[ip].C12)*stgO[0].phi[ip];
      stiffness_c[X].C44 += (stiffness_phase_n[ip].C44)*stgO[0].phi[ip];
    }
    stiffness_c[Y] = stiffness_c[X];
    stiffness_c[Z] = stiffness_c[X];


    stgO_lf[0] = gridinfoO[ (y  -1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_lf[1] = gridinfoO[ (y-1-1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_lf[2] = gridinfoO[ (y+1-1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_lf[3] = gridinfoO[ (y  -1) + ny*( (x  )*nz + (z-1) ) ];
    stgO_lf[4] = gridinfoO[ (y  -1) + ny*( (x  )*nz + (z+1) ) ];
    stgO_lf[5] = gridinfoO[ (y  -1) + ny*( (x-1)*nz + (z  ) ) ];
    stgO_lf[6] = gridinfoO[ (y  -1) + ny*( (x+1)*nz + (z  ) ) ];

    st_it_gO_lf[0] = it_gridinfoO[ (y  -1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_lf[1] = it_gridinfoO[ (y-1-1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_lf[2] = it_gridinfoO[ (y+1-1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_lf[3] = it_gridinfoO[ (y  -1) + ny*( (x  )*nz + (z-1) ) ];
    st_it_gO_lf[4] = it_gridinfoO[ (y  -1) + ny*( (x  )*nz + (z+1) ) ];
    st_it_gO_lf[5] = it_gridinfoO[ (y  -1) + ny*( (x-1)*nz + (z  ) ) ];
    st_it_gO_lf[6] = it_gridinfoO[ (y  -1) + ny*( (x+1)*nz + (z  ) ) ];


    stgO_rt[0] = gridinfoO[ (y  +1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_rt[1] = gridinfoO[ (y-1+1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_rt[2] = gridinfoO[ (y+1+1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_rt[3] = gridinfoO[ (y  +1) + ny*( (x  )*nz + (z-1) ) ];
    stgO_rt[4] = gridinfoO[ (y  +1) + ny*( (x  )*nz + (z+1) ) ];
    stgO_rt[5] = gridinfoO[ (y  +1) + ny*( (x-1)*nz + (z  ) ) ];
    stgO_rt[6] = gridinfoO[ (y  +1) + ny*( (x+1)*nz + (z  ) ) ];

    st_it_gO_rt[0] = it_gridinfoO[ (y  +1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_rt[1] = it_gridinfoO[ (y-1+1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_rt[2] = it_gridinfoO[ (y+1+1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_rt[3] = it_gridinfoO[ (y  +1) + ny*( (x  )*nz + (z-1) ) ];
    st_it_gO_rt[4] = it_gridinfoO[ (y  +1) + ny*( (x  )*nz + (z+1) ) ];
    st_it_gO_rt[5] = it_gridinfoO[ (y  +1) + ny*( (x-1)*nz + (z  ) ) ];
    st_it_gO_rt[6] = it_gridinfoO[ (y  +1) + ny*( (x+1)*nz + (z  ) ) ];


    stgO_bt[0] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  -1) ) ];
    stgO_bt[1] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  -1) ) ];
    stgO_bt[2] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  -1) ) ];
    stgO_bt[3] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1-1) ) ];
    stgO_bt[4] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1-1) ) ];
    stgO_bt[5] = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  -1) ) ];
    stgO_bt[6] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  -1) ) ];

    st_it_gO_bt[0] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  -1) ) ];
    st_it_gO_bt[1] = it_gridinfoO[ (y-1) + ny*( (x  )*nz + (z  -1) ) ];
    st_it_gO_bt[2] = it_gridinfoO[ (y+1) + ny*( (x  )*nz + (z  -1) ) ];
    st_it_gO_bt[3] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1-1) ) ];
    st_it_gO_bt[4] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1-1) ) ];
    st_it_gO_bt[5] = it_gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  -1) ) ];
    st_it_gO_bt[6] = it_gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  -1) ) ];


    stgO_tp[0] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  +1) ) ];
    stgO_tp[1] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  +1) ) ];
    stgO_tp[2] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  +1) ) ];
    stgO_tp[3] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1+1) ) ];
    stgO_tp[4] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1+1) ) ];
    stgO_tp[5] = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  +1) ) ];
    stgO_tp[6] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  +1) ) ];

    st_it_gO_tp[0] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  +1) ) ];
    st_it_gO_tp[1] = it_gridinfoO[ (y-1) + ny*( (x  )*nz + (z  +1) ) ];
    st_it_gO_tp[2] = it_gridinfoO[ (y+1) + ny*( (x  )*nz + (z  +1) ) ];
    st_it_gO_tp[3] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z-1+1) ) ];
    st_it_gO_tp[4] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z+1+1) ) ];
    st_it_gO_tp[5] = it_gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  +1) ) ];
    st_it_gO_tp[6] = it_gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  +1) ) ];


    stgO_bk[0] = gridinfoO[ (y  ) + ny*( (x  -1)*nz + (z  ) ) ];
    stgO_bk[1] = gridinfoO[ (y-1) + ny*( (x  -1)*nz + (z  ) ) ];
    stgO_bk[2] = gridinfoO[ (y+1) + ny*( (x  -1)*nz + (z  ) ) ];
    stgO_bk[3] = gridinfoO[ (y  ) + ny*( (x  -1)*nz + (z-1) ) ];
    stgO_bk[4] = gridinfoO[ (y  ) + ny*( (x  -1)*nz + (z+1) ) ];
    stgO_bk[5] = gridinfoO[ (y  ) + ny*( (x-1-1)*nz + (z  ) ) ];
    stgO_bk[6] = gridinfoO[ (y  ) + ny*( (x+1-1)*nz + (z  ) ) ];

    st_it_gO_bk[0] = it_gridinfoO[ (y  ) + ny*( (x  -1)*nz + (z  ) ) ];
    st_it_gO_bk[1] = it_gridinfoO[ (y-1) + ny*( (x  -1)*nz + (z  ) ) ];
    st_it_gO_bk[2] = it_gridinfoO[ (y+1) + ny*( (x  -1)*nz + (z  ) ) ];
    st_it_gO_bk[3] = it_gridinfoO[ (y  ) + ny*( (x  -1)*nz + (z-1) ) ];
    st_it_gO_bk[4] = it_gridinfoO[ (y  ) + ny*( (x  -1)*nz + (z+1) ) ];
    st_it_gO_bk[5] = it_gridinfoO[ (y  ) + ny*( (x-1-1)*nz + (z  ) ) ];
    st_it_gO_bk[6] = it_gridinfoO[ (y  ) + ny*( (x+1-1)*nz + (z  ) ) ];


    stgO_ft[0] = gridinfoO[ (y  ) + ny*( (x  +1)*nz + (z  ) ) ];
    stgO_ft[1] = gridinfoO[ (y-1) + ny*( (x  +1)*nz + (z  ) ) ];
    stgO_ft[2] = gridinfoO[ (y+1) + ny*( (x  +1)*nz + (z  ) ) ];
    stgO_ft[3] = gridinfoO[ (y  ) + ny*( (x  +1)*nz + (z-1) ) ];
    stgO_ft[4] = gridinfoO[ (y  ) + ny*( (x  +1)*nz + (z+1) ) ];
    stgO_ft[5] = gridinfoO[ (y  ) + ny*( (x-1+1)*nz + (z  ) ) ];
    stgO_ft[6] = gridinfoO[ (y  ) + ny*( (x+1+1)*nz + (z  ) ) ];

    st_it_gO_ft[0] = it_gridinfoO[ (y  ) + ny*( (x  +1)*nz + (z  ) ) ];
    st_it_gO_ft[1] = it_gridinfoO[ (y-1) + ny*( (x  +1)*nz + (z  ) ) ];
    st_it_gO_ft[2] = it_gridinfoO[ (y+1) + ny*( (x  +1)*nz + (z  ) ) ];
    st_it_gO_ft[3] = it_gridinfoO[ (y  ) + ny*( (x  +1)*nz + (z-1) ) ];
    st_it_gO_ft[4] = it_gridinfoO[ (y  ) + ny*( (x  +1)*nz + (z+1) ) ];
    st_it_gO_ft[5] = it_gridinfoO[ (y  ) + ny*( (x-1+1)*nz + (z  ) ) ];
    st_it_gO_ft[6] = it_gridinfoO[ (y  ) + ny*( (x+1+1)*nz + (z  ) ) ];

    eigen_strain_lf[0].xx = 0.0;
    eigen_strain_lf[0].yy = 0.0;
    eigen_strain_lf[0].zz = 0.0;
    eigen_strain_lf[0].yz = 0.0;
    eigen_strain_lf[0].xz = 0.0;
    eigen_strain_lf[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_lf[X].xx += eigen_strain_phase[ip].xx*stgO_lf[0].phi[ip];
      eigen_strain_lf[X].yy += eigen_strain_phase[ip].yy*stgO_lf[0].phi[ip];
      eigen_strain_lf[X].zz += eigen_strain_phase[ip].zz*stgO_lf[0].phi[ip];
      eigen_strain_lf[X].yz += eigen_strain_phase[ip].yz*stgO_lf[0].phi[ip];
      eigen_strain_lf[X].xz += eigen_strain_phase[ip].xz*stgO_lf[0].phi[ip];
      eigen_strain_lf[X].xy += eigen_strain_phase[ip].xy*stgO_lf[0].phi[ip];
    }
    eigen_strain_lf[Y] = eigen_strain_lf[X];
    eigen_strain_lf[Z] = eigen_strain_lf[X];

    eigen_strain_rt[0].xx = 0.0;
    eigen_strain_rt[0].yy = 0.0;
    eigen_strain_rt[0].zz = 0.0;
    eigen_strain_rt[0].yz = 0.0;
    eigen_strain_rt[0].xz = 0.0;
    eigen_strain_rt[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_rt[X].xx += eigen_strain_phase[ip].xx*stgO_rt[0].phi[ip];
      eigen_strain_rt[X].yy += eigen_strain_phase[ip].yy*stgO_rt[0].phi[ip];
      eigen_strain_rt[X].zz += eigen_strain_phase[ip].zz*stgO_rt[0].phi[ip];
      eigen_strain_rt[X].yz += eigen_strain_phase[ip].yz*stgO_rt[0].phi[ip];
      eigen_strain_rt[X].xz += eigen_strain_phase[ip].xz*stgO_rt[0].phi[ip];
      eigen_strain_rt[X].xy += eigen_strain_phase[ip].xy*stgO_rt[0].phi[ip];
    }
    eigen_strain_rt[Y] = eigen_strain_rt[X];
    eigen_strain_rt[Z] = eigen_strain_rt[X];

    eigen_strain_bt[0].xx = 0.0;
    eigen_strain_bt[0].yy = 0.0;
    eigen_strain_bt[0].zz = 0.0;
    eigen_strain_bt[0].yz = 0.0;
    eigen_strain_bt[0].xz = 0.0;
    eigen_strain_bt[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_bt[X].xx += eigen_strain_phase[ip].xx*stgO_bt[0].phi[ip];
      eigen_strain_bt[X].yy += eigen_strain_phase[ip].yy*stgO_bt[0].phi[ip];
      eigen_strain_bt[X].zz += eigen_strain_phase[ip].zz*stgO_bt[0].phi[ip];
      eigen_strain_bt[X].yz += eigen_strain_phase[ip].yz*stgO_bt[0].phi[ip];
      eigen_strain_bt[X].xz += eigen_strain_phase[ip].xz*stgO_bt[0].phi[ip];
      eigen_strain_bt[X].xy += eigen_strain_phase[ip].xy*stgO_bt[0].phi[ip];
    }
    eigen_strain_bt[Y] = eigen_strain_bt[X];
    eigen_strain_bt[Z] = eigen_strain_bt[X];

    eigen_strain_tp[0].xx = 0.0;
    eigen_strain_tp[0].yy = 0.0;
    eigen_strain_tp[0].zz = 0.0;
    eigen_strain_tp[0].yz = 0.0;
    eigen_strain_tp[0].xz = 0.0;
    eigen_strain_tp[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_tp[X].xx += eigen_strain_phase[ip].xx*stgO_tp[0].phi[ip];
      eigen_strain_tp[X].yy += eigen_strain_phase[ip].yy*stgO_tp[0].phi[ip];
      eigen_strain_tp[X].zz += eigen_strain_phase[ip].zz*stgO_tp[0].phi[ip];
      eigen_strain_tp[X].yz += eigen_strain_phase[ip].yz*stgO_tp[0].phi[ip];
      eigen_strain_tp[X].xz += eigen_strain_phase[ip].xz*stgO_tp[0].phi[ip];
      eigen_strain_tp[X].xy += eigen_strain_phase[ip].xy*stgO_tp[0].phi[ip];
    }
    eigen_strain_tp[Y] = eigen_strain_tp[X];
    eigen_strain_tp[Z] = eigen_strain_tp[X];

    eigen_strain_bk[0].xx = 0.0;
    eigen_strain_bk[0].yy = 0.0;
    eigen_strain_bk[0].zz = 0.0;
    eigen_strain_bk[0].yz = 0.0;
    eigen_strain_bk[0].xz = 0.0;
    eigen_strain_bk[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_bk[X].xx += eigen_strain_phase[ip].xx*stgO_bk[0].phi[ip];
      eigen_strain_bk[X].yy += eigen_strain_phase[ip].yy*stgO_bk[0].phi[ip];
      eigen_strain_bk[X].zz += eigen_strain_phase[ip].zz*stgO_bk[0].phi[ip];
      eigen_strain_bk[X].yz += eigen_strain_phase[ip].yz*stgO_bk[0].phi[ip];
      eigen_strain_bk[X].xz += eigen_strain_phase[ip].xz*stgO_bk[0].phi[ip];
      eigen_strain_bk[X].xy += eigen_strain_phase[ip].xy*stgO_bk[0].phi[ip];
    }
    eigen_strain_bk[Y] = eigen_strain_bk[X];
    eigen_strain_bk[Z] = eigen_strain_bk[X];

    eigen_strain_ft[0].xx = 0.0;
    eigen_strain_ft[0].yy = 0.0;
    eigen_strain_ft[0].zz = 0.0;
    eigen_strain_ft[0].yz = 0.0;
    eigen_strain_ft[0].xz = 0.0;
    eigen_strain_ft[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_ft[X].xx += eigen_strain_phase[ip].xx*stgO_ft[0].phi[ip];
      eigen_strain_ft[X].yy += eigen_strain_phase[ip].yy*stgO_ft[0].phi[ip];
      eigen_strain_ft[X].zz += eigen_strain_phase[ip].zz*stgO_ft[0].phi[ip];
      eigen_strain_ft[X].yz += eigen_strain_phase[ip].yz*stgO_ft[0].phi[ip];
      eigen_strain_ft[X].xz += eigen_strain_phase[ip].xz*stgO_ft[0].phi[ip];
      eigen_strain_ft[X].xy += eigen_strain_phase[ip].xy*stgO_ft[0].phi[ip];
    }
    eigen_strain_ft[Y] = eigen_strain_ft[X];
    eigen_strain_ft[Z] = eigen_strain_ft[X];

    strain_lf[X].xx = 0.5*(st_it_gO_lf[6].disp[X][2] - st_it_gO_lf[5].disp[X][2]) - eigen_strain_lf[X].xx;
    strain_lf[X].yy = 0.5*(st_it_gO_lf[2].disp[Y][2] - st_it_gO_lf[1].disp[Y][2]) - eigen_strain_lf[X].yy;
    strain_lf[X].zz = 0.5*(st_it_gO_lf[4].disp[Z][2] - st_it_gO_lf[3].disp[Z][2]) - eigen_strain_lf[X].zz;
    strain_lf[Y].xx = strain_lf[X].xx;
    strain_lf[Y].yy = strain_lf[X].yy;
    strain_lf[Y].zz = strain_lf[X].zz;
    strain_lf[Z].yy = strain_lf[X].yy;
    strain_lf[Z].xx = strain_lf[X].xx;
    strain_lf[Z].zz = strain_lf[X].zz;

    strain_rt[X].xx = 0.5*(st_it_gO_rt[6].disp[X][2] - st_it_gO_rt[5].disp[X][2]) - eigen_strain_rt[X].xx;
    strain_rt[X].yy = 0.5*(st_it_gO_rt[2].disp[Y][2] - st_it_gO_rt[1].disp[Y][2]) - eigen_strain_rt[X].yy;
    strain_rt[X].zz = 0.5*(st_it_gO_rt[4].disp[Z][2] - st_it_gO_rt[3].disp[Z][2]) - eigen_strain_rt[X].zz;
    strain_rt[Y].xx = strain_rt[X].xx;
    strain_rt[Y].yy = strain_rt[X].yy;
    strain_rt[Y].zz = strain_rt[X].zz;
    strain_rt[Z].yy = strain_rt[X].yy;
    strain_rt[Z].xx = strain_rt[X].xx;
    strain_rt[Z].zz = strain_rt[X].zz;

    strain_bt[X].xx = 0.5*(st_it_gO_bt[6].disp[X][2] - st_it_gO_bt[5].disp[X][2]) - eigen_strain_bt[X].xx;
    strain_bt[X].yy = 0.5*(st_it_gO_bt[2].disp[Y][2] - st_it_gO_bt[1].disp[Y][2]) - eigen_strain_bt[X].yy;
    strain_bt[X].zz = 0.5*(st_it_gO_bt[4].disp[Z][2] - st_it_gO_bt[3].disp[Z][2]) - eigen_strain_bt[X].zz;
    strain_bt[Y].xx = strain_bt[X].xx;
    strain_bt[Y].yy = strain_bt[X].yy;
    strain_bt[Y].zz = strain_bt[X].zz;
    strain_bt[Z].yy = strain_bt[X].yy;
    strain_bt[Z].xx = strain_bt[X].xx;
    strain_bt[Z].zz = strain_bt[X].zz;

    strain_tp[X].xx = 0.5*(st_it_gO_tp[6].disp[X][2] - st_it_gO_tp[5].disp[X][2]) - eigen_strain_tp[X].xx;
    strain_tp[X].yy = 0.5*(st_it_gO_tp[2].disp[Y][2] - st_it_gO_tp[1].disp[Y][2]) - eigen_strain_tp[X].yy;
    strain_tp[X].zz = 0.5*(st_it_gO_tp[4].disp[Z][2] - st_it_gO_tp[3].disp[Z][2]) - eigen_strain_tp[X].zz;
    strain_tp[Y].xx = strain_tp[X].xx;
    strain_tp[Y].yy = strain_tp[X].yy;
    strain_tp[Y].zz = strain_tp[X].zz;
    strain_tp[Z].yy = strain_tp[X].yy;
    strain_tp[Z].xx = strain_tp[X].xx;
    strain_tp[Z].zz = strain_tp[X].zz;

    strain_bk[X].xx = 0.5*(st_it_gO_bk[6].disp[X][2] - st_it_gO_bk[5].disp[X][2]) - eigen_strain_bk[X].xx;
    strain_bk[X].yy = 0.5*(st_it_gO_bk[2].disp[Y][2] - st_it_gO_bk[1].disp[Y][2]) - eigen_strain_bk[X].yy;
    strain_bk[X].zz = 0.5*(st_it_gO_bk[4].disp[Z][2] - st_it_gO_bk[3].disp[Z][2]) - eigen_strain_bk[X].zz;
    strain_bk[Y].xx = strain_bk[X].xx;
    strain_bk[Y].yy = strain_bk[X].yy;
    strain_bk[Y].zz = strain_bk[X].zz;
    strain_bk[Z].yy = strain_bk[X].yy;
    strain_bk[Z].xx = strain_bk[X].xx;
    strain_bk[Z].zz = strain_bk[X].zz;

    strain_ft[X].xx = 0.5*(st_it_gO_ft[6].disp[X][2] - st_it_gO_ft[5].disp[X][2]) - eigen_strain_ft[X].xx;
    strain_ft[X].yy = 0.5*(st_it_gO_ft[2].disp[Y][2] - st_it_gO_ft[1].disp[Y][2]) - eigen_strain_ft[X].yy;
    strain_ft[X].zz = 0.5*(st_it_gO_ft[4].disp[Z][2] - st_it_gO_ft[3].disp[Z][2]) - eigen_strain_ft[X].zz;
    strain_ft[Y].xx = strain_ft[X].xx;
    strain_ft[Y].yy = strain_ft[X].yy;
    strain_ft[Y].zz = strain_ft[X].zz;
    strain_ft[Z].yy = strain_ft[X].yy;
    strain_ft[Z].xx = strain_ft[X].xx;
    strain_ft[Z].zz = strain_ft[X].zz;

    strain_lf[X].xy = 0.25*((st_it_gO_lf[2].disp[X][2] - st_it_gO_lf[1].disp[X][2]) + (st_it_gO_lf[6].disp[Y][2] - st_it_gO_lf[5].disp[Y][2]));
    strain_lf[X].xz = 0.25*((st_it_gO_lf[4].disp[X][2] - st_it_gO_lf[3].disp[X][2]) + (st_it_gO_lf[6].disp[Z][2] - st_it_gO_lf[5].disp[Z][2]));
    strain_lf[X].yz = 0.25*((st_it_gO_lf[2].disp[Z][2] - st_it_gO_lf[1].disp[Z][2]) + (st_it_gO_lf[4].disp[Y][2] - st_it_gO_lf[3].disp[Y][2]));
    strain_lf[Y].xy = strain_lf[X].xy;
    strain_lf[Y].xz = strain_lf[X].xz;
    strain_lf[Y].yz = strain_lf[X].yz;
    strain_lf[Z].xy = strain_lf[X].xy;
    strain_lf[Z].xz = strain_lf[X].xz;
    strain_lf[Z].yz = strain_lf[X].yz;

    strain_rt[X].xy = 0.25*((st_it_gO_rt[2].disp[X][2] - st_it_gO_rt[1].disp[X][2]) + (st_it_gO_rt[6].disp[Y][2] - st_it_gO_rt[5].disp[Y][2]));
    strain_rt[X].xz = 0.25*((st_it_gO_rt[4].disp[X][2] - st_it_gO_rt[3].disp[X][2]) + (st_it_gO_rt[6].disp[Z][2] - st_it_gO_rt[5].disp[Z][2]));
    strain_rt[X].yz = 0.25*((st_it_gO_rt[2].disp[Z][2] - st_it_gO_rt[1].disp[Z][2]) + (st_it_gO_rt[4].disp[Y][2] - st_it_gO_rt[3].disp[Y][2]));
    strain_rt[Y].xy = strain_rt[X].xy;
    strain_rt[Y].xz = strain_rt[X].xz;
    strain_rt[Y].yz = strain_rt[X].yz;
    strain_rt[Z].xy = strain_rt[X].xy;
    strain_rt[Z].xz = strain_rt[X].xz;
    strain_rt[Z].yz = strain_rt[X].yz;

    strain_bt[X].xy = 0.25*((st_it_gO_bt[2].disp[X][2] - st_it_gO_bt[1].disp[X][2]) + (st_it_gO_bt[6].disp[Y][2] - st_it_gO_bt[5].disp[Y][2]));
    strain_bt[X].xz = 0.25*((st_it_gO_bt[4].disp[X][2] - st_it_gO_bt[3].disp[X][2]) + (st_it_gO_bt[6].disp[Z][2] - st_it_gO_bt[5].disp[Z][2]));
    strain_bt[X].yz = 0.25*((st_it_gO_bt[2].disp[Z][2] - st_it_gO_bt[1].disp[Z][2]) + (st_it_gO_bt[4].disp[Y][2] - st_it_gO_bt[3].disp[Y][2]));
    strain_bt[Y].xy = strain_bt[X].xy;
    strain_bt[Y].xz = strain_bt[X].xz;
    strain_bt[Y].yz = strain_bt[X].yz;
    strain_bt[Z].xy = strain_bt[X].xy;
    strain_bt[Z].xz = strain_bt[X].xz;
    strain_bt[Z].yz = strain_bt[X].yz;

    strain_tp[X].xy = 0.25*((st_it_gO_tp[2].disp[X][2] - st_it_gO_tp[1].disp[X][2]) + (st_it_gO_tp[6].disp[Y][2] - st_it_gO_tp[5].disp[Y][2]));
    strain_tp[X].xz = 0.25*((st_it_gO_tp[4].disp[X][2] - st_it_gO_tp[3].disp[X][2]) + (st_it_gO_tp[6].disp[Z][2] - st_it_gO_tp[5].disp[Z][2]));
    strain_tp[X].yz = 0.25*((st_it_gO_tp[2].disp[Z][2] - st_it_gO_tp[1].disp[Z][2]) + (st_it_gO_tp[4].disp[Y][2] - st_it_gO_tp[3].disp[Y][2]));
    strain_tp[Y].xy = strain_tp[X].xy;
    strain_tp[Y].xz = strain_tp[X].xz;
    strain_tp[Y].yz = strain_tp[X].yz;
    strain_tp[Z].xy = strain_tp[X].xy;
    strain_tp[Z].xz = strain_tp[X].xz;
    strain_tp[Z].yz = strain_tp[X].yz;

    strain_bk[X].xy = 0.25*((st_it_gO_bk[2].disp[X][2] - st_it_gO_bk[1].disp[X][2]) + (st_it_gO_bk[6].disp[Y][2] - st_it_gO_bk[5].disp[Y][2]));
    strain_bk[X].xz = 0.25*((st_it_gO_bk[4].disp[X][2] - st_it_gO_bk[3].disp[X][2]) + (st_it_gO_bk[6].disp[Z][2] - st_it_gO_bk[5].disp[Z][2]));
    strain_bk[X].yz = 0.25*((st_it_gO_bk[2].disp[Z][2] - st_it_gO_bk[1].disp[Z][2]) + (st_it_gO_bk[4].disp[Y][2] - st_it_gO_bk[3].disp[Y][2]));
    strain_bk[Y].xy = strain_bk[X].xy;
    strain_bk[Y].xz = strain_bk[X].xz;
    strain_bk[Y].yz = strain_bk[X].yz;
    strain_bk[Z].xy = strain_bk[X].xy;
    strain_bk[Z].xz = strain_bk[X].xz;
    strain_bk[Z].yz = strain_bk[X].yz;

    strain_ft[X].xy = 0.25*((st_it_gO_ft[2].disp[X][2] - st_it_gO_ft[1].disp[X][2]) + (st_it_gO_ft[6].disp[Y][2] - st_it_gO_ft[5].disp[Y][2]));
    strain_ft[X].xz = 0.25*((st_it_gO_ft[4].disp[X][2] - st_it_gO_ft[3].disp[X][2]) + (st_it_gO_ft[6].disp[Z][2] - st_it_gO_ft[5].disp[Z][2]));
    strain_ft[X].yz = 0.25*((st_it_gO_ft[2].disp[Z][2] - st_it_gO_ft[1].disp[Z][2]) + (st_it_gO_ft[4].disp[Y][2] - st_it_gO_ft[3].disp[Y][2]));
    strain_ft[Y].xy = strain_ft[X].xy;
    strain_ft[Y].xz = strain_ft[X].xz;
    strain_ft[Y].yz = strain_ft[X].yz;
    strain_ft[Z].xy = strain_ft[X].xy;
    strain_ft[Z].xz = strain_ft[X].xz;
    strain_ft[Z].yz = strain_ft[X].yz;

    stiffness_c_lf[0].C11 = 0.0;
    stiffness_c_lf[0].C12 = 0.0;
    stiffness_c_lf[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_lf[X].C11 += (stiffness_phase_n[ip].C11)*stgO_lf[0].phi[ip];
      stiffness_c_lf[X].C12 += (stiffness_phase_n[ip].C12)*stgO_lf[0].phi[ip];
      stiffness_c_lf[X].C44 += (stiffness_phase_n[ip].C44)*stgO_lf[0].phi[ip];
    }
    stiffness_c_lf[Y] = stiffness_c_lf[X];
    stiffness_c_lf[Z] = stiffness_c_lf[X];

    stiffness_c_rt[0].C11 = 0.0;
    stiffness_c_rt[0].C12 = 0.0;
    stiffness_c_rt[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_rt[X].C11 += (stiffness_phase_n[ip].C11)*stgO_rt[0].phi[ip];
      stiffness_c_rt[X].C12 += (stiffness_phase_n[ip].C12)*stgO_rt[0].phi[ip];
      stiffness_c_rt[X].C44 += (stiffness_phase_n[ip].C44)*stgO_rt[0].phi[ip];
    }
    stiffness_c_rt[Y] = stiffness_c_rt[X];
    stiffness_c_rt[Z] = stiffness_c_rt[X];

    stiffness_c_bt[0].C11 = 0.0;
    stiffness_c_bt[0].C12 = 0.0;
    stiffness_c_bt[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_bt[X].C11 += (stiffness_phase_n[ip].C11)*stgO_bt[0].phi[ip];
      stiffness_c_bt[X].C12 += (stiffness_phase_n[ip].C12)*stgO_bt[0].phi[ip];
      stiffness_c_bt[X].C44 += (stiffness_phase_n[ip].C44)*stgO_bt[0].phi[ip];
    }
    stiffness_c_bt[Y] = stiffness_c_bt[X];
    stiffness_c_bt[Z] = stiffness_c_bt[X];

    stiffness_c_tp[0].C11 = 0.0;
    stiffness_c_tp[0].C12 = 0.0;
    stiffness_c_tp[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_tp[X].C11 += (stiffness_phase_n[ip].C11)*stgO_tp[0].phi[ip];
      stiffness_c_tp[X].C12 += (stiffness_phase_n[ip].C12)*stgO_tp[0].phi[ip];
      stiffness_c_tp[X].C44 += (stiffness_phase_n[ip].C44)*stgO_tp[0].phi[ip];
    }
    stiffness_c_tp[Y] = stiffness_c_tp[X];
    stiffness_c_tp[Z] = stiffness_c_tp[X];

    stiffness_c_bk[0].C11 = 0.0;
    stiffness_c_bk[0].C12 = 0.0;
    stiffness_c_bk[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_bk[X].C11 += (stiffness_phase_n[ip].C11)*stgO_bk[0].phi[ip];
      stiffness_c_bk[X].C12 += (stiffness_phase_n[ip].C12)*stgO_bk[0].phi[ip];
      stiffness_c_bk[X].C44 += (stiffness_phase_n[ip].C44)*stgO_bk[0].phi[ip];
    }
    stiffness_c_bk[Y] = stiffness_c_bk[X];
    stiffness_c_bk[Z] = stiffness_c_bk[X];

    stiffness_c_ft[0].C11 = 0.0;
    stiffness_c_ft[0].C12 = 0.0;
    stiffness_c_ft[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_ft[X].C11 += (stiffness_phase_n[ip].C11)*stgO_ft[0].phi[ip];
      stiffness_c_ft[X].C12 += (stiffness_phase_n[ip].C12)*stgO_ft[0].phi[ip];
      stiffness_c_ft[X].C44 += (stiffness_phase_n[ip].C44)*stgO_ft[0].phi[ip];
    }
    stiffness_c_ft[Y] = stiffness_c_ft[X];
    stiffness_c_ft[Z] = stiffness_c_ft[X];

    trace_strain_front     = strain_ft[X].xx + strain_ft[X].yy + strain_ft[X].zz;
    trace_strain_back      = strain_bk[X].xx + strain_bk[X].yy + strain_bk[X].zz;
    trace_strain_right     = strain_rt[Y].xx + strain_rt[Y].yy + strain_rt[Y].zz;
    trace_strain_left      = strain_lf[Y].xx + strain_lf[Y].yy + strain_lf[Y].zz;
    trace_strain_top       = strain_tp[Z].xx + strain_tp[Z].yy + strain_tp[Z].zz;
    trace_strain_bottom    = strain_bt[Z].xx + strain_bt[Z].yy + strain_bt[Z].zz;


    lambda_front    = stiffness_c_ft[X].C12;  //C12
    lambda_back     = stiffness_c_bk[X].C12;   //C12

    mu_front        = stiffness_c_ft[X].C44;  //C44
    mu_back         = stiffness_c_bk[X].C44;   //C44

    mu_prime_front  = stiffness_c_ft[X].C11 - stiffness_c_ft[X].C12 - 2.0*stiffness_c_ft[X].C44;
    mu_prime_back   = stiffness_c_bk[X].C11 - stiffness_c_bk[X].C12 - 2.0*stiffness_c_bk[X].C44;

    lambda_right    = stiffness_c_rt[Y].C12;  //C12
    lambda_left     = stiffness_c_lf[Y].C12;   //C12

    mu_right        = stiffness_c_rt[Y].C44;  //C44
    mu_left         = stiffness_c_lf[Y].C44;   //C44

    mu_prime_right  = stiffness_c_rt[Y].C11 - stiffness_c_rt[Y].C12 - 2.0*stiffness_c_rt[Y].C44;
    mu_prime_left   = stiffness_c_lf[Y].C11 - stiffness_c_lf[Y].C12 - 2.0*stiffness_c_lf[Y].C44;

    lambda_top      = stiffness_c_tp[Z].C12;  //C12
    lambda_bottom   = stiffness_c_bt[Z].C12;   //C12

    mu_top          = stiffness_c_tp[Z].C44;  //C44
    mu_bottom       = stiffness_c_bt[Z].C44;   //C44

    mu_prime_top    = stiffness_c_tp[Z].C11 - stiffness_c_tp[Z].C12 - 2.0*stiffness_c_tp[Z].C44;
    mu_prime_bottom = stiffness_c_bt[Z].C11 - stiffness_c_bt[Z].C12 - 2.0*stiffness_c_bt[Z].C44;

    sigma_front.xx  = lambda_front*trace_strain_front  + 2.0*mu_front*strain_ft[X].xx		//Cubic
                      + mu_prime_front*strain_ft[X].xx;
    sigma_back.xx   = lambda_back*trace_strain_back  + 2.0*mu_back*strain_bk[X].xx
                      + mu_prime_back*strain_bk[X].xx;

    sigma_right.xy  = 2.0*mu_right*strain_rt[X].xy;
    sigma_left.xy   = 2.0*mu_left*strain_lf[X].xy;

    sigma_right.yz  = 2.0*mu_right*strain_rt[Z].yz;
    sigma_left.yz   = 2.0*mu_left*strain_lf[Z].yz;

    sigma_top.xz    = 2.0*mu_top*strain_tp[X].xz;
    sigma_bottom.xz = 2.0*mu_bottom*strain_bt[X].xz;

    sigma_top.yz    = 2.0*mu_top*strain_tp[Y].yz;
    sigma_bottom.yz = 2.0*mu_bottom*strain_bt[Y].yz;

    sigma_front.xy  = 2.0*mu_front*strain_ft[Y].xy;
    sigma_back.xy   = 2.0*mu_back*strain_bk[Y].xy;

    sigma_front.xz  = 2.0*mu_front*strain_ft[Z].xz;
    sigma_back.xz   = 2.0*mu_back*strain_bk[Z].xz;

    sigma_right.yy  = lambda_right*trace_strain_right   + 2.0*mu_right*strain_rt[Y].yy + mu_prime_right*strain_rt[Y].yy;

    sigma_left.yy   = lambda_left*trace_strain_left     + 2.0*mu_left*strain_lf[Y].yy + mu_prime_left*strain_lf[Y].yy;

    sigma_top.zz    = lambda_top*trace_strain_top       + 2.0*mu_top*strain_tp[Z].zz + mu_prime_top*strain_tp[Z].zz;

    sigma_bottom.zz = lambda_bottom*trace_strain_bottom + 2.0*mu_bottom*strain_bt[Z].zz + mu_prime_bottom*strain_bt[Z].zz;

    forceX          = (sigma_front.xx - sigma_back.xx)  + (sigma_right.xy - sigma_left.xy) + (sigma_top.xz   - sigma_bottom.xz);
    forceX         *= 0.5;

    forceY          = (sigma_front.xy - sigma_back.xy)  + (sigma_top.yz - sigma_bottom.yz) + (sigma_right.yy - sigma_left.yy);
    forceY         *= 0.5;

    forceZ          = (sigma_front.xz - sigma_back.xz)  + (sigma_right.yz - sigma_left.yz) + (sigma_top.zz   - sigma_bottom.zz);
    forceZ         *= 0.5;

    //printf("%d, %d, %d, %d, %le, %le, %le, %le, %le, %le, %le\n", x, y, z , index, sigma_front.xx, sigma_back.xx, sigma_right.xy, sigma_left.xy, sigma_top.xz, sigma_bottom.xz, forceX);

    it_gridinfoO[center].disp[X][0] = it_gridinfoO[center].disp[X][1];
    it_gridinfoO[center].disp[Y][0] = it_gridinfoO[center].disp[Y][1];
    it_gridinfoO[center].disp[Z][0] = it_gridinfoO[center].disp[Z][1];
    it_gridinfoO[center].disp[X][1] = it_gridinfoO[center].disp[X][2];
    it_gridinfoO[center].disp[Y][1] = it_gridinfoO[center].disp[Y][2];
    it_gridinfoO[center].disp[Z][1] = it_gridinfoO[center].disp[Z][2];

    it_gridinfoO[center].disp[X][2] = (((deltat_e*deltat_e)/rho)*forceX - (1 - damping_factor*deltat_e)*it_gridinfoO[center].disp[X][0] + 2*it_gridinfoO[center].disp[X][1])/(1.0 + damping_factor*deltat_e);
    it_gridinfoO[center].disp[Y][2] = (((deltat_e*deltat_e)/rho)*forceY - (1 - damping_factor*deltat_e)*it_gridinfoO[center].disp[Y][0] + 2*it_gridinfoO[center].disp[Y][1])/(1.0 + damping_factor*deltat_e);
    it_gridinfoO[center].disp[Z][2] = (((deltat_e*deltat_e)/rho)*forceZ - (1 - damping_factor*deltat_e)*it_gridinfoO[center].disp[Z][0] + 2*it_gridinfoO[center].disp[Z][1])/(1.0 + damping_factor*deltat_e);

    //printf("%d, %d, %d, %d, %le, %le, %le\n", x, y, z , index, damping_factor, rho, forceZ);

    //printf("%d, %d, %d, %d, %d, %le, %le, %le\n", tstep[0], x, y, z, index, it_gridinfoO[center].disp[X][2], it_gridinfoO[center].disp[Y][2], it_gridinfoO[center].disp[Z][2]);




  }
  }
  else {

  if ( ( x > 1 ) && ( x < (nx-2) ) && ( y > 1 ) && ( y < (ny-2) ) ) {
  //if ( x > 1 && x < nx-2 && y > 1 && y < ny-2 && z > 1 && z < nz-2 ) {

    index = y + ny*(z + x*nz);
    center = index;

    //printf("#%d, %d, %d, %d, %le, %le, %le\n", x, y, z , index, it_gridinfoO[center].disp[X][2], it_gridinfoO[center].disp[Y][2], it_gridinfoO[center].disp[Z][2]);

    stgO[0] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    stgO[1] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    stgO[2] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    stgO[3] = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    stgO[4] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];

    st_it_gO[0] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO[1] = it_gridinfoO[ (y-1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO[2] = it_gridinfoO[ (y+1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO[3] = it_gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  ) ) ];
    st_it_gO[4] = it_gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  ) ) ];

    eigen_strain[0].xx = 0.0;
    eigen_strain[0].yy = 0.0;
    eigen_strain[0].xy = 0.0;

    for (ip = 0; ip < npha; ip++) {

      eigen_strain[X].xx += eigen_strain_phase[ip].xx*stgO[0].phi[ip];
      eigen_strain[X].yy += eigen_strain_phase[ip].yy*stgO[0].phi[ip];
      eigen_strain[X].xy += eigen_strain_phase[ip].xy*stgO[0].phi[ip];

    }

    for (i1 = 1; i1 < 3; i1++) {

      eigen_strain[i1].xx = eigen_strain[X].xx;
      eigen_strain[i1].yy = eigen_strain[X].yy;
      eigen_strain[i1].xy = eigen_strain[X].xy;

    }
    // or
    eigen_strain[Y] = eigen_strain[X];

    strain[X].xx = 0.5*(st_it_gO[4].disp[X][2] - st_it_gO[3].disp[X][2]) - eigen_strain[X].xx;
    strain[X].yy = 0.5*(st_it_gO[2].disp[Y][2] - st_it_gO[1].disp[Y][2]) - eigen_strain[X].yy;

    strain[Y].xx = strain[X].xx;
    strain[Y].yy = strain[X].yy;

    strain[X].xy = 0.25*((st_it_gO[2].disp[X][2] - st_it_gO[1].disp[X][2]) + (st_it_gO[4].disp[Y][2] - st_it_gO[3].disp[Y][2]));

    strain[Y].xy = strain[X].xy;

    stiffness_c[0].C11 = 0.0;
    stiffness_c[0].C12 = 0.0;
    stiffness_c[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c[X].C11 += (stiffness_phase_n[ip].C11)*stgO[0].phi[ip];
      stiffness_c[X].C12 += (stiffness_phase_n[ip].C12)*stgO[0].phi[ip];
      stiffness_c[X].C44 += (stiffness_phase_n[ip].C44)*stgO[0].phi[ip];
    }
    stiffness_c[Y] = stiffness_c[X];


    stgO_lf[0] = gridinfoO[ (y  -1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_lf[1] = gridinfoO[ (y-1-1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_lf[2] = gridinfoO[ (y+1-1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_lf[3] = gridinfoO[ (y  -1) + ny*( (x-1)*nz + (z  ) ) ];
    stgO_lf[4] = gridinfoO[ (y  -1) + ny*( (x+1)*nz + (z  ) ) ];

    st_it_gO_lf[0] = it_gridinfoO[ (y  -1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_lf[1] = it_gridinfoO[ (y-1-1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_lf[2] = it_gridinfoO[ (y+1-1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_lf[3] = it_gridinfoO[ (y  -1) + ny*( (x-1)*nz + (z  ) ) ];
    st_it_gO_lf[4] = it_gridinfoO[ (y  -1) + ny*( (x+1)*nz + (z  ) ) ];




    stgO_rt[0] = gridinfoO[ (y  +1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_rt[1] = gridinfoO[ (y-1+1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_rt[2] = gridinfoO[ (y+1+1) + ny*( (x  )*nz + (z  ) ) ];
    stgO_rt[3] = gridinfoO[ (y  +1) + ny*( (x-1)*nz + (z  ) ) ];
    stgO_rt[4] = gridinfoO[ (y  +1) + ny*( (x+1)*nz + (z  ) ) ];



    st_it_gO_rt[0] = it_gridinfoO[ (y  +1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_rt[1] = it_gridinfoO[ (y-1+1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_rt[2] = it_gridinfoO[ (y+1+1) + ny*( (x  )*nz + (z  ) ) ];
    st_it_gO_rt[3] = it_gridinfoO[ (y  +1) + ny*( (x-1)*nz + (z  ) ) ];
    st_it_gO_rt[4] = it_gridinfoO[ (y  +1) + ny*( (x+1)*nz + (z  ) ) ];




    stgO_bt[0] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  -1) ) ];
    stgO_bt[1] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  -1) ) ];
    stgO_bt[2] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  -1) ) ];
    stgO_bt[3] = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  -1) ) ];
    stgO_bt[4] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  -1) ) ];



    st_it_gO_bt[0] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  -1) ) ];
    st_it_gO_bt[1] = it_gridinfoO[ (y-1) + ny*( (x  )*nz + (z  -1) ) ];
    st_it_gO_bt[2] = it_gridinfoO[ (y+1) + ny*( (x  )*nz + (z  -1) ) ];
    st_it_gO_bt[3] = it_gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  -1) ) ];
    st_it_gO_bt[4] = it_gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  -1) ) ];




    stgO_tp[0] = gridinfoO[ (y  ) + ny*( (x  )*nz + (z  +1) ) ];
    stgO_tp[1] = gridinfoO[ (y-1) + ny*( (x  )*nz + (z  +1) ) ];
    stgO_tp[2] = gridinfoO[ (y+1) + ny*( (x  )*nz + (z  +1) ) ];
    stgO_tp[3] = gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  +1) ) ];
    stgO_tp[4] = gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  +1) ) ];



    st_it_gO_tp[0] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  +1) ) ];
    st_it_gO_tp[1] = it_gridinfoO[ (y-1) + ny*( (x  )*nz + (z  +1) ) ];
    st_it_gO_tp[2] = it_gridinfoO[ (y+1) + ny*( (x  )*nz + (z  +1) ) ];
    st_it_gO_tp[3] = it_gridinfoO[ (y  ) + ny*( (x-1)*nz + (z  +1) ) ];
    st_it_gO_tp[4] = it_gridinfoO[ (y  ) + ny*( (x+1)*nz + (z  +1) ) ];




    stgO_bk[0] = gridinfoO[ (y  ) + ny*( (x  -1)*nz + (z  ) ) ];
    stgO_bk[1] = gridinfoO[ (y-1) + ny*( (x  -1)*nz + (z  ) ) ];
    stgO_bk[2] = gridinfoO[ (y+1) + ny*( (x  -1)*nz + (z  ) ) ];
    stgO_bk[3] = gridinfoO[ (y  ) + ny*( (x-1-1)*nz + (z  ) ) ];
    stgO_bk[4] = gridinfoO[ (y  ) + ny*( (x+1-1)*nz + (z  ) ) ];



    st_it_gO_bk[0] = it_gridinfoO[ (y  ) + ny*( (x  -1)*nz + (z  ) ) ];
    st_it_gO_bk[1] = it_gridinfoO[ (y-1) + ny*( (x  -1)*nz + (z  ) ) ];
    st_it_gO_bk[2] = it_gridinfoO[ (y+1) + ny*( (x  -1)*nz + (z  ) ) ];
    st_it_gO_bk[3] = it_gridinfoO[ (y  ) + ny*( (x-1-1)*nz + (z  ) ) ];
    st_it_gO_bk[4] = it_gridinfoO[ (y  ) + ny*( (x+1-1)*nz + (z  ) ) ];




    stgO_ft[0] = gridinfoO[ (y  ) + ny*( (x  +1)*nz + (z  ) ) ];
    stgO_ft[1] = gridinfoO[ (y-1) + ny*( (x  +1)*nz + (z  ) ) ];
    stgO_ft[2] = gridinfoO[ (y+1) + ny*( (x  +1)*nz + (z  ) ) ];
    stgO_ft[3] = gridinfoO[ (y  ) + ny*( (x-1+1)*nz + (z  ) ) ];
    stgO_ft[4] = gridinfoO[ (y  ) + ny*( (x+1+1)*nz + (z  ) ) ];



    st_it_gO_ft[0] = it_gridinfoO[ (y  ) + ny*( (x  +1)*nz + (z  ) ) ];
    st_it_gO_ft[1] = it_gridinfoO[ (y-1) + ny*( (x  +1)*nz + (z  ) ) ];
    st_it_gO_ft[2] = it_gridinfoO[ (y+1) + ny*( (x  +1)*nz + (z  ) ) ];
    st_it_gO_ft[3] = it_gridinfoO[ (y  ) + ny*( (x-1+1)*nz + (z  ) ) ];
    st_it_gO_ft[4] = it_gridinfoO[ (y  ) + ny*( (x+1+1)*nz + (z  ) ) ];



    eigen_strain_lf[0].xx = 0.0;
    eigen_strain_lf[0].yy = 0.0;
    eigen_strain_lf[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_lf[X].xx += eigen_strain_phase[ip].xx*stgO_lf[0].phi[ip];
      eigen_strain_lf[X].yy += eigen_strain_phase[ip].yy*stgO_lf[0].phi[ip];
      eigen_strain_lf[X].xy += eigen_strain_phase[ip].xy*stgO_lf[0].phi[ip];
    }
    eigen_strain_lf[Y] = eigen_strain_lf[X];

    eigen_strain_rt[0].xx = 0.0;
    eigen_strain_rt[0].yy = 0.0;
    eigen_strain_rt[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_rt[X].xx += eigen_strain_phase[ip].xx*stgO_rt[0].phi[ip];
      eigen_strain_rt[X].yy += eigen_strain_phase[ip].yy*stgO_rt[0].phi[ip];
      eigen_strain_rt[X].xy += eigen_strain_phase[ip].xy*stgO_rt[0].phi[ip];
    }
    eigen_strain_rt[Y] = eigen_strain_rt[X];

    eigen_strain_bt[0].xx = 0.0;
    eigen_strain_bt[0].yy = 0.0;
    eigen_strain_bt[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_bt[X].xx += eigen_strain_phase[ip].xx*stgO_bt[0].phi[ip];
      eigen_strain_bt[X].yy += eigen_strain_phase[ip].yy*stgO_bt[0].phi[ip];
      eigen_strain_bt[X].xy += eigen_strain_phase[ip].xy*stgO_bt[0].phi[ip];
    }
    eigen_strain_bt[Y] = eigen_strain_bt[X];

    eigen_strain_tp[0].xx = 0.0;
    eigen_strain_tp[0].yy = 0.0;
    eigen_strain_tp[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_tp[X].xx += eigen_strain_phase[ip].xx*stgO_tp[0].phi[ip];
      eigen_strain_tp[X].yy += eigen_strain_phase[ip].yy*stgO_tp[0].phi[ip];
      eigen_strain_tp[X].xy += eigen_strain_phase[ip].xy*stgO_tp[0].phi[ip];
    }
    eigen_strain_tp[Y] = eigen_strain_tp[X];

    eigen_strain_bk[0].xx = 0.0;
    eigen_strain_bk[0].yy = 0.0;
    eigen_strain_bk[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_bk[X].xx += eigen_strain_phase[ip].xx*stgO_bk[0].phi[ip];
      eigen_strain_bk[X].yy += eigen_strain_phase[ip].yy*stgO_bk[0].phi[ip];
      eigen_strain_bk[X].xy += eigen_strain_phase[ip].xy*stgO_bk[0].phi[ip];
    }
    eigen_strain_bk[Y] = eigen_strain_bk[X];

    eigen_strain_ft[0].xx = 0.0;
    eigen_strain_ft[0].yy = 0.0;
    eigen_strain_ft[0].xy = 0.0;
    for (ip = 0; ip < npha; ip++) {
      eigen_strain_ft[X].xx += eigen_strain_phase[ip].xx*stgO_ft[0].phi[ip];
      eigen_strain_ft[X].yy += eigen_strain_phase[ip].yy*stgO_ft[0].phi[ip];
      eigen_strain_ft[X].xy += eigen_strain_phase[ip].xy*stgO_ft[0].phi[ip];
    }
    eigen_strain_ft[Y] = eigen_strain_ft[X];

    strain_lf[X].xx = 0.5*(st_it_gO_lf[4].disp[X][2] - st_it_gO_lf[3].disp[X][2]) - eigen_strain_lf[X].xx;
    strain_lf[X].yy = 0.5*(st_it_gO_lf[2].disp[Y][2] - st_it_gO_lf[1].disp[Y][2]) - eigen_strain_lf[X].yy;
    strain_lf[Y].xx = strain_lf[X].xx;
    strain_lf[Y].yy = strain_lf[X].yy;

    strain_rt[X].xx = 0.5*(st_it_gO_rt[4].disp[X][2] - st_it_gO_rt[3].disp[X][2]) - eigen_strain_rt[X].xx;
    strain_rt[X].yy = 0.5*(st_it_gO_rt[2].disp[Y][2] - st_it_gO_rt[1].disp[Y][2]) - eigen_strain_rt[X].yy;
    strain_rt[Y].xx = strain_rt[X].xx;
    strain_rt[Y].yy = strain_rt[X].yy;

    strain_bt[X].xx = 0.5*(st_it_gO_bt[4].disp[X][2] - st_it_gO_bt[3].disp[X][2]) - eigen_strain_bt[X].xx;
    strain_bt[X].yy = 0.5*(st_it_gO_bt[2].disp[Y][2] - st_it_gO_bt[1].disp[Y][2]) - eigen_strain_bt[X].yy;
    strain_bt[Y].xx = strain_bt[X].xx;
    strain_bt[Y].yy = strain_bt[X].yy;

    strain_tp[X].xx = 0.5*(st_it_gO_tp[4].disp[X][2] - st_it_gO_tp[3].disp[X][2]) - eigen_strain_tp[X].xx;
    strain_tp[X].yy = 0.5*(st_it_gO_tp[2].disp[Y][2] - st_it_gO_tp[1].disp[Y][2]) - eigen_strain_tp[X].yy;
    strain_tp[Y].xx = strain_tp[X].xx;
    strain_tp[Y].yy = strain_tp[X].yy;

    strain_bk[X].xx = 0.5*(st_it_gO_bk[4].disp[X][2] - st_it_gO_bk[3].disp[X][2]) - eigen_strain_bk[X].xx;
    strain_bk[X].yy = 0.5*(st_it_gO_bk[2].disp[Y][2] - st_it_gO_bk[1].disp[Y][2]) - eigen_strain_bk[X].yy;
    strain_bk[Y].xx = strain_bk[X].xx;
    strain_bk[Y].yy = strain_bk[X].yy;

    strain_ft[X].xx = 0.5*(st_it_gO_ft[4].disp[X][2] - st_it_gO_ft[3].disp[X][2]) - eigen_strain_ft[X].xx;
    strain_ft[X].yy = 0.5*(st_it_gO_ft[2].disp[Y][2] - st_it_gO_ft[1].disp[Y][2]) - eigen_strain_ft[X].yy;
    strain_ft[Y].xx = strain_ft[X].xx;
    strain_ft[Y].yy = strain_ft[X].yy;

    strain_lf[X].xy = 0.25*((st_it_gO_lf[2].disp[X][2] - st_it_gO_lf[1].disp[X][2]) + (st_it_gO_lf[4].disp[Y][2] - st_it_gO_lf[3].disp[Y][2]));
    strain_lf[Y].xy = strain_lf[X].xy;

    strain_rt[X].xy = 0.25*((st_it_gO_rt[2].disp[X][2] - st_it_gO_rt[1].disp[X][2]) + (st_it_gO_rt[4].disp[Y][2] - st_it_gO_rt[3].disp[Y][2]));
    strain_rt[Y].xy = strain_rt[X].xy;

    strain_bt[X].xy = 0.25*((st_it_gO_bt[2].disp[X][2] - st_it_gO_bt[1].disp[X][2]) + (st_it_gO_bt[4].disp[Y][2] - st_it_gO_bt[3].disp[Y][2]));
    strain_bt[Y].xy = strain_bt[X].xy;

    strain_tp[X].xy = 0.25*((st_it_gO_tp[2].disp[X][2] - st_it_gO_tp[1].disp[X][2]) + (st_it_gO_tp[4].disp[Y][2] - st_it_gO_tp[3].disp[Y][2]));
    strain_tp[Y].xy = strain_tp[X].xy;

    strain_bk[X].xy = 0.25*((st_it_gO_bk[2].disp[X][2] - st_it_gO_bk[1].disp[X][2]) + (st_it_gO_bk[4].disp[Y][2] - st_it_gO_bk[3].disp[Y][2]));
    strain_bk[Y].xy = strain_bk[X].xy;

    strain_ft[X].xy = 0.25*((st_it_gO_ft[2].disp[X][2] - st_it_gO_ft[1].disp[X][2]) + (st_it_gO_ft[4].disp[Y][2] - st_it_gO_ft[3].disp[Y][2]));
    strain_ft[Y].xy = strain_ft[X].xy;
    strain_ft[Z].xy = strain_ft[X].xy;

    stiffness_c_lf[0].C11 = 0.0;
    stiffness_c_lf[0].C12 = 0.0;
    stiffness_c_lf[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_lf[X].C11 += (stiffness_phase_n[ip].C11)*stgO_lf[0].phi[ip];
      stiffness_c_lf[X].C12 += (stiffness_phase_n[ip].C12)*stgO_lf[0].phi[ip];
      stiffness_c_lf[X].C44 += (stiffness_phase_n[ip].C44)*stgO_lf[0].phi[ip];
    }
    stiffness_c_lf[Y] = stiffness_c_lf[X];

    stiffness_c_rt[0].C11 = 0.0;
    stiffness_c_rt[0].C12 = 0.0;
    stiffness_c_rt[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_rt[X].C11 += (stiffness_phase_n[ip].C11)*stgO_rt[0].phi[ip];
      stiffness_c_rt[X].C12 += (stiffness_phase_n[ip].C12)*stgO_rt[0].phi[ip];
      stiffness_c_rt[X].C44 += (stiffness_phase_n[ip].C44)*stgO_rt[0].phi[ip];
    }
    stiffness_c_rt[Y] = stiffness_c_rt[X];

    stiffness_c_bt[0].C11 = 0.0;
    stiffness_c_bt[0].C12 = 0.0;
    stiffness_c_bt[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_bt[X].C11 += (stiffness_phase_n[ip].C11)*stgO_bt[0].phi[ip];
      stiffness_c_bt[X].C12 += (stiffness_phase_n[ip].C12)*stgO_bt[0].phi[ip];
      stiffness_c_bt[X].C44 += (stiffness_phase_n[ip].C44)*stgO_bt[0].phi[ip];
    }
    stiffness_c_bt[Y] = stiffness_c_bt[X];

    stiffness_c_tp[0].C11 = 0.0;
    stiffness_c_tp[0].C12 = 0.0;
    stiffness_c_tp[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_tp[X].C11 += (stiffness_phase_n[ip].C11)*stgO_tp[0].phi[ip];
      stiffness_c_tp[X].C12 += (stiffness_phase_n[ip].C12)*stgO_tp[0].phi[ip];
      stiffness_c_tp[X].C44 += (stiffness_phase_n[ip].C44)*stgO_tp[0].phi[ip];
    }
    stiffness_c_tp[Y] = stiffness_c_tp[X];

    stiffness_c_bk[0].C11 = 0.0;
    stiffness_c_bk[0].C12 = 0.0;
    stiffness_c_bk[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_bk[X].C11 += (stiffness_phase_n[ip].C11)*stgO_bk[0].phi[ip];
      stiffness_c_bk[X].C12 += (stiffness_phase_n[ip].C12)*stgO_bk[0].phi[ip];
      stiffness_c_bk[X].C44 += (stiffness_phase_n[ip].C44)*stgO_bk[0].phi[ip];
    }
    stiffness_c_bk[Y] = stiffness_c_bk[X];

    stiffness_c_ft[0].C11 = 0.0;
    stiffness_c_ft[0].C12 = 0.0;
    stiffness_c_ft[0].C44 = 0.0;
    for (ip = 0; ip < npha; ip++) {
      stiffness_c_ft[X].C11 += (stiffness_phase_n[ip].C11)*stgO_ft[0].phi[ip];
      stiffness_c_ft[X].C12 += (stiffness_phase_n[ip].C12)*stgO_ft[0].phi[ip];
      stiffness_c_ft[X].C44 += (stiffness_phase_n[ip].C44)*stgO_ft[0].phi[ip];
    }
    stiffness_c_ft[Y] = stiffness_c_ft[X];

    trace_strain_front     = strain_ft[X].xx + strain_ft[X].yy;// + strain_ft[X].zz;
    trace_strain_back      = strain_bk[X].xx + strain_bk[X].yy;// + strain_bk[X].zz;
    trace_strain_right     = strain_rt[Y].xx + strain_rt[Y].yy;// + strain_rt[Y].zz;
    trace_strain_left      = strain_lf[Y].xx + strain_lf[Y].yy;// + strain_lf[Y].zz;
    trace_strain_top       = 0.0;//strain_tp[Z].xx + strain_tp[Z].yy;// + strain_tp[Z].zz;
    trace_strain_bottom    = 0.0;//strain_bt[Z].xx + strain_bt[Z].yy;// + strain_bt[Z].zz;


    lambda_front    = stiffness_c_ft[X].C12;  //C12
    lambda_back     = stiffness_c_bk[X].C12;   //C12

    mu_front        = stiffness_c_ft[X].C44;  //C44
    mu_back         = stiffness_c_bk[X].C44;   //C44

    mu_prime_front  = stiffness_c_ft[X].C11 - stiffness_c_ft[X].C12 - 2.0*stiffness_c_ft[X].C44;
    mu_prime_back   = stiffness_c_bk[X].C11 - stiffness_c_bk[X].C12 - 2.0*stiffness_c_bk[X].C44;

    lambda_right    = stiffness_c_rt[Y].C12;  //C12
    lambda_left     = stiffness_c_lf[Y].C12;   //C12

    mu_right        = stiffness_c_rt[Y].C44;  //C44
    mu_left         = stiffness_c_lf[Y].C44;   //C44

    mu_prime_right  = stiffness_c_rt[Y].C11 - stiffness_c_rt[Y].C12 - 2.0*stiffness_c_rt[Y].C44;
    mu_prime_left   = stiffness_c_lf[Y].C11 - stiffness_c_lf[Y].C12 - 2.0*stiffness_c_lf[Y].C44;

    lambda_top      = 0.0;//stiffness_c_tp[Z].C12;  //C12
    lambda_bottom   = 0.0;//stiffness_c_bt[Z].C12;   //C12

    mu_top          = 0.0;//stiffness_c_tp[Z].C44;  //C44
    mu_bottom       = 0.0;//stiffness_c_bt[Z].C44;   //C44

    mu_prime_top    = 0.0;//stiffness_c_tp[Z].C11 - stiffness_c_tp[Z].C12 - 2.0*stiffness_c_tp[Z].C44;
    mu_prime_bottom = 0.0;//stiffness_c_bt[Z].C11 - stiffness_c_bt[Z].C12 - 2.0*stiffness_c_bt[Z].C44;

    sigma_front.xx  = lambda_front*trace_strain_front  + 2.0*mu_front*strain_ft[X].xx		//Cubic
                      + mu_prime_front*strain_ft[X].xx;
    sigma_back.xx   = lambda_back*trace_strain_back  + 2.0*mu_back*strain_bk[X].xx
                      + mu_prime_back*strain_bk[X].xx;

    sigma_right.xy  = 2.0*mu_right*strain_rt[X].xy;
    sigma_left.xy   = 2.0*mu_left*strain_lf[X].xy;

    sigma_right.yz  = 0.0;//2.0*mu_right*strain_rt[Z].yz;
    sigma_left.yz   = 0.0;//2.0*mu_left*strain_lf[Z].yz;

    sigma_top.xz    = 0.0;//2.0*mu_top*strain_tp[X].xz;
    sigma_bottom.xz = 0.0;//2.0*mu_bottom*strain_bt[X].xz;

    sigma_top.yz    = 0.0;//2.0*mu_top*strain_tp[Y].yz;
    sigma_bottom.yz = 0.0;//2.0*mu_bottom*strain_bt[Y].yz;

    sigma_front.xy  = 2.0*mu_front*strain_ft[Y].xy;
    sigma_back.xy   = 2.0*mu_back*strain_bk[Y].xy;

    sigma_front.xz  = 0.0;//2.0*mu_front*strain_ft[Z].xz;
    sigma_back.xz   = 0.0;//2.0*mu_back*strain_bk[Z].xz;

    sigma_right.yy  = lambda_right*trace_strain_right   + 2.0*mu_right*strain_rt[Y].yy + mu_prime_right*strain_rt[Y].yy;

    sigma_left.yy   = lambda_left*trace_strain_left     + 2.0*mu_left*strain_lf[Y].yy + mu_prime_left*strain_lf[Y].yy;

    sigma_top.zz    = 0.0;//lambda_top*trace_strain_top       + 2.0*mu_top*strain_tp[Z].zz + mu_prime_top*strain_tp[Z].zz;

    sigma_bottom.zz = 0.0;//lambda_bottom*trace_strain_bottom + 2.0*mu_bottom*strain_bt[Z].zz + mu_prime_bottom*strain_bt[Z].zz;

    forceX          = (sigma_front.xx - sigma_back.xx)  + (sigma_right.xy - sigma_left.xy) + (sigma_top.xz   - sigma_bottom.xz);
    forceX         *= 0.5;

    forceY          = (sigma_front.xy - sigma_back.xy)  + (sigma_top.yz - sigma_bottom.yz) + (sigma_right.yy - sigma_left.yy);
    forceY         *= 0.5;

    forceZ          = 0.0;//(sigma_front.xz - sigma_back.xz)  + (sigma_right.yz - sigma_left.yz) + (sigma_top.zz   - sigma_bottom.zz);
    forceZ         *= 0.0;//0.5;

    //printf("%d, %d, %d, %d, %le, %le, %le, %le, %le, %le, %le\n", x, y, z , index, sigma_front.xx, sigma_back.xx, sigma_right.xy, sigma_left.xy, sigma_top.xz, sigma_bottom.xz, forceX);

    it_gridinfoO[center].disp[X][0] = it_gridinfoO[center].disp[X][1];
    it_gridinfoO[center].disp[Y][0] = it_gridinfoO[center].disp[Y][1];
    it_gridinfoO[center].disp[Z][0] = 0.0;//it_gridinfoO[center].disp[Z][1];
    it_gridinfoO[center].disp[X][1] = it_gridinfoO[center].disp[X][2];
    it_gridinfoO[center].disp[Y][1] = it_gridinfoO[center].disp[Y][2];
    it_gridinfoO[center].disp[Z][1] = 0.0;//it_gridinfoO[center].disp[Z][2];

    it_gridinfoO[center].disp[X][2] = (((deltat_e*deltat_e)/rho)*forceX - (1 - damping_factor*deltat_e)*it_gridinfoO[center].disp[X][0] + 2*it_gridinfoO[center].disp[X][1])/(1.0 + damping_factor*deltat_e);
    it_gridinfoO[center].disp[Y][2] = (((deltat_e*deltat_e)/rho)*forceY - (1 - damping_factor*deltat_e)*it_gridinfoO[center].disp[Y][0] + 2*it_gridinfoO[center].disp[Y][1])/(1.0 + damping_factor*deltat_e);
    it_gridinfoO[center].disp[Z][2] = 0.0;//(((deltat_e*deltat_e)/rho)*forceZ - (1 - damping_factor*deltat_e)*it_gridinfoO[center].disp[Z][0] + 2*it_gridinfoO[center].disp[Z][1])/(1.0 + damping_factor*deltat_e);

    //printf("%d, %d, %d, %d, %le, %le, %le\n", x, y, z , index, damping_factor, rho, forceZ);

    //printf("%d, %d, %d, %d, %d, %le, %le, %le\n", tstep[0], x, y, z, index, it_gridinfoO[center].disp[X][2], it_gridinfoO[center].disp[Y][2], it_gridinfoO[center].disp[Z][2]);




  }
  }


}



__kernel void apply_BC_phi_y0_noflux(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index =  (y  ) + ny*( (x  )*nz + (z  ) );
  
  if ( y == 0 ) {
    for ( ip = 0; ip < npha; ip++ ) {
      gridinfo[ (y  ) + ny*( (x  )*nz + (z  ) ) ].phi[ip] = gridinfo[ (y+1) + ny*( (x  )*nz + (z  ) ) ].phi[ip];
    }
  }
    
}

__kernel void apply_BC_phi_yn_noflux(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index =  (y  ) + ny*( (x  )*nz + (z  ) );
  
  if ( y == (ny-1) ) {
    for ( ip = 0; ip < npha; ip++ ) {
      gridinfo[ (y  ) + ny*( (x  )*nz + (z  ) ) ].phi[ip] = gridinfo[ (y-1) + ny*( (x  )*nz + (z  ) ) ].phi[ip];
    }
  }
    
}

__kernel void apply_BC_phi_y0_periodic(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index =  (y  ) + ny*( (x  )*nz + (z  ) );
  
  if ( y == 0 ) {
    for ( ip = 0; ip < npha; ip++ ) {
      gridinfo[ (y  ) + ny*( (x  )*nz + (z  ) ) ].phi[ip] = gridinfo[ (ny-2) + ny*( (x  )*nz + (z  ) ) ].phi[ip];
    }
  }
    
}

__kernel void apply_BC_phi_yn_periodic(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index =  (y  ) + ny*( (x  )*nz + (z  ) );
  
  if ( y == (ny-1) ) {
    for ( ip = 0; ip < npha; ip++ ) {
      gridinfo[ (y  ) + ny*( (x  )*nz + (z  ) ) ].phi[ip] = gridinfo[ (1 ) + ny*( (x  )*nz + (z  ) ) ].phi[ip];
    }
  }
    
}

__kernel void apply_BC_phi_z0_noflux(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index = y + ny*(x*nz + z);
  
  if ( z == 0 ) {
    for ( ip = 0; ip < npha; ip++ ) {
      gridinfo[ y + ny*(x*nz + z) ].phi[ip] = gridinfo[ y + ny*(x*nz + (z+1)) ].phi[ip];
    }
  }
    
}

__kernel void apply_BC_phi_zn_noflux(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index = y + ny*(x*nz + z);
  
  if ( z == (nz-1) ) {
    for ( ip = 0; ip < npha; ip++ ) {
      gridinfo[ y + ny*(x*nz + z) ].phi[ip] = gridinfo[ y + ny*(x*nz + (z-1)) ].phi[ip];
    }
  }
    
}

__kernel void apply_BC_phi_z0_periodic(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index = y + ny*(x*nz + z);
  
  if ( z == 0 ) {
    for ( ip = 0; ip < npha; ip++ ) {
      gridinfo[ y + ny*(x*nz + z) ].phi[ip] = gridinfo[ y + ny*(x*nz + (nz-2)) ].phi[ip];
    }
  }
  
}

__kernel void apply_BC_phi_zn_periodic(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index = y + ny*(x*nz + z);
  
  if ( z == (nz-1) ) {
    for ( ip = 0; ip < npha; ip++ ) {
      gridinfo[ y + ny*(x*nz + z) ].phi[ip] = gridinfo[ y + ny*(x*nz + 1) ].phi[ip];
    }
  }
    
}

__kernel void apply_BC_phi_x0_noflux(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index =  (y  ) + ny*( (x  )*nz + (z  ) );
  
  if ( x == 0 ) {
    for ( ip = 0; ip < npha; ip++ ) {
      gridinfo[ (y  ) + ny*( (x  )*nz + (z  ) ) ].phi[ip] = gridinfo[ (y ) + ny*( (x+1)*nz + (z  ) ) ].phi[ip];
    }

    for (is = 0; is < nsol; is++) {
      gridinfo[ (y  ) + ny*( (x  )*nz + (z  ) ) ].com[is] = gridinfo[ (y ) + ny*( (x+1)*nz + (z  ) ) ].com[is];
      gridinfo[ (y  ) + ny*( (x  )*nz + (z  ) ) ].mu[is] =  gridinfo[ (y ) + ny*( (x+1)*nz + (z  ) ) ].mu[is];
    }
    cscl[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = cscl[ (y ) + ny*( (x+1)*nz + (z  ) ) ];

  }
    
}

__kernel void apply_BC_phi_xn_noflux(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index =  (y  ) + ny*( (x  )*nz + (z  ) );
  
  if ( x == (nx-1) ) {
    for ( ip = 0; ip < npha; ip++ ) {
      gridinfo[ (y  ) + ny*( (x  )*nz + (z  ) ) ].phi[ip] = gridinfo[ (y ) + ny*( (x-1)*nz + (z  ) ) ].phi[ip];
    }

    for (is = 0; is < nsol; is++) {
      gridinfo[ (y  ) + ny*( (x  )*nz + (z  ) ) ].com[is] = gridinfo[ (y ) + ny*( (x-1)*nz + (z  ) ) ].com[is];
      gridinfo[ (y  ) + ny*( (x  )*nz + (z  ) ) ].mu[is] =  gridinfo[ (y ) + ny*( (x-1)*nz + (z  ) ) ].mu[is];
    }
    cscl[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = cscl[ (y ) + ny*( (x-1)*nz + (z  ) ) ];

  }
    
}

__kernel void apply_BC_phi_x0_periodic(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index =  (y  ) + ny*( (x  )*nz + (z  ) );
  
  if ( x == 0 ) {
    for ( ip = 0; ip < npha; ip++ ) {
      gridinfo[ (y  ) + ny*( (x  )*nz + (z  ) ) ].phi[ip] = gridinfo[ (y ) + ny*( (nx-2)*nz + (z  ) ) ].phi[ip];
    }
  }
    
}

__kernel void apply_BC_phi_xn_periodic(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index =  (y  ) + ny*( (x  )*nz + (z  ) );
  
  if ( x == (nx-1) ) {
    for ( ip = 0; ip < npha; ip++ ) {
      gridinfo[ (y  ) + ny*( (x  )*nz + (z  ) ) ].phi[ip] = gridinfo[ (y  ) + ny*( (1  )*nz + (z  ) ) ].phi[ip];
    }
  }
    
}

__kernel void apply_BC_com_y0_noflux(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index = y + ny*(x*nz + z);
  
  if ( y == 0 ) {
    for (is = 0; is < nsol; is++) { 
      gridinfo[ y + ny*(x*nz + z) ].com[is] = gridinfo[ (y+1) + ny*(x*nz + z) ].com[is];
      gridinfo[ y + ny*(x*nz + z) ].mu[is]  = gridinfo[ (y+1) + ny*(x*nz + z) ].mu[is];
    }
    cscl[ y + ny*(x*nz + z) ] = cscl[ (y+1) + ny*(x*nz + z) ];
  }
    
}

__kernel void apply_BC_com_yn_noflux(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index = y + ny*(x*nz + z);
  
  if ( y == (ny-1) ) {
    for (is = 0; is < nsol; is++) { 
      gridinfo[ y + ny*(x*nz + z) ].com[is] = gridinfo[ (y-1) + ny*(x*nz + z) ].com[is];
      gridinfo[ y + ny*(x*nz + z) ].mu[is]  = gridinfo[ (y-1) + ny*(x*nz + z) ].mu[is];
    }
    cscl[ y + ny*(x*nz + z) ] = cscl[ (y-1) + ny*(x*nz + z) ];
  }
    
}

__kernel void apply_BC_com_y0_periodic(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index = y + ny*(x*nz + z);
  
  if ( y == 0 ) {
    for (is = 0; is < nsol; is++) { 
      gridinfo[ y + ny*(x*nz + z) ].com[is] = gridinfo[ (ny-2) + ny*(x*nz + z) ].com[is];
      gridinfo[ y + ny*(x*nz + z) ].mu[is]  = gridinfo[ (ny-2) + ny*(x*nz + z) ].mu[is];
    }
    cscl[ y + ny*(x*nz + z) ] = cscl[ (ny-2) + ny*(x*nz + z) ];
  }
    
}

__kernel void apply_BC_com_yn_periodic(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index = y + ny*(x*nz + z);
  
  if ( y == (ny-1) ) {
    for (is = 0; is < nsol; is++) { 
      gridinfo[ y + ny*(x*nz + z) ].com[is] = gridinfo[ (1) + ny*(x*nz + z) ].com[is];
      gridinfo[ y + ny*(x*nz + z) ].mu[is]  = gridinfo[ (1) + ny*(x*nz + z) ].mu[is];
    }
    cscl[ y + ny*(x*nz + z) ] = cscl[ (1) + ny*(x*nz + z) ];
  }
    
}

__kernel void apply_BC_com_z0_noflux(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index = y + ny*(x*nz + z);
  
  if ( z == 0 ) {
    for (is = 0; is < nsol; is++) { 
      gridinfo[ y + ny*(x*nz + z) ].com[is] = gridinfo[ y + ny*(x*nz + (z+1)) ].com[is];
      gridinfo[ y + ny*(x*nz + z) ].mu[is]  = gridinfo[ y + ny*(x*nz + (z+1)) ].mu[is];
    }
    cscl[ y + ny*(x*nz + z) ] = cscl[ y + ny*(x*nz + (z+1)) ];
  }
    
}

__kernel void apply_BC_com_zn_noflux(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index = y + ny*(x*nz + z);
  
  if ( z == (nz-1) ) {
    for (is = 0; is < nsol; is++) { 
      gridinfo[ y + ny*(x*nz + z) ].com[is] = gridinfo[ y + ny*(x*nz + (z-1)) ].com[is];
      gridinfo[ y + ny*(x*nz + z) ].mu[is]  = gridinfo[ y + ny*(x*nz + (z-1)) ].mu[is];
    }
    cscl[ y + ny*(x*nz + z) ] = cscl[ y + ny*(x*nz + (z-1)) ];
  }
  
    
}

__kernel void apply_BC_com_z0_periodic(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index = y + ny*(x*nz + z);
  
  if ( z == 0 ) {
    for (is = 0; is < nsol; is++) { 
      gridinfo[ y + ny*(x*nz + z) ].com[is] = gridinfo[ y + ny*(x*nz + (nz-2)) ].com[is];
      gridinfo[ y + ny*(x*nz + z) ].mu[is]  = gridinfo[ y + ny*(x*nz + (nz-2)) ].mu[is];
    }
    cscl[ y + ny*(x*nz + z) ] = cscl[ y + ny*(x*nz + (nz-2)) ];
  }
    
}

__kernel void apply_BC_com_zn_periodic(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index = y + ny*(x*nz + z);
  
  if ( z == (nz-1) ) {
    for (is = 0; is < nsol; is++) { 
      gridinfo[ y + ny*(x*nz + z) ].com[is] = gridinfo[ y + ny*(x*nz + 1) ].com[is];
      gridinfo[ y + ny*(x*nz + z) ].mu[is]  = gridinfo[ y + ny*(x*nz + 1) ].mu[is];
    }
    cscl[ y + ny*(x*nz + z) ] = cscl[ y + ny*(x*nz + 1) ];
  }
    
}

__kernel void apply_BC_com_x0_noflux(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index = y + ny*(x*nz + z);
  
  if ( x == 0 ) {

    for ( ip = 0; ip < npha; ip++ ) {
      gridinfo[ y + ny*(x*nz + z) ].phi[ip] = gridinfo[ (y) + ny*((x+1)*nz + z) ].phi[ip];
    }

    for (is = 0; is < nsol; is++) { 
      gridinfo[ y + ny*(x*nz + z) ].com[is] = gridinfo[ (y) + ny*((x+1)*nz + z) ].com[is];
      gridinfo[ y + ny*(x*nz + z) ].mu[is]  = gridinfo[ (y) + ny*((x+1)*nz + z) ].mu[is];
    }
    cscl[ y + ny*(x*nz + z) ] = cscl[ (y) + ny*((x+1)*nz + z) ];
  }
    
}

__kernel void apply_BC_com_xn_noflux(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index = y + ny*(x*nz + z);
  
  if ( x == (nx-1) ) {

    for ( ip = 0; ip < npha; ip++ ) {
      gridinfo[ y + ny*(x*nz + z) ].phi[ip] = gridinfo[ y + ny*((x-1)*nz + z) ].phi[ip];
    }


    for (is = 0; is < nsol; is++) { 
      gridinfo[ y + ny*(x*nz + z) ].com[is] = gridinfo[ y + ny*((x-1)*nz + z) ].com[is];
      gridinfo[ y + ny*(x*nz + z) ].mu[is]  = gridinfo[ y + ny*((x-1)*nz + z) ].mu[is];
    }
    cscl[ y + ny*(x*nz + z) ] = cscl[ y + ny*((x-1)*nz + z) ];
  }
    
}

__kernel void apply_BC_com_x0_periodic(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index = y + ny*(x*nz + z);
  
  if ( x == 0 ) {
    for (is = 0; is < nsol; is++) { 
      gridinfo[ y + ny*(x*nz + z) ].com[is] = gridinfo[ y + ny*((nx-2)*nz + z) ].com[is];
      gridinfo[ y + ny*(x*nz + z) ].mu[is]  = gridinfo[ y + ny*((nx-2)*nz + z) ].mu[is];
    }
    cscl[ y + ny*(x*nz + z) ] = cscl[ y + ny*((nx-2)*nz + z) ];
  }
    
}

__kernel void apply_BC_com_xn_periodic(__global struct fields *gridinfo, __constant struct pfmval *pfmdat, __global struct csle *cscl) {
  
  int x, y, z, nx, ny, nz, index;
  int is, ip;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index = y + ny*(x*nz + z);
  
  if ( x == (nx-1) ) {
    for (is = 0; is < nsol; is++) { 
      gridinfo[ y + ny*(x*nz + z) ].com[is] = gridinfo[ y + ny*(1*nz + z) ].com[is];
      gridinfo[ y + ny*(x*nz + z) ].mu[is]  = gridinfo[ y + ny*(1*nz + z) ].mu[is];
    }
    cscl[ y + ny*(x*nz + z) ] = cscl[ y + ny*(1*nz + z) ];
  }
    
}

__kernel void addNoise(__global struct fields *gridinfo, __constant struct pfmval *pfmdat) {
  
  int i;
  int j;
  int k;
  int nx;
  int ny;
  int nz;
  int index;
  int is, js, ip;
  
  j = get_global_id(1);
  i = get_global_id(0);
  
  nx = pfmdat->Nx;
  ny = pfmdat->Ny;
  index = (i*ny + j);

  double minr = -1.0;
  double maxr = 1.0;
  double rt, noise;
  
  for ( ip = 0; ip < npha; ip++ ) { 

    if ( ( gridinfo[index].phi[ip] > 0.1 ) && ( gridinfo[index].phi[ip] < 0.5 ) ) { 
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
      
      //gridinfo[index].phi[ip] = gridinfo[index].phi[ip] - ( pfmdat->NoiseFac * noise );
      
      gridinfo[index].phi[ip] = gridinfo[index].phi[ip] - ( pfmdat->NoiseFac * noise ) * gridinfo[index].phi[ip] * ( 1.0 - gridinfo[index].phi[ip] ) * ip / (npha*npha);
      
    }
    
    if (gridinfo[index].phi[ip] < 0.0) {
      gridinfo[index].phi[ip] = 0.0;
    }
    if (gridinfo[index].phi[ip] > 1.0) {
      gridinfo[index].phi[ip] = 1.0;
    }

  }

}

__kernel void copy_New_To_Old(__global struct fields *gridinfoO, __global struct fields *gridinfo, __constant struct pfmval *pfmdat) {
  
  int x, y, z;
  int nx, ny, nz;
  int index;
  
  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);
  
  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);
  
  index = y + ny*(z + x*nz);
  
  gridinfoO[index] = gridinfo[index];

}

__kernel void update_temp_UC(__global double *temp, __global int *tstep, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar) {

    int i;
    

}

__kernel void update_temp_DS(__global double *temp, __global int *tstep, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar) {

    int i;

}
  
  
__kernel void update_temp_CR(__global double *temp, __global int *tstep, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar) { 
    
    int i;

}
  
__kernel void apply_BC_temp_it_noflux(__global double *temp, __constant struct pfmval *pfmdat) {

    int i;

}
  
__kernel void apply_BC_temp_ib_noflux(__global double *temp, __constant struct pfmval *pfmdat) {

    int i;

}


__kernel void apply_BC_ela_y0_noflux(__global struct iter_variables *it_gridinfoO, __constant struct pfmval *pfmdat, __global struct csle *cscl) {

  int x, y, z, nx, ny, nz, index;
  int is, ip;

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  index =  (y  ) + ny*( (x  )*nz + (z  ) );

  if ( y == 0 ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (3  ) + ny*( (x  )*nz + (z  ) ) ];
  }

  if ( y == 1 ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (2  ) + ny*( (x  )*nz + (z  ) ) ];
  }

}

__kernel void apply_BC_ela_yn_noflux(__global struct iter_variables *it_gridinfoO, __constant struct pfmval *pfmdat, __global struct csle *cscl) {

  int x, y, z, nx, ny, nz, index;
  int is, ip;

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  index =  (y  ) + ny*( (x  )*nz + (z  ) );

  if ( y == (ny-1) ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (ny-4  ) + ny*( (x  )*nz + (z  ) ) ];
  }

  if ( y == (ny-2) ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (ny-3  ) + ny*( (x  )*nz + (z  ) ) ];
  }

}

__kernel void apply_BC_ela_y0_periodic(__global struct iter_variables *it_gridinfoO, __constant struct pfmval *pfmdat, __global struct csle *cscl) {

  int x, y, z, nx, ny, nz, index;
  int is, ip;

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  index =  (y  ) + ny*( (x  )*nz + (z  ) );

  if ( y == 0 ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (ny-4  ) + ny*( (x  )*nz + (z  ) ) ];
  }

  if ( y == 1 ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (ny-3  ) + ny*( (x  )*nz + (z  ) ) ];
  }

}

__kernel void apply_BC_ela_yn_periodic(__global struct iter_variables *it_gridinfoO, __constant struct pfmval *pfmdat, __global struct csle *cscl) {

  int x, y, z, nx, ny, nz, index;
  int is, ip;

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  index =  (y  ) + ny*( (x  )*nz + (z  ) );

  if ( y == (ny-1) ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (3  ) + ny*( (x  )*nz + (z  ) ) ];
  }

  if ( y == (ny-2) ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (2  ) + ny*( (x  )*nz + (z  ) ) ];
  }

}

__kernel void apply_BC_ela_z0_noflux(__global struct iter_variables *it_gridinfoO, __constant struct pfmval *pfmdat, __global struct csle *cscl) {

  int x, y, z, nx, ny, nz, index;
  int is, ip;

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  index = y + ny*(x*nz + z);

  if ( z == 0 ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (3  ) ) ];
  }

  if ( z == 1 ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (2  ) ) ];
  }

}

__kernel void apply_BC_ela_zn_noflux(__global struct iter_variables *it_gridinfoO, __constant struct pfmval *pfmdat, __global struct csle *cscl) {

  int x, y, z, nx, ny, nz, index;
  int is, ip;

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  index = y + ny*(x*nz + z);

  if ( z == (nz-1) ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (nz-4  ) ) ];
  }

  if ( z == (nz-2) ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (nz-3  ) ) ];
  }

}

__kernel void apply_BC_ela_z0_periodic(__global struct iter_variables *it_gridinfoO, __constant struct pfmval *pfmdat, __global struct csle *cscl) {

  int x, y, z, nx, ny, nz, index;
  int is, ip;

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  index = y + ny*(x*nz + z);

  if ( z == 0 ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (nz-4  ) ) ];
  }

  if ( z == 1 ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (nz-3  ) ) ];
  }

}

__kernel void apply_BC_ela_zn_periodic(__global struct iter_variables *it_gridinfoO, __constant struct pfmval *pfmdat, __global struct csle *cscl) {

  int x, y, z, nx, ny, nz, index;
  int is, ip;

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  index = y + ny*(x*nz + z);

  if ( z == (nz-1) ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (3  ) ) ];
  }

  if ( z == (nz-2) ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (y  ) + ny*( (x  )*nz + (2  ) ) ];
  }

}

__kernel void apply_BC_ela_x0_noflux(__global struct iter_variables *it_gridinfoO, __constant struct pfmval *pfmdat, __global struct csle *cscl) {

  int x, y, z, nx, ny, nz, index;
  int is, ip;

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  index =  (y  ) + ny*( (x  )*nz + (z  ) );

  if ( x == 0 ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (y  ) + ny*( (3  )*nz + (z  ) ) ];
  }

  if ( x == 1 ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (y  ) + ny*( (2  )*nz + (z  ) ) ];
  }

}

__kernel void apply_BC_ela_xn_noflux(__global struct iter_variables *it_gridinfoO, __constant struct pfmval *pfmdat, __global struct csle *cscl) {

  int x, y, z, nx, ny, nz, index;
  int is, ip;

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  index =  (y  ) + ny*( (x  )*nz + (z  ) );

  if ( x == (nx-1) ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (nx  ) + ny*( (nx-4  )*nz + (z  ) ) ];
  }

  if ( x == (nx-2) ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (nx  ) + ny*( (nx-3  )*nz + (z  ) ) ];
  }

}

__kernel void apply_BC_ela_x0_periodic(__global struct iter_variables *it_gridinfoO, __constant struct pfmval *pfmdat, __global struct csle *cscl) {

  int x, y, z, nx, ny, nz, index;
  int is, ip;

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  index =  (y  ) + ny*( (x  )*nz + (z  ) );

  if ( x == 0 ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (y  ) + ny*( (nx-4  )*nz + (z  ) ) ];
  }

  if ( x == 1 ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (y  ) + ny*( (nx-3  )*nz + (z  ) ) ];
  }

}

__kernel void apply_BC_ela_xn_periodic(__global struct iter_variables *it_gridinfoO, __constant struct pfmval *pfmdat, __global struct csle *cscl) {

  int x, y, z, nx, ny, nz, index;
  int is, ip;

  y = get_global_id(0);
  z = get_global_id(1);
  x = get_global_id(2);

  ny = get_global_size(0);
  nz = get_global_size(1);
  nx = get_global_size(2);

  index =  (y  ) + ny*( (x  )*nz + (z  ) );

  if ( x == (nx-1) ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (nx  ) + ny*( (3  )*nz + (z  ) ) ];
  }

  if ( x == (nx-2) ) {
    it_gridinfoO[ (y  ) + ny*( (x  )*nz + (z  ) ) ] = it_gridinfoO[ (nx  ) + ny*( (2  )*nz + (z  ) ) ];
  }

}
