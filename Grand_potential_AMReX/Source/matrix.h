#ifndef MATRIX_H
#define MATRIX_H

#include <AMReX_Utility.H>
#include <Variables.h>
using namespace amrex;

/////////////////////////////////////////////////////////////////////////////////////////////////////////

                                 //This Function is not in use

/////////////////////////////////////////////////////////////////////////////////////////////////////////

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE
// void mat_inv(Array2D <Real,0,compcount-2,0,compcount-2> &mat, Array2D<Real,0,compcount-2,0,compcount-2> &matinv, int sz, int a){
    
//     if(sz==2){
//     Real det = mat(0,0)*mat(1,1) - mat(0,1)*mat(1,0);
//     matinv(0,0) = mat(1,1)/det;
//     matinv(0,1) = -mat(0,1)/det;
//     matinv(1,0) = -mat(1,0)/det;
//     matinv(1,1) = mat(0,0)/det;
//     }

//     else if(sz==1){
//         matinv(0,0) = 1/(mat(0,0));
//     }
     
// }


// void mat_inv(Array2D <Real,0,compcount-2,0,compcount-2> &mat, Vector<Vector<Vector<Real>>> &matinv, int sz, int a){
    
//     if(sz==2){
//     Real det = mat(0,0)*mat(1,1) - mat(0,1)*mat(1,0);
//     matinv[a][0][0] = mat(1,1)/det;
//     matinv[a][0][1] = -mat(0,1)/det;
//     matinv[a][1][0] = -mat(1,0)/det;
//     matinv[a][1][1] = mat(0,0)/det;
//     }

//     else if(sz==1){
//         matinv[a][0][0] = 1/(mat(0,0));
//     }
     
// }


// void mat_inv(Vector<Vector<Real>> mat, Vector<Vector<Real>> matinv, int sz){
    
//     if(sz==2){
//     Real det = mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
//     matinv[0][0] = mat[1][1]/det;
//     matinv[0][1] = -mat[0][1]/det;
//     matinv[1][0] = -mat[1][0]/det;
//     matinv[1][1] = mat[0][0]/det;
//     }

//     else if(sz==1){
//         matinv[0][0] = 1/(mat[0][0]);
//     }
     
// }

// void multiply_2D(Vector<Vector<Vector<Real>>> dif,Vector<Vector<Vector<Real>>> dcdmu,Vector<Vector<Real>> prod, int sz){
    
//     if(sz==3){
//         prod[0][0] = dcdmu[nump-1][0][0]*dif[nump-1][0][0] + dcdmu[nump-1][0][1]*dif[nump-1][1][0];
//         prod[0][1] = dcdmu[nump-1][0][1]*dif[nump-1][0][0] + dcdmu[nump-1][1][1]*dif[nump-1][0][1];
//         prod[1][0] = dcdmu[nump-1][0][0]*dif[nump-1][1][0] + dcdmu[nump-1][1][0]*dif[nump-1][1][1];
//         prod[1][1] = dcdmu[nump-1][1][0]*dif[nump-1][0][1] + dcdmu[nump-1][1][1]*dif[nump-1][1][1];
//     }

//     else if(sz==2){
//         prod[0][0] = dcdmu[nump-1][0][0]*dif[nump-1][0][0];
//     }
// }

// void multiply(Vector<Vector<Real>> inv_dcdmu,Vector<Real> deltac, Vector<Real> deltamu,int sz){

//     if(sz==2){
//         deltamu[0] = inv_dcdmu[0][0]*deltac[0] + inv_dcdmu[0][1]*deltac[1];
//         deltamu[1] = inv_dcdmu[1][0]*deltac[0] + inv_dcdmu[1][1]*deltac[1]; 
//     }

//     else if(sz==1){
//         deltamu[0] = inv_dcdmu[0][0]*deltac[0]; 
//     }
// }

#endif