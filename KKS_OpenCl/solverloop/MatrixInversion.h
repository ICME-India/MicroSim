#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#include "solverloop/defines.h"
#include "solverloop/defines1.h"

void MatrixInversion_njac(double *Mat, double *InvMat); 

void MatrixInversion_nsol(double *Mat, double *InvMat); 

void MatrixInversion_njac(double *Mat, double *InvMat)  { 
    int is, js;
    int is1, il1, i1, is2, il2, i2, j1, il, jl, ill, jll, ii, kl;
    int ip;
    int ilmax, ks, jl1;
    int indx[njac];
    int n;
    int isingular;
    
    double big, sum, dum;
    double vv[njac], d;
    double z[njac];
    double tmpp;
    
          d = 1.0;
          
          for ( il = 0; il < njac; il++ ) { 
            big = 0.0;
            for ( jl = 0; jl < njac; jl++ ) { 
              tmpp=fabs(Mat[il*njac+jl]);
              if (tmpp > big) {
                big=tmpp;
              }
            }
            if (big == 0.0) {
              //nrerror("Singular matrix in routine ludcmp");
              isingular = 1;
              //printf("Singular1\n");
              //exit(1);
              //return LUdata;
            }
            vv[il]=1.0/big;
          }
          
          for (jl = 0; jl < njac; jl++) {
            for (il = 0; il < jl; il++) { 
              sum = Mat[il*njac+jl];
              for ( ks = 0; ks < il; ks++ ) {
                sum -= Mat[il*njac+ks]*Mat[ks*njac+jl];
              }
              Mat[il*njac+jl] = sum;
            }
            
            big=0.0;
            
            for (il = jl; il < njac; il++) {
              sum = Mat[il*njac+jl];
              for (ks = 0; ks < jl; ks++) {
                sum -= Mat[il*njac+ks]*Mat[ks*njac+jl];
              }
              Mat[il*njac+jl] = sum;
              dum = vv[il]*fabs(sum);
              if ( dum >= big ) {
                big = dum;
                ilmax = il;
              }
            }
            
            if (jl != ilmax) {
              for (ks = 0; ks < njac; ks++) {
                dum = Mat[ilmax*njac+ks];
                Mat[ilmax*njac+ks] = Mat[jl*njac+ks];
                Mat[jl*njac+ks] = dum;
              }
              d = -1.0*d;
              vv[ilmax] = vv[jl];
            }
          indx[jl] = ilmax;
          
          //if (Mat[j][j] == 0.0) Mat[j][j]=TINY;
          if ( Mat[jl*njac+jl] == 0.0 ) { 
            //Mat[jl*njac+jl] = TINY;
            isingular = 1;
            //printf("Singular2\n");
            //exit(1);
            //return LUdata;
          }
          
          if (jl != (njac-1)) {
            dum = 1.0/(Mat[jl*njac+jl]);
            for ( il = jl+1; il < njac; il++) { 
              Mat[il*njac+jl] *= dum;
            }
          }
        }
        
        
        for ( jl1 = 0; jl1 < njac; jl1++ ) { 
          for ( il1 = 0; il1 < njac; il1++ ) {
            z[il1] = 0.0;
          }
          z[jl1] = 1.0;
          ii=-1;
          for (il = 0;il < njac; il++) {
            ip = indx[il];
            sum = z[ip];
            z[ip] = z[il];
            if (ii != -1) {
              for (jl = ii; jl < il; jl++) { 
                sum -= Mat[il*njac+jl]*z[jl];
              }
            }
            else if (sum) {
              ii = il;
            }
            z[il] = sum;
          }
          
          for ( il = njac-1; il >= 0; il--) {
            sum = z[il];
            for ( jl = il+1; jl < njac; jl++) {
              sum -= Mat[il*njac+jl]*z[jl];
            }
            z[il] = sum/Mat[il*njac+il];
          }
          for ( il1 = 0; il1 < njac; il1++ ) { 
            InvMat[il1*njac+jl1] = z[il1];
          }
        }

}


void MatrixInversion_nsol(double *Mat, double *InvMat)  { 
    int is, js;
    int is1, il1, i1, is2, il2, i2, j1, il, jl, ill, jll, ii, kl;
    int ip;
    int ilmax, ks, jl1;
    int indx[nsol];
    int n;
    int isingular;
    
    double big, sum, dum;
    double vv[nsol], d;
    double z[nsol];
    double tmpp;
    
          d = 1.0;
          
          for ( il = 0; il < nsol; il++ ) { 
            big = 0.0;
            for ( jl = 0; jl < nsol; jl++ ) { 
              tmpp=fabs(Mat[il*nsol+jl]);
              if (tmpp > big) {
                big=tmpp;
              }
            }
            if (big == 0.0) {
              //nrerror("Singular matrix in routine ludcmp");
              isingular = 1;
              //printf("Singular1\n");
              //exit(1);
              //return LUdata;
            }
            vv[il]=1.0/big;
          }
          
          for (jl = 0; jl < nsol; jl++) {
            for (il = 0; il < jl; il++) { 
              sum = Mat[il*nsol+jl];
              for ( ks = 0; ks < il; ks++ ) {
                sum -= Mat[il*nsol+ks]*Mat[ks*nsol+jl];
              }
              Mat[il*nsol+jl] = sum;
            }
            
            big=0.0;
            
            for (il = jl; il < nsol; il++) {
              sum = Mat[il*nsol+jl];
              for (ks = 0; ks < jl; ks++) {
                sum -= Mat[il*nsol+ks]*Mat[ks*nsol+jl];
              }
              Mat[il*nsol+jl] = sum;
              dum = vv[il]*fabs(sum);
              if ( dum >= big ) {
                big = dum;
                ilmax = il;
              }
            }
            
            if (jl != ilmax) {
              for (ks = 0; ks < nsol; ks++) {
                dum = Mat[ilmax*nsol+ks];
                Mat[ilmax*nsol+ks] = Mat[jl*nsol+ks];
                Mat[jl*nsol+ks] = dum;
              }
              d = -1.0*d;
              vv[ilmax] = vv[jl];
            }
          indx[jl] = ilmax;
          
          //if (Mat[j][j] == 0.0) Mat[j][j]=TINY;
          if ( Mat[jl*nsol+jl] == 0.0 ) { 
            //Mat[jl*nsol+jl] = TINY;
            isingular = 1;
            //printf("Singular2\n");
            //exit(1);
            //return LUdata;
          }
          
          if (jl != (nsol-1)) {
            dum = 1.0/(Mat[jl*nsol+jl]);
            for ( il = jl+1; il < nsol; il++) { 
              Mat[il*nsol+jl] *= dum;
            }
          }
        }
        
        
        for ( jl1 = 0; jl1 < nsol; jl1++ ) { 
          for ( il1 = 0; il1 < nsol; il1++ ) {
            z[il1] = 0.0;
          }
          z[jl1] = 1.0;
          ii=-1;
          for (il = 0;il < nsol; il++) {
            ip = indx[il];
            sum = z[ip];
            z[ip] = z[il];
            if (ii != -1) {
              for (jl = ii; jl < il; jl++) { 
                sum -= Mat[il*nsol+jl]*z[jl];
              }
            }
            else if (sum) {
              ii = il;
            }
            z[il] = sum;
          }
          
          for ( il = nsol-1; il >= 0; il--) {
            sum = z[il];
            for ( jl = il+1; jl < nsol; jl++) {
              sum -= Mat[il*nsol+jl]*z[jl];
            }
            z[il] = sum/Mat[il*nsol+il];
          }
          for ( il1 = 0; il1 < nsol; il1++ ) { 
            InvMat[il1*nsol+jl1] = z[il1];
          }
        }

}