#pragma OPENCL EXTENSION cl_khr_fp64 : enable
//#include "solverloop/GibbsEnergyFunctions/Thermo.h"
//#include "solverloop/GibbsEnergyFunctions/Thermo.c"


void c_mu(double *mu, double *cp, double T, int phas, double *c_gues, int t, int iq, int jq);

void c_mu(double *mu, double *cp, double T, int phas, double *c_gues, int t, int iq, int jq) { 

  double fun[nsol], JacInv[nsol*nsol], cn[nsol], co[nsol], cas[npha*(nsol+1)], cie[npha][nsol];
  double retmu[npha*nsol], retdmu[npha*nsol*nsol], ddgedcdc[npha][nsol*nsol], dgedc[npha][nsol], tmp0, norm; 
  int count, ip, is, is1, is2;
  
  for ( ip = 0; ip < npha; ip++ ) { 
    for ( is = 0; is < nsol; is++ ) { 
      cie[ip][is] = c_gues[is];
      co[is] = c_gues[is];
      cn[is] = c_gues[is];
    }
  }

  //printf("t = %d, i = %d, j = %d, mu = %le, co = %le\n", t, iq, jq, mu[0], c_gues[0]);

  count = 0;
  do { 
    count++;

    for ( is = 0; is < nsol; is++ ) { 
      co[is] = cn[is];
    }
    
    for ( ip = 0; ip < npha; ip++ ) { 
      for ( is = 0; is < nsol; is++ ) { 
        cie[ip][is] = 0.0;
      }
      for ( is = 0; is < nsol+1; is++ ) { 
        cas[ip*(nsol+1)+is] = 0.0;
      }
    }
    
    tmp0 = 0.0;
    for ( is = 0; is < nsol; is++ ) { 
      cas[phas*(nsol+1)+is] = co[is];
      tmp0 += co[is];
    }
    cas[phas*(nsol+1)+nsol] = 1.0 - tmp0;
    
    mus(T, cas, retmu);
    for ( ip = 0; ip < npha; ip++ ) { 
      for ( is = 0; is < nsol; is++ ) { 
        dgedc[ip][is] = retmu[ip*nsol+is];
      }
    }

    for ( is = 0; is < nsol; is++ ) { 
      fun[is] = ( dgedc[phas][is] - mu[is] ) / ( 8.314 * T ); 
    }
    
    dmudcs(T, cas, retdmu);
    for ( ip = 0; ip < npha; ip++ ) { 
      for ( is1 = 0; is1 < nsol; is1++ ) { 
        for ( is2 = 0; is2 < nsol; is2++ ) { 
          ddgedcdc[ip][is1*nsol+is2] = retdmu[(ip*nsol+is1)*nsol+is2];
        }
      }
    }
    
    //ip = phas;
    //for ( ip = 0; ip < npha; ip++ ) { 
      //MatrixInversion_nsol(ddgedcdc[ip], dc_dmu[ip]);
      MatrixInversion_nsol(ddgedcdc[phas], JacInv);
    //}
    
    //ip = phas;
    for ( is1 = 0; is1 < nsol; is1++ ) { 
      tmp0 = 0.0; 
      for ( is2 = 0; is2 < nsol; is2++ ) { 
        tmp0 += JacInv[is1*nsol+is2] * fun[is2];
      }
      cn[is1] = co[is1] - tmp0;
    }

    // for ( is = 0; is < nsol; is++ ) { 
    //   cie[phas][is] = cn[is];
    //   //co[is] = cn[is];
    // }
    
    norm = 0.0;
    for ( is = 0; is < nsol; is++ ) { 
      norm += ( cn[is] - co[is] ) * ( cn[is] - co[is] ); 
    }

    if ( count > 1000000 ) { 
      printf("Too many iterations > at t = %d, i = %d, j = %d, count = %d, cg = %le, co = %le, cn = %le, init_mu = %le\n", t, iq, jq, count, c_gues[0], co[0], cn[0], mu[0]);
      break;
    }

  } while (fabs(cn[0]-co[0]) > 1e-8);
  //} while (norm > 1e-12);
  
  for ( is = 0; is < nsol; is++ ) { 
    cp[is] = cn[is];
  }
  
  //printf("fi = %d, t = %d, i = %d, j = %d, fas = %d,  count = %d, mu = %le, cg = %le, co = %le, cn = %le\n", phicker, t, iq, jq, phas, count, mu[0], c_gues[0], cx, cn);

}
