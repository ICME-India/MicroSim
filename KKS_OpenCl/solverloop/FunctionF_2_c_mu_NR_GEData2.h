#pragma OPENCL EXTENSION cl_khr_fp64 : enable

void c_mu(double *mu, double *cp, double T, int phas, double *c_gues, int t, int yq, int zq, int xq, int tphas);

void c_mu(double *mu, double *cp, double T, int phas, double *c_gues, int t, int yq, int zq, int xq, int tphas) {

  double fun[nsol], JacInv[nsol*nsol], cn[nsol], co[nsol], cas[npha*(nsol+1)], cie[npha][nsol];
  double retmu[npha*nsol], retdmu[npha*nsol*nsol], ddgedcdc[npha][nsol*nsol], dgedc[npha][nsol], tmp0, norm; 
  double retmuphase[nsol], retdmuphase[nsol*nsol], y[nsol+1], mu0[nsol]; 
  int count, ip, is, is1, is2;
  
  tmp0 = 0.0;
  for ( is = 0; is < nsol; is++ ) { 
    cn[is] = c_gues[is];
    co[is] = c_gues[is];
    y[is] = c_gues[is];
    tmp0 += y[is];
  }
  y[nsol] = 1.0 -tmp0;

  //printf("t = %d, i = %d, j = %d, mu = %le, co = %le\n", t, iq, jq, mu[0], c_gues[0]);

  count = 0;
  do { 
    count++;

    tmp0 = 0.0;
    for ( is = 0; is < nsol; is++ ) { 
      co[is] = cn[is];
      y[is] = co[is];
      tmp0 += co[is];
    }
    y[nsol] = 1.0 - tmp0;
    
    Mu(T, y, mu0, tphas);

    for ( is = 0; is < nsol; is++ ) { 
      fun[is] = ( mu0[is] - mu[is] );// / ( 8.3142 * T ); 
    }
    
    dMudc(T, y, retdmuphase, tphas);
    MatrixInversion_nsol(retdmuphase, JacInv);
    
    for ( is1 = 0; is1 < nsol; is1++ ) { 
      tmp0 = 0.0; 
      for ( is2 = 0; is2 < nsol; is2++ ) { 
        tmp0 += JacInv[is1*nsol+is2] * fun[is2];// * ( 8.3142 * T );
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

    if ( count > 1000 ) { 
      printf("Too many iterations > at t = %d, y = %d, z = %d, x = %d, count = %d, cg = %le, co = %le, cn = %le, init_mu = %le\n", t, yq, zq, xq, count, c_gues[0], co[0], cn[0], mu[0]);
      break;
    }

    //printf("t = %d, i = %d, j = %d, fas = %d,  count = %d, mu = %le, cg = %le, co = %le, cn = %le, cdiff = %le, mu0 = %le, dmu = %le\n", t, iq, jq, tphas, count, mu[0], c_gues[0], co[0], cn[0], fabs(co[0]-cn[0]), mu0[0], retdmuphase[0]);

  //} while (fabs(cn[0]-co[0]) > 1e-6);
  } while (norm > 1e-12);
  
  for ( is = 0; is < nsol; is++ ) { 
    cp[is] = cn[is];
    if (cn[is] > 1 || cn[is] < 0) {
      printf("Out of bounds > at t = %d, y = %d, z = %d, x = %d, count = %d, cg = %le, co = %le, cn = %le, init_mu = %le\n", t, yq, zq, xq, count, c_gues[0], co[0], cn[0], mu[0]);
    }
  }
  
  // printf("t = %d, i = %d, j = %d, fas = %d,  count = %d, mu = %le, cg = %le, co = %le, cn = %le\n", t, iq, jq, tphas, count, mu[0], c_gues[0], co[0], cn[0]);

}
