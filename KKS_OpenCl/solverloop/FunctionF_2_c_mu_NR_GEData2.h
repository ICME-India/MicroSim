#pragma OPENCL EXTENSION cl_khr_fp64 : enable

void c_mu(double *mu, double *cp, double T, int phas, double *c_gues, int t, int iq, int jq, int tphas);

void c_mu(double *mu, double *cp, double T, int phas, double *c_gues, int t, int iq, int jq, int tphas) { 

  double fun[nsol], JacInv[nsol*nsol], cn[nsol], co[nsol], cas[npha*(nsol+1)], cie[npha][nsol];
  double tmp0, norm; 
  double retmuphase[nsol], retdmuphase[nsol*nsol], y[nsol+1], mu0[nsol]; 
  int count, ip, is, is1, is2;
  

    cn[0] = c_gues[0];
    co[0] = c_gues[0];
     y[0] = c_gues[0];
    
    y[1] = 1.0 - y[0];

  //printf("t = %d, i = %d, j = %d, mu = %le, co = %le\n", t, iq, jq, mu[0], c_gues[0]);

  count = 0;
  do { 
    count++;
    
    co[0] = cn[0];
    y[0] = co[0];
    y[1] = 1.0 - y[0];
    
    Mu(T, y, mu0, tphas);
    
    fun[0] = ( mu0[0] - mu[0] );
      
    dMudc(T, y, retdmuphase, tphas);
    
    cn[0] = co[0] - fun[0] / retdmuphase[0];

    if ( count > 1000 ) { 
      //printf("Too many iterations > at t = %d, i = %d, j = %d, count = %d, cg = %le, co = %le, cn = %le, init_mu = %le, %le, %le, %d, %d\n", t, iq, jq, count, c_gues[0], co[0], cn[0], mu[0], mu0[0], retdmuphase[0], phas, tphas);
      break;
    }

  } while (fabs(cn[0]-co[0]) > 1e-8);
  
  cp[0] = cn[0];
  
  
  // printf("t = %d, i = %d, j = %d, fas = %d,  count = %d, mu = %le, cg = %le, co = %le, cn = %le\n", t, iq, jq, tphas, count, mu[0], c_gues[0], co[0], cn[0]);

}
