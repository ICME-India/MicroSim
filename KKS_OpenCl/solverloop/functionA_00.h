
double functionA_00(double *phi, double *divphi, double *WH,  double *Gamma_abc, int a);

double functionA_00(double *phi, double *divphi, double *WH,  double *Gamma_abc, int a) {
  int b,c;
  double sum=0.0;
  double phibphic;
  
  for ( ip = 0; ip < npha; ip++ ) { 
    gradx_phi[ip][0] = ( stgridO[5].phi[ip] - stgridO[3].phi[ip] ) / ( 2.0*(pfmvar->deltax) );
    gradx_phi[ip][1] = ( stgridO[8].phi[ip] - stgridO[6].phi[ip] ) / ( 2.0*(pfmvar->deltax) );
    gradx_phi[ip][2] = ( stgridO[2].phi[ip] - stgridO[0].phi[ip] ) / ( 2.0*(pfmvar->deltax) );
    
    grady_phi[ip][0] = ( stgridO[7].phi[ip] - stgridO[1].phi[ip] ) / ( 2.0*(pfmvar->deltay) );
    grady_phi[ip][1] = ( stgridO[8].phi[ip] - stgridO[2].phi[ip] ) / ( 2.0*(pfmvar->deltay) );
    grady_phi[ip][2] = ( stgridO[6].phi[ip] - stgridO[0].phi[ip] ) / ( 2.0*(pfmvar->deltay) );
    
    gradz_phi[ip][0] = 0.0;
    gradz_phi[ip][1] = 0.0;
    gradz_phi[ip][2] = 0.0;
    
  }
  
  for ( ip = 0; ip < npha; ip++ ) { 
    phix[ip][0] = gradx_phi[ip][0];
    phix[ip][1] = ( stgridO[5].phi[ip] - stgridO[4].phi[ip] ) / ( (pfmvar->deltax) );
    phix[ip][2] = ( stgridO[4].phi[ip] - stgridO[3].phi[ip] ) / ( (pfmvar->deltax) );
    phix[ip][3] = ( gradx_phi[ip][0] + gradx_phi[ip][1] ) / ( 2.0 );
    phix[ip][4] = ( gradx_phi[ip][0] + gradx_phi[ip][2] ) / ( 2.0 );
    
    phiy[ip][0] = grady_phi[ip][0];
    phiy[ip][1] = ( grady_phi[ip][0] + grady_phi[ip][1] ) / ( 2.0 );
    phiy[ip][2] = ( grady_phi[ip][0] + grady_phi[ip][2] ) / ( 2.0 );
    phiy[ip][3] = ( stgridO[7].phi[ip] - stgridO[4].phi[ip] ) / ( (pfmvar->deltay) );
    phiy[ip][4] = ( stgridO[4].phi[ip] - stgridO[1].phi[ip] ) / ( (pfmvar->deltay) );
    
    phiz[ip][0] = 0.0; 
    phiz[ip][1] = 0.0; 
    phiz[ip][2] = 0.0; 
    phiz[ip][3] = 0.0; 
    phiz[ip][4] = 0.0; 
    
  }
      
      
  //The well potential
  for (b=0; b < npha; b++) {
    if (b!=a) {
      if (fabs(divphi[b]) > 0.0) {
        sum += 2.0*WH[a*npha+b]*(phi[b]*phi[b]*phi[a]);
      }
    }
  }
  //sum *= 9.0;
  //The third order potential
  for (b=0; b < npha; b++) {
    for (c=0; c < npha; c++) {
      if (b!=a && c!=a && b < c) {
        if (fabs(divphi[b]) > 0.0 && fabs(divphi[c]) > 0.0) {
          phibphic = phi[b]*phi[c];
            sum += Gamma_abc[(a*npha+b)*npha+c]*phibphic;
        }
      }
    }
  }
  return sum;
}
