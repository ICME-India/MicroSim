void propf3Hostupdate(struct propmatf3 *propf3) { 
  int ip, ip1, ip2, is, is1, is2;

  printf(" In propf3Hosrupdate \n");
  
  for ( ip1 = 0; ip1 < npha; ip1 ++ ) { 
    for ( ip2 = 0; ip2 < npha; ip2++ ) { 
      for ( is = 0; is < nsol; is++ ) { 
        propf3->ceq[ip1][ip2][is]   = ceq[ip1][ip2][is]; 
        propf3->cfill[ip1][ip2][is] = cfill[ip1][ip2][is]; 
        propf3->slopes[ip1][ip2][is] = slopes[ip1][ip2][is];
        propf3->dcbdT[ip1][ip2][is] = dcbdT[ip1][ip2][is];
      }
      propf3->DELTA_T[ip1][ip2] = DELTA_T[ip1][ip2];
    }
  }
  for ( ip = 0; ip < npha; ip++ ) { 
    for ( is = 0; is < nsol; is++ ) { 
      propf3->DELTA_C[ip][is] = DELTA_C[ip][is]; 
      propf3->dcbdT_phase[ip][is] = dcbdT_phase[ip][is]; 
      propf3->B[ip][is] = B[ip][is]; 
      printf("%d, %d, B = %le\n", ip, is, propf3->B[ip][is]);
      propf3->Beq[ip][is] = Beq[ip][is]; 
      printf("%d, %d, Beq = %le\n", ip, is, propf3->Beq[ip][is]);
      propf3->dBbdT[ip][is] = dBbdT[ip][is]; 
      printf("%d, %d, dBbdT = %le\n", ip, is, propf3->dBbdT[ip][is]);
    }
    propf3->C[ip] = C[ip]; 
    printf("%d, C = %le\n", ip, propf3->C[ip]);
  }
  
  for ( ip = 0; ip < npha; ip++ ) { 
    for ( is1 = 0; is1 < nsol; is1++ ) { 
      for ( is2 = 0; is2 < nsol; is2++ ) { 
        propf3->A[ip][is1][is2] = A[ip][is1][is2];
        propf3->cmu[ip][is1][is2] = cmu[ip][is1][is2];
        propf3->muc[ip][is1][is2] = muc[ip][is1][is2];
        printf("%d, %d, %d,   A = %le\n", ip, is1, is2, propf3->A[ip][is1][is2]);
        printf("%d, %d, %d, cmu = %le\n", ip, is1, is2, propf3->cmu[ip][is1][is2]);
        printf("%d, %d, %d, muc = %le\n", ip, is1, is2, propf3->muc[ip][is1][is2]);
      }
    }
  }
}
