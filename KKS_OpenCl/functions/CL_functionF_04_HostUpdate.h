void propf4Hostupdate(struct propmatf4 *propf4) { 
  int ip, ip1, ip2, is, is1, is2;

  printf(" In propf4Hosrupdate \n");
  
  for ( ip1 = 0; ip1 < NUMPHASES; ip1 ++ ) { 
    for ( ip2 = 0; ip2 < NUMPHASES; ip2++ ) { 
      for ( is = 0; is < (NUMCOMPONENTS-1); is++ ) { 
        propf4->ceq[ip1][ip2][is]   = ceq[ip1][ip2][is]; 
        propf4->cfill[ip1][ip2][is] = cfill[ip1][ip2][is]; 
        propf4->slopes[ip1][ip2][is] = slopes[ip1][ip2][is];
        propf4->dcbdT[ip1][ip2][is] = dcbdT[ip1][ip2][is];
      }
      propf4->DELTA_T[ip1][ip2] = DELTA_T[ip1][ip2];
    }
  }
  for ( ip = 0; ip < NUMPHASES; ip++ ) { 
    for ( is = 0; is < (NUMCOMPONENTS-1); is++ ) { 
      propf4->DELTA_C[ip][is] = DELTA_C[ip][is]; 
      propf4->dcbdT_phase[ip][is] = dcbdT_phase[ip][is]; 
      propf4->B[ip][is] = B[ip][is]; 
      printf("%d, %d, B = %le\n", ip, is, propf4->B[ip][is]);
      propf4->Beq[ip][is] = Beq[ip][is]; 
      //printf("%d, %d, Beq = %le\n", ip, is, propf4->Beq[ip][is]);
      propf4->dBbdT[ip][is] = dBbdT[ip][is]; 
      //printf("%d, %d, dBbdT = %le\n", ip, is, propf4->dBbdT[ip][is]);
    }
    propf4->C[ip] = C[ip]; 
    printf("%d, C = %le\n", ip, propf4->C[ip]);
  }
  
  for ( ip = 0; ip < NUMPHASES; ip++ ) { 
    for ( is1 = 0; is1 < (NUMCOMPONENTS-1); is1++ ) { 
      for ( is2 = 0; is2 < (NUMCOMPONENTS-1); is2++ ) { 
        propf4->A[ip][is1][is2] = A[ip][is1][is2];
        propf4->cmu[ip][is1][is2] = cmu[ip][is1][is2];
        propf4->muc[ip][is1][is2] = muc[ip][is1][is2];
        printf("%d, %d, %d,   A = %le\n", ip, is1, is2, propf4->A[ip][is1][is2]);
        printf("%d, %d, %d, cmu = %le\n", ip, is1, is2, propf4->cmu[ip][is1][is2]);
        printf("%d, %d, %d, muc = %le\n", ip, is1, is2, propf4->muc[ip][is1][is2]);
      }
    }
  }

  printf(" Exit propf4Hosrupdate \n");

}
