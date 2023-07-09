void CL_Initialize_domain() {

  //int restart=0;
  int x, y, z;
  long index;
  int ip, is;

  for (x=0; x<mpiparam.rows_x; x++) {
    for(z=0; z<mpiparam.rows_z; z++) {
      for (y=0; y<mpiparam.rows_y; y++) {

        index = y + mpiparam.rows_y*(z + mpiparam.rows_z*x );

        for ( ip = 0; ip < NUMPHASES; ip++ ) { 
          for ( is = 0; is < (NUMCOMPONENTS-1); is++ ) { 
            cscl[index].comie[ip][is] = ceq[ip][ip][is]; //gridinfo[index].composition[is]; //c_guess[ip][ip][is];
          }
        }

        //printf("%le ", gridinfo[index].compi[0]);

        //printf("%le, %le\n", cscl[index].comie[0][0], cscl[index].comie[1][0]);

      }
      //printf("\n");
    }
  }


}
