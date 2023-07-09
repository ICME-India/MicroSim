double calcAnisotropy_01(double *phi_1,
                         double *dab, double *eps_ab,
                         int phase,
                         double DELTA_X, double DELTA_Y, double DELTA_Z, double *Rot_mat, double *Inv_Rot_mat);


double calcAnisotropy_01(double *phi_1,
                         double *dab, double *eps_ab,
                         int phase,
                         double DELTA_X, double DELTA_Y, double DELTA_Z, double *Rot_mat, double *Inv_Rot_mat)
{

  double sum1, sum2, sum3;

  int i, j, k, X, ix, iy, iz, index, i1, i2;

  int ip1, ip2;

  double gradPhiMid[npha][7][3];
  double phiMid[npha][7];

  double qMid[npha][7][3], Rotated_qMid[npha][7][3];

  double ac[npha][7];
  double dadq[npha][7][3], Rotated_dadq[npha][7][3];

  double dqdphi[npha][1][3], Rotated_dqdphi[npha][1][3];

  double phi[npha][3][3][3];

  double Rot_Matrix[npha][npha][9], Inv_Rot_Matrix[npha][npha][9];

  //printf("%d,-----------\n", npha);

  ix = 0;
  iy = 1;
  iz = 2;

  for ( ip1 = 0; ip1 < npha; ip1++ ) {
    for ( ip2 = 0; ip2 < npha; ip2++ ) {
      for ( i1 = 0; i1 < 3; i1++ ) {
        for ( i2 = 0; i2 < 3; i2++ ) {
          Rot_Matrix[ip1][ip2][i1*3+i2] = Rot_mat[((ip1*npha+ip2)*3+i1)*3+i2];
          Inv_Rot_Matrix[ip1][ip2][i1*3+i2] = Inv_Rot_mat[((ip1*npha+ip2)*3+i1)*3+i2];
        }
      }
    }
  }

  for (ip1 = 0; ip1 < npha; ip1++) {
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        for (k = 0; k < 3; k++) {
          index = ((ip1*3+i)*3+j)*3+k;
          phi[ip1][i][j][k] = phi_1[index];
        }
      }
    }
  }

  for (ip1 = 0; ip1 < npha; ip1++) {
    /*
    * 0 -> i,j,k
    * 1 -> i-1/2,j,k
    * 2 -> i+1/2,j,k
    * 3 -> i,j-1/2,k
    * 4 -> i,j+1/2,k
    * 5 -> i,j,k-1/2
    * 6 -> i,j,k+1/2
    */

    phiMid[ip1][0] = phi[ip1][1][1][1];
    phiMid[ip1][1] = (phi[ip1][0][1] + phi[ip1][1][1])/2.0;
    phiMid[ip1][2] = (phi[ip1][2][1] + phi[ip1][1][1])/2.0;
    phiMid[ip1][3] = (phi[ip1][1][0] + phi[ip1][1][1])/2.0;
    phiMid[ip1][4] = (phi[ip1][1][2] + phi[ip1][1][1])/2.0;

  }

  for (ip1 = 0; ip1 < npha; ip1++) {
    /*
     * Second index:
     * 0 -> i,j,k
     * 1 -> i-1/2, j, k
     * 2 -> i+1/2, j, k
     * 3 -> i, j-1/2, k
     * 4 -> i, j+1/2, k
     * 5 -> i, j, k-1/2
     * 6 -> i, j, k+1/2
     *
     * Third index:
     * 0 -> x
     * 1 -> y
     * 2 -> z
     *
    */

    gradPhiMid[ip1][0][ix] = (phi[ip1][2][1] - phi[ip1][0][1])/(2.0*DELTA_X);
    gradPhiMid[ip1][1][ix] = (phi[ip1][1][1] - phi[ip1][0][1])/(DELTA_X);
    gradPhiMid[ip1][2][ix] = (phi[ip1][2][1] - phi[ip1][1][1])/(DELTA_X);
    gradPhiMid[ip1][3][ix] = ((phi[ip1][2][0] - phi[ip1][0][0])/(2.0*DELTA_X) + gradPhiMid[ip1][0][ix])/2.0;
    gradPhiMid[ip1][4][ix] = ((phi[ip1][2][2] - phi[ip1][0][2])/(2.0*DELTA_X) + gradPhiMid[ip1][0][ix])/2.0;

    gradPhiMid[ip1][0][iy] =  (phi[ip1][1][2] - phi[ip1][1][0])/(2.0*DELTA_Y);
    gradPhiMid[ip1][1][iy] = ((phi[ip1][0][2] - phi[ip1][0][0])/(2.0*DELTA_Y) + gradPhiMid[ip1][0][iy])/2.0;
    gradPhiMid[ip1][2][iy] = ((phi[ip1][2][2] - phi[ip1][2][0])/(2.0*DELTA_Y) + gradPhiMid[ip1][0][iy])/2.0;
    gradPhiMid[ip1][3][iy] =  (phi[ip1][1][1] - phi[ip1][1][0])/(DELTA_Y);
    gradPhiMid[ip1][4][iy] =  (phi[ip1][1][2] - phi[ip1][1][1])/(DELTA_Y);

    gradPhiMid[ip1][0][iz] = 0.0;
    gradPhiMid[ip1][1][iz] = 0.0;
    gradPhiMid[ip1][2][iz] = 0.0;
    gradPhiMid[ip1][3][iz] = 0.0;
    gradPhiMid[ip1][4][iz] = 0.0;

  }

  // for (ip1 = 0; ip1 < npha; ip1++) {
  //     for (i = 0; i < 5; i++) {
  //         for (j = 0; j < 3; j++ ) {
  //             printf("ip=%d, id=%d, ix=%d, grad=%le\n", ip1, i, j, gradPhiMid[ip1][i][j]);
  //         }
  //     }
  // }

  for (ip1 = 0; ip1 < npha; ip1++) {
    /*
     * Second index:
     * 0 -> i,j,k
     * 1 -> i-1/2, j, k
     * 2 -> i+1/2, j, k
     * 3 -> i, j-1/2, k
     * 4 -> i, j+1/2, k
     * 5 -> i, j, k-1/2
     * 6 -> i, j, k+1/2
     *
     * Third index:
     * 0 -> x
     * 1 -> y
     * 2 -> z
     *
    */

    /*
     * q_x_i-1/2, q_x_i+1/2, q_x_j-1/2, q_x_j+1/2, q_x_k-1/2, q_x_k+1/2, ... , q_z_k+1/2
     *
     * q^{\alpha\beta}_{x}_{i-1/2,j,k} = \phi_{\alpha}_{i-1/2,j,k}\nabla\phi_{\beta}_{i-1/2,j,k} - \phi_{\beta}_{i-1/2,j,k}\nabla\phi_{\alpha}_{i-1/2,j,k}
     *
    */

    if (ip1 == phase) {
      for (i = 0; i < 5; i++) {
        for (X = 0; X < 3; X++) {
          qMid[ip1][i][X] = 0.0;
        }
      }
    }
    else {
      for (i = 0; i < 5; i++) {
        for (X = 0; X < 3; X++) {
          qMid[ip1][i][X] = phiMid[phase][i]*gradPhiMid[ip1][i][X] - phiMid[ip1][i]*gradPhiMid[phase][i][X];
          qMid[ip1][i][2] = 0.0;
          //printf("alpa=%d, beta = %d, i = %d, X = %d, qMid = %le\n", phase, ip1, i, X, qMid[ip1][i][X]);
        }
      }
    }
  }

  for (ip1 = 0; ip1 < npha; ip1++) {
    if ( ip1 != phase ) {
      for (i = 0; i < 5; i++) {
        multiply_3x1(Rot_Matrix[phase][ip1], qMid[ip1][i], Rotated_qMid[ip1][i]);
      }
    }
  }

  for (ip1 = 0; ip1 < npha; ip1++) {
    for (i = 0; i < 5; i++) {
      for (X = 0; X < 3; X++) {
        qMid[ip1][i][X] = Rotated_qMid[ip1][i][X];
      }
    }
  }

  for (ip1 = 0; ip1 < npha; ip1++) {
    if (ip1 != phase) {
      for (i = 0; i < 5; i++) {
        ac[ip1][i] = anisotropy_01_function_ac(qMid[ip1][i], phase, ip1, dab);
        anisotropy_01_dAdq(qMid[ip1][i], dadq[ip1][i], phase, ip1, dab);
        //printf("%le, %le, %le, %le\n", ac[ip1][i], dadq[ip1][i][0], dadq[ip1][i][1], dadq[ip1][i][2]);
      }
    }
  }

  for (ip1 = 0; ip1 < npha; ip1++) {
    if ( ip1 != phase ) {
      for (i = 0; i < 5; i++) {
        multiply_3x1(Inv_Rot_Matrix[phase][ip1], dadq[ip1][i], Rotated_dadq[ip1][i]);
      }
    }
  }

  for (ip1 = 0; ip1 < npha; ip1++) {
    for (X = 0; X < 3; X++) {
      dqdphi[ip1][0][X] = gradPhiMid[ip1][0][X];
    }
  }

  for (ip1 = 0; ip1 < npha; ip1++) {
    if ( ip1 != phase ) {
      multiply_3x1(Rot_Matrix[phase][ip1], dqdphi[ip1][0], Rotated_dqdphi[ip1][0]);
    }
  }

  for (ip1 = 0; ip1 < npha; ip1++) {
    for (X = 0; X < 3; X++) {
      dqdphi[ip1][0][X] = Rotated_dqdphi[ip1][0][X];
    }
  }

  sum1 = 0.0;
  sum2 = 0.0;
  sum3 = 0.0;

  for (ip1 = 0; ip1 < npha; ip1++) {
    if (ip1 != phase) {

      sum1 += (2.0*eps_ab[phase*npha + ip1]*ac[ip1][2]*(Rotated_dadq[ip1][2][ix])*(-phiMid[ip1][2])*(gradPhiMid[phase][2][ix]*gradPhiMid[ip1][2][ix] + gradPhiMid[phase][2][iy]*gradPhiMid[ip1][2][iy]))/DELTA_X;
      sum1 -= (2.0*eps_ab[phase*npha + ip1]*ac[ip1][1]*(Rotated_dadq[ip1][1][ix])*(-phiMid[ip1][1])*(gradPhiMid[phase][1][ix]*gradPhiMid[ip1][1][ix] + gradPhiMid[phase][1][iy]*gradPhiMid[ip1][1][iy]))/DELTA_X;

      sum1 += (2.0*eps_ab[phase*npha + ip1]*ac[ip1][4]*(Rotated_dadq[ip1][4][iy])*(-phiMid[ip1][4])*(gradPhiMid[phase][4][ix]*gradPhiMid[ip1][4][ix] + gradPhiMid[phase][4][iy]*gradPhiMid[ip1][4][iy]))/DELTA_Y;
      sum1 -= (2.0*eps_ab[phase*npha + ip1]*ac[ip1][3]*(Rotated_dadq[ip1][3][iy])*(-phiMid[ip1][3])*(gradPhiMid[phase][3][ix]*gradPhiMid[ip1][3][ix] + gradPhiMid[phase][3][iy]*gradPhiMid[ip1][3][iy]))/DELTA_Y;

      //sum1 += (2.0*eps_ab[phase*npha + ip1]*ac[ip1][6]*(Rotated_dadq[ip1][6][iz])*(-phiMid[ip1][6])*(gradPhiMid[phase][6][ix]*gradPhiMid[ip1][6][ix] + gradPhiMid[phase][6][iy]*gradPhiMid[ip1][6][iy] + gradPhiMid[phase][6][iz]*gradPhiMid[ip1][6][iz]))/DELTA_Z;
      //sum1 -= (2.0*eps_ab[phase*npha + ip1]*ac[ip1][5]*(Rotated_dadq[ip1][5][iz])*(-phiMid[ip1][5])*(gradPhiMid[phase][5][ix]*gradPhiMid[ip1][5][ix] + gradPhiMid[phase][5][iy]*gradPhiMid[ip1][5][iy] + gradPhiMid[phase][5][iz]*gradPhiMid[ip1][5][iz]))/DELTA_Z;

      sum2 += eps_ab[phase*npha + ip1]*ac[ip1][2]*ac[ip1][2]*(gradPhiMid[ip1][2][ix])/DELTA_X;
      sum2 -= eps_ab[phase*npha + ip1]*ac[ip1][1]*ac[ip1][1]*(gradPhiMid[ip1][1][ix])/DELTA_X;

      sum2 += eps_ab[phase*npha + ip1]*ac[ip1][4]*ac[ip1][4]*(gradPhiMid[ip1][4][iy])/DELTA_Y;
      sum2 -= eps_ab[phase*npha + ip1]*ac[ip1][3]*ac[ip1][3]*(gradPhiMid[ip1][3][iy])/DELTA_Y;

      //sum2 += eps_ab[phase*npha + ip1]*ac[ip1][6]*ac[ip1][6]*(gradPhiMid[ip1][6][iz])/DELTA_Z;
      //sum2 -= eps_ab[phase*npha + ip1]*ac[ip1][5]*ac[ip1][5]*(gradPhiMid[ip1][5][iz])/DELTA_Z;

      sum3 += -2.0*eps_ab[phase*npha + ip1]*ac[ip1][0]*(dadq[ip1][0][ix]*dqdphi[ip1][0][ix] + dadq[ip1][0][iy]*dqdphi[ip1][0][iy])*(gradPhiMid[phase][0][ix]*gradPhiMid[ip1][0][ix] + gradPhiMid[phase][0][iy]*gradPhiMid[ip1][0][iy]);
      //printf("\n%le, %le, %le, %le, %le, %le, %le, %le, %le, %le\n", eps_ab[phase*npha + ip1], ac[ip1][0], dadq[ip1][0][ix], gradPhiMid[ip1][0][ix], dadq[ip1][0][iy], gradPhiMid[ip1][0][iy], gradPhiMid[phase][0][ix], gradPhiMid[ip1][0][ix], gradPhiMid[phase][0][iy], gradPhiMid[ip1][0][iy]);
        
    }

  }

  //printf("aniso=%le\n", sum1);
  return (sum1+sum2 + sum3);

}
