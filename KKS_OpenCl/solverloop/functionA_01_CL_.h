double calcAnisotropy_01(double *phi_1,
                         double *dab, double *eps_ab,
                         long phase,
                         double DELTA_X, double DELTA_Y, double DELTA_Z);


double calcAnisotropy_01(double *phi_1,
                         double *dab, double *eps_ab,
                         long phase,
                         double DELTA_X, double DELTA_Y, double DELTA_Z)
{
    double sum1, sum2, sum3;

    long i, j, k, X, ix, iy, iz, index;

    long ip1;

    double gradPhiMid[npha][7][3];
    double phiMid[npha][7];

    double qMid[npha][7][3];

    double ac[npha][7];
    double dadq[npha][7][3];

    double phi[npha][3][3][3];

    ix = 0; 
    iy = 1; 
    iz = 2;

    for (ip1 = 0; ip1 < npha; ip1++) { 
        for (i = 0; i < 3; i++) { 
            for (j = 0; j < 3; j++) { 
                for (k = 0; k < 3; k++) { 
                    index = ip1*npha*3*3 + i*3*3 + j*3 + k; 
                    phi[ip1][i][j][k] = phi_1[index];
                }
            }
        }
    }

    for (ip1 = 0; ip1 < npha; ip1++)
    {
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
        phiMid[ip1][1] = (phi[ip1][0][1][1] + phi[ip1][1][1][1])/2.0;
        phiMid[ip1][2] = (phi[ip1][2][1][1] + phi[ip1][1][1][1])/2.0;
        phiMid[ip1][3] = (phi[ip1][1][0][1] + phi[ip1][1][1][1])/2.0;
        phiMid[ip1][4] = (phi[ip1][1][2][1] + phi[ip1][1][1][1])/2.0;
        phiMid[ip1][5] = 0.0;//(phi[ip1][1][1][0] + phi[ip1][1][1][1])/2.0;
        phiMid[ip1][6] = 0.0;//(phi[ip1][1][1][2] + phi[ip1][1][1][1])/2.0;
    }

    for (ip1 = 0; ip1 < npha; ip1++)
    {
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

        gradPhiMid[ip1][0][ix] = (phi[ip1][2][1][1] - phi[ip1][0][1][1])/(2.0*DELTA_X);
        gradPhiMid[ip1][1][ix] = (phi[ip1][1][1][1] - phi[ip1][0][1][1])/(DELTA_X);
        gradPhiMid[ip1][2][ix] = (phi[ip1][2][1][1] - phi[ip1][1][1][1])/(DELTA_X);
        gradPhiMid[ip1][3][ix] = ((phi[ip1][2][0][1] - phi[ip1][0][0][1])/(2.0*DELTA_X) + gradPhiMid[ip1][0][ix])/2.0;
        gradPhiMid[ip1][4][ix] = ((phi[ip1][2][2][1] - phi[ip1][0][2][1])/(2.0*DELTA_X) + gradPhiMid[ip1][0][ix])/2.0;
        gradPhiMid[ip1][5][ix] = 0.0;//((phi[ip1][2][1][0] - phi[ip1][0][1][0])/(2.0*DELTA_X) + gradPhiMid[ip1][0][ix])/2.0;
        gradPhiMid[ip1][6][ix] = 0.0;//((phi[ip1][2][1][2] - phi[ip1][0][1][2])/(2.0*DELTA_X) + gradPhiMid[ip1][0][ix])/2.0;

        gradPhiMid[ip1][0][iy] = (phi[ip1][1][2][1] - phi[ip1][1][0][1])/(2.0*DELTA_Y);
        gradPhiMid[ip1][1][iy] = ((phi[ip1][0][2][1] - phi[ip1][0][0][1])/(2.0*DELTA_Y) + gradPhiMid[ip1][0][iy])/2.0;
        gradPhiMid[ip1][2][iy] = ((phi[ip1][2][2][1] - phi[ip1][2][0][1])/(2.0*DELTA_Y) + gradPhiMid[ip1][0][iy])/2.0;
        gradPhiMid[ip1][3][iy] = (phi[ip1][1][1][1] - phi[ip1][1][0][1])/(DELTA_Y);
        gradPhiMid[ip1][4][iy] = (phi[ip1][1][2][1] - phi[ip1][1][1][1])/(DELTA_Y);
        gradPhiMid[ip1][5][iy] = 0.0;//((phi[ip1][1][2][0] - phi[ip1][1][0][0])/(2.0*DELTA_Y) + gradPhiMid[ip1][0][iy])/2.0;
        gradPhiMid[ip1][6][iy] = 0.0;//((phi[ip1][1][2][2] - phi[ip1][1][0][2])/(2.0*DELTA_Y) + gradPhiMid[ip1][0][iy])/2.0;

        gradPhiMid[ip1][0][iz] = 0.0;//(phi[ip1][1][1][2] - phi[ip1][1][1][0])/(2.0*DELTA_Z);
        gradPhiMid[ip1][1][iz] = 0.0;//((phi[ip1][0][1][2] - phi[ip1][0][1][0])/(2.0*DELTA_Z) + gradPhiMid[ip1][0][iz])/2.0;
        gradPhiMid[ip1][2][iz] = 0.0;//((phi[ip1][2][1][2] - phi[ip1][2][1][0])/(2.0*DELTA_Z) + gradPhiMid[ip1][0][iz])/2.0;
        gradPhiMid[ip1][3][iz] = 0.0;//((phi[ip1][1][0][2] - phi[ip1][1][0][0])/(2.0*DELTA_Z) + gradPhiMid[ip1][0][iz])/2.0;
        gradPhiMid[ip1][4][iz] = 0.0;//((phi[ip1][1][2][2] - phi[ip1][1][2][0])/(2.0*DELTA_Z) + gradPhiMid[ip1][0][iz])/2.0;
        gradPhiMid[ip1][5][iz] = 0.0;//(phi[ip1][1][1][1] - phi[ip1][1][1][0])/(DELTA_Z);
        gradPhiMid[ip1][6][iz] = 0.0;//(phi[ip1][1][1][2] - phi[ip1][1][1][1])/(DELTA_Z);
    }

    for (ip1 = 0; ip1 < npha; ip1++)
    {
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

        if (ip1 == phase)
        {
            for (i = 0; i < 7; i++)
            {
                for (X = 0; X < 3; X++)
                {
                    qMid[ip1][i][X] = 0.0;
                }
            }
        }
        else
        {
            for (i = 0; i < 7; i++)
            {
                for (X = 0; X < 3; X++)
                {
                    qMid[ip1][i][X] = phiMid[phase][i]*gradPhiMid[ip1][i][X] - phiMid[ip1][i]*gradPhiMid[phase][i][X];
                }
            }
        }
    }

    for (ip1 = 0; ip1 < npha; ip1++)
    {
        if (ip1 == phase)
            continue;

        for (i = 0; i < 7; i++)
        {
            ac[ip1][i] = anisotropy_01_function_ac(qMid[ip1][i], phase, ip1, dab);
            anisotropy_01_dAdq(qMid[ip1][i], dadq[ip1][i], phase, ip1, dab);
        }
    }

    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;

    for (ip1 = 0; ip1 < npha; ip1++)
    {
        if (ip1 == phase)
            continue;

        sum1 += (2.0*eps_ab[phase*npha + ip1]*ac[ip1][2]*(dadq[ip1][2][ix])*(-phiMid[ip1][2])*(gradPhiMid[phase][2][ix]*gradPhiMid[ip1][2][ix] + gradPhiMid[phase][2][iy]*gradPhiMid[ip1][2][iy] + gradPhiMid[phase][2][iz]*gradPhiMid[ip1][2][iz]))/DELTA_X;
        sum1 -= (2.0*eps_ab[phase*npha + ip1]*ac[ip1][1]*(dadq[ip1][1][ix])*(-phiMid[ip1][1])*(gradPhiMid[phase][1][ix]*gradPhiMid[ip1][1][ix] + gradPhiMid[phase][1][iy]*gradPhiMid[ip1][1][iy] + gradPhiMid[phase][1][iz]*gradPhiMid[ip1][1][iz]))/DELTA_X;

        sum1 += (2.0*eps_ab[phase*npha + ip1]*ac[ip1][4]*(dadq[ip1][4][iy])*(-phiMid[ip1][4])*(gradPhiMid[phase][4][ix]*gradPhiMid[ip1][4][ix] + gradPhiMid[phase][4][iy]*gradPhiMid[ip1][4][iy] + gradPhiMid[phase][4][iz]*gradPhiMid[ip1][4][iz]))/DELTA_Y;
        sum1 -= (2.0*eps_ab[phase*npha + ip1]*ac[ip1][3]*(dadq[ip1][3][iy])*(-phiMid[ip1][3])*(gradPhiMid[phase][3][ix]*gradPhiMid[ip1][3][ix] + gradPhiMid[phase][3][iy]*gradPhiMid[ip1][3][iy] + gradPhiMid[phase][3][iz]*gradPhiMid[ip1][3][iz]))/DELTA_Y;

        //sum1 += (2.0*eps_ab[phase*npha + ip1]*ac[ip1][6]*(dadq[ip1][6][iz])*(-phiMid[ip1][6])*(gradPhiMid[phase][6][ix]*gradPhiMid[ip1][6][ix] + gradPhiMid[phase][6][iy]*gradPhiMid[ip1][6][iy] + gradPhiMid[phase][6][iz]*gradPhiMid[ip1][6][iz]))/DELTA_Z;
        //sum1 -= (2.0*eps_ab[phase*npha + ip1]*ac[ip1][5]*(dadq[ip1][5][iz])*(-phiMid[ip1][5])*(gradPhiMid[phase][5][ix]*gradPhiMid[ip1][5][ix] + gradPhiMid[phase][5][iy]*gradPhiMid[ip1][5][iy] + gradPhiMid[phase][5][iz]*gradPhiMid[ip1][5][iz]))/DELTA_Z;

        sum2 += eps_ab[phase*npha + ip1]*ac[ip1][2]*ac[ip1][2]*(gradPhiMid[ip1][2][ix])/DELTA_X;
        sum2 -= eps_ab[phase*npha + ip1]*ac[ip1][1]*ac[ip1][1]*(gradPhiMid[ip1][1][ix])/DELTA_X;

        sum2 += eps_ab[phase*npha + ip1]*ac[ip1][4]*ac[ip1][4]*(gradPhiMid[ip1][4][iy])/DELTA_Y;
        sum2 -= eps_ab[phase*npha + ip1]*ac[ip1][3]*ac[ip1][3]*(gradPhiMid[ip1][3][iy])/DELTA_Y;

        //sum2 += eps_ab[phase*npha + ip1]*ac[ip1][6]*ac[ip1][6]*(gradPhiMid[ip1][6][iz])/DELTA_Z;
        //sum2 -= eps_ab[phase*npha + ip1]*ac[ip1][5]*ac[ip1][5]*(gradPhiMid[ip1][5][iz])/DELTA_Z;

        sum3 += -2.0*eps_ab[phase*npha + ip1]*ac[ip1][0]*(dadq[ip1][0][ix]*gradPhiMid[ip1][0][ix] + dadq[ip1][0][iy]*gradPhiMid[ip1][0][iy] + dadq[ip1][0][iz]*gradPhiMid[ip1][0][iz])*(gradPhiMid[phase][0][ix]*gradPhiMid[ip1][0][ix] + gradPhiMid[phase][0][iy]*gradPhiMid[ip1][0][iy] + gradPhiMid[phase][0][iz]*gradPhiMid[ip1][0][iz]);

    }

    return sum1+sum2 + sum3;
}
