/***********************************************************************
 *       Eigenstrain
 ***********************************************************************/

void calc_strn(double eigen_strn[3][3], double eps11, double eps22, double eps33, double eps12, double eps13, double eps23)
{
    eigen_strn[0][0]=eps11;
    eigen_strn[1][1]=eps22;
    eigen_strn[2][2]=eps33;

    eigen_strn[0][1] = eigen_strn[1][0] = eps12;
    eigen_strn[0][2] = eigen_strn[2][0] = eps13;
    eigen_strn[1][2] = eigen_strn[2][1] = eps23;
}

/***********************************************************************
 *        (A, mubar, nubar)  ---> (C11, C12, C44)
 ***********************************************************************/
// void elast_const(double cub_el[3])
// {
//     cub_el[0] = mubar*(2.0*(2+Aniso)/(1+Aniso) - (1-4*nubar)/(1-2*nubar));
//     cub_el[1] = mubar*(2.0*Aniso/(1+Aniso) - (1-4*nubar)/(1-2*nubar));
//     cub_el[2] = mubar*(2.0*Aniso/(1+Aniso));
// }

/*************************************************************
 *       Function to translate C[6][6] to lambda[3][3][3][3] *
 *************************************************************/
void eval_lambda(double C[6][6], double lambda[3][3][3][3])
{
    int i,j,k,l,m,n;

    for(i = 0; i < 3; ++i)
        for(j = 0; j < 3; ++j)
            for(k = 0; k < 3; ++k)
                for(l = 0; l < 3; ++l)
                    lambda[i][j][k][l] = 0.0;
                m=0;
            n=0;
        for(i = 0; i < 3; ++i)
            for(j = 0; j < 3; ++j)
                for(k = 0; k < 3; ++k)
                    for(l = 0; l < 3; ++l)
                    {
                        if(i == j)
                            m = i;
                        if(k == l)
                            n = k;
                        if((i == 1 && j == 2) || (i == 2 && j == 1))
                            m = 3;
                        if((k == 1 && l == 2) || (k == 2 && l == 1))
                            n = 3;
                        if((i == 0 && j == 2) || (i == 2 && j == 0))
                            m = 4;
                        if((k == 0 && l == 2) || (k == 2 && l == 0))
                            n = 4;
                        if((i == 0 && j == 1) || (i == 1 && j == 0))
                            m = 5;
                        if((k == 0 && l == 1) || (k == 1 && l == 0))
                            n = 5;
                        lambda[i][j][k][l] = C[m][n];
                    }
}

/***********************************************************************
 * Function to calculate inverse of omega, which is normalized G_inverse
 * Khachaturyan calls G as Fourier transform of Green function of anisotropic
 * elasticity
 ***********************************************************************/
void calc_inv_omega(double kx, double ky, double lambda[3][3][3][3], double inv_omega[3][3])
{
    int i,j,k,l;
    double n[3];

    /*
     *  if(mod_k > 1.0e-7)
     *  {
     *    n[0] = kx/mod_k;
     *    n[1] = ky/mod_k;
     *    n[2] = 0.0;
}
*/
    n[0] = kx;
    n[1] = ky;
    n[2] = 0.0;

    for(i = 0; i < 3; ++i)
        for(j = 0; j < 3; ++j)
            inv_omega[i][j] = 0.0;

        for(i = 0; i < 3; ++i)
            for(j = 0; j < 3; ++j)
                for(k = 0; k < 3; ++k)
                    for(l = 0; l < 3; ++l)
                        inv_omega[i][j] += lambda[i][k][l][j]*n[k]*n[l];
}

/*******************************************************************
 * This function just inverts a three by three matrix
 *******************************************************************/
void calc_inverse(double a[3][3], double ainv[3][3])
{
    int i,j;
    double det;

    for(i = 0; i < 3; ++i)
        for(j = 0; j < 3; ++j)
            ainv[i][j] = 0.0;

    det = a[0][0]*(a[1][1]*a[2][2] - a[1][2]*a[2][1]) -
        a[0][1]*(a[1][0]*a[2][2] - a[1][2]*a[2][0]) +
        a[0][2]*(a[1][0]*a[2][1] - a[1][1]*a[2][0]);

    if(fabs(det) > 1.0e-7)
    {
        ainv[0][0] = (a[1][1]*a[2][2] - a[1][2]*a[2][1])/det;
        ainv[0][1] = (a[2][1]*a[0][2] - a[0][1]*a[2][2])/det;
        ainv[0][2] = (a[0][1]*a[1][2] - a[0][2]*a[1][1])/det;
        ainv[1][0] = (a[2][0]*a[1][2] - a[1][0]*a[2][2])/det;
        ainv[1][1] = (a[0][0]*a[2][2] - a[2][0]*a[0][2])/det;
        ainv[1][2] = (a[1][0]*a[0][2] - a[0][0]*a[1][2])/det;
        ainv[2][0] = (a[1][0]*a[2][1] - a[2][0]*a[1][1])/det;
        ainv[2][1] = (a[2][0]*a[0][1] - a[0][0]*a[2][1])/det;
        ainv[2][2] = (a[0][0]*a[1][1] - a[1][0]*a[0][1])/det;
    }
}

/************************************************************
 * This is linear elastic Hooke's law. Here, we are
 * interested in calculating the stress-like quantity calculated
 * from the eigen strain.
 ************************************************************/
void calc_stress(double eigen_strn[3][3], double lambda[3][3][3][3], double stress[3][3])
{
    int i,j,k,l;

    for(i = 0; i < 3; ++i)
        for(j = 0; j < 3; ++j)
            stress[i][j] = 0.0;

        for(i = 0; i < 3; ++i)
            for(j = 0; j < 3; ++j)
                for(k = 0; k < 3; ++k)
                    for(l = 0; l < 3; ++l)
                        stress[i][j] += lambda[i][j][k][l]*eigen_strn[k][l];
}

/****************************************************************
 * Given the components of the Fourier vector, the following function
 * evaluates B, which is the Fourier transform of the strain induced
 * energy
 ****************************************************************/
double elast(double kx, double ky, double omega[3][3], double
eigen_strn1[3][3], double eigen_strn2[3][3], double stress1[3][3], 
double stress2[3][3], double lambda[3][3][3][3])
{
    int i,j,k,l;
    double Bn, n[3];

    n[0] = kx;
    n[1] = ky;
    n[2] = 0.0;


    Bn = 0.0;

    for(i = 0; i < 3; ++i)
        for(j = 0; j < 3; ++j)
            for(k = 0; k < 3; ++k)
                for(l = 0; l < 3; ++l)
                    Bn += lambda[i][j][k][l]*eigen_strn1[i][j]*eigen_strn2[k][l]  -
                    n[i]*stress1[i][j]*omega[j][k]*stress2[k][l]*n[l];
    return(Bn);
}

void calculate_Bn(double *B, double eigen_strn1[3][3],double eigen_strn2[3][3],
                  double stress1[3][3], double stress2[3][3])
{
    int i, j;
    double omega[3][3], inv_omega[3][3], kx, ky, del_kx, del_ky;
    double cub_el[3],C[6][6],lambda[3][3][3][3];

    int nx = MESH_X, ny = MESH_Y;
    double dx = DELTA_X, dy = DELTA_Y;

    cub_el[0] = Stiffness_c[1].C11;
    cub_el[1] = Stiffness_c[1].C12;
    cub_el[2] = Stiffness_c[1].C44;

    for(i = 0; i < 6; ++i)
        for(j = 0; j < 6; ++j)
            C[i][j] = 0.0;

        /****************************************************************************
         *    The following relations are valid only for cases with atmost cubic symmetry
         ****************************************************************************/
    C[0][0] = cub_el[0];
    C[1][1] = cub_el[0];
    C[2][2] = cub_el[0];

    C[0][1] = cub_el[1];
    C[1][0] = cub_el[1];
    C[0][2] = cub_el[1];
    C[2][0] = cub_el[1];
    C[1][2] = cub_el[1];
    C[2][1] = cub_el[1];

    C[3][3] = cub_el[2];
    C[4][4] = cub_el[2];
    C[5][5] = cub_el[2];

    eval_lambda(C,lambda);

    calc_strn(eigen_strn1, eps_b11, eps_b22, eps_b33, eps_b12, eps_b13, eps_b23);
    calc_strn(eigen_strn2, eps_b11, eps_b22, eps_b33, eps_b12, eps_b13, eps_b23);

    calc_stress(eigen_strn1,lambda,stress1);

    calc_stress(eigen_strn2,lambda,stress2);

    del_kx = 2.0 * PI / ((double) nx * dx);
    del_ky = 2.0 * PI / ((double) ny * dy);

    for (i = 0; i < nx; ++i)
    {
        if (i <= nx/2)
            kx = (double) i * del_kx;
        else
            kx = del_kx * (double)(i - nx);

        for (j = 0; j < ny; ++j)
        {
            if (j <= ny/2)
                ky = (double)j * del_ky;
            else
                ky = del_ky * (double)(j - ny);
            calc_inv_omega(kx,ky,lambda,inv_omega);
            calc_inverse(inv_omega,omega);
            B[j + i*ny] = elast(kx,ky,omega,eigen_strn1,eigen_strn2,stress1,stress2,lambda);
        }
    }
}
