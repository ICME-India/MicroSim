#include "matrix.cuh"

/*
 *  LU Decomposition
 */
int LUPDecompose(double **A, int N, double Tol, int *P)
{

    int i, j, k, imax;
    double maxA, *ptr, absA;

    for (i = 0; i <= N; i++)
        P[i] = i; //Unit permutation matrix, P[N] initialized with N

        for (i = 0; i < N; i++) {
            maxA = 0.0;
            imax = i;

            for (k = i; k < N; k++)
                if ((absA = fabs(A[k][i])) > maxA)
                {
                    maxA = absA;
                    imax = k;
                }

                if (maxA < Tol) return 0; //failure, matrix is degenerate

                if (imax != i)
                {
                    //pivoting P
                    j = P[i];
                    P[i] = P[imax];
                    P[imax] = j;

                    //pivoting rows of A
                    ptr = A[i];
                    A[i] = A[imax];
                    A[imax] = ptr;

                    //counting pivots starting from N (for determinant)
                    P[N]++;
                }

                for (j = i + 1; j < N; j++)
                {
                    A[j][i] /= A[i][i];

                    for (k = i + 1; k < N; k++)
                        A[j][k] -= A[j][i] * A[i][k];
                }
        }

        return 1;  //decomposition done
}

/*
 *  Inversion after LU-decomposition
 */
void LUPInvert(double **A, int *P, int N, double **IA)
{
    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < N; i++)
        {
            IA[i][j] = P[i] == j ? 1.0 : 0.0;

            for (int k = 0; k < i; k++)
                IA[i][j] -= A[i][k] * IA[k][j];
        }

        for (int i = N - 1; i >= 0; i--)
        {
            for (int k = i + 1; k < N; k++)
                IA[i][j] -= A[i][k] * IA[k][j];

            IA[i][j] /= A[i][i];
        }
    }
}

/*
 *  AB = C
 *  N is the dimension of all 3 matrices
 */
void matrixMultiply(double **A, double **B, double **C, int N)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            C[i][j] = 0.0;

            for (int k = 0; k < N; k++)
            {
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
}

extern __device__
int LUPDecomposeC1(double A[][MAX_NUM_COMP], long N, double Tol, int *P)
{
    long i, j, k, imax;
    double maxA, ptr, absA;

    for (i = 0; i <= N; i++)
    {
        P[i] = i; //Unit permutation matrix, P[N] initialized with N
    }

    for (i = 0; i < N; i++)
    {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A[k][i])) > maxA)
            {
                maxA = absA;
                imax = k;
            }

            if (maxA < Tol) return 0; //failure, matrix is degenerate

            if (imax != i)
            {
                //pivoting P
                j = P[i];
                P[i] = P[imax];
                P[imax] = j;

                //pivoting rows of A
                for (j = 0; j < N; j++)
                {
                    ptr = A[i][j];
                    A[i][j] = A[imax][j];
                    A[imax][j] = ptr;
                }

                //counting pivots starting from N (for determinant)
                P[N]++;
            }

            for (j = i + 1; j < N; j++)
            {
                A[j][i] /= A[i][i];

                for (k = i + 1; k < N; k++)
                    A[j][k] -= A[j][i] * A[i][k];
            }
    }

    return 1;  //decomposition done
}

extern __device__
int LUPDecomposeC2(double A[(MAX_NUM_COMP)*(MAX_NUM_COMP)], long N, double Tol, int *P)
{
    long i, j, k, imax;
    double maxA, ptr, absA;

    for (i = 0; i <= N; i++)
    {
        P[i] = i; //Unit permutation matrix, P[N] initialized with N
    }

    for (i = 0; i < N; i++)
    {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A[k*N+i])) > maxA)
            {
                maxA = absA;
                imax = k;
            }

            if (maxA < Tol) return 0; //failure, matrix is degenerate

            if (imax != i)
            {
                //pivoting P
                j = P[i];
                P[i] = P[imax];
                P[imax] = j;

                //pivoting rows of A
                for (j = 0; j < N; j++)
                {
                    ptr = A[i*N+j];
                    A[i*N+j] = A[imax*N+j];
                    A[imax*N+j] = ptr;
                }

                //counting pivots starting from N (for determinant)
                P[N]++;
            }

            for (j = i + 1; j < N; j++)
            {
                A[j*N+i] /= A[i*N+i];

                for (k = i + 1; k < N; k++)
                    A[j*N+k] -= A[j*N+i] * A[i*N+k];
            }
    }

    return 1;  //decomposition done
}

extern __device__
void LUPSolveC1(double A[][MAX_NUM_COMP], int *P, double *b, long N, double *x)
{
    for (long i = 0; i < N; i++) {
        x[i] = b[P[i]];

        for (long k = 0; k < i; k++)
            x[i] -= A[i][k] * x[k];
    }

    for (long i = N - 1; i >= 0; i--) {
        for (long k = i + 1; k < N; k++)
            x[i] -= A[i][k] * x[k];

        x[i] /= A[i][i];
    }
}

extern __device__
void LUPSolveC2(double A[(MAX_NUM_COMP)*(MAX_NUM_COMP)], int *P, double *b, long N, double *x)
{
    for (long i = 0; i < N; i++) {
        x[i] = b[P[i]];

        for (long k = 0; k < i; k++)
            x[i] -= A[i*N+k] * x[k];
    }

    for (long i = N - 1; i >= 0; i--) {
        for (long k = i + 1; k < N; k++)
            x[i] -= A[i*N+k] * x[k];

        x[i] /= A[i*N+i];
    }
}

extern __device__
void LUPInvertC1(double A[][MAX_NUM_COMP], int *P, long N, double IA[][MAX_NUM_COMP])
{
    for (long j = 0; j < N; j++) {
        for (long i = 0; i < N; i++) {
            IA[i][j] = P[i] == j ? 1.0 : 0.0;

            for (long k = 0; k < i; k++)
                IA[i][j] -= A[i][k] * IA[k][j];
        }

        for (long i = N - 1; i >= 0; i--) {
            for (long k = i + 1; k < N; k++)
                IA[i][j] -= A[i][k] * IA[k][j];

            IA[i][j] /= A[i][i];
        }
    }
}

extern __device__
void LUPInvertC2(double A[(MAX_NUM_COMP)*(MAX_NUM_COMP)], int *P, long N, double IA[][MAX_NUM_COMP])
{
    for (long j = 0; j < N; j++) {
        for (long i = 0; i < N; i++) {
            IA[i][j] = P[i] == j ? 1.0 : 0.0;

            for (long k = 0; k < i; k++)
                IA[i][j] -= A[i*N+k] * IA[k][j];
        }

        for (long i = N - 1; i >= 0; i--) {
            for (long k = i + 1; k < N; k++)
                IA[i][j] -= A[i*N+k] * IA[k][j];

            IA[i][j] /= A[i*N+i];
        }
    }
}

extern __device__
int LUPDecomposePC1(double A[][MAX_NUM_PHASE_COMP], long N, double Tol, int *P)
{
    long i, j, k, imax;
    double maxA, ptr, absA;

    for (i = 0; i <= N; i++)
    {
        P[i] = i; //Unit permutation matrix, P[N] initialized with N
    }

    for (i = 0; i < N; i++)
    {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A[k][i])) > maxA)
            {
                maxA = absA;
                imax = k;
            }

            if (maxA < Tol) return 0; //failure, matrix is degenerate

            if (imax != i)
            {
                //pivoting P
                j = P[i];
                P[i] = P[imax];
                P[imax] = j;

                //pivoting rows of A
                for (j = 0; j < N; j++)
                {
                    ptr = A[i][j];
                    A[i][j] = A[imax][j];
                    A[imax][j] = ptr;
                }

                //counting pivots starting from N (for determinant)
                P[N]++;
            }

            for (j = i + 1; j < N; j++)
            {
                A[j][i] /= A[i][i];

                for (k = i + 1; k < N; k++)
                    A[j][k] -= A[j][i] * A[i][k];
            }
    }

    return 1;  //decomposition done
}

extern __device__
int LUPDecomposePC2(double A[(MAX_NUM_PHASE_COMP)*(MAX_NUM_PHASE_COMP)], long N, double Tol, int *P)
{
    long i, j, k, imax;
    double maxA, ptr, absA;

    for (i = 0; i <= N; i++)
    {
        P[i] = i; //Unit permutation matrix, P[N] initialized with N
    }

    for (i = 0; i < N; i++)
    {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A[k*N+i])) > maxA)
            {
                maxA = absA;
                imax = k;
            }

            if (maxA < Tol) return 0; //failure, matrix is degenerate

            if (imax != i)
            {
                //pivoting P
                j = P[i];
                P[i] = P[imax];
                P[imax] = j;

                //pivoting rows of A
                for (j = 0; j < N; j++)
                {
                    ptr = A[i*N+j];
                    A[i*N+j] = A[imax*N+j];
                    A[imax*N+j] = ptr;
                }

                //counting pivots starting from N (for determinant)
                P[N]++;
            }

            for (j = i + 1; j < N; j++)
            {
                A[j*N+i] /= A[i*N+i];

                for (k = i + 1; k < N; k++)
                    A[j*N+k] -= A[j*N+i] * A[i*N+k];
            }
    }

    return 1;  //decomposition done
}

extern __device__
void LUPSolvePC1(double A[][MAX_NUM_PHASE_COMP], int *P, double *b, long N, double *x)
{
    for (long i = 0; i < N; i++) {
        x[i] = b[P[i]];

        for (long k = 0; k < i; k++)
            x[i] -= A[i][k] * x[k];
    }

    for (long i = N - 1; i >= 0; i--) {
        for (long k = i + 1; k < N; k++)
            x[i] -= A[i][k] * x[k];

        x[i] /= A[i][i];
    }
}

extern __device__
void LUPSolvePC2(double A[(MAX_NUM_PHASE_COMP)*(MAX_NUM_PHASE_COMP)], int *P, double *b, long N, double *x)
{
    for (long i = 0; i < N; i++) {
        x[i] = b[P[i]];

        for (long k = 0; k < i; k++)
            x[i] -= A[i*N+k] * x[k];
    }

    for (long i = N - 1; i >= 0; i--) {
        for (long k = i + 1; k < N; k++)
            x[i] -= A[i*N+k] * x[k];

        x[i] /= A[i*N+i];
    }
}

extern __device__
void LUPInvertPC1(double A[][MAX_NUM_PHASE_COMP], int *P, long N, double IA[][MAX_NUM_PHASE_COMP])
{
    for (long j = 0; j < N; j++) {
        for (long i = 0; i < N; i++) {
            IA[i][j] = P[i] == j ? 1.0 : 0.0;

            for (long k = 0; k < i; k++)
                IA[i][j] -= A[i][k] * IA[k][j];
        }

        for (long i = N - 1; i >= 0; i--) {
            for (long k = i + 1; k < N; k++)
                IA[i][j] -= A[i][k] * IA[k][j];

            IA[i][j] /= A[i][i];
        }
    }
}

extern __device__
void LUPInvertPC2(double A[(MAX_NUM_PHASE_COMP)*(MAX_NUM_PHASE_COMP)], int *P, long N, double IA[][MAX_NUM_PHASE_COMP])
{
    for (long j = 0; j < N; j++) {
        for (long i = 0; i < N; i++) {
            IA[i][j] = P[i] == j ? 1.0 : 0.0;

            for (long k = 0; k < i; k++)
                IA[i][j] -= A[i*N+k] * IA[k][j];
        }

        for (long i = N - 1; i >= 0; i--) {
            for (long k = i + 1; k < N; k++)
                IA[i][j] -= A[i*N+k] * IA[k][j];

            IA[i][j] /= A[i*N+i];
        }
    }
}
