#include "matrix.cuh"

//Sum of two vectors
void vectorsum(double *y1,double *y2, double *sum, long size) {
	int j;
	for(j=0;j < (size);j++)
               sum[j]=y1[j]+y2[j];
}
//Multipliction of matrix and vector
void multiply(double **inv,double *y,double *prod,long size) {
	int i,j;
	double sum;
	for(i=0;i<size;i++)
	{
		sum=0;
		for(j=0;j<size;j++)
			sum=sum+inv[i][j]*y[j];
		prod[i]=sum;
	}

}
//Multiplication of two matrices
void multiply2d(double **m1,double **m2,double **prod,long size) {
	int i,j,k;
	double sum;
	for(k=0;k<size;k++)
	{
		for(i=0;i<size;i++)
		{
			sum=0;
			for(j=0;j<size;j++)
				sum=sum+m1[k][j]*m2[j][i];
			prod[k][i]=sum;
		}
	}
}
//Matrix inversion using LU decomposition
void matinvnew(double **coeffmatrix,double **inv,long size) {
	int i,j,k,tag[size];
	double **factor,**iden, **inv1, **prod;
	double *vec1,*vec,fact;
	factor = MallocM(size,size);
	inv1   = MallocM(size,size);
	iden   = MallocM(size,size);
	prod   = MallocM(size,size);
	vec1   = MallocV(size);
	vec    = MallocV(size);
	//Making the Upper Triangular Matrix.
	for(k=0; k <size; k++)
		tag[k]=k;
	for(k=0; k < size; k++) {
		pivot(coeffmatrix,factor,k,tag,size);
		for(i=k+1; i < size;i++)
		{
			fact=-coeffmatrix[i][k]/coeffmatrix[k][k];
			factor[i][k]=-fact;
			for(j=k;j<=(size-1);j++)
				coeffmatrix[i][j]=fact*coeffmatrix[k][j]+coeffmatrix[i][j];
		}
	}
	for(i=0;i < size; i++) {
		for(j=0;j<size;j++)
		{
			if(i==j)
				factor[i][j]=1;
			if(j>i)
				factor[i][j]=0;
		}
	}
	/*************************************************************************/
	//The Identity Matrix.
	for(i=0;i<(size);i++)
	{
		for(j=0;j<(size);j++)
		{
			if(i==j)
				iden[i][j]=1;
			else
				iden[i][j]=0;
		}
	}
	/*************************************************************************/
	//Forward and backward substitution to get the final identity matrix.
	for(i=0;i<(size);i++)
	{
		substitutef(factor,iden,i,vec1,size);
		substituteb(coeffmatrix,vec1,vec,size);
		for(j=0;j<(size);j++)
			inv1[j][i]=vec[j];
	}
	/**************************************************************************/
	colswap(inv1,inv,tag,size);
	multiply2d(factor,coeffmatrix,prod,size);
	rowswap(prod,coeffmatrix,tag,size);

	FreeM(factor,size);
	FreeM(iden,size);
	FreeM(inv1,size);
	FreeM(prod,size);
	free(vec1);
	free(vec);
}
/***************************************************************************************/
//Back Substitution.
void substituteb(double **fac,double *y,double *vec,long size)
{
	int i,j;
	double sum;
	vec[size-1]=y[size-1]*pow(fac[size-1][size-1],-1);
	for(i=(size-2);i>=0;i--)
	{
		sum=0;
		for(j=i+1;j<(size);j++)
			sum=sum-fac[i][j]*vec[j];
		vec[i]=(y[i]+sum)*pow(fac[i][i],-1);
	}
}
/*********************************************************************************************/
//Forward Substitution.
void substitutef(double **fac,double **y1,int index,double *vec,long size)
{
	int i,j;
	double d[size],sum;

	for(i=0;i<size;i++)
		d[i]=y1[i][index];
  	vec[0]=d[0];
	for(i=1;i<size;i++)
	{
		sum=0;
		for(j=0;j<i;j++)
			sum=sum-fac[i][j]*vec[j];
		vec[i]=d[i]+sum;
	}
// 	for(i=0;i<(size);i++)
// 		newmat[i][index]=vec[i];
}
/********************************************************************************************/
//Modulus operator
double mod(double k)
{
	if(k<0)
		return(-k);
	else
		return(k);
}
/********************************************************************************************/
//Pivoting.
void pivot(double **coeffmatrix,double **factor,int k,int *tag,long size)
{
	double swap,big;
	int tagswap,i,tag1;
	big=mod(coeffmatrix[k][k]);
	tag1=k;
	for(i=k+1;i<(size);i++)
	{
		if(mod(coeffmatrix[i][k])>big)
		{
			tag1=i;
			big=coeffmatrix[i][k];
		}
	}
	tagswap=tag[k];
	tag[k]=tag[tag1];
	tag[tag1]=tagswap;

	for(i=0;i<(size);i++)
	{
		swap=coeffmatrix[k][i];
		coeffmatrix[k][i]=coeffmatrix[tag1][i];
		coeffmatrix[tag1][i]=swap;
	}
	for(i=0;i<k;i++)
	{
		swap=factor[k][i];
		factor[k][i]=factor[tag1][i];
		factor[tag1][i]=swap;
	}
}
/*******************************************************************************************/			//Swapping Coloumns To get the final identity matrix because of the initial swappping for pivoting.
void colswap(double **m1,double **m2,int *tag, long size)
{
	int j,k,p;
	for(k=0;k < size;k++)
	{
		for(j=0;j<size;j++)
		{
			for(p=0;p <size; p++)
				m2[p][tag[j]]=m1[p][j];
		}
	}
}
/********************************************************************************************/
//Switching rows
void rowswap(double **m1,double **m2,int *tag,long size)
{
	int j,k,p;
	for(k=0;k<(size);k++)
	{
		for(j=0;j< size ;j++)
		{
			for(p=0; p < size;p++)
				m2[tag[j]][p]=m1[j][p];
		}
	}
}
/********************************************************************************************/

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

/*
 * N -> Dimension
 * AB = C
 * A is an N*N*3*3 matrix, and B is an N-element vector
 * C is of length N
 */
extern __device__
void multiply(double *A, double *B, double *C, long ip1, long ip2, long NUMPHASES, int N)
{
    int i, j;

    double sum;

    for (i = 0; i < N; i++)
    {
        sum = 0.0;
        for (j = 0; j < N; j++)
        {
            sum += A[((ip1*NUMPHASES + ip2)*3 + i)*3 + j]*B[j];
        }
        C[i] = sum;
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
