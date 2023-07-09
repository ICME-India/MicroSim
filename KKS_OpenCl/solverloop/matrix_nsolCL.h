//Matrix operations
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
// #include "solverloop/defines.h"
// #include "solverloop/defines1.h"
// #include <stdio.h>
// #include <math.h>
// #include "functions1.h"
// #include "open.h"
void matinvnew_nsol(double *coeffmatrix,double *inv);
void multiply_nsol(double *inv,double *y,double *prod);
void multiply_3x1(double *inv,double *y,double *prod);
void multiply2d_nsol(double *m1,double *m2,double *prod);
void multiply2d_3x3(double *m1,double *m2,double *prod);
void vectorsum_nsol(double *y1,double *y2, double *sum);
void substituteb_nsol(double *fac,double *y,double *vec);
void substitutef_nsol(double *fac,double *y1,int index,double *vec);
void pivot_nsol(double *coeffmatrix,double *factor,int k,int *tag);
void colswap_nsol(double *m1,double *m2,int *tag);
void rowswap_nsol(double *m1,double *m2,int *tag);
//Sum of two vectors
void vectorsum_nsol(double *y1,double *y2, double *sum) {
	int j;
	for(j=0;j < (nsol);j++)
               sum[j]=y1[j]+y2[j];
}
//Multipliction of matrix and vector
void multiply_nsol(double *inv,double *y,double *prod) {
	int i,j,k;
	double sum;
	for(i=0;i<nsol;i++)
	{
		sum=0;
		for(j=0;j<nsol;j++)
			sum=sum+inv[i*nsol+j]*y[j];
		prod[i]=sum;
	}

}
//Multipliction of matrix and vector
void multiply_3x1(double *inv,double *y,double *prod) {
	int i,j,k;
	double sum;
	int vsize=3;
	for(i=0;i<vsize;i++)
	{
		sum=0;
		for(j=0;j<vsize;j++)
			sum=sum+inv[i*vsize+j]*y[j];
		prod[i]=sum;
	}

}
//Multiplication of two matrices
void multiply2d_nsol(double *m1,double *m2,double *prod) {
	int i,j,k;
	double sum;
	for(k=0;k<nsol;k++)
	{
		for(i=0;i<nsol;i++)
		{
			sum=0;
			for(j=0;j<nsol;j++)
				sum=sum+m1[k*nsol+j]*m2[j*nsol+i];
			prod[k*nsol+i]=sum;
		}
	}
}
//Multiplication of two matrices
void multiply2d_3x3(double *m1,double *m2,double *prod) {
	int i,j,k;
	double sum;
	int matsize=3;
	for(k=0;k<matsize;k++)
	{
		for(i=0;i<matsize;i++)
		{
			sum=0;
			for(j=0;j<matsize;j++)
				sum=sum+m1[k*matsize+j]*m2[j*matsize+i];
			prod[k*matsize+i]=sum;
		}
	}
}
//Matrix inversion using LU decomposition
void matinvnew_nsol(double *coeffmatrix,double *inv) {
	int i,j,k,p,q,tag[nsol];
	double factor[nsol*nsol],iden[nsol*nsol], inv1[nsol*nsol], prod[nsol*nsol];
	double vec1[nsol],vec[nsol],fact;
	
	//Making the Upper Triangular Matrix.
	for(k=0; k <nsol; k++)
		tag[k]=k;
	for(k=0; k < nsol; k++) {
		pivot_nsol(coeffmatrix,factor,k,tag);
		for(i=k+1; i < nsol;i++)
		{
			fact=-coeffmatrix[i*nsol+k]/coeffmatrix[k*nsol+k];
			factor[i*nsol+k]=-fact;
			for(j=k;j<=(nsol-1);j++)
				coeffmatrix[i*nsol+j]=fact*coeffmatrix[k*nsol+j]+coeffmatrix[i*nsol+j];
		}
	}
	for(i=0;i < nsol; i++) {
		for(j=0;j<nsol;j++)
		{
			if(i==j)
				factor[i*nsol+j]=1;
			if(j>i)
				factor[i*nsol+j]=0;
		}
	}
	/*************************************************************************/
	//The Identity Matrix.
	for(i=0;i<(nsol);i++)
	{
		for(j=0;j<(nsol);j++)
		{
			if(i==j)
				iden[i*nsol+j]=1;
			else
				iden[i*nsol+j]=0;
		}
	}
	/*************************************************************************/
	//Forward and backward substitution to get the final identity matrix. 
	for(i=0;i<(nsol);i++)
	{
		substitutef_nsol(factor,iden,i,vec1);
		substituteb_nsol(coeffmatrix,vec1,vec);
		for(j=0;j<(nsol);j++)
			inv1[j*nsol+i]=vec[j];
	}
	/**************************************************************************/
	colswap_nsol(inv1,inv,tag);
	multiply2d_nsol(factor,coeffmatrix,prod);
	rowswap_nsol(prod,coeffmatrix,tag);
	
}
/***************************************************************************************/
//Back Substitution.
void substituteb_nsol(double *fac,double *y,double *vec)
{
	int i,j;
	double sum;
	vec[nsol-1]=y[nsol-1]*pow(fac[(nsol-1)*nsol+(nsol-1)],-1);
	for(i=(nsol-2);i>=0;i--)
	{
		sum=0;
		for(j=i+1;j<(nsol);j++)
			sum=sum-fac[i*nsol+j]*vec[j];
		vec[i]=(y[i]+sum)*pow(fac[i*nsol+i],-1);
	}
}
/*********************************************************************************************/
//Forward Substitution.
void substitutef_nsol(double *fac,double *y1,int index,double *vec)
{
	int i,j;
	double d[nsol],sum;
	
	for(i=0;i<nsol;i++)
		d[i]=y1[i*nsol+index];
  	vec[0]=d[0];
	for(i=1;i<nsol;i++)
	{
		sum=0;
		for(j=0;j<i;j++)
			sum=sum-fac[i*nsol+j]*vec[j];
		vec[i]=d[i]+sum;
	}
// 	for(i=0;i<(nsol);i++)
// 		newmat[i*nsol+index]=vec[i];
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
void pivot_nsol(double *coeffmatrix,double *factor,int k,int *tag)
{
	double swap,big;
	int tagswap,i,p,q,tag1;
	big=mod(coeffmatrix[k*nsol+k]);
	tag1=k;
	for(i=k+1;i<(nsol);i++)
	{
		if(mod(coeffmatrix[i*nsol+k])>big)
		{
			tag1=i;
			big=coeffmatrix[i*nsol+k];
		}
	}
	tagswap=tag[k];
	tag[k]=tag[tag1];
	tag[tag1]=tagswap;
	
	for(i=0;i<(nsol);i++)
	{
		swap=coeffmatrix[k*nsol+i];
		coeffmatrix[k*nsol+i]=coeffmatrix[tag1*nsol+i];
		coeffmatrix[tag1*nsol+i]=swap;
	}
	for(i=0;i<k;i++)
	{
		swap=factor[k*nsol+i];
		factor[k*nsol+i]=factor[tag1*nsol+i];
		factor[tag1*nsol+i]=swap;
	}
}
/*******************************************************************************************/			//Swapping Coloumns To get the final identity matrix because of the initial swappping for pivoting.
void colswap_nsol(double *m1,double *m2,int *tag)
{	
	int i,j,k,p;
	for(k=0;k < nsol;k++)
	{
		for(j=0;j<nsol;j++)
		{
			for(p=0;p <nsol; p++)
				m2[p*nsol+tag[j]]=m1[p*nsol+j];
		}
	}
}
/********************************************************************************************/	 
//Switching rows
void rowswap_nsol(double *m1,double *m2,int *tag)
{	
	int i,j,k,p;
	for(k=0;k<(nsol);k++)
	{
		for(j=0;j< nsol ;j++)
		{
			for(p=0; p < nsol;p++)
				m2[tag[j]*nsol+p]=m1[j*nsol+p];
		}
	}
}
/********************************************************************************************/
