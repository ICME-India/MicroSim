/******************************************************************************
 *                       Code generated with sympy 1.8                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                       This file is part of 'project'                       *
 ******************************************************************************/


#ifndef THERMO_CUH_
#define THERMO_CUH_

extern __device__ __host__ void GE_0(double T, double *y, double *Ge);
extern __device__ __host__ void GE_1(double T, double *y, double *Ge);
extern __device__ __host__ void Mu_0(double T, double *y, double *Mu);
extern __device__ __host__ void Mu_1(double T, double *y, double *Mu);
extern __device__ __host__ void dmudc_0(double T, double *y, double *Dmudc);
extern __device__ __host__ void dmudc_1(double T, double *y, double *Dmudc);

static void(*free_energy_tdb[])(double T, double *y, double *Ge) = {GE_0,GE_1};
static __device__ void(*free_energy_tdb_dev[])(double T, double *y, double *Ge) = {GE_0,GE_1};

static void(*Mu_tdb[])(double T, double *y, double *Mu) = {Mu_0,Mu_1};
static __device__ void(*Mu_tdb_dev[])(double T, double *y, double *Mu) = {Mu_0,Mu_1};

static void(*dmudc_tdb[])(double T, double *y, double *Dmudc) = {dmudc_0,dmudc_1};
static __device__ void(*dmudc_tdb_dev[])(double T, double *y, double *Dmudc) = {dmudc_0,dmudc_1};

#endif   

