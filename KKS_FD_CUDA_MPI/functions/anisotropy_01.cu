#include "anisotropy_01.cuh"

extern __device__
void anisotropy_01_dAdq(double *qab, double *dadq, long a, long b, double *dab, long NUMPHASES)
{
    long X = 0, Y = 1, Z = 2;

    double qx3 = qab[X]*qab[X]*qab[X];
    double qy3 = qab[Y]*qab[Y]*qab[Y];
    double qz3 = qab[Z]*qab[Z]*qab[Z];

    double q2  = qab[X]*qab[X] + qab[Y]*qab[Y] + qab[Z]*qab[Z];
    double q22 = q2*q2;
    double q23 = q2*q2*q2;
    double q4  = qx3*qab[X] + qy3*qab[Y] + qz3*qab[Z];

    if (fabs(q2) > 1.0e-15)
    {
        dadq[X] = 16.0*dab[a*NUMPHASES + b]*(qx3/q22 - qab[X]*q4/q23);
        dadq[Y] = 16.0*dab[a*NUMPHASES + b]*(qy3/q22 - qab[Y]*q4/q23);
        dadq[Z] = 16.0*dab[a*NUMPHASES + b]*(qz3/q22 - qab[Z]*q4/q23);
    }
    else
    {
        dadq[X] = 0.0;
        dadq[Y] = 0.0;
        dadq[Z] = 0.0;
    }
}

extern __device__
double anisotropy_01_function_ac(double *qab, long a, long b, double *dab, long NUMPHASES)
{
    long X = 0, Y = 1, Z = 2;

    double qx2 = qab[X]*qab[X];
    double qx4 = qx2*qx2;

    double qy2 = qab[Y]*qab[Y];
    double qy4 = qy2*qy2;

    double qz2 = qab[Z]*qab[Z];
    double qz4 = qz2*qz2;

    double q2  = qx2 + qy2 + qz2;
    double ac;

    if (fabs(q2) > 1.0e-15)
    {
        ac = 1.0 - dab[a*NUMPHASES + b]*(3.0 - 4.0*(qx4 + qy4 + qz4)/(q2*q2));
    }
    else
    {
        ac = 1.0;
    }

    return ac;
}
