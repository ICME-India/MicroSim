#include "functionH.cuh"

extern __device__
double calcInterp3rd(double **phi, long a, long idx, long NUMPHASES)
{
    if (NUMPHASES < 2)
        return 0.0;

    double phiValue = phi[a][idx];
    double ans = 3.0*phiValue*phiValue - 2.0*phiValue*phiValue*phiValue;

    double sum1 = 0.0;

    for (long b = 0; b < NUMPHASES; b++)
    {
        for (long c = 0; c < NUMPHASES; c++)
        {
            if (b != a && c != a && b < c)
            {
                sum1 += phi[b][idx]*phi[c][idx];
            }
        }
    }

    sum1 *= 2.0*phiValue;

    return ans + sum1;
}

extern __device__
double calcInterp3rdDiff(double **phi, long a, long b, long idx, long NUMPHASES)
{
    long e, f;
    double sum = 0.0;

  if (a == b)
  {
    sum = 6.0*phi[b][idx]*(1.0-phi[b][idx]);
    for (e=0; e < NUMPHASES; e++) {
      for (f=0; f < NUMPHASES; f++) {
        if (e!=a && f!=a && e<f) {
          sum += 2.0*phi[e][idx]*phi[f][idx];
        }
      }
    }
  } else {
    for (e=0; e < NUMPHASES; e++) {
      if (e!=b && e!=a) {
        sum += 2.0*phi[e][idx];
      }
    }
    sum *= phi[a][idx];
  }

  return sum;
}

extern __device__
double calcInterp5th(double **phi, long a, long idx, long NUMPHASES)
{
    if (NUMPHASES < 2)
        return 0.0;

    double ans = 0.0, temp = 0.0;
    double phiValue = phi[a][idx];
    double const1 = 7.5*((double)NUMPHASES-2.0)/((double)NUMPHASES-1.0);

    ans  = pow(phiValue, 5)*(6.0 - const1);
    ans += pow(phiValue, 4)*(-15.0 + 3.0*const1);

    for (long i = 1; i < NUMPHASES; i++)
    {
        if (i != a)
        {
            for (long j = 0; j < i; j++)
            {
                if (j != a)
                {
                    temp += (phi[j][idx] - phi[i][idx])*(phi[j][idx] - phi[i][idx]);
                }
            }
        }
    }
    temp *= 7.5*(phiValue - 1.0)/((double)NUMPHASES-1.0);
    temp += phiValue*(10.0 - 3.0*const1) + const1;
    ans  += pow(phiValue, 2)*temp;

    temp  = 0.0;
    if (NUMPHASES > 3)
    {
        for (long i = 2; i < NUMPHASES; i++)
        {
            if (i != a)
            {
                for (long j = 1; j < i; j++)
                {
                    if (j != a)
                    {
                        for (long k = 0; k < j; k++)
                        {
                            if (k != a)
                            {
                                temp += phi[i][idx]*phi[j][idx]*phi[k][idx];
                            }
                        }
                    }
                }
            }
        }
        ans += 15.0*phiValue*phiValue*temp;
        temp = 0.0;
    }

    if (NUMPHASES > 4)
    {
        for (long i = 3; i < NUMPHASES; i++)
        {
            if (i != a)
            {
                for (long j = 2; j < i; j++)
                {
                    if (j != a)
                    {
                        for (long k = 1; k < j; k++)
                        {
                            if (k != a)
                            {
                                for (long l = 0; l < k; l++)
                                {
                                    if (l != a)
                                    {
                                        temp += phi[i][idx]*phi[j][idx]*phi[k][idx]*phi[l][idx];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        ans += 24.0*phiValue*temp;
    }

    return ans;
}

extern __device__
double calcInterp5thDiff(double **phi, long a, long b, long idx, long NUMPHASES)
{
    if (NUMPHASES < 2)
        return 0.0;

    double ans = 0.0;
    double temp = 0.0;

    double phiValue = phi[a][idx];
    double const1 = 7.5*((double)NUMPHASES-2.0)/((double)NUMPHASES-1.0);

    if (a == b)
    {
        ans  = 5.0*pow(phiValue, 4)*(6.0 - const1);
        ans += 4.0*pow(phiValue, 3)*(-15.0 + 3.0*const1);

        for (long i = 1; i < NUMPHASES; i++)
        {
            if (i != a)
            {
                for (long j = 0; j < i; j++)
                {
                    if (j != a)
                    {
                        temp += (phi[j][idx] - phi[i][idx])*(phi[j][idx] - phi[i][idx]);
                    }
                }
            }
        }

        temp *= 7.5*(3.0*phiValue - 2.0)/((double)NUMPHASES-1.0);
        temp += 3.0*phiValue*(10.0 - 3.0*const1) + 2.0*const1;
        ans  += phiValue*temp;

        temp  = 0.0;
        if (NUMPHASES > 3)
        {
            for (long i = 2; i < NUMPHASES; i++)
            {
                if (i != a)
                {
                    for (long j = 1; j < i; j++)
                    {
                        if (j != a)
                        {
                            for (long k = 0; k < j; k++)
                            {
                                if (k != a)
                                {
                                    temp += phi[i][idx]*phi[j][idx]*phi[k][idx];
                                }
                            }
                        }
                    }
                }
            }
            ans += 30.0*phiValue*temp;
            temp = 0.0;
        }

        if (NUMPHASES > 4)
        {
            for (long i = 3; i < NUMPHASES; i++)
            {
                if (i != a)
                {
                    for (long j = 2; j < i; j++)
                    {
                        if (j != a)
                        {
                            for (long k = 1; k < j; k++)
                            {
                                if (k != a)
                                {
                                    for (long l = 0; l < k; l++)
                                    {
                                        if (l != a)
                                        {
                                            temp += phi[i][idx]*phi[j][idx]*phi[k][idx]*phi[l][idx];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            ans += 24.0*temp;
        }
    }
    else
    {
        if (NUMPHASES == 2)
        {
            ans = -30.0*phiValue*phiValue*(1.0-phiValue)*(1.0-phiValue);
        }

        if (NUMPHASES > 2)
        {
            for (long i = 1; i < NUMPHASES; i++)
            {
                if (i != a)
                {
                    for (long j = 0; j < i; j++)
                    {
                        if (j != a)
                        {
                            if (i == b)
                                temp += -2.0*(phi[j][idx] - phi[i][idx]);
                            else if (j == b)
                                temp += 2.0*(phi[j][idx] - phi[i][idx]);
                        }
                    }
                }
            }
            ans = (phiValue-1.0)*phiValue*phiValue*7.5*temp/((double)NUMPHASES-1.0);
            temp = 0.0;
        }
        if (NUMPHASES > 3)
        {
            for (long i = 2; i < NUMPHASES; i++)
            {
                if (i != a)
                {
                    for (long j = 1; j < i; j++)
                    {
                        if (j != a)
                        {
                            for (long k = 0; k < j; k++)
                            {
                                if (k != a)
                                {
                                    if (i == b)
                                        temp += phi[j][idx]*phi[k][idx];
                                    else if (j == b)
                                        temp += phi[i][idx]*phi[k][idx];
                                    else if (k == b)
                                        temp += phi[i][idx]*phi[j][idx];
                                }
                            }
                        }
                    }
                }
            }
            ans += 15*phiValue*phiValue*temp;
            temp = 0.0;
        }

        if (NUMPHASES > 4)
        {
            for (long i = 3; i < NUMPHASES; i++)
            {
                if (i != a)
                {
                    for (long j = 2; j < i; j++)
                    {
                        if (j != a)
                        {
                            for (long k = 1; k < j; k++)
                            {
                                if (k != a)
                                {
                                    for (long l = 0; l < k; l++)
                                    {
                                        if (l != a)
                                        {
                                            if (i == b)
                                                temp += phi[j][idx]*phi[k][idx]*phi[l][idx];
                                            else if (j == b)
                                                temp += phi[i][idx]*phi[k][idx]*phi[l][idx];
                                            else if (k == b)
                                                temp += phi[i][idx]*phi[j][idx]*phi[l][idx];
                                            else if (l == b)
                                                temp += phi[i][idx]*phi[j][idx]*phi[k][idx];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            ans += 24.0*phiValue*temp;
        }
    }

    return ans;
}
