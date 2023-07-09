
double hfhi(double *phi, int a);
double dhfhi(double *phi, int a, int b);

double dhfhi(double *phi, int a, int b) {
  int e, f, i, j, k, l;
  double sum=0.0;
  
  double ans = 0.0;
    double temp = 0.0;

    double phiValue = phi[a];
    double const1 = 7.5*((double)npha-2.0)/((double)npha-1.0);

    if (a == b)
    {
        ans  = 5.0*pow(phiValue, 4)*(6.0 - const1);
        ans += 4.0*pow(phiValue, 3)*(-15.0 + 3.0*const1);

        for (long i = 1; i < npha; i++)
        {
            if (i != a)
            {
                for (long j = 0; j < i; j++)
                {
                    if (j != a)
                    {
                        temp += (phi[j] - phi[i])*(phi[j] - phi[i]);
                    }
                }
            }
        }

        temp *= 7.5*(3.0*phiValue - 2.0)/((double)npha-1.0);
        temp += 3.0*phiValue*(10.0 - 3.0*const1) + 2.0*const1;
        ans  += phiValue*temp;

        temp  = 0.0;
        if (npha > 3)
        {
            for (long i = 2; i < npha; i++)
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
                                    temp += phi[i]*phi[j]*phi[k];
                                }
                            }
                        }
                    }
                }
            }
            ans += 30.0*phiValue*temp;
            temp = 0.0;
        }

        if (npha > 4)
        {
            for (long i = 3; i < npha; i++)
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
                                            temp += phi[i]*phi[j]*phi[k]*phi[l];
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
        if (npha == 2)
        {
            ans = 0.0;//-30.0*phiValue*phiValue*(1.0-phiValue)*(1.0-phiValue);
        }

        if (npha > 2)
        {
            for (long i = 1; i < npha; i++)
            {
                if (i != a)
                {
                    for (long j = 0; j < i; j++)
                    {
                        if (j != a)
                        {
                            if (i == b)
                                temp += -2.0*(phi[j] - phi[i]);
                            else if (j == b)
                                temp += 2.0*(phi[j] - phi[i]);
                        }
                    }
                }
            }
            ans = (phiValue-1.0)*phiValue*phiValue*7.5*temp/((double)npha-1.0);
            temp = 0.0;
        }
        if (npha > 3)
        {
            for (long i = 2; i < npha; i++)
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
                                        temp += phi[j]*phi[k];
                                    else if (j == b)
                                        temp += phi[i]*phi[k];
                                    else if (k == b)
                                        temp += phi[i]*phi[j];
                                }
                            }
                        }
                    }
                }
            }
            ans += 15*phiValue*phiValue*temp;
            temp = 0.0;
        }

        if (npha > 4)
        {
            for (long i = 3; i < npha; i++)
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
                                                temp += phi[j]*phi[k]*phi[l];
                                            else if (j == b)
                                                temp += phi[i]*phi[k]*phi[l];
                                            else if (k == b)
                                                temp += phi[i]*phi[j]*phi[l];
                                            else if (l == b)
                                                temp += phi[i]*phi[j]*phi[k];
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
double hfhi(double *phi, int a) {
  int b,c, i, j, k, l;
  double sum1 = 0.0;
  double sum; 
  
  double ans = 0.0, temp = 0.0;
    double phiValue = phi[a];
    double const1 = 7.5*((double)npha-2.0)/((double)npha-1.0);

    ans  = pow(phiValue, 5)*(6.0 - const1);
    ans += pow(phiValue, 4)*(-15.0 + 3.0*const1);

    for (long i = 1; i < npha; i++)
    {
        if (i != a)
        {
            for (long j = 0; j < i; j++)
            {
                if (j != a)
                {
                    temp += (phi[j] - phi[i])*(phi[j] - phi[i]);
                }
            }
        }
    }
    temp *= 7.5*(phiValue - 1.0)/((double)npha-1.0);
    temp += phiValue*(10.0 - 3.0*const1) + const1;
    ans  += pow(phiValue, 2)*temp;

    temp  = 0.0;
    if (npha > 3)
    {
        for (long i = 2; i < npha; i++)
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
                                temp += phi[i]*phi[j]*phi[k];
                            }
                        }
                    }
                }
            }
        }
        ans += 15.0*phiValue*phiValue*temp;
        temp = 0.0;
    }

    if (npha > 4)
    {
        for (long i = 3; i < npha; i++)
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
                                        temp += phi[i]*phi[j]*phi[k]*phi[l];
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


