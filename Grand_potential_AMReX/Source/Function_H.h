#ifndef FUNCTION_H_H_
#define FUNCTION_H_H_

#include <AMReX_Utility.H>
#include "Function_F4.h"

using namespace amrex;

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
double hphi(int i, int j, int k, amrex::Array4<amrex::Real const> const& phi, long a, int numphase)
{
    if (numphase < 2){
        return 0.0;
       }

    double ans = 0.0, temp = 0.0;
    double const1 = 7.5*((double)numphase-2.0)/((double)numphase-1.0);

    ans = pow(phi(i,j,k,a),4)*(phi(i,j,k,a)*(6.0 - const1) - 15.0 + 3.0*const1);

    for (long m = 1; m < numphase; m++)
    {
        if (m!= a)
        {
            for (long n = 0; n < m; n++)
            {
                if (n!= a)
                {
                    temp += (phi(i,j,k,n) - phi(i,j,k,m))*(phi(i,j,k,n) - phi(i,j,k,m));
                }
            }
        }
    }
    temp *= 7.5*(phi(i,j,k,a) - 1.0)/((double)numphase-1.0);
    temp += phi(i,j,k,a)*(10.0 - 3.0*const1) + const1;
    ans  += phi(i,j,k,a)*phi(i,j,k,a)*temp;

    if (numphase == 2)
        return ans;

    temp  = 0.0;

    if (numphase > 3)
    {
        for (long m = 2; m < numphase; m++)
        {
            if (m != a)
            {
                for (long n = 1; n < m; n++)
                {
                    if (n != a)
                    {
                        for (long p = 0; p < n; p++)
                        {
                            if (p != a)
                            {
                                temp += phi(i,j,k,m)*phi(i,j,k,n)*phi(i,j,k,p);
                            }
                        }
                    }
                }
            }
        }
        ans += 15.0*phi(i,j,k,a)*phi(i,j,k,a)*temp;
        temp = 0.0;
    }

    if (numphase > 4)
    {
        for (long m = 3; m < numphase; m++)
        {
            if (m != a)
            {
                for (long n = 2; n < m; n++)
                {
                    if (n != a)
                    {
                        for (long p = 1; p < n; p++)
                        {
                            if (p != a)
                            {
                                for (long q = 0; q < p; q++)
                                {
                                    if (q != a)
                                    {
                                        temp += phi(i,j,k,m)*phi(i,j,k,n)*phi(i,j,k,p)*phi(i,j,k,q);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        ans += 24.0*phi(i,j,k,a)*temp;
    }

    return ans;
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
double dhphi(int i, int j, int k, amrex::Array4<amrex::Real const> const& phi, long a, long b,int numphase){
    
    if (numphase < 2){
        return 0.0;
     }

    double ans = 0.0;
    double temp = 0.0;
    double const1 = 7.5*((double)numphase-2.0)/((double)numphase-1.0);

    if (a == b)
    {
        ans = pow(phi(i,j,k,a),3)*(5.0*phi(i,j,k,a)*(6.0 - const1) + 4.0*(3.0*const1 - 15.0));

        for (long m = 1; m < numphase; m++)
        {
            if (m!= a)
            {
                for (long n = 0; n < m; n++)
                {
                    if (n != a)
                    {
                        temp += (phi(i,j,k,n) - phi(i,j,k,m))*(phi(i,j,k,n) - phi(i,j,k,m));
                    }
                }
            }
        }

        temp *= 7.5*(3.0*phi(i,j,k,a) - 2.0)/((double)numphase-1.0);
        temp += 3.0*phi(i,j,k,a)*(10.0 - 3.0*const1) + 2.0*const1;
        ans  += phi(i,j,k,a)*temp;

        if (numphase == 2){
            return ans;
         }

        temp  = 0.0;
        if (numphase > 3)
        {
            for (long m = 2; m < numphase; m++)
            {
                if (m != a)
                {
                    for (long n = 1; n < m; n++)
                    {
                        if (n != a)
                        {
                            for (long p = 0; p < n; p++)
                            {
                                if (p != a)
                                {
                                    temp += phi(i,j,k,m)*phi(i,j,k,n)*phi(i,j,k,p);
                                }
                            }
                        }
                    }
                }
            }
            ans += 30.0*phi(i,j,k,a)*temp;
            temp = 0.0;
        }

        if (numphase > 4)
        {
            for (long m = 3; m < numphase; m++)
            {
                if (m != a)
                {
                    for (long n = 2; n < m; n++)
                    {
                        if (n != a)
                        {
                            for (long p = 1; p < n; p++)
                            {
                                if (p != a)
                                {
                                    for (long q = 0; q < p; q++)
                                    {
                                        if (q != a)
                                        {
                                            temp += phi(i,j,k,m)*phi(i,j,k,n)*phi(i,j,k,p)*phi(i,j,k,q);
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
        if (numphase == 2)
        {
            ans = 0.0;
            //ans = -30.0*phi(i,j,k,a)*phi(i,j,k,a)*(1.0-phi(i,j,k,a))*(1.0-phi(i,j,k,a));
            return ans;
        }

        if (numphase > 2)
        {
            for (long m = 1; m < numphase; m++)
            {
                if (m != a)
                {
                    for (long n = 0; n < m; n++)
                    {
                        if (n != a)
                        {
                            if (m == b)
                                temp += -2.0*(phi(i,j,k,n) - phi(i,j,k,m));
                            else if (n == b)
                                temp += 2.0*(phi(i,j,k,n) - phi(i,j,k,m));
                        }
                    }
                }
            }
            ans = (phi(i,j,k,a)-1.0)*phi(i,j,k,a)*phi(i,j,k,a)*7.5*temp/((double)numphase-1.0);
            temp = 0.0;
        }
        if (numphase > 3)
        {
            for (long m = 2; m < numphase; m++)
            {
                if (m!= a)
                {
                    for (long n = 1; n < m; n++)
                    {
                        if (n!= a)
                        {
                            for (long p = 0; p < n; p++)
                            {
                                if (p!= a)
                                {
                                    if (m == b)
                                        temp += phi(i,j,k,n)*phi(i,j,k,p);
                                    else if (n == b)
                                        temp += phi(i,j,k,m)*phi(i,j,k,p);
                                    else if (p == b)
                                        temp += phi(i,j,k,m)*phi(i,j,k,n);
                                }
                            }
                        }
                    }
                }
            }
            ans += 15*phi(i,j,k,a)*phi(i,j,k,a)*temp;
            temp = 0.0;
        }

        if (numphase > 4)
        {
            for (long m = 3; m < numphase; m++)
            {
                if (m!= a)
                {
                    for (long n = 2; n < m; n++)
                    {
                        if (n!= a)
                        {
                            for (long p= 1; p < n; p++)
                            {
                                if (p!= a)
                                {
                                    for (long q = 0; q < p; q++)
                                    {
                                        if (q!= a)
                                        {
                                            if (m == b)
                                                temp += phi(i,j,k,n)*phi(i,j,k,p)*phi(i,j,k,q);
                                            else if (n == b)
                                                temp += phi(i,j,k,m)*phi(i,j,k,p)*phi(i,j,k,q);
                                            else if (p == b)
                                                temp += phi(i,j,k,m)*phi(i,j,k,n)*phi(i,j,k,q);
                                            else if (q == b)
                                                temp += phi(i,j,k,m)*phi(i,j,k,n)*phi(i,j,k,p);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            ans += 24.0*phi(i,j,k,a)*temp;
        }
    }

    return ans;
}


#endif
