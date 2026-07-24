#ifndef FUNCTION_F_H_
#define FUNCTION_F_H_

using namespace amrex;
using namespace std;


//####################################################################################################################

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void c_mu(int i, int j, int k, 
        amrex::Array4<Real const> const& mu, 
        Array1D<Real,0,compcount-2> &c, 
        Array3D <Real,0,phasecount-1,0,compcount-2,0,compcount-2, Order::C> dcdmu_a, 
        Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> BB, 
        Array3D<Real,0,phasecount-1,0,compcount-2,0,compcount-2, 
        Order::C> AA, int numcomp, int a)
{
                
    for(int l=0; l < numcomp-1; l++) 
    {
        double sum = 0.0;
        for (int m=0; m < numcomp-1; m++) 
        {
            sum += dcdmu_a(a,l,m)*(mu(i,j,k,m)-BB(a,m)); 
        }
        c(l) = sum;
                    
    }
}

//####################################################################################################################

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void free_energy(int i, int j, int k, int numphase, 
        Array3D<Real,0,phasecount-1,0,compcount-2,0,compcount-2, Order::C> AA, 
        Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> BB, 
        Array1D <Real,0,phasecount-1> CC, 
        Real &fe ,
        Array1D<Real,0,compcount-2> &c, 
        int numcomp, int a)
{
    double sum=0.0;
        for (int l=0;l<numcomp-1;l++) 
        {
            for (int m=0;m<numcomp-1;m++) 
            {
                if (l<=m) 
                {
                    sum += AA(a,l,m)*c(l)*c(m);
                }
            }
            sum += BB(a,l)*c(l);
        }
    sum += CC(a);
    fe = sum;
}
//####################################################################################################################
#endif
