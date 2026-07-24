#ifndef GP_UTILITY_H
#define GP_UTILITY_H

#include <AMReX_Utility.H>
using namespace amrex;


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Vec_rot(int i, int j, int k, int a, int b, Array2D<Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> &vec, Array2D<Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> &r_vec, Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_x,Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_y,Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_z){

    
#if(AMREX_SPACEDIM == 2)
    
        for(int p = 0; p < vec.xlen(); p++){
            r_vec(p,X) = vec(p,X)*cos(matrot_z(a,b))-vec(p,Y)*sin(matrot_z(a,b));
            r_vec(p,Y) = vec(p,X)*sin(matrot_z(a,b))+vec(p,Y)*cos(matrot_z(a,b));
        }

#endif

#if(AMREX_SPACEDIM>2)

    for(int p = 0; p < vec.xlen(); p++){
            r_vec(p,X) = vec(p,X)*(cos(matrot_y(a,b))*cos(matrot_z(a,b)))-vec(p,Y)*(cos(matrot_y(a,b))*sin(matrot_z(a,b)))-vec(p,Z)*sin(matrot_y(a,b));
            r_vec(p,Y) = vec(p,X)*(-cos(matrot_z(a,b))*sin(matrot_x(a,b))*sin(matrot_y(a,b))+cos(matrot_x(a,b))*sin(matrot_z(a,b)))+vec(p,Y)*(sin(matrot_x(a,b))*sin(matrot_y(a,b))*sin(matrot_z(a,b))+cos(matrot_x(a,b))*cos(matrot_z(a,b)))-vec(p,Z)*(sin(matrot_x(a,b))*cos(matrot_y(a,b)));
            r_vec(p,Z) = vec(p,X)*(cos(matrot_x(a,b))*sin(matrot_y(a,b))*cos(matrot_z(a,b))+sin(matrot_z(a,b))*sin(matrot_x(a,b)))+vec(p,Y)*(-cos(matrot_x(a,b))*sin(matrot_y(a,b))*sin(matrot_z(a,b))+cos(matrot_z(a,b))*sin(matrot_x(a,b)))+vec(p,Z)*(cos(matrot_x(a,b))*cos(matrot_y(a,b)));
        }

#endif

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Vec_rot(int i, int j, int k, int a, int b, Array1D <Real, 0, AMREX_SPACEDIM-1> &vec, Array1D <Real, 0, AMREX_SPACEDIM-1> &r_vec,Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_x,Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_y,Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_z){

#if(AMREX_SPACEDIM == 2)
    
            r_vec(X) = vec(X)*cos(matrot_z(a,b))-vec(Y)*sin(matrot_z(a,b));
            r_vec(Y) = vec(X)*sin(matrot_z(a,b))+vec(Y)*cos(matrot_z(a,b));


#endif

#if(AMREX_SPACEDIM>2)

            r_vec(X) = vec(X)*(cos(matrot_y(a,b))*cos(matrot_z(a,b)))-vec(Y)*(cos(matrot_y(a,b))*sin(matrot_z(a,b)))-vec(Z)*sin(matrot_y(a,b));
            r_vec(Y) = vec(X)*(-cos(matrot_z(a,b))*sin(matrot_x(a,b))*sin(matrot_y(a,b))+cos(matrot_x(a,b))*sin(matrot_z(a,b)))+vec(Y)*(sin(matrot_x(a,b))*sin(matrot_y(a,b))*sin(matrot_z(a,b))+cos(matrot_x(a,b))*cos(matrot_z(a,b)))-vec(Z)*(sin(matrot_x(a,b))*cos(matrot_y(a,b)));
            r_vec(Z) = vec(X)*(cos(matrot_x(a,b))*sin(matrot_y(a,b))*cos(matrot_z(a,b))+sin(matrot_z(a,b))*sin(matrot_x(a,b)))+vec(Y)*(-cos(matrot_x(a,b))*sin(matrot_y(a,b))*sin(matrot_z(a,b))+cos(matrot_z(a,b))*sin(matrot_x(a,b)))+vec(Z)*(cos(matrot_x(a,b))*cos(matrot_y(a,b)));
        

#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Inv_vec_rot(int i, int j, int k, int a, int b, Array2D <Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> &vec, Array2D <Real,0,AMREX_SPACEDIM*2,0,AMREX_SPACEDIM-1,Order::C> &r_vec, Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_x, Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_y, Array2D <Real,0,phasecount-1,0,phasecount-1,Order::C> matrot_z){
    
#if(AMREX_SPACEDIM == 2)
    
        for(int p = 0; p < vec.xlen(); p++){
            r_vec(p,X) = vec(p,X)*cos(matrot_z(a,b))+vec(p,Y)*sin(matrot_z(a,b));
            r_vec(p,Y) = -vec(p,X)*sin(matrot_z(a,b))+vec(p,Y)*cos(matrot_z(a,b));
        }

#endif

#if(AMREX_SPACEDIM>2)

    for(int p = 0; p < vec.xlen(); p++){
            r_vec(p,X) = vec(p,X)*(cos(matrot_y(a,b))*cos(matrot_z(a,b)))+vec(p,Y)*(cos(matrot_y(a,b))*sin(matrot_z(a,b)))+vec(p,Z)*sin(matrot_y(a,b));
            r_vec(p,Y) = vec(p,X)*(-cos(matrot_z(a,b))*sin(matrot_x(a,b))*sin(matrot_y(a,b))-cos(matrot_x(a,b))*sin(matrot_z(a,b)))+vec(p,Y)*(-sin(matrot_x(a,b))*sin(matrot_y(a,b))*sin(matrot_z(a,b))+cos(matrot_x(a,b))*cos(matrot_z(a,b)))+vec(p,Z)*(sin(matrot_x(a,b))*cos(matrot_y(a,b)));
            r_vec(p,Z) = vec(p,X)*(-cos(matrot_x(a,b))*sin(matrot_y(a,b))*cos(matrot_z(a,b))+sin(matrot_z(a,b))*sin(matrot_x(a,b)))+vec(p,Y)*(-cos(matrot_x(a,b))*sin(matrot_y(a,b))*sin(matrot_z(a,b))-cos(matrot_z(a,b))*sin(matrot_x(a,b)))+vec(p,Z)*(cos(matrot_x(a,b))*cos(matrot_y(a,b)));
        }

#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mat_inv(Array2D<Real,0,compcount-2,0,compcount-2, Order::C> &denom, Array2D<Real,0,compcount-2,0,compcount-2, Order::C> &inv_denom,int numcomp){

    int size = numcomp-1;
    Real fact{};
    Array2D<Real,0,compcount-2,0,compcount-2, Order::C> factor{};
    Array2D<Real,0,compcount-2,0,compcount-2, Order::C> iden{};
    Array1D<Real,0,compcount-2> vec1{};
    Array1D<Real,0,compcount-2> vec{};
    //making upper traingular matrix
    for(int m=0; m<size; m++){
        //pivot
        Real big = amrex::Math::abs(denom(m,m));
        int tag1 = m;
        Real swap{0.0};
        for(int n=m+1; n<size;n++){
            if(amrex::Math::abs(denom(n,m))>big){
                big=denom(n,m);
            }
        }

        for(int n=0;n<size;n++)
	    {
            swap = denom(m,n);
            denom(m,n)=denom(tag1,n);
            denom(tag1,n)=swap;
        }
        for(int n=0; n<m; n++)
        {
            swap=factor(m,n);
            factor(m,n)=factor(tag1,n);
            factor(tag1,n)=swap;
        }

        //back to loop
        for(int n=m+1; n<size; n++){
            fact = -denom(n,m)/denom(m,m);
            factor(n,m) = -fact;
            for(int p=m; (p<=size-1); p++){
                denom(n,p) = fact*denom(m,p)+denom(n,p);
            }
        }
    }

    for(int m=0; m<size; m++) {
		for(int n=0; n<size; n++)
		{
			if(m==n)
				factor(m,n)=1;
			if(n>m)
				factor(m,n)=0;
		}
	}
    
    //Identity matrix
    for(int m=0; m<(size);m++)
	{
		for(int n=0;n<(size);n++)
		{
			if(m==n)
				iden(m,n)=1;
			else
				iden(m,n)=0;
		}
	}

    //Forward and backward substitution

    for(int m=0; m<size; m++){
      
   

    //forward substitution 
      double sum;
      Array1D<Real,0,compcount-2> d;
    
      for(int n=0;n<size;n++){
        d(n)=iden(n,m);
      }
        vec1(0)=d(0);
      for(int n=1; n<size; n++)
      {
        sum=0;
        for(int l=0; l<n; l++){
          sum=sum-factor(n,l)*vec1(l);
        }
        vec1(n)=d(n)+sum;
      }



    //backward substitution
      vec(size-1)=vec1(size-1)*pow(denom(size-1,size-1),-1);
      for(int n=(size-2); n>=0; n--)
      {
        sum=0;
        for(int l=n+1; l<(size); l++){
          sum=sum-denom(n,l)*vec(l);
        }

        vec(n)=(vec1(n)+sum)*pow(denom(n,n),-1);
      }


    //copying to inv
      for(int n=0; n<(size); n++){
			  inv_denom(n,m)=vec(n);
      }


    }

    

}

void mat_mul2D(Vector<Vector<Vector<Real>>> &dc_dmu,Vector<Vector<Vector<double>>> &diffu,Vector<Vector<Real>> &prod, int sz) {
	
    Real sum;
    for(int l=0; l<sz; l++)
    {
        for(int m=0; m<sz; m++)
        {
          sum=0;
            for(int n=0; n<sz; n++){
              sum=sum+dc_dmu[nump-1][l][n]*diffu[nump-1][n][m];
            }
          prod[l][m]=sum;
        }
    }
}

void mat_mul1D(Vector<Vector<Real>> &inv_dcdmu, Vector<Real> &deltac, Vector<Real> &deltamu, int sz){
  
   Real sum;
    for(int i=0; i<sz; i++)
    {
        sum = 0;
        for(int j=0; j<sz; j++){
            sum = sum + inv_dcdmu[i][j]*deltac[j];
        }
        deltamu[i] = sum;
    }
}
#endif