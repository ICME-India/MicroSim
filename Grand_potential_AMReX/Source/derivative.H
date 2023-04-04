#ifndef _DERIVATIVE_H_
#define _DERIVATIVE_H_

using namespace amrex;

//0 is iph
//1 is imh
//2 is jph
//3 is jmh

void compute_der(MultiFab& phi_old, MultiFab& derivx, MultiFab& derivy, Geometry const& geom)
{   
	BL_PROFILE("compute_der()");  

  GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& gbx = mfi.validbox();
        Array4<Real> const& phiOld = phi_old.array(mfi);
        Array4<Real> const& derivxx = derivx.array(mfi);
        Array4<Real> const& derivyy = derivy.array(mfi);
        
        amrex::ParallelFor(gbx, 
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            derivxx(i,j,k,0)=(phiOld(i+1,j,k)-phiOld(i,j,k))/dx[0];
            derivxx(i,j,k,1)=(phiOld(i,j,k)-phiOld(i-1,j,k))/dx[0];
            derivxx(i,j,k,2)=(phiOld(i+1,j+1,k)-phiOld(i-1,j+1,k)+phiOld(i+1,j,k)-phiOld(i-1,j,k))/(4*dx[0]);
            derivxx(i,j,k,3)=(phiOld(i+1,j,k)-phiOld(i-1,j,k)+phiOld(i+1,j-1,k)-phiOld(i-1,j-1,k))/(4*dx[0]);

            derivyy(i,j,k,0)=(phiOld(i+1,j+1,k)-phiOld(i+1,j-1,k)+phiOld(i,j+1,k)-phiOld(i,j-1,k))/(4*dx[1]);
            derivyy(i,j,k,1)=(phiOld(i,j+1,k)-phiOld(i,j-1,k)+phiOld(i-1,j+1,k)-phiOld(i-1,j-1,k))/(4*dx[1]);
            derivyy(i,j,k,2)=(phiOld(i,j+1,k)-phiOld(i,j,k))/dx[1];
            derivyy(i,j,k,3)=(phiOld(i,j,k)-phiOld(i,j-1,k))/dx[1];

        });
}
}


#endif
