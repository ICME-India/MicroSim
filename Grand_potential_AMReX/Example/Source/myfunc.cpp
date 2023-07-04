#include "myfunc.H"
#include "mykernel.H"
#include <AMReX_BCRec.H>
#include <AMReX_BCUtil.H>
#include <AMReX_BC_TYPES.H>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <AMReX_ParallelDescriptor.H>


void advance (MultiFab& phi_old,
              MultiFab& phi_new,
	      MultiFab& lap_temp,
	      MultiFab& dfdcon,
	      Real& K,
	      Real& M,              
              Real& dt,
              Geometry const& geom)
{
	phi_old.FillBoundary(geom.periodicity());
    // Fill the ghost cells of each grid from the other grids
    // includes periodic domain boundaries
   // phi_old.FillBoundary(geom.periodicity());
    // There are physical boundaries to fill.
   //FillDomainBoundary(phi_old, geom, BoundaryCondition);

    //const BCRec& bc = BoundaryCondition[0];

    //
    // Note that this simple example is not optimized.
    // The following two MFIter loops could be merged
    // and we do not have to use flux MultiFab.
    //
    // =======================================================

    // This example supports both 2D and 3D.  Otherwise,
    // we would not need to use AMREX_D_TERM.
//--------------------------------------------------------------   
 	

//	laplacian(phi_old, lap_temp, geom);	

//	lap_temp.mult(-K,2);

//	MultiFab::Add(lap_temp, dfdcon, 0, 0, 1, 2);
	
//	laplacian(phi_old, lap_temp, geom);

//	lap_temp.mult(M,2);
//------------------------------------------------------------


    // Advance the solution one grid at a time
    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();
        
        auto const& phiOld = phi_old.array(mfi);
        auto const& phiNew = phi_new.array(mfi);
	auto const& laptemp = lap_temp.array(mfi);
	auto const& dc = dfdcon.array(mfi);

        amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
	    derivative(i, j, k, phiOld, dc);
	    laplacian(i, j, k, phiOld, laptemp, geom);      
        });
    }

    lap_temp.mult(-K);
    MultiFab::Add(lap_temp, dfdcon, 0, 0, 1, 0);
    lap_temp.FillBoundary(geom.periodicity());

    

    for ( MFIter mfi(lap_temp); mfi.isValid(); ++mfi )
    {
        const Box& validbx = mfi.validbox();
        
        auto const& phiOLD = phi_old.array(mfi);
        auto const& phiNEW = phi_new.array(mfi);
	auto const& lapTEMP = lap_temp.array(mfi);
	auto const& df = dfdcon.array(mfi);

        amrex::ParallelFor(validbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {

	    laplacian(i, j, k, lapTEMP, df, geom);

            update_phi(i,j,k,phiOLD,phiNEW,
                       df,M,
                       dt);
        });
    }

	
	


}

void init_phi(amrex::MultiFab& phi_new, amrex::Geometry const& geom, amrex::Real c0){

    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();
    GpuArray<Real,AMREX_SPACEDIM> prob_lo = geom.ProbLoArray();
    srand(time(0));
    // =======================================
    // Initialize phi_new by calling a Fortran routine.
    // MFIter = MultiFab Iterator
    for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phiNew = phi_new.array(mfi);
        amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            init_phi(i,j,k,phiNew,c0);
        });
    }
}


/*
void compute_energy(MultiFab const& phiold, Real& energ, Real const& K, Geometry const& geom)
{
	

	for (MFIter mfi(phiold); mfi.isValid(); ++mfi)
	{
		const Box& vvbx = mfi.validbox();
		auto const& ma = phiold.array(mfi);

		amrex::ParallelFor(vvbx,
		[=] AMREX_GPU_DEVICE (int i, int j, int k)
		{
			energy1cal(i,j,k,phiold,energ);
			energy2cal(i,j,k,phiold,energ,K,geom);
			
		
		});

	
	}	

	
	
	std::cout<<energ<<"\n";
	
	
	
}

*/






