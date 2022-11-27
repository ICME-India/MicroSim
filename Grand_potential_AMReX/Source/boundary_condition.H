#ifndef BOUNDARY_CONDITIONS_H_
#define BOUNDARY_CONDITIONS_H_

#include <AMReX_BCRec.H>
//#include <AMReX_BCUtil.H>
#include "Variables.H"

using namespace amrex;

void bound_cond(MultiFab& phi_old, Vector<BCRec>& bc, Vector<int> bc_hi, Vector<int> bc_lo){

    //bc_lo is the boundary condition at the lower edge(eg. X-, Y-, Z-)
   // Vector<int> bc_lo(AMREX_SPACEDIM,0);
	//bc_hi is the boundary condition at the upper edge(eg. X+, Y+, Z+)
   // Vector<int> bc_hi(AMREX_SPACEDIM,0);

	//Reading the boundary condition from the infile into bc_lo and bc_hi
	//bc_hi[0] = stod(bound[0][1]);
	//bc_hi[1] = stod(bound[0][3]);
	//bc_lo[0] = stod(bound[0][2]);
	//bc_lo[1] = stod(bound[0][4]);

    
    for (int n = 0; n < phi_old.nComp(); ++n)
    {
        for(int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            //3: Periodic
            if (bc_lo[idim] == 3) {
                bc[n].setLo(idim, BCType::int_dir);
            }
            //1: Neumann
            else if (bc_lo[idim] == 1) {
                bc[n].setLo(idim, BCType::foextrap);
            }
            //2:Dirichlet
            else if(bc_lo[idim] == 2) {
                bc[n].setLo(idim, BCType::ext_dir);
            }
            else {
                amrex::Abort("Invalid bc_lo");
            }

            //3: Periodic
            if (bc_hi[idim] == 3) {
                bc[n].setHi(idim, BCType::int_dir);
            }
            //First Order Extrapolation for Neumann boundary conditions or bc_lo, bc_hi = 2
            else if (bc_hi[idim] == 1) {
                bc[n].setHi(idim, BCType::foextrap);
            }
            //External Dirichlet Boundary Condition, or bc_lo, bc_hi = 3
            else if(bc_hi[idim] == 2) {
                bc[n].setHi(idim, BCType::ext_dir);
            }
            else {
                amrex::Abort("Invalid bc_hi");
            }
        }
    }

}
#endif
