#ifndef BOUNDARY_CONDITIONS_H_
#define BOUNDARY_CONDITIONS_H_

#include <AMReX_BCRec.H>
#include <AMReX_BCUtil.H>

using namespace amrex;

//####################################################################################################################
void bound_cond(Vector<BCRec>& bc, Vector<int>& bc_hi, Vector<int>& bc_lo, int comp)
{

	BL_PROFILE("bound_cond()");
    for (int n = 0; n < comp; ++n)
    {   
        for(int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            //1: Neumann
            if (bc_lo[idim] == 1) 
            {
                bc[n].setLo(idim, BCType::foextrap);
            }

            //2:Dirichlet
            else if(bc_lo[idim] == 2) 
            {
                bc[n].setLo(idim, BCType::ext_dir);
            }

            //3: Periodic
            else if (bc_lo[idim] == 3) 
            {
                bc[n].setLo(idim, BCType::int_dir);
            }

            else 
            {
                amrex::Abort("Invalid bc_lo");
            }

            //1: Neumann 
            if (bc_hi[idim] == 1) 
            {
                bc[n].setHi(idim, BCType::foextrap);
            }

            //2: Dirichlet 
            else if(bc_hi[idim] == 2) 
            {
                bc[n].setHi(idim, BCType::ext_dir);
            }

             //3: Periodic
            else if (bc_hi[idim] == 3) 
            {
                bc[n].setHi(idim, BCType::int_dir);
            }
            else 
            {
                amrex::Abort("Invalid bc_hi");
            }
        }
    }

}

//####################################################################################################################
#endif
