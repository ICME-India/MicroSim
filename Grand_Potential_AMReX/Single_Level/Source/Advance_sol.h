#ifndef _ADVANCE_SOL_H_
#define _ADVANCE_SOL_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <iostream>
#include <fstream>
#include <string>
#include <AMReX_BCRec.H>
#include <AMReX_BCUtil.H>

#include "Variables.h"
#include "Filling_fv.h"
#include "Adv_phi.h"
#include "Function_W.h"
#include "Function_dpsi.h"
#include "Function_A.h"
#include "Adv_Chem_Pot.h"
#include "Functionf_elast.h"

using namespace amrex;
			

void advance(	MultiFab& phi_old, 
		MultiFab& phi_new,
		MultiFab& mu_old, 
		MultiFab& mu_new,
		MultiFab& comp_old,
		MultiFab& comp_new,
		MultiFab& term1, 
		MultiFab& term2,
		MultiFab& term3,
		MultiFab& term4,
		MultiFab& dispx,
		MultiFab& dispy,
		MultiFab& dispz,
		MultiFab& psi,
		MultiFab& lambad,
		Vector<BCRec> bc_phi,
		Vector<BCRec> bc_mu,
		Vector<BCRec> bc_comp,
		Vector<BCRec> bc_disp,
		Geometry const& geom)
{
	//Fill the ghost cells
	phi_old.FillBoundary(geom.periodicity());
	if(GG!=1){
		mu_old.FillBoundary(geom.periodicity());
		comp_old.FillBoundary(geom.periodicity());
	}
	//mu_old.FillBoundary(geom.periodicity());
	//comp_old.FillBoundary(geom.periodicity());
	if(ELASTICITY==1){
		dispx.FillBoundary(geom.periodicity());
		dispy.FillBoundary(geom.periodicity());
		#if(AMREX_SPACEDIM>2)
			dispz.FillBoundary(geom.periodicity());
		#endif
	}


	//Filling the boundary cells for phi_old and mu_old----------------------------------------------------------------------------------
	FillDomainBoundary(phi_old, geom, bc_phi);
	if(GG!=1){
	 FillDomainBoundary(mu_old, geom, bc_mu);
	 FillDomainBoundary(comp_old, geom, bc_comp);
	}
		
	if(ELASTICITY == 1){
		FillDomainBoundary(dispx, geom, bc_disp);
		FillDomainBoundary(dispy, geom, bc_disp);
		#if(AMREX_SPACEDIM>2)
			FillDomainBoundary(dispz, geom, bc_disp);
		#endif
	}

	//Print()<<"3\n";
	//Computing the anisotropy term(term1) in the phi evolution equation (Refer to Function_A.h for the formulation)
	aniso_term(term1, phi_old, geom);

	//Print()<<"4\n";
	//Computing the Double well Potential(term2) in the phi evolution equation (Refer to Function_W.h for the formulation of double well potential calculation)
	dwdphi(term2, phi_old, geom);

	//Print()<<"5\n";
	//Computing the Psi equation(term3) in the phi evolution equation (Refer to Function_dpsi.h for the formulation of psi calculation)
	//dpsi(mu_old, term3, phi_old, psi, geom);
	if(funcf!=5){
		dpsi(mu_old, term3, phi_old, psi, geom);
	}
	if(funcf==5){
		func_dpsi(phi_old, term3);
	}

	if(ELASTICITY==1){
	#if(AMREX_SPACEDIM==2)	
		df_elast_2D(phi_old, dispx, dispy, term4);
	#elif(AMREX_SPACEDIM>2)
		df_elast_3D(phi_old, dispx, dispy, dispz, term4);
	#endif
	
	}

	//Print()<<"6\n";
	//Now we have all the terms for terms for phi evolution, we simply add them (Refer to Adv_phi.h for the formulation of update_phi function)
	update_phi(phi_new, phi_old, term1,term2,term3,term4,lambad,geom);

	//Print()<<"7\n";
	//Fill ghost cells of phi_new along with periodic boundaries
	// if(PERIODIC_BC){
		
	//}
	//else{
		//FillDomainBoundary(phi_new, geom, bc_phi);
	//}
	//Print()<<"8\n";
	//Phi is already updated now here we update mu (Refer to dmudt fucntion in adv_Chem_Pot.h for formulation)
	if(GG!=1){
		phi_new.FillBoundary(geom.periodicity());
		Chem_pot(mu_new, mu_old, phi_new, phi_old, comp_new, comp_old,geom);
	}

}

#endif
