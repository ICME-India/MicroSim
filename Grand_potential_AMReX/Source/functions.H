#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <iostream>
#include <fstream>
#include <string>
#include <AMReX_BCRec.H>
//#include <AMReX_BCUtil.H>

#include "Variables.H"
#include "Initialize.H"
#include "calc.H"
#include "FunctionW.H"
#include "Thcoupling.H"
#include "anisotropy.H"
#include "interp.H"

using namespace amrex;
			

void advance(	MultiFab& phi_old, 
		MultiFab& phi_new,
		MultiFab& mu_old, 
		MultiFab& mu_new,
		MultiFab& term1, 
		MultiFab& term2,
		MultiFab& term3,
		MultiFab& derivx,
		MultiFab& derivy,
		MultiFab& h_phi,
		MultiFab& dh_dphi,
		MultiFab& psi,
		Geometry const& geom,
		Vector<BCRec> const& BoundaryCondition)
{
	//Fill the ghost cells
	phi_old.FillBoundary(geom.periodicity());		
	mu_old.FillBoundary(geom.periodicity());

//	FillDomainBoundary(phi_old,geom,BoundaryCondition);

//	FillDomainBoundary(mu_old,geom,BoundaryCondition);

	//Computing the interpolation function (Refer to interp.H)
	Func_H(h_phi, dh_dphi,phi_old);	

	h_phi.FillBoundary(geom.periodicity());
	dh_dphi.FillBoundary(geom.periodicity());

	//Computing the anisotropy term(term1) in the phi evolution equation (Refer to anisotropy.H for the formulation)
	computeterm1(term1, phi_old, derivx, derivy, geom);

	//Fill ghost cells for term1 along with periodic boundaries
	term1.FillBoundary(geom.periodicity());

	//Computing the Double well Potential(term2) in the phi evolution equation (Refer to FunctionW.H for the formulation of double well potential calculation)
	computeterm2(term2, phi_old);

	//Fill ghost cells for term2 along with periodic boundaries
	term2.FillBoundary(geom.periodicity());

	//Computing the Psi equation(term3) in the phi evolution equation (Refer to Thcoupling.H for the formulation of psi calculation)
	computeterm3(mu_old, term3, phi_old,dh_dphi,psi);

	//Fill ghost cells for term3 along with periodic boundaries
	term3.FillBoundary(geom.periodicity());

	//Now we have all the terms for terms for phi evolution, we simply add them (Refer to calc.H for the formulation of update_phi function)
	update_phi(phi_new, phi_old, term1,term2,term3);

	//Fill ghost cells of phi_new along with periodic boundaries
	phi_new.FillBoundary(geom.periodicity());

	//Phi is already updated now here we update mu (Refer to dmudt fucntion in calc.H for formulation)
	dmudt(mu_new, mu_old, phi_new, phi_old, derivx, derivy, h_phi, dh_dphi, geom);
	
}

#endif
