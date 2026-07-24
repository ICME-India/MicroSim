#ifndef MYFUNC_H_
#define MYFUNC_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>

using namespace amrex;

void main_main ();

void advance (amrex::MultiFab& phi_old,
              amrex::MultiFab& phi_new,
	      amrex::MultiFab& lap_temp,
	      amrex::MultiFab& dfdcon,
              amrex::Real& K,
              amrex::Real& M,
              amrex::Real& dt,
              amrex::Geometry const& geom);

void init_phi (amrex::MultiFab& phi_new, amrex::Geometry const& geom, amrex::Real c0);

//void compute_energy(amrex::MultiFab const& phiold, Real const& energ, amrex::Real const& K, amrex::Geometry const& geom);

//void derivative(amrex::MultiFab& phi_old,
//		amrex::MultiFab& dfdcon);

//void laplacian(amrex::MultiFab& phi_old,
//		amrex::MultiFab& lap_temp, amrex::Geometry const& geom);

#endif
