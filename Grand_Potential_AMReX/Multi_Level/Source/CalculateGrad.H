#ifndef Gradient_H_
#define Gradient_H_

using namespace amrex;

//####################################################################################################################
void Gradient_x(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& phiOld, 
        amrex::Array4<amrex::Real> const& grad_x)
    {    
        grad_x(i,j,k) = amrex::Math::abs(phiOld(i+1,j,k) - phiOld(i,j,k));
    }    

//####################################################################################################################
void Gradient_y(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& phiOld, 
        amrex::Array4<amrex::Real> const& grad_y)
    {
        grad_y(i,j,k) = amrex::Math::abs(phiOld(i,j+1,k) - phiOld(i,j,k));
    } 

//####################################################################################################################
void Gradient_z(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& phiOld, 
        amrex::Array4<amrex::Real> const& grad_z)
    {
        grad_z(i,j,k) = amrex::Math::abs(phiOld(i,j,k+1) - phiOld(i,j,k));
    } 

//####################################################################################################################
#endif
