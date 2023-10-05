/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    phaseFieldDynamic

Description
    Solves the two-phase two-component coupled phase-field and chemical potential equations based on grand-potential formulation.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include <fstream>
#include <sstream>
//#include "scalarMatrices.H"
//#include "LUscalarMatrix.H"
//#include "RectangularMatrix.H"
//#include "SquareMatrix.H"
//#include "Random.H"
//! Using GSL for interpolating thermodynamic data
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
// Using dynamicFvMesh for dynamic interface refinement
//#include "dynamicFvMesh.H"

gsl_spline *spline_hess1[4];
gsl_spline *spline_hess2[4];
gsl_spline *spline_hess3[4];
gsl_spline *spline_comp1[4];
gsl_spline *spline_comp2[4];
gsl_spline *spline_L;
gsl_spline *spline_C;

gsl_interp_accel *acc_hess1[4];
gsl_interp_accel *acc_hess2[4];
gsl_interp_accel *acc_hess3[4];
gsl_interp_accel *acc_comp1[4];
gsl_interp_accel *acc_comp2[4];
gsl_interp_accel *acc10;
gsl_interp_accel *acc_L;
gsl_interp_accel *acc_C;

double calculate_melting_point(double cl, double T_guess, double T_high, double T_low) {
  //int iter=0;
  double tol=1e-6;
  dimensionedScalar T = T_guess;
  dimensionedScalar fval, dfdT;  
  
  do {
    dimensionedScalar A_Sol = gsl_spline_eval (spline_hess1[0], T.value(), acc_hess1[0]);
    dimensionedScalar A_Liq = gsl_spline_eval (spline_hess1[1], T.value(), acc_hess1[1]);
  
    dimensionedScalar A_SoldT = gsl_spline_eval_deriv (spline_hess1[0], T.value(), acc_hess1[0]);
    dimensionedScalar A_LiqdT = gsl_spline_eval_deriv (spline_hess1[1], T.value(), acc_hess1[1]);
  
    dimensionedScalar c_Sol = gsl_spline_eval (spline_comp1[0], T.value(), acc_comp1[0]);
    dimensionedScalar c_Liq = gsl_spline_eval (spline_comp1[1], T.value(), acc_comp1[1]);

    dimensionedScalar c_SoldT = gsl_spline_eval_deriv (spline_comp1[0], T.value(), acc_comp1[0]);
    dimensionedScalar c_LiqdT = gsl_spline_eval_deriv (spline_comp1[1], T.value(), acc_comp1[1]);
    
    dimensionedScalar B_Sol = 2.0*(A_Liq*c_Liq - A_Sol*c_Sol);
    dimensionedScalar dBdT = 2.0*(A_LiqdT*c_Liq + A_Liq*c_LiqdT - A_SoldT*c_Sol - A_Sol*c_SoldT);
    
    dimensionedScalar B_Liq = 0.0; //("B_Liq", 0.*mu);

    dimensionedScalar D_Sol = -A_Liq*c_Liq*c_Liq + A_Sol*c_Sol*c_Sol;
    dimensionedScalar D_Liq = 0.0; //("D_Liq", 0.*mu);
    
    dimensionedScalar D_Sol_dT = -(A_LiqdT)*c_Liq*c_Liq + (A_SoldT)*c_Sol*c_Sol -2.0*A_Liq*c_Liq*c_LiqdT + 2.0*A_Sol*c_Sol*c_SoldT;
    
    
  //Info << "Melting_point: "  << T.value() << endl;
    
    fval  = -0.25*(2.0*A_Liq*cl - B_Sol)*(2.0*A_Liq*cl - B_Sol)/A_Sol + A_Liq*cl*cl + D_Sol;
    dfdT  = -0.5*((2.0*A_Liq*cl - B_Sol)/A_Sol)*(2.0*A_LiqdT*cl - dBdT) +  D_Sol_dT +  (0.25*(2.0*A_Liq*cl - B_Sol)*(2.0*A_Liq*cl - B_Sol)/(A_Sol*A_Sol))*A_SoldT;
    dfdT += A_LiqdT*cl*cl;

    /*if ((T.value() > 870) || (T.value() < 793)) {
      T.value() = 870;
      break;
    }*/
    
    T.value() = T.value() - fval.value()/dfdT.value();

    if ((T.value() > T_high) || (T.value() < T_low)) {
      T.value() = T_high;
      break;
    }

  } while(fabs(fval.value()) > tol);
  return T.value();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Overloading function calculate_melting_point for ternary

double calculate_melting_point(scalarRectangularMatrix cl, double T_guess, double T_high, double T_low) {
  //int iter=0;
  double tol=1e-6;
  dimensionedScalar T = T_guess;
  dimensionedScalar fval, dfdT;

   scalarRectangularMatrix dmudc_a(2,2,0);
   scalarRectangularMatrix dmudc_l(2,2,0);
   scalarRectangularMatrix dmudc_adT(2,2,0);
   scalarRectangularMatrix dmudc_ldT(2,2,0);
   scalarRectangularMatrix dcdmu_a(2,2,0);
   scalarRectangularMatrix dcdmu_l(2,2,0);
   scalarRectangularMatrix dcdmu_adT(2,2,0);
   scalarRectangularMatrix dcdmu_ldT(2,2,0);
   scalarRectangularMatrix ceq_a(2,1,0);
   scalarRectangularMatrix ceq_l(2,1,0);
   scalarRectangularMatrix ceq_adT(2,1,0);
   scalarRectangularMatrix ceq_ldT(2,1,0);
  
  do {
    dmudc_a[0][0] = gsl_spline_eval (spline_hess1[0], T.value(), acc_hess1[0]);
    dmudc_a[0][1] = gsl_spline_eval (spline_hess2[0], T.value(), acc_hess2[0]);
    dmudc_a[1][1] = gsl_spline_eval (spline_hess3[0], T.value(), acc_hess3[0]);

    dmudc_adT[0][0] = gsl_spline_eval_deriv (spline_hess1[0], T.value(), acc_hess1[0]);
    dmudc_adT[0][1] = gsl_spline_eval_deriv (spline_hess2[0], T.value(), acc_hess2[0]);
    dmudc_adT[1][1] = gsl_spline_eval_deriv (spline_hess3[0], T.value(), acc_hess3[0]);
    
    dmudc_l[0][0] = gsl_spline_eval (spline_hess1[1], T.value(), acc_hess1[1]);
    dmudc_l[0][1] = gsl_spline_eval (spline_hess2[1], T.value(), acc_hess2[1]);
    dmudc_l[1][1] = gsl_spline_eval (spline_hess3[1], T.value(), acc_hess3[1]);

    dmudc_ldT[0][0] = gsl_spline_eval_deriv (spline_hess1[1], T.value(), acc_hess1[1]);
    dmudc_ldT[0][1] = gsl_spline_eval_deriv (spline_hess2[1], T.value(), acc_hess2[1]);
    dmudc_ldT[1][1] = gsl_spline_eval_deriv (spline_hess3[1], T.value(), acc_hess3[1]);
    
    ceq_a[0][0] = gsl_spline_eval (spline_comp1[0], T.value(), acc_comp1[0]);
    ceq_a[1][0] = gsl_spline_eval (spline_comp2[0], T.value(), acc_comp2[0]);

    ceq_adT[0][0] = gsl_spline_eval_deriv (spline_comp1[0], T.value(), acc_comp1[0]);
    ceq_adT[1][0] = gsl_spline_eval_deriv (spline_comp2[0], T.value(), acc_comp2[0]);
    
    ceq_l[0][0] = gsl_spline_eval (spline_comp1[1], T.value(), acc_comp1[1]);
    ceq_l[1][0] = gsl_spline_eval (spline_comp2[1], T.value(), acc_comp2[1]);

    ceq_ldT[0][0] = gsl_spline_eval_deriv (spline_comp1[1], T.value(), acc_comp1[1]);
    ceq_ldT[1][0] = gsl_spline_eval_deriv (spline_comp2[1], T.value(), acc_comp2[1]);

    dmudc_a[1][0] = dmudc_a[0][1];
    dmudc_l[1][0] = dmudc_l[0][1];

    dmudc_adT[1][0] = dmudc_adT[0][1];
    dmudc_ldT[1][0] = dmudc_ldT[0][1];

    //Info << T.value() << " " << ceq_l[0][0] << " " << ceq_l[1][0] << " " << Cp << endl;

    scalar detdmudc_a = (dmudc_a[0][0]*dmudc_a[1][1]-dmudc_a[0][1]*dmudc_a[0][1]);

    dcdmu_adT[0][0] = -dmudc_a[1][1]*(dmudc_adT[0][0]*dmudc_a[1][1] + dmudc_a[0][0]*dmudc_adT[1][1] - 2*dmudc_adT[0][1]*dmudc_a[0][1])/(detdmudc_a*detdmudc_a) + dmudc_adT[1][1]/detdmudc_a;
    dcdmu_adT[0][1] = dmudc_a[0][1]*(dmudc_adT[0][0]*dmudc_a[1][1] + dmudc_a[0][0]*dmudc_adT[1][1] - 2*dmudc_adT[0][1]*dmudc_a[0][1])/(detdmudc_a*detdmudc_a) - dmudc_adT[0][1]/detdmudc_a;
    dcdmu_adT[1][1] = -dmudc_a[0][0]*(dmudc_adT[0][0]*dmudc_a[1][1] + dmudc_a[0][0]*dmudc_adT[1][1] - 2*dmudc_adT[0][1]*dmudc_a[0][1])/(detdmudc_a*detdmudc_a) + dmudc_adT[0][0]/detdmudc_a;

    scalar detdmudc_l = (dmudc_l[0][0]*dmudc_l[1][1]-dmudc_l[0][1]*dmudc_l[0][1]);

    dcdmu_ldT[0][0] = -dmudc_l[1][1]*(dmudc_ldT[0][0]*dmudc_l[1][1] + dmudc_l[0][0]*dmudc_ldT[1][1] - 2*dmudc_ldT[0][1]*dmudc_l[0][1])/(detdmudc_l*detdmudc_l) + dmudc_ldT[1][1]/detdmudc_l;
    dcdmu_ldT[0][1] = dmudc_l[0][1]*(dmudc_ldT[0][0]*dmudc_l[1][1] + dmudc_l[0][0]*dmudc_ldT[1][1] - 2*dmudc_ldT[0][1]*dmudc_l[0][1])/(detdmudc_l*detdmudc_l) - dmudc_ldT[0][1]/detdmudc_l;
    dcdmu_ldT[1][1] = -dmudc_l[0][0]*(dmudc_ldT[0][0]*dmudc_l[1][1] + dmudc_l[0][0]*dmudc_ldT[1][1] - 2*dmudc_ldT[0][1]*dmudc_l[0][1])/(detdmudc_l*detdmudc_l) + dmudc_ldT[0][0]/detdmudc_l;

    dcdmu_a = SVDinv(dmudc_a);
    dcdmu_l = SVDinv(dmudc_l);
    
    
    dimensionedScalar B_a1 = (dmudc_l[0][0]*ceq_l[0][0] - dmudc_a[0][0]*ceq_a[0][0]) + (dmudc_l[0][1]*ceq_l[1][0] - dmudc_a[0][1]*ceq_a[1][0]);
   
   dimensionedScalar dB_a1dT = (dmudc_ldT[0][0]*ceq_l[0][0] + dmudc_l[0][0]*ceq_ldT[0][0] - dmudc_adT[0][0]*ceq_a[0][0] - dmudc_a[0][0]*ceq_adT[0][0]) + (dmudc_ldT[0][1]*ceq_l[1][0] + dmudc_l[0][1]*ceq_ldT[1][0] - dmudc_adT[0][1]*ceq_a[1][0] - dmudc_a[0][1]*ceq_adT[1][0]);
   
   dimensionedScalar B_a2 = (dmudc_l[1][1]*ceq_l[1][0] - dmudc_a[1][1]*ceq_a[1][0]) + (dmudc_l[0][1]*ceq_l[0][0] - dmudc_a[0][1]*ceq_a[0][0]);
   
   dimensionedScalar dB_a2dT = (dmudc_ldT[1][1]*ceq_l[1][0] + dmudc_l[1][1]*ceq_ldT[1][0] - dmudc_adT[1][1]*ceq_a[1][0] - dmudc_a[1][1]*ceq_adT[1][0]) + (dmudc_ldT[0][1]*ceq_l[0][0] + dmudc_l[0][1]*ceq_ldT[0][0] - dmudc_adT[0][1]*ceq_a[0][0] - dmudc_a[0][1]*ceq_adT[0][0]);
   
   //dimensionedScalar B_Liq = 0.0; //("B_Liq", 0.*mu);
   
   dimensionedScalar DD_a = -(0.5*dmudc_l[0][0]*ceq_l[0][0]*ceq_l[0][0] + 0.5*dmudc_l[1][1]*ceq_l[1][0]*ceq_l[1][0] + dmudc_l[0][1]*ceq_l[0][0]*ceq_l[1][0]) + (0.5*dmudc_a[0][0]*ceq_a[0][0]*ceq_a[0][0] + 0.5*dmudc_a[1][1]*ceq_a[1][0]*ceq_a[1][0] + dmudc_a[0][1]*ceq_a[0][0]*ceq_a[1][0]);
   
   //dimensionedScalar D_Liq = 0.0; //("D_Liq", 0.*mu);
   
   dimensionedScalar DD_a_dT = -((0.5*dmudc_ldT[0][0]*ceq_l[0][0]*ceq_l[0][0] + dmudc_l[0][0]*ceq_ldT[0][0]*ceq_l[0][0]) + (0.5*dmudc_ldT[1][1]*ceq_l[1][0]*ceq_l[1][0] + dmudc_l[1][1]*ceq_ldT[1][0]*ceq_l[1][0]) + (dmudc_ldT[0][1]*ceq_l[0][0]*ceq_l[1][0] + dmudc_l[0][1]*ceq_ldT[0][0]*ceq_l[1][0] + dmudc_l[0][1]*ceq_l[0][0]*ceq_ldT[1][0])) + ((0.5*dmudc_adT[0][0]*ceq_a[0][0]*ceq_a[0][0] + dmudc_a[0][0]*ceq_adT[0][0]*ceq_a[0][0]) + (0.5*dmudc_adT[1][1]*ceq_a[1][0]*ceq_a[1][0] + 0.5*dmudc_a[1][1]*ceq_adT[1][0]*ceq_a[1][0]) + (dmudc_adT[0][1]*ceq_a[0][0]*ceq_a[1][0] + dmudc_a[0][1]*ceq_adT[0][0]*ceq_a[1][0] + dmudc_a[0][1]*ceq_a[0][0]*ceq_adT[1][0]));
    
    dimensionedScalar mu1_trial = (dmudc_l[0][0]*cl[0][0] + dmudc_l[0][1]*cl[1][0] - B_a1);
    dimensionedScalar mu2_trial = (dmudc_l[0][1]*cl[0][0] + dmudc_l[1][1]*cl[1][0] - B_a2);

    dimensionedScalar mu1_trialdT = (dmudc_ldT[0][0]*cl[0][0] + dmudc_ldT[0][1]*cl[1][0] - dB_a1dT);
    dimensionedScalar mu2_trialdT = (dmudc_ldT[0][1]*cl[0][0] + dmudc_ldT[1][1]*cl[1][0] - dB_a2dT);

   fval  = -(0.5*dmudc_a[0][0]*(dcdmu_a[0][0]*mu1_trial + dcdmu_a[0][1]*mu2_trial)*(dcdmu_a[0][0]*mu1_trial + dcdmu_a[0][1]*mu2_trial)
    + 0.5*dmudc_a[1][1]*(dcdmu_a[0][1]*mu1_trial + dcdmu_a[1][1]*mu2_trial)*(dcdmu_a[0][1]*mu1_trial + dcdmu_a[1][1]*mu2_trial)
    + dmudc_a[0][1]*(dcdmu_a[0][0]*mu1_trial + dcdmu_a[0][1]*mu2_trial)*(dcdmu_a[0][1]*mu1_trial + dcdmu_a[1][1]*mu2_trial))
    + (0.5*dmudc_l[0][0]*cl[0][0]*cl[0][0] + 0.5*dmudc_l[1][1]*cl[1][0]*cl[1][0] + dmudc_l[0][1]*cl[0][0]*cl[1][0]) + DD_a;

    dfdT  = -(0.5*dmudc_adT[0][0]*(dcdmu_a[0][0]*mu1_trial + dcdmu_a[0][1]*mu2_trial)*(dcdmu_a[0][0]*mu1_trial + dcdmu_a[0][1]*mu2_trial)
    + 0.5*dmudc_adT[1][1]*(dcdmu_a[0][1]*mu1_trial + dcdmu_a[1][1]*mu2_trial)*(dcdmu_a[0][1]*mu1_trial + dcdmu_a[1][1]*mu2_trial)
    + dmudc_adT[0][1]*(dcdmu_a[0][0]*mu1_trial + dcdmu_a[0][1]*mu2_trial)*(dcdmu_a[0][1]*mu1_trial + dcdmu_a[1][1]*mu2_trial)
    + dmudc_a[0][0]*(dcdmu_a[0][0]*mu1_trial + dcdmu_a[0][1]*mu2_trial)*(dcdmu_adT[0][0]*mu1_trial + dcdmu_adT[0][1]*mu2_trial + dcdmu_a[0][0]*mu1_trialdT + dcdmu_a[0][1]*mu2_trialdT)
    + dmudc_a[1][1]*(dcdmu_a[0][1]*mu1_trial + dcdmu_a[1][1]*mu2_trial)*(dcdmu_adT[0][1]*mu1_trial + dcdmu_adT[1][1]*mu2_trial + dcdmu_a[0][1]*mu1_trialdT + dcdmu_a[1][1]*mu2_trialdT)
    + dmudc_a[0][1]*((dcdmu_adT[0][0]*mu1_trial + dcdmu_adT[0][1]*mu2_trial + dcdmu_a[0][0]*mu1_trialdT + dcdmu_a[0][1]*mu2_trialdT)*(dcdmu_a[0][1]*mu1_trial + dcdmu_a[1][1]*mu2_trial) + (dcdmu_a[0][0]*mu1_trial + dcdmu_a[0][1]*mu2_trial)*(dcdmu_adT[0][1]*mu1_trial + dcdmu_adT[1][1]*mu2_trial + dcdmu_a[0][1]*mu1_trialdT + dcdmu_a[1][1]*mu2_trialdT))
    );

    dfdT += DD_a_dT + (0.5*dmudc_ldT[0][0]*cl[0][0]*cl[0][0] + 0.5*dmudc_ldT[1][1]*cl[1][0]*cl[1][0] + dmudc_ldT[0][1]*cl[0][0]*cl[1][0]);

    /*if ((T.value() > 870) || (T.value() < 793)) {
      T.value() = 870;
      break;
    }*/
    
    T.value() = T.value() - fval.value()/dfdT.value();

    if ((T.value() > T_high) || (T.value() < T_low)) {
      T.value() = T_high;
      break;
    }

  } while(fabs(fval.value()) > tol);
  return T.value();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    //#include "createDynamicFvMesh.H"
    #include "createMesh.H"
    #include "createControls.H"
    
    simpleControl simple(mesh);

    #include "createFields.H"
    #include "initializeFields.H"
    #include "createTol.H"

    int iter_num = 0;
    int nucleation_event =1;


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating phase-field and chemical potential distributions\n" << endl;

   //! The scalar and vector fields initialized below
   volVectorField q=dimx*fvc::grad(phi_1);
   volVectorField q2=dimx*fvc::grad(phi_2);
   volVectorField q3=dimx*fvc::grad(phi_3);
   volVectorField q4=dimx*fvc::grad(phi_4);
   volScalarField magq = 0.0*phi_1;
   //volScalarField T = 0.0*phi + initial;
   //volScalarField sumLHS = 0.0*phi;
   //volVectorField q_4 = 0.0*q;
   //volVectorField q_6 = 0.0*q;
   volScalarField ac_01 =0.0*phi_1;
   volScalarField ac_02 =0.0*phi_2;
   volScalarField ac_03 =0.0*phi_3;
   volScalarField ac_04 =0.0*phi_4;
   volVectorField dAdq01= phi_1*vector(0,0,0);
   volVectorField dadgradPhi=q*0.0;
   volScalarField grad_qt_sqr = 0.0*phi_1;
    //volScalarField A_Sol = 0.0*phi;
    //volScalarField A_Liq = 0.0*phi;
    //volScalarField c_Sol = 0.0*phi;
    //volScalarField c_Liq = 0.0*phi;
    //volScalarField dcdmu_Liq = 0.0*phi;
    //volScalarField omega = 0.0*phi;

//! The unit normal vector to the interface with a small number in denominator to prevent solution from diverging
//n=dimx*fvc::grad(phi_1)/(1E-20+mag(dimx*fvc::grad(phi_1)));
volVectorField n=dimx*fvc::grad(phi_1);
volVectorField n2=dimx*fvc::grad(phi_2);
volVectorField n3=dimx*fvc::grad(phi_3);
volVectorField n4=dimx*fvc::grad(phi_4);

volScalarField hphi = 0.0*phi_1; //-phi*phi*phi*(10*phi*phi*phi - 36*phi*phi + 45*phi - 20);
volScalarField hphi2 = 0.0*phi_2;
volScalarField hphi3 = 0.0*phi_3;
volScalarField hphi4 = 0.0*phi_4;
    
    //! Done reading thermodynamic database
    
    //! Allocation for GSL acceleration and interpolation
    
    //if (phases == 2)
    //{
    double *temp_array = new double[np];
    
    if (components == 2)
    {
    for (i_phase = 0; i_phase < phases; i_phase++){
    for (int i = 0; i < np; i++){
    temp_array[i] = A_arr[i][i_phase];
    }
    spline_hess1[i_phase] = gsl_spline_alloc (gsl_interp_cspline, np);
    acc_hess1[i_phase] = gsl_interp_accel_alloc ();
    gsl_spline_init (spline_hess1[i_phase], T1, temp_array, np);
    
    for (int i = 0; i < np; i++){
    temp_array[i] = c_arr[i][i_phase];
    }
    spline_comp1[i_phase] = gsl_spline_alloc (gsl_interp_cspline, np);
    acc_comp1[i_phase] = gsl_interp_accel_alloc ();
    gsl_spline_init (spline_comp1[i_phase], T1, temp_array, np);
    }

    spline_L = gsl_spline_alloc (gsl_interp_cspline, np1);
    acc_L = gsl_interp_accel_alloc ();
    gsl_spline_init (spline_L, T2, Lf_arr, np1);

    spline_C = gsl_spline_alloc (gsl_interp_cspline, np1);
    acc_C = gsl_interp_accel_alloc ();
    gsl_spline_init (spline_C, T2, Cp_arr, np1);
    }
    else if (components == 3)
    {
    for (i_phase = 0; i_phase < phases; i_phase++){
    for (int i = 0; i < np; i++){
    temp_array[i] = H11_arr[i][i_phase];
    }
    spline_hess1[i_phase] = gsl_spline_alloc (gsl_interp_cspline, np);
    acc_hess1[i_phase] = gsl_interp_accel_alloc ();
    gsl_spline_init (spline_hess1[i_phase], T1, temp_array, np);

    for (int i = 0; i < np; i++){
    temp_array[i] = H12_arr[i][i_phase];
    }
    spline_hess2[i_phase] = gsl_spline_alloc (gsl_interp_cspline, np);
    acc_hess2[i_phase] = gsl_interp_accel_alloc ();
    gsl_spline_init (spline_hess2[i_phase], T1, temp_array, np);

    for (int i = 0; i < np; i++){
    temp_array[i] = H22_arr[i][i_phase];
    }
    spline_hess3[i_phase] = gsl_spline_alloc (gsl_interp_cspline, np);
    acc_hess3[i_phase] = gsl_interp_accel_alloc ();
    gsl_spline_init (spline_hess3[i_phase], T1, temp_array, np);
    
    for (int i = 0; i < np; i++){
    temp_array[i] = c1_arr[i][i_phase];
    }
    spline_comp1[i_phase] = gsl_spline_alloc (gsl_interp_cspline, np);
    acc_comp1[i_phase] = gsl_interp_accel_alloc ();
    gsl_spline_init (spline_comp1[i_phase], T1, temp_array, np);
    
    for (int i = 0; i < np; i++){
    temp_array[i] = c2_arr[i][i_phase];
    }
    spline_comp2[i_phase] = gsl_spline_alloc (gsl_interp_cspline, np);
    acc_comp2[i_phase] = gsl_interp_accel_alloc ();
    gsl_spline_init (spline_comp2[i_phase], T1, temp_array, np);
    }
    
    spline_L = gsl_spline_alloc (gsl_interp_cspline, np1);
    acc_L = gsl_interp_accel_alloc ();
    gsl_spline_init (spline_L, T2, Lf_arr, np1);

    spline_C = gsl_spline_alloc (gsl_interp_cspline, np1);
    acc_C = gsl_interp_accel_alloc ();
    gsl_spline_init (spline_C, T2, Cp_arr, np1);
    }
    //}
    
    delete [] temp_array;
    
    //! Current temperature
    dimensionedScalar T = initial;
    dimensionedScalar Tm = T0;
    
    //! File to write temperature evolution
    const fileName pathToTemp = "tempvstime.dat";
    ofstream outpft(pathToTemp);
    if (!outpft.is_open())
    {
        Info << "No output file found" << endl;
    }
    
   scalar InitialResidual_0 = 0.0;
   scalar InitialResidual_10 = 0.0;
   scalar InitialResidual_11 = 0.0;
   scalar InitialResidual_20 = 0.0;
   scalar InitialResidual_21 = 0.0;
   scalar InitialResidual_22 = 0.0;
   scalar InitialResidual_23 = 0.0;
   scalar InitialResidual_3 = 0.0;
   
   scalar InitialDeltaT = runTime.deltaTValue();
   
   int iCorr = 0;
   
   //! For binary
   double A[int(phases)];
   double A_dT[int(phases)];
   double c_bin[int(phases)];
   double c_bin_dT[int(phases)];
   double B_bin[int(phases)];
   double B_bin_dT[int(phases)];
   double D_bin[int(phases)];
   
   dimensionedScalar DeltaC = 0.0;
        
        dimensionedScalar dcdmu_Liq = 0.0;
            
   dimensionedScalar omega = 0.0;

        
        dimensionedScalar avg_liq_conc = 0.0;
   
    //if (phases == 2)
    //{
    if (components == 2)
    {
    dimensionedScalar c_a = gsl_spline_eval (spline_comp1[0], Tm.value(), acc_comp1[0]);
    dimensionedScalar c_l = gsl_spline_eval (spline_comp1[int(phases)-1], Tm.value(), acc_comp1[int(phases)-1]);
    
    //Info << T << "AS " << A_Sol << "AL " << A_Liq << "cS " << c_Sol << "cL " << c_Liq << endl;
    
    DeltaC = c_l-c_a;
    }
    //}
    
    //! For ternary
   double H11[int(phases)];
   double H12[int(phases)];
   double H22[int(phases)];
   double H11_dT[int(phases)];
   double H12_dT[int(phases)];
   double H22_dT[int(phases)];
   double c_ter1[int(phases)];
   double c_ter2[int(phases)];
   double c_ter1_dT[int(phases)];
   double c_ter2_dT[int(phases)];
   double B_ter1[int(phases)-1];
   double B_ter2[int(phases)-1];
   double B_ter1_dT[int(phases)-1];
   double B_ter2_dT[int(phases)-1];
   double D_ter[int(phases)-1];
   
   scalarRectangularMatrix dmudc(2*phases,2,0);
   scalarRectangularMatrix dmudc_dT(2*phases,2,0);
   scalarRectangularMatrix dcdmu(2*phases,2,0);
   scalarRectangularMatrix dcdmu_dT(2*phases,2,0);

   scalarRectangularMatrix Mob(2*phases,2,0);
   
   scalarSquareMatrix dmudc_temp(2,0);
   scalarSquareMatrix dmudc_dT_temp(2,0);
   scalarRectangularMatrix dcdmu_temp(2,0);
   scalarSquareMatrix dcdmu_dT_temp(2,0);
   scalarRectangularMatrix diff_temp(2,2,0);
   scalarRectangularMatrix M_temp(2,2,0);

   scalarRectangularMatrix c_at(2,1,0);
   scalarRectangularMatrix c_lt(2,1,0);
   scalarRectangularMatrix DeltaCt(2,1,0);
   scalarRectangularMatrix avg_liq_conct(2,1,0);
   scalarRectangularMatrix testMatrix(2,2,0);
   scalarRectangularMatrix testMatrix_1(2,2,0);

        scalarRectangularMatrix M_l(2,2,0);

    
    //if (phases == 2)
    //{
    if (components == 3)
    {
    c_at[0][0] = gsl_spline_eval (spline_comp1[0], Tm.value(), acc_comp1[0]);
    c_at[1][0] = gsl_spline_eval (spline_comp2[0], Tm.value(), acc_comp2[0]);
    
    c_lt[0][0] = gsl_spline_eval (spline_comp1[int(phases)-1], Tm.value(), acc_comp1[int(phases)-1]);
    c_lt[1][0] = gsl_spline_eval (spline_comp2[int(phases)-1], Tm.value(), acc_comp2[int(phases)-1]);
   
   DeltaCt = c_lt-c_at;

   }
   //}

    dimensionedScalar totalVol = fvc::domainIntegrate(1-0.0*phi_1);
   
    //! Total volume of solidus
   dimensionedScalar solidVolFracOld = fvc::domainIntegrate(phi_1)/totalVol;
   dimensionedScalar solidVolFrac = fvc::domainIntegrate(phi_1)/totalVol;
   
   dimensionedScalar dsolidVolFrac = (solidVolFrac - solidVolFracOld);
   dimensionedScalar dsolidVolFracOld = dsolidVolFrac;
    
    
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        #include "variabledeltaT.H"

        #include "readSolidDisplacementFoamControls.H"
        
            // The imposed temperature field as a function of thermal gradient
            // in the x direction, G and pulling velocity, v
            //T = G*( (1/dimx)* mesh.C().component(vector::X) - v*runTime.value()) + initial;
        
        //! Computing temperature evolution
    
        #include "TEqn.H"

        //! Computing the thermodynamic parameters in phase diagram
    //if (phases == 2)
    //{
        if (components == 2)
    	{
    	for (i_phase = 0; i_phase < phases; i_phase++){
        A[i_phase] = gsl_spline_eval (spline_hess1[i_phase], T.value(), acc_hess1[i_phase]);
        A_dT[i_phase] = gsl_spline_eval_deriv (spline_hess1[i_phase], T.value(), acc_hess1[i_phase]);
        c_bin[i_phase] = gsl_spline_eval (spline_comp1[i_phase], T.value(), acc_comp1[i_phase]);
        c_bin_dT[i_phase] = gsl_spline_eval_deriv (spline_comp1[i_phase], T.value(), acc_comp1[i_phase]);
        }
        
        
        dcdmu_Liq = 0.5/A[int(phases)-1];
            
        omega = epsilon*0.18*DeltaC*DeltaC/(dcdmu_Liq*diff_Liq*Vm);

	//Info<< "omega= " << omega.value() << nl << endl;
        
        for (i_phase = 0; i_phase < int(phases)-1; i_phase++){
        B_bin[i_phase] = 2.0*(A[int(phases)-1]*c_bin[int(phases)-1] - A[i_phase]*c_bin[i_phase]);
        B_bin_dT[i_phase] = 2.0*(A_dT[int(phases)-1]*c_bin[int(phases)-1] + A[int(phases)-1]*c_bin_dT[int(phases)-1] - A_dT[i_phase]*c_bin[i_phase] - A[i_phase]*c_bin_dT[i_phase]);
        
        D_bin[i_phase] = -A[int(phases)-1]*c_bin[int(phases)-1]*c_bin[int(phases)-1] + A[i_phase]*c_bin[i_phase]*c_bin[i_phase];
        }
        
        B_bin[int(phases)-1] = 0.0; //("B_Liq", 0.*mu);
        
        D_bin[int(phases)-1] = 0.0; //("D_Liq", 0.*mu);
        //////////////////////

	}
        
        if (components == 3)
    	{
    	//i_phase = 0;
    	for (i_phase = 0; i_phase < phases; i_phase++){
    H11[i_phase] = gsl_spline_eval (spline_hess1[i_phase], T.value(), acc_hess1[i_phase]);
    H12[i_phase] = gsl_spline_eval (spline_hess2[i_phase], T.value(), acc_hess2[i_phase]);
    H22[i_phase] = gsl_spline_eval (spline_hess3[i_phase], T.value(), acc_hess3[i_phase]);

    H11_dT[i_phase] = gsl_spline_eval_deriv (spline_hess1[i_phase], T.value(), acc_hess1[i_phase]);
    H12_dT[i_phase] = gsl_spline_eval_deriv (spline_hess2[i_phase], T.value(), acc_hess2[i_phase]);
    H22_dT[i_phase] = gsl_spline_eval_deriv (spline_hess3[i_phase], T.value(), acc_hess3[i_phase]);
    
    c_ter1[i_phase] = gsl_spline_eval (spline_comp1[i_phase], T.value(), acc_comp1[i_phase]);
    c_ter2[i_phase] = gsl_spline_eval (spline_comp2[i_phase], T.value(), acc_comp2[i_phase]);

    c_ter1_dT[i_phase] = gsl_spline_eval_deriv (spline_comp1[i_phase], T.value(), acc_comp1[i_phase]);
    c_ter2_dT[i_phase] = gsl_spline_eval_deriv (spline_comp2[i_phase], T.value(), acc_comp2[i_phase]);
    
    //Info << T.value() << " " << c_ter1[i_phase] << " " << c_ter2[i_phase] << " " << H11[i_phase] << " " << H22[i_phase] << endl;
    }

    
    for (i_phase = 0; i_phase < phases; i_phase++){
    
    dmudc_temp[0][0] = H11[i_phase];
    dmudc_temp[0][1] = H12[i_phase];
    dmudc_temp[1][0] = H12[i_phase];
    dmudc_temp[1][1] = H22[i_phase];
    
    //Info << T.value() << " " << dmudc_temp[0][0] << " " << dmudc_temp[0][1] << " " << dmudc_temp[1][0] << " " << dmudc_temp[1][1] << endl;
    
    scalar detdmudc = (H11[i_phase]*H22[i_phase]-H12[i_phase]*H12[i_phase]);
    
    dcdmu_dT[2*i_phase][0] = -H22[i_phase]*(H11_dT[i_phase]*H22[i_phase] + H11[i_phase]*H22_dT[i_phase] - 2*H12_dT[i_phase]*H12[i_phase])/(detdmudc*detdmudc) + H22_dT[i_phase]/detdmudc;
    dcdmu_dT[2*i_phase][1] = H12[i_phase]*(H11_dT[i_phase]*H22[i_phase] + H11[i_phase]*H22_dT[i_phase] - 2*H12_dT[i_phase]*H12[i_phase])/(detdmudc*detdmudc) - H12_dT[i_phase]/detdmudc;
    dcdmu_dT[2*i_phase+1][0] = dcdmu_dT[2*i_phase][1];
    dcdmu_dT[2*i_phase+1][1] = -H11[i_phase]*(H11_dT[i_phase]*H22[i_phase] + H11[i_phase]*H22_dT[i_phase] - 2*H12_dT[i_phase]*H12[i_phase])/(detdmudc*detdmudc) + H11_dT[i_phase]/detdmudc;
    
    dcdmu_temp = SVDinv(dmudc_temp);
    
    dcdmu[2*i_phase][0] = dcdmu_temp[0][0];
    dcdmu[2*i_phase][1] = dcdmu_temp[0][1];
    dcdmu[2*i_phase+1][0] = dcdmu_temp[1][0];
    dcdmu[2*i_phase+1][1] = dcdmu_temp[1][1];
    
    //Info << T.value() << " " << dcdmu_dT[2*i_phase][0] << " " << dcdmu_dT[2*i_phase][1] << " " << dcdmu_dT[2*i_phase+1][0] << " " << dcdmu_dT[2*i_phase+1][1] << endl;
    //Info << T.value() << " " << dcdmu_temp[0][0] << " " << dcdmu_temp[0][1] << " " << dcdmu_temp[1][0] << " " << dcdmu_temp[1][1] << endl;
    
    diff_temp[0][0] = diff_ter[2*i_phase][0];
    diff_temp[0][1] = diff_ter[2*i_phase][1];
    diff_temp[1][0] = diff_ter[2*i_phase+1][0];
    diff_temp[1][1] = diff_ter[2*i_phase+1][1];
    
    //Info << T.value() << " " << diff_temp[0][0] << " " << diff_temp[0][1] << " " << diff_temp[1][0] << " " << diff_temp[1][1] << endl;
    
    multiply(M_temp, diff_temp, dcdmu_temp);
    
    Mob[2*i_phase][0] = M_temp[0][0];
    Mob[2*i_phase][1] = M_temp[0][1];
    Mob[2*i_phase+1][0] = M_temp[1][0];
    Mob[2*i_phase+1][1] = M_temp[1][1];
    }
	
   M_l[0][0]   = Mob[2*(int(phases)-1)][0];      
   M_l[0][1]   = Mob[2*(int(phases)-1)][1];
   M_l[1][0]   = Mob[2*(int(phases)-1)+1][0];
   M_l[1][1]   = Mob[2*(int(phases)-1)+1][1];
   
      testMatrix = SVDinv(M_l);
   
   multiply(testMatrix_1,DeltaCt.T(), testMatrix);
   
   multiply(testMatrix, testMatrix_1, DeltaCt);
   
//    transportProperties.set("omega",epsilon*0.18*testMatrix[0][0]);
   
   omega = epsilon*0.18*testMatrix[0][0]/Vm;
   
   //Info<< "omega= " << omega.value() << nl << endl;
      
   for (i_phase = 0; i_phase < int(phases)-1; i_phase++){
   
   B_ter1[i_phase] = (H11[int(phases)-1]*c_ter1[int(phases)-1] - H11[i_phase]*c_ter1[i_phase]) + (H12[int(phases)-1]*c_ter2[int(phases)-1] - H12[i_phase]*c_ter2[i_phase]);
   
   B_ter1_dT[i_phase] = (H11_dT[int(phases)-1]*c_ter1[int(phases)-1] + H11[int(phases)-1]*c_ter1_dT[int(phases)-1] - H11_dT[i_phase]*c_ter1[i_phase] - H11[i_phase]*c_ter1_dT[i_phase]) + (H12_dT[int(phases)-1]*c_ter2[int(phases)-1] + H12[int(phases)-1]*c_ter2_dT[int(phases)-1] - H12_dT[i_phase]*c_ter2[i_phase] - H12[i_phase]*c_ter2_dT[i_phase]);
   
   B_ter2[i_phase] = (H22[int(phases)-1]*c_ter2[int(phases)-1] - H22[i_phase]*c_ter2[i_phase]) + (H12[int(phases)-1]*c_ter1[int(phases)-1] - H12[i_phase]*c_ter1[i_phase]);
   
   B_ter2_dT[i_phase] = (H22_dT[int(phases)-1]*c_ter2[int(phases)-1] + H22[int(phases)-1]*c_ter2_dT[int(phases)-1] - H22_dT[i_phase]*c_ter2[i_phase] - H22[i_phase]*c_ter2_dT[i_phase]) + (H12_dT[int(phases)-1]*c_ter1[int(phases)-1] + H12[int(phases)-1]*c_ter1_dT[int(phases)-1] - H12_dT[i_phase]*c_ter1[i_phase] - H12[i_phase]*c_ter1_dT[i_phase]);
   
   D_ter[i_phase] = -(0.5*H11[int(phases)-1]*c_ter1[int(phases)-1]*c_ter1[int(phases)-1] + 0.5*H22[int(phases)-1]*c_ter2[int(phases)-1]*c_ter2[int(phases)-1] + H12[int(phases)-1]*c_ter1[int(phases)-1]*c_ter2[int(phases)-1]) + (0.5*H11[i_phase]*c_ter1[i_phase]*c_ter1[i_phase] + 0.5*H22[i_phase]*c_ter2[i_phase]*c_ter2[i_phase] + H12[i_phase]*c_ter1[i_phase]*c_ter2[i_phase]);
   
   //Info << T.value() << " " << B_ter1[i_phase] << " " << B_ter2[i_phase] << " " << B_ter1_dT[i_phase] << " " << B_ter2_dT[i_phase] << " " << D_ter[i_phase] << endl;
   }
   
    	}
    	//}


          //! For orthogonal correction of the finite volume mesh
          while (simple.correctNonOrthogonal())
        {


          /*label curTimeIndex = mesh.time().timeIndex();

	for (int i=0; i<30; i++)
        {
		if(curTimeIndex == 0)
				{
                 //! The interpolation equation is used for smoothening of the
                 //! phase-field variable
                                    #include "smoothTheta.H"
                                }

	}*/

                //mesh.update();
                //! Solving the phase-field and chemical potential equations after
                //! updating the mesh
                #include "alphaEqnp2c2.H"
                #include "alphaEqnp2c3.H"
                #include "alphaEqnp3c2.H"
                #include "alphaEqnp3c3.H"
                #include "alphaEqnp4c2.H"
                #include "alphaEqnp4c3.H"
                
                #include "muEqnp2c2.H"
                #include "muEqnp2c3.H"
                #include "muEqnp3c2.H"
                #include "muEqnp3c3.H"
                #include "muEqnp4c2.H"
                #include "muEqnp4c3.H"

                #include "thetaEqn.H"



	}
    
	#include "DEqnp2.H"
	#include "DEqnp3.H"
	#include "DEqnp4.H"

    //#include "nucleateFields.H"
    //for (i_phase = 0; i_phase < phases; i_phase++){
    if (components == 2)
    {
    c_1 = 0.5*(mu_1-B_bin[0])/A[0];
    }
    if (components == 3)
    {      
      c_1 = dcdmu[0][0]*(mu_1 - B_ter1[0]) + dcdmu[0][1]*(mu_2 - B_ter2[0]);
      
      c_2 = dcdmu[1][0]*(mu_1 - B_ter1[0]) + dcdmu[1][1]*(mu_2 - B_ter2[0]);
    }
    //}
    

    if (phases == 2)
    {
    if (components == 2)
    {
    avg_liq_conc = fvc::domainIntegrate((0.5*(mu_1-B_bin[int(phases)-1])/A[int(phases)-1])*(1-phi_1)).value()/fvc::domainIntegrate((1-phi_1)).value();
    }
    
    if (components == 3)
    {
    avg_liq_conct[0][0] = fvc::domainIntegrate((dcdmu[2*(int(phases)-1)][0]*(mu_1) + dcdmu[2*(int(phases)-1)][1]*(mu_2))*(1-phi_1)).value()/fvc::domainIntegrate((1-phi_1)).value();

      avg_liq_conct[1][0] = fvc::domainIntegrate((dcdmu[2*(int(phases)-1)+1][0]*(mu_1) + dcdmu[2*(int(phases)-1)+1][1]*(mu_2))*(1-phi_1)).value()/fvc::domainIntegrate((1-phi_1)).value();
    }
    }
    
    //Info << "avg_liq_conc: " << avg_liq_conc.value() << endl;
    
    //! Initial conditions for cooling simulation
if ((swcool == 1)&&(swch == 1))
{
    if (phases == 2)
    {
    if (components == 2)
    {
    Tm.value() = calculate_melting_point(avg_liq_conc.value(), Tm.value(), T1[np-1], T1[0]);
      
       Info << "Melting_point: "  << Tm.value() << " " << avg_liq_conc.value() << endl;
    }
    
    if (components == 3)
    {
      Tm.value() = calculate_melting_point(avg_liq_conct, Tm.value(), T1[np-1], T1[0]);
      
       Info << "Melting_point: "  << Tm.value() << " " << avg_liq_conct[0][0] << " " << avg_liq_conct[1][0] << endl;
    }
    }
       
    if ((T.value() < Tm.value())&&(nucleation_event < numSeeds)) {
        scalar deltaTemp = Tm.value() - T.value();
        scalar nprob;
        if (deltaTemp > 0.0) {
          nprob = (4.5/((dt_s)*Foam::pow((2*M_PI),0.5)))*(1/(deltaTemp))*Foam::exp(-0.5*Foam::pow(((Foam::log(deltaTemp)-Foam::log(dt_0))/dt_s),2));
        } else {
          nprob = 0.0;
        }
        scalar random = randNumber.globalScalar01();
        if (random < nprob) {
          if ((int(runTime.value()/runTime.deltaTValue()))%(int(nucleation_interval))==0) {
              //#include "nucleateSeeds.H"
              //theta_fill = 2.0*M_PI*randNumber.scalar01();
    	      //r          = (droplet_radius)*randNumber.scalar01();

              Info<< "creating a new seed" << endl;
            xCenter = randNumber.globalScalar01()*xMax;
            yCenter = randNumber.globalScalar01()*yMax;
            //xCenter = 0.5*xMax;
            //yCenter = 0.5*yMax;
            if (phases == 2)
            {
            if (dimensions == 2)
            {
            Info<< "xCenter, yCenter: " << xCenter << ", " << yCenter << endl;
            randTheta[2] = randNumber.globalScalar01()*(pi.value()/2);
            Info<< "random theta: " << randTheta[2] << endl;

            Info<< "Filling phi and theta fields in seeds" << endl;
            volScalarField gaussianSeed = (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius[0][0]*seedRadius[0][0]));

            theta += randTheta[2]*gaussianSeed*vector(0,0,1);
            phi_1 += gaussianSeed;
            }

            if (dimensions == 3)
            {
            zCenter = randNumber.globalScalar01()*zMax;
            //xCenter = 0.5*xMax;
            //yCenter = 0.5*yMax;
            Info<< "xCenter, yCenter, zCenter: " << xCenter << ", " << yCenter << ", " << zCenter << endl;
            randTheta[2] = randNumber.globalScalar01()*(pi.value()/4);
            Info<< "random thetaz: " << randTheta[2] << endl;
            randTheta[0] = randNumber.globalScalar01()*(pi.value()/4);
            randTheta[1] = randNumber.globalScalar01()*(pi.value()/4);
            Info<< "random thetax: " << randTheta[0] << endl;
            Info<< "random thetay: " << randTheta[1] << endl;

            Info<< "Filling phi and theta fields in seeds" << endl;
            volScalarField gaussianSeed = (1-phi_1)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(seedRadius[0][0]*seedRadius[0][0]));


            theta += gaussianSeed*(randTheta[0]*vector(1,0,0) + randTheta[1]*vector(0,1,0) + randTheta[2]*vector(0,0,1));
            phi_1 += gaussianSeed;
            }
            }

              nucleation_event++;
          }
        }
      }
}

    //! Writing the results according to keywords in controlDict
    //if (runTime.writeTime())
    //{
    //  runTime.write();
    //}
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        runTime.write();
        
    }

    Info<< "End\n" << endl;
    
    //! releasing reserved memory from GSL interpolation
    for (i_phase = 0; i_phase < 4; i_phase++){
    gsl_spline_free (spline_hess1[i_phase]);
    gsl_spline_free (spline_hess2[i_phase]);
    gsl_spline_free (spline_hess3[i_phase]);
    gsl_spline_free (spline_comp1[i_phase]);
    gsl_spline_free (spline_comp2[i_phase]);
    gsl_interp_accel_free (acc_hess1[i_phase]);
    gsl_interp_accel_free (acc_hess2[i_phase]);
    gsl_interp_accel_free (acc_hess3[i_phase]);
    gsl_interp_accel_free (acc_comp1[i_phase]);
    gsl_interp_accel_free (acc_comp2[i_phase]);    
    }
    
    gsl_spline_free (spline_L);
    gsl_spline_free (spline_C);
    
    gsl_interp_accel_free (acc_L);
    gsl_interp_accel_free (acc_C);


    outpft.close();
    
    return 0;
}
