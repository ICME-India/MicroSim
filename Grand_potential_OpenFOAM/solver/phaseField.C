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

gsl_spline *spline1;
gsl_spline *spline2;
gsl_spline *spline3;
gsl_spline *spline4;
gsl_spline *spline5;
gsl_spline *spline6;
gsl_spline *spline7;
gsl_spline *spline8;
gsl_spline *spline9;
gsl_spline *spline10;
gsl_spline *spline_L;
gsl_spline *spline_C;

gsl_interp_accel *acc1;
gsl_interp_accel *acc2;
gsl_interp_accel *acc3;
gsl_interp_accel *acc4;
gsl_interp_accel *acc5;
gsl_interp_accel *acc6;
gsl_interp_accel *acc7;
gsl_interp_accel *acc8;
gsl_interp_accel *acc9;
gsl_interp_accel *acc10;
gsl_interp_accel *acc_L;
gsl_interp_accel *acc_C;

double calculate_melting_point(double cl, double T_guess, double T_high, double T_low) {
  //int iter=0;
  double tol=1e-6;
  dimensionedScalar T = T_guess;
  dimensionedScalar fval, dfdT;  
  
  do {
    dimensionedScalar A_Sol = gsl_spline_eval (spline1, T.value(), acc1);
    dimensionedScalar A_Liq = gsl_spline_eval (spline2, T.value(), acc2);
  
    dimensionedScalar A_SoldT = gsl_spline_eval_deriv (spline1, T.value(), acc1);
    dimensionedScalar A_LiqdT = gsl_spline_eval_deriv (spline2, T.value(), acc2);
  
    dimensionedScalar c_Sol = gsl_spline_eval (spline3, T.value(), acc3);
    dimensionedScalar c_Liq = gsl_spline_eval (spline4, T.value(), acc4);

    dimensionedScalar c_SoldT = gsl_spline_eval_deriv (spline3, T.value(), acc3);
    dimensionedScalar c_LiqdT = gsl_spline_eval_deriv (spline4, T.value(), acc4);
    
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
    dmudc_a[0][0] = gsl_spline_eval (spline1, T.value(), acc1);
    dmudc_a[0][1] = gsl_spline_eval (spline2, T.value(), acc2);
    dmudc_a[1][1] = gsl_spline_eval (spline3, T.value(), acc3);

    dmudc_adT[0][0] = gsl_spline_eval_deriv (spline1, T.value(), acc1);
    dmudc_adT[0][1] = gsl_spline_eval_deriv (spline2, T.value(), acc2);
    dmudc_adT[1][1] = gsl_spline_eval_deriv (spline3, T.value(), acc3);
    
    dmudc_l[0][0] = gsl_spline_eval (spline4, T.value(), acc4);
    dmudc_l[0][1] = gsl_spline_eval (spline5, T.value(), acc5);
    dmudc_l[1][1] = gsl_spline_eval (spline6, T.value(), acc6);

    dmudc_ldT[0][0] = gsl_spline_eval_deriv (spline4, T.value(), acc4);
    dmudc_ldT[0][1] = gsl_spline_eval_deriv (spline5, T.value(), acc5);
    dmudc_ldT[1][1] = gsl_spline_eval_deriv (spline6, T.value(), acc6);
    
    ceq_a[0][0] = gsl_spline_eval (spline7, T.value(), acc7);
    ceq_a[1][0] = gsl_spline_eval (spline8, T.value(), acc8);

    ceq_adT[0][0] = gsl_spline_eval_deriv (spline7, T.value(), acc7);
    ceq_adT[1][0] = gsl_spline_eval_deriv (spline8, T.value(), acc8);
    
    ceq_l[0][0] = gsl_spline_eval (spline9, T.value(), acc9);
    ceq_l[1][0] = gsl_spline_eval (spline10, T.value(), acc10);

    ceq_ldT[0][0] = gsl_spline_eval_deriv (spline9, T.value(), acc9);
    ceq_ldT[1][0] = gsl_spline_eval_deriv (spline10, T.value(), acc10);

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
   volVectorField q=dimx*fvc::grad(phi);
   volScalarField magq = 0.0*phi;
   //volScalarField T = 0.0*phi + initial;
   //volScalarField sumLHS = 0.0*phi;
   //volVectorField q_4 = 0.0*q;
   //volVectorField q_6 = 0.0*q;
   volScalarField ac_01 =0.0*phi;
   volVectorField dAdq01= phi*vector(0,0,0);
   volVectorField dadgradPhi=q*0.0;
   volScalarField grad_qt_sqr = 0.0*phi;
    //volScalarField A_Sol = 0.0*phi;
    //volScalarField A_Liq = 0.0*phi;
    //volScalarField c_Sol = 0.0*phi;
    //volScalarField c_Liq = 0.0*phi;
    //volScalarField dcdmu_Liq = 0.0*phi;
    //volScalarField omega = 0.0*phi;

    
    //! Done reading thermodynamic database
    
    //! Allocation for GSL acceleration and interpolation
    
    if (components == 2)
    {
    spline1 = gsl_spline_alloc (gsl_interp_cspline, np);
    acc1 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline1, T1, ASol, np);
    
    spline2 = gsl_spline_alloc (gsl_interp_cspline, np);
    acc2 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline2, T1, ALiq, np);
    
    spline3 = gsl_spline_alloc (gsl_interp_cspline, np);
    acc3 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline3, T1, cSol, np);
    
    spline4 = gsl_spline_alloc (gsl_interp_cspline, np);
    acc4 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline4, T1, cLiq, np);

    spline_L = gsl_spline_alloc (gsl_interp_cspline, np1);
    acc_L = gsl_interp_accel_alloc ();
    gsl_spline_init (spline_L, T2, Lf_arr, np1);

    spline_C = gsl_spline_alloc (gsl_interp_cspline, np1);
    acc_C = gsl_interp_accel_alloc ();
    gsl_spline_init (spline_C, T2, Cp_arr, np1);
    }
    else if (components == 3)
    {
    spline1 = gsl_spline_alloc (gsl_interp_cspline, np);
    acc1 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline1, T1, HSol11, np);

    spline2 = gsl_spline_alloc (gsl_interp_cspline, np);
    acc2 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline2, T1, HSol12, np);

    spline3 = gsl_spline_alloc (gsl_interp_cspline, np);
    acc3 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline3, T1, HSol22, np);
    
    spline4 = gsl_spline_alloc (gsl_interp_cspline, np);
    acc4 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline4, T1, HLiq11, np);
    
    spline5 = gsl_spline_alloc (gsl_interp_cspline, np);
    acc5 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline5, T1, HLiq12, np);
    
    spline6 = gsl_spline_alloc (gsl_interp_cspline, np);
    acc6 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline6, T1, HLiq22, np);
    
    spline7 = gsl_spline_alloc (gsl_interp_cspline, np);
    acc7 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline7, T1, cSol1, np);
    
    spline8 = gsl_spline_alloc (gsl_interp_cspline, np);
    acc8 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline8, T1, cSol2, np);
    
    spline9 = gsl_spline_alloc (gsl_interp_cspline, np);
    acc9 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline9, T1, cLiq1, np);

    spline10 = gsl_spline_alloc (gsl_interp_cspline, np);
    acc10 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline10, T1, cLiq2, np);
    
    spline_L = gsl_spline_alloc (gsl_interp_cspline, np1);
    acc_L = gsl_interp_accel_alloc ();
    gsl_spline_init (spline_L, T2, Lf_arr, np1);

    spline_C = gsl_spline_alloc (gsl_interp_cspline, np1);
    acc_C = gsl_interp_accel_alloc ();
    gsl_spline_init (spline_C, T2, Cp_arr, np1);
    }
    
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
   dimensionedScalar DeltaC = 0.0;

   dimensionedScalar A_Sol = 0.0;
   dimensionedScalar A_Liq = 0.0;
        dimensionedScalar A_SoldT = 0.0;
        dimensionedScalar A_LiqdT = 0.0;
        
        dimensionedScalar c_Sol = 0.0;
        dimensionedScalar c_Liq = 0.0;

        dimensionedScalar c_SoldT = 0.0;
        dimensionedScalar c_LiqdT = 0.0;
        
        
        dimensionedScalar dcdmu_Liq = 0.0;
            
   dimensionedScalar omega = 0.0;
    
        dimensionedScalar B_Sol = 0.0;
        dimensionedScalar dBdT = 0.0;
        
        dimensionedScalar B_Liq = 0.0; //("B_Liq", 0.*mu);

        dimensionedScalar D_Sol = 0.0;
        dimensionedScalar D_Liq = 0.0;
        
        dimensionedScalar avg_liq_conc = 0.0;
   
    
    if (components == 2)
    {
    dimensionedScalar c_a = gsl_spline_eval (spline3, Tm.value(), acc3);
    dimensionedScalar c_l = gsl_spline_eval (spline4, Tm.value(), acc4);
    
    //Info << T << "AS " << A_Sol << "AL " << A_Liq << "cS " << c_Sol << "cL " << c_Liq << endl;
    
    DeltaC = c_l-c_a;
    }
    
    //! For ternary
   scalarRectangularMatrix c_at(2,1,0);
   scalarRectangularMatrix c_lt(2,1,0);
   scalarRectangularMatrix DeltaCt(2,1,0);
   scalarRectangularMatrix avg_liq_conct(2,1,0);
   scalarRectangularMatrix testMatrix(2,2,0);
   scalarRectangularMatrix testMatrix_1(2,2,0);

   scalarRectangularMatrix dmudc_adT(2,2,0);
   scalarRectangularMatrix dmudc_ldT(2,2,0);
   scalarRectangularMatrix dcdmu_adT(2,2,0);
   scalarRectangularMatrix dcdmu_ldT(2,2,0);
   scalarRectangularMatrix ceq_adT(2,1,0);
   scalarRectangularMatrix ceq_ldT(2,1,0);
   
   dimensionedScalar B_a1 = 0.0;
   dimensionedScalar dB_a1dT = 0.0;
   dimensionedScalar B_a2 = 0.0;
   dimensionedScalar dB_a2dT = 0.0;
   dimensionedScalar DD_a = 0.0;
    
    if (components == 3)
    {
    c_at[0][0] = gsl_spline_eval (spline7, Tm.value(), acc7);
    c_at[1][0] = gsl_spline_eval (spline8, Tm.value(), acc8);
    
    c_lt[0][0] = gsl_spline_eval (spline9, Tm.value(), acc9);
    c_lt[1][0] = gsl_spline_eval (spline10, Tm.value(), acc10);
   
   DeltaCt = c_lt-c_at;
   }

    dimensionedScalar totalVol = fvc::domainIntegrate(1-0.0*phi);
   
    //! Total volume of solidus
   dimensionedScalar solidVolFracOld = fvc::domainIntegrate(phi)/totalVol;
   dimensionedScalar solidVolFrac = fvc::domainIntegrate(phi)/totalVol;
   
   dimensionedScalar dsolidVolFrac = (solidVolFrac - solidVolFracOld);
   dimensionedScalar dsolidVolFracOld = dsolidVolFrac;
    
    
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        //! Timestep increment for accelerating simulation.
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//! Checking monotonicity for stability of phi equation
	//! Change of average phi showing no oscillation
	if (sign(dsolidVolFrac.value()) == sign(dsolidVolFracOld.value()))
	//if (dsolidVolFrac.value() >= 0)
	{
	iter_num += 1;

	if (iter_num > 100 && max(max(InitialResidual_0,max(InitialResidual_10,InitialResidual_11)),max(max(max(InitialResidual_20,InitialResidual_21),max(InitialResidual_22,InitialResidual_23)),InitialResidual_3)) < Tol) //.value() && runTime.value() < 1000*InitialDeltaT
	{   
    	runTime.setDeltaT
        	(
            	runTime.deltaTValue() + dtf*InitialDeltaT
        	);
    	// set the iter_num = 0
    	iter_num = 0;
    	Info<< "deltaT increased to " <<  runTime.deltaTValue() << endl;
    	
    
	} else if (iter_num > 100 && max(max(InitialResidual_0,max(InitialResidual_10,InitialResidual_11)),max(max(max(InitialResidual_20,InitialResidual_21),max(InitialResidual_22,InitialResidual_23)),InitialResidual_3)) > Tol)
	{
    
	runTime.setDeltaT
        	(
            	runTime.deltaTValue() - dtf*InitialDeltaT 
        	);
	if (runTime.deltaTValue() < InitialDeltaT) //|| mag(Tdot) > )
	  {
          	runTime.setDeltaT(InitialDeltaT);  
          }
    	iter_num =0;
    	Info<< "deltaT decreased to " <<  runTime.deltaTValue() << endl;

	}

	Info<< "iter_num = " <<  iter_num << endl;
	}
	//! Change of average phi showing oscillation
	if (sign(dsolidVolFrac.value()) != sign(dsolidVolFracOld.value()))
	//if (dsolidVolFrac.value() < 0)
	{
	runTime.setDeltaT(InitialDeltaT);
    	Info<< "deltaT decreased to " <<  runTime.deltaTValue() << endl;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

        #include "readSolidDisplacementFoamControls.H"
        
            // The imposed temperature field as a function of thermal gradient
            // in the x direction, G and pulling velocity, v
            //T = G*( (1/dimx)* mesh.C().component(vector::X) - v*runTime.value()) + initial;
        
        //! Computing temperature evolution
    
        dimensionedScalar T_old = T;
        
    solidVolFrac = fvc::domainIntegrate(phi)/totalVol;
    
    dsolidVolFracOld = dsolidVolFrac;
    dsolidVolFrac = (solidVolFrac - solidVolFracOld);
    
        if (swcool == 1)
        {
            dimensionedScalar Lf = gsl_spline_eval (spline_L, T.value(), acc_L);
            dimensionedScalar Cp = gsl_spline_eval (spline_C, T.value(), acc_C);
            T = T + (Lf/Cp)*dsolidVolFrac - (qdot/Cp)*runTime.deltaTValue();
        
        }

        if (T.value() > Tm.value()) {
        T.value() = Tm.value();
    	}
    	else if (T.value() < T1[0]) {
        T.value() = T1[0];
        break;
    	}
    
        solidVolFracOld = solidVolFrac;

        outpft << "time " << runTime.timeName() << " " << T.value() << " " << solidVolFrac.value() << nl << endl;
        Info << "Temperature: "  << T.value() << " K" << endl;
        
        dimensionedScalar Tdot = (T-T_old)/runTime.deltaTValue();

        //! Computing the thermodynamic parameters in phase diagram
        if (components == 2)
    	{
        A_Sol = gsl_spline_eval (spline1, T.value(), acc1);
        A_Liq = gsl_spline_eval (spline2, T.value(), acc2);
        
        A_SoldT = gsl_spline_eval_deriv (spline1, T.value(), acc1);
        A_LiqdT = gsl_spline_eval_deriv (spline2, T.value(), acc2);
        
        c_Sol = gsl_spline_eval (spline3, T.value(), acc3);
        c_Liq = gsl_spline_eval (spline4, T.value(), acc4);

        c_SoldT = gsl_spline_eval_deriv (spline3, T.value(), acc3);
        c_LiqdT = gsl_spline_eval_deriv (spline4, T.value(), acc4);
        
        dcdmu_Liq = 0.5/A_Liq;
            
        omega = epsilon*0.18*DeltaC*DeltaC/(dcdmu_Liq*diff_Liq*Vm);
    
        B_Sol = 2.0*(A_Liq*c_Liq - A_Sol*c_Sol);
        dBdT = 2.0*(A_LiqdT*c_Liq + A_Liq*c_LiqdT - A_SoldT*c_Sol - A_Sol*c_SoldT);
        
        B_Liq = 0.0; //("B_Liq", 0.*mu);

        D_Sol = -A_Liq*c_Liq*c_Liq + A_Sol*c_Sol*c_Sol;
        D_Liq = 0.0; //("D_Liq", 0.*mu);
	}
        
        if (components == 3)
    	{
    dmudc_a[0][0] = gsl_spline_eval (spline1, T.value(), acc1);
    dmudc_a[0][1] = gsl_spline_eval (spline2, T.value(), acc2);
    dmudc_a[1][1] = gsl_spline_eval (spline3, T.value(), acc3);

    dmudc_adT[0][0] = gsl_spline_eval_deriv (spline1, T.value(), acc1);
    dmudc_adT[0][1] = gsl_spline_eval_deriv (spline2, T.value(), acc2);
    dmudc_adT[1][1] = gsl_spline_eval_deriv (spline3, T.value(), acc3);
    
    dmudc_l[0][0] = gsl_spline_eval (spline4, T.value(), acc4);
    dmudc_l[0][1] = gsl_spline_eval (spline5, T.value(), acc5);
    dmudc_l[1][1] = gsl_spline_eval (spline6, T.value(), acc6);

    dmudc_ldT[0][0] = gsl_spline_eval_deriv (spline4, T.value(), acc4);
    dmudc_ldT[0][1] = gsl_spline_eval_deriv (spline5, T.value(), acc5);
    dmudc_ldT[1][1] = gsl_spline_eval_deriv (spline6, T.value(), acc6);
    
    ceq_a[0][0] = gsl_spline_eval (spline7, T.value(), acc7);
    ceq_a[1][0] = gsl_spline_eval (spline8, T.value(), acc8);

    ceq_adT[0][0] = gsl_spline_eval_deriv (spline7, T.value(), acc7);
    ceq_adT[1][0] = gsl_spline_eval_deriv (spline8, T.value(), acc8);
    
    ceq_l[0][0] = gsl_spline_eval (spline9, T.value(), acc9);
    ceq_l[1][0] = gsl_spline_eval (spline10, T.value(), acc10);

    ceq_ldT[0][0] = gsl_spline_eval_deriv (spline9, T.value(), acc9);
    ceq_ldT[1][0] = gsl_spline_eval_deriv (spline10, T.value(), acc10);
    
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
    
   multiply(M_a, D_a, dcdmu_a);
   multiply(M_l, D_l, dcdmu_l);
   
      testMatrix = SVDinv(M_l);
   
   multiply(testMatrix_1,DeltaCt.T(), testMatrix);
   
   multiply(testMatrix, testMatrix_1, DeltaCt);
   
//    transportProperties.set("omega",epsilon*0.18*testMatrix[0][0]);
   
   omega = epsilon*0.18*testMatrix[0][0]/Vm;
   
   Info<< "omega= " << omega.value() << nl << endl;
      
   B_a1 = (dmudc_l[0][0]*ceq_l[0][0] - dmudc_a[0][0]*ceq_a[0][0]) + (dmudc_l[0][1]*ceq_l[1][0] - dmudc_a[0][1]*ceq_a[1][0]);
   
   dB_a1dT = (dmudc_ldT[0][0]*ceq_l[0][0] + dmudc_l[0][0]*ceq_ldT[0][0] - dmudc_adT[0][0]*ceq_a[0][0] - dmudc_a[0][0]*ceq_adT[0][0]) + (dmudc_ldT[0][1]*ceq_l[1][0] + dmudc_l[0][1]*ceq_ldT[1][0] - dmudc_adT[0][1]*ceq_a[1][0] - dmudc_a[0][1]*ceq_adT[1][0]);
   
   B_a2 = (dmudc_l[1][1]*ceq_l[1][0] - dmudc_a[1][1]*ceq_a[1][0]) + (dmudc_l[0][1]*ceq_l[0][0] - dmudc_a[0][1]*ceq_a[0][0]);
   
   dB_a2dT = (dmudc_ldT[1][1]*ceq_l[1][0] + dmudc_l[1][1]*ceq_ldT[1][0] - dmudc_adT[1][1]*ceq_a[1][0] - dmudc_a[1][1]*ceq_adT[1][0]) + (dmudc_ldT[0][1]*ceq_l[0][0] + dmudc_l[0][1]*ceq_ldT[0][0] - dmudc_adT[0][1]*ceq_a[0][0] - dmudc_a[0][1]*ceq_adT[0][0]);
   
   DD_a = -(0.5*dmudc_l[0][0]*ceq_l[0][0]*ceq_l[0][0] + 0.5*dmudc_l[1][1]*ceq_l[1][0]*ceq_l[1][0] + dmudc_l[0][1]*ceq_l[0][0]*ceq_l[1][0]) + (0.5*dmudc_a[0][0]*ceq_a[0][0]*ceq_a[0][0] + 0.5*dmudc_a[1][1]*ceq_a[1][0]*ceq_a[1][0] + dmudc_a[0][1]*ceq_a[0][0]*ceq_a[1][0]);
    	}

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
                #include "alphaEqn.H"
                #include "muEqn.H"
                #include "thetaEqn.H"



	}
    
	#include "DEqn.H"

    //#include "nucleateFields.H"
    if (components == 2)
    {
    avg_liq_conc = fvc::domainIntegrate((0.5*(mu_1-B_Liq)/A_Liq)*(1-phi)).value()/fvc::domainIntegrate((1-phi)).value();
    }
    
    if (components == 3)
    {
    avg_liq_conct[0][0] = fvc::domainIntegrate((dcdmu_l[0][0]*(mu_1) + dcdmu_l[0][1]*(mu_2))*(1-phi)).value()/fvc::domainIntegrate((1-phi)).value();

      avg_liq_conct[1][0] = fvc::domainIntegrate((dcdmu_l[1][0]*(mu_1) + dcdmu_l[1][1]*(mu_2))*(1-phi)).value()/fvc::domainIntegrate((1-phi)).value();
    }
    
    //Info << "avg_liq_conc: " << avg_liq_conc.value() << endl;
    
    //! Initial conditions for cooling simulation
if ((swcool == 1)&&(swch == 1))
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
          if ((int)(runTime.value()/runTime.deltaTValue())%(int)(nucleation_interval)==0) {
              //#include "nucleateSeeds.H"
              //theta_fill = 2.0*M_PI*randNumber.scalar01();
    	      //r          = (droplet_radius)*randNumber.scalar01();

              Info<< "creating a new seed" << endl;
            xCenter = randNumber.globalScalar01()*xMax;
            yCenter = randNumber.globalScalar01()*yMax;
            //xCenter = 0.5*xMax;
            //yCenter = 0.5*yMax;
            if (dimensions == 2)
            {
            Info<< "xCenter, yCenter: " << xCenter << ", " << yCenter << endl;
            randTheta[2] = randNumber.globalScalar01()*(pi.value()/2);
            Info<< "random theta: " << randTheta[2] << endl;

            Info<< "Filling phi and theta fields in seeds" << endl;
            volScalarField gaussianSeed = (1-phi)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius*seedRadius));

            theta = theta + randTheta[2]*gaussianSeed*vector(0,0,1);
            phi = phi + gaussianSeed;
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
            volScalarField gaussianSeed = (1-phi)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter) + (mesh.C().component(vector::Z)/dimx-zCenter)*(mesh.C().component(vector::Z)/dimx-zCenter))/(seedRadius*seedRadius));


            theta = theta + gaussianSeed*(randTheta[0]*vector(1,0,0) + randTheta[1]*vector(0,1,0) + randTheta[2]*vector(0,0,1));
            phi = phi + gaussianSeed;
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
    gsl_spline_free (spline1);
    gsl_spline_free (spline2);
    gsl_spline_free (spline3);
    gsl_spline_free (spline4);
    gsl_spline_free (spline5);
    gsl_spline_free (spline6);
    gsl_spline_free (spline7);
    gsl_spline_free (spline8);
    gsl_spline_free (spline9);
    gsl_spline_free (spline10);
    gsl_spline_free (spline_L);
    gsl_spline_free (spline_C);

    
    gsl_interp_accel_free (acc1);
    gsl_interp_accel_free (acc2);
    gsl_interp_accel_free (acc3);
    gsl_interp_accel_free (acc4);    
    gsl_interp_accel_free (acc5);
    gsl_interp_accel_free (acc6);
    gsl_interp_accel_free (acc7);
    gsl_interp_accel_free (acc8);    
    gsl_interp_accel_free (acc9);
    gsl_interp_accel_free (acc10);
    gsl_interp_accel_free (acc_L);
    gsl_interp_accel_free (acc_C);


    outpft.close();
    
    return 0;
}
