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
    laplacianFoam

Description
    Solves a simple Laplace equation, e.g. for thermal diffusion in a solid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include <fstream>
#include <sstream>
#include "scalarMatrices.H"
#include "LUscalarMatrix.H"
#include "RectangularMatrix.H"
#include "SquareMatrix.H"
//#include "Random.H"
//! Using GSL for interpolating thermodynamic data
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
// #include "dynamicFvMesh.H"

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
gsl_spline *spline11;
gsl_spline *spline12;

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
gsl_interp_accel *acc11;
gsl_interp_accel *acc12;

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

int main(int argc, char *argv[]) {
  
    //#include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControls.H"
    #include "createFields.H"
    #include "initializeFields.H"
    //#include "initializeFields3D.H"
    #include "createTol.H"

    int iter_num = 0;
    int nucleation_event =0;

//     #include "createDynamicFvMesh.H"

    simpleControl simple(mesh);


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating phase-field and chemical potential distributions\n" << endl;
  
    volVectorField q=dimx*fvc::grad(phi);
//    volScalarField magq = 0.0*phi;
//    volScalarField T = 0.0*phi + initial;
//    volScalarField sumLHS = 0.0*phi;
//    volVectorField q_4 = 0.0*q;
//    volVectorField q_6 = 0.0*q;
    volScalarField ac_01 =0.0*phi;
    volVectorField dAdq01= phi*vector(0,0,0);
    volVectorField dadgradPhi=q*0.0;
   

    // Done reading thermodynamic database
    
    //! Allocation for GSL acceleration and interpolation
    
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
    
    spline11 = gsl_spline_alloc (gsl_interp_cspline, np1);
    acc11 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline11, T2, Lf_arr, np1);

    spline12 = gsl_spline_alloc (gsl_interp_cspline, np1);
    acc12 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline12, T2, Cp_arr, np1);

    
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
    
   
   scalar InitialResidual_p = 0.0;
   scalar InitialResidual_t = 0.0;
   scalar InitialResidual_1 = 0.0;
   scalar InitialResidual_2 = 0.0;
   scalar InitialResidual_3 = 0.0;
   
   scalar InitialDeltaT = runTime.deltaTValue();
   
   int iCorr = 0;
   
   scalarRectangularMatrix c_a(2,1,0);
   scalarRectangularMatrix c_l(2,1,0);
   scalarRectangularMatrix DeltaC(2,1,0);
   scalarRectangularMatrix avg_liq_conc(2,1,0);
   scalarRectangularMatrix testMatrix(2,2,0);
   scalarRectangularMatrix testMatrix_1(2,2,0);

   scalarRectangularMatrix dmudc_adT(2,2,0);
   scalarRectangularMatrix dmudc_ldT(2,2,0);
   scalarRectangularMatrix dcdmu_adT(2,2,0);
   scalarRectangularMatrix dcdmu_ldT(2,2,0);
   scalarRectangularMatrix ceq_adT(2,1,0);
   scalarRectangularMatrix ceq_ldT(2,1,0);
   //    scalarRectangularMatrix suscept(2,2,0);
   
    c_a[0][0] = gsl_spline_eval (spline7, Tm.value(), acc7);
    c_a[1][0] = gsl_spline_eval (spline8, Tm.value(), acc8);
    
    c_l[0][0] = gsl_spline_eval (spline9, Tm.value(), acc9);
    c_l[1][0] = gsl_spline_eval (spline10, Tm.value(), acc10);
   
   DeltaC = c_l-c_a;
   
   dimensionedScalar totalVol = fvc::domainIntegrate(1-0.0*phi);
   
    //! Total volume of solidus
   dimensionedScalar solidVolFracOld = fvc::domainIntegrate(phi)/totalVol;
   dimensionedScalar solidVolFrac = fvc::domainIntegrate(phi)/totalVol;
   
   dimensionedScalar dsolidVolFrac = (solidVolFrac - solidVolFracOld);
   dimensionedScalar dsolidVolFracOld = dsolidVolFrac;

   
  while (runTime.loop()) {
    Info<< "Time = " << runTime.timeName() << nl << endl;
    
       
        //! Timestep increment for accelerating simulation.
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//! Checking monotonicity for stability of phi equation
	//! Change of average phi showing no oscillation
	if (sign(dsolidVolFrac.value()) == sign(dsolidVolFracOld.value()))
	//if (dsolidVolFrac.value() >= 0)
	{
	iter_num += 1;

	if (iter_num > 100 && max(max(InitialResidual_p,InitialResidual_t),max(max(InitialResidual_1,InitialResidual_2),InitialResidual_3)) < Tol) //.value() && runTime.value() < 1000*InitialDeltaT
	{   
    	runTime.setDeltaT
        	(
            	runTime.deltaTValue() + dtf*InitialDeltaT
        	);
    	// set the iter_num = 0
    	iter_num = 0;
    	Info<< "deltaT increased to " <<  runTime.deltaTValue() << endl;
    	
    
	} else if (iter_num > 100 && max(max(InitialResidual_p,InitialResidual_t),max(max(InitialResidual_1,InitialResidual_2),InitialResidual_3)) > Tol)
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

    //GSL accepts only type double
    //! Computing temperature evolution
    
    dimensionedScalar T_old = T;
    
    
    solidVolFrac = fvc::domainIntegrate(phi)/totalVol;
    
    dsolidVolFracOld = dsolidVolFrac;
    dsolidVolFrac = (solidVolFrac - solidVolFracOld);
    
    if (swcool == 1)
    {
        dimensionedScalar Lf = gsl_spline_eval (spline11, T.value(), acc11);
        dimensionedScalar Cp = gsl_spline_eval (spline12, T.value(), acc12);
        T = T + (Lf/Cp)*dsolidVolFrac - (qdot/Cp)*runTime.deltaTValue();
        
    }
    
    	if (T.value() > Tm.value()) {
        T.value() = Tm.value();
    	}
    	else if (T.value() < T1[0]) {
        T.value() = T1[0];
        break;
    	}
    
    //dimensionedScalar solidVolFrac = solidVol/totalVol;
    
    solidVolFracOld = solidVolFrac;
    
    //if (runTime.writeTime()) //runTime == 
    //{
    outpft << "time " << runTime.timeName() << " " << T.value() << " " << solidVolFrac.value() << nl << endl;
    //}
    Info << "Temperature: "  << T.value() << " K" << endl;
    //Info << T.value() << " " << Lf << " " << Cp << " " << solidVol << " " << totalVol << endl;
    
    dimensionedScalar Tdot = (T-T_old)/runTime.deltaTValue();
    
    //! Computing the thermodynamic parameters in phase diagram
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
   
   multiply(testMatrix_1,DeltaC.T(), testMatrix);
   
   multiply(testMatrix, testMatrix_1, DeltaC);
   
//    transportProperties.set("omega",epsilon*0.18*testMatrix[0][0]);
   
   dimensionedScalar omega = epsilon*0.18*testMatrix[0][0]/Vm;
   
   Info<< "omega= " << omega.value() << nl << endl;
      
   dimensionedScalar B_a1 = (dmudc_l[0][0]*ceq_l[0][0] - dmudc_a[0][0]*ceq_a[0][0]) + (dmudc_l[0][1]*ceq_l[1][0] - dmudc_a[0][1]*ceq_a[1][0]);
   
   dimensionedScalar dB_a1dT = (dmudc_ldT[0][0]*ceq_l[0][0] + dmudc_l[0][0]*ceq_ldT[0][0] - dmudc_adT[0][0]*ceq_a[0][0] - dmudc_a[0][0]*ceq_adT[0][0]) + (dmudc_ldT[0][1]*ceq_l[1][0] + dmudc_l[0][1]*ceq_ldT[1][0] - dmudc_adT[0][1]*ceq_a[1][0] - dmudc_a[0][1]*ceq_adT[1][0]);
   
   dimensionedScalar B_a2 = (dmudc_l[1][1]*ceq_l[1][0] - dmudc_a[1][1]*ceq_a[1][0]) + (dmudc_l[0][1]*ceq_l[0][0] - dmudc_a[0][1]*ceq_a[0][0]);
   
   dimensionedScalar dB_a2dT = (dmudc_ldT[1][1]*ceq_l[1][0] + dmudc_l[1][1]*ceq_ldT[1][0] - dmudc_adT[1][1]*ceq_a[1][0] - dmudc_a[1][1]*ceq_adT[1][0]) + (dmudc_ldT[0][1]*ceq_l[0][0] + dmudc_l[0][1]*ceq_ldT[0][0] - dmudc_adT[0][1]*ceq_a[0][0] - dmudc_a[0][1]*ceq_adT[0][0]);
   
   dimensionedScalar DD_a = -(0.5*dmudc_l[0][0]*ceq_l[0][0]*ceq_l[0][0] + 0.5*dmudc_l[1][1]*ceq_l[1][0]*ceq_l[1][0] + dmudc_l[0][1]*ceq_l[0][0]*ceq_l[1][0]) + (0.5*dmudc_a[0][0]*ceq_a[0][0]*ceq_a[0][0] + 0.5*dmudc_a[1][1]*ceq_a[1][0]*ceq_a[1][0] + dmudc_a[0][1]*ceq_a[0][0]*ceq_a[1][0]);
   
    
   while (simple.correctNonOrthogonal()){
    //! Solving the phase-field and chemical potential equations
    #include "alphaEqn.H"
    
    volVectorField n=dimx*fvc::grad(phi)/(1E-20+mag(dimx*fvc::grad(phi)));

    volScalarField hphi = phi*phi*phi*(6*phi*phi - 15*phi + 10);
    
    /*
     c_a1 = (dcdmu_1[0][0]*(mu_1 - B_a1) + dcdmu_1[0][1]*(mu_2 - B_a2));
     c_l1 = (dcdmu_2[0][0]*mu_1 + dcdmu_2[0][1]*mu_2);
     c_a2 = (dcdmu_1[1][0]*(mu_1 - B_a1) + dcdmu_1[1][1]*(mu_2 - B_a2));
     c_l2 = (dcdmu_2[1][0]*mu_1 + dcdmu_2[1][1]*mu_2);
     */
    
    do {
      
      fvScalarMatrix muEqn_1 (
        (dcdmu_a[0][0]*hphi + dcdmu_l[0][0]*(1-hphi))*dimt*fvm::ddt(mu_1) +  (dcdmu_a[0][1]*hphi + dcdmu_l[0][1]*(1-hphi))*dimt*fvc::ddt(mu_2) 
        == dimx*dimx*fvm::laplacian(M_a[0][0]*phi + M_l[0][0]*(1-phi), mu_1) + dimx*dimx*fvc::laplacian(M_a[0][1]*phi + M_l[0][1]*(1-phi), mu_2)
        -((dcdmu_a[0][0]*(mu_1 - B_a1) + dcdmu_a[0][1]*(mu_2 - B_a2)) - (dcdmu_l[0][0]*mu_1 + dcdmu_l[0][1]*mu_2))* dimt*fvc::ddt(phi)*30.0*phi*phi*(1.0-phi)*(1.0-phi) - Tdot*((dcdmu_a[0][0]*(-dB_a1dT) + dcdmu_adT[0][0]*(mu_1 - B_a1) + dcdmu_a[0][1]*(-dB_a2dT) + dcdmu_adT[0][1]*(mu_2 - B_a2))*hphi + (dcdmu_ldT[0][0]*mu_1 + dcdmu_ldT[0][1]*mu_2)*(1-hphi)) - anti_trap*epsilon*((dcdmu_a[0][0]*(mu_1 - B_a1) + dcdmu_a[0][1]*(mu_2 - B_a2)) - (dcdmu_l[0][0]*mu_1 + dcdmu_l[0][1]*mu_2))*dimx*fvc::div((n*dimt*fvc::ddt(phi)))
      );
      
      InitialResidual_1 = muEqn_1.solve().max().initialResidual();
      
      fvScalarMatrix muEqn_2 (
        (dcdmu_a[1][0]*hphi + dcdmu_l[1][0]*(1-hphi))*dimt*fvc::ddt(mu_1) +  (dcdmu_a[1][1]*hphi + dcdmu_l[1][1]*(1-hphi))*dimt*fvm::ddt(mu_2) == 
        dimx*dimx*fvc::laplacian(M_a[1][0]*phi + M_l[1][0]*(1-phi), mu_1) + dimx*dimx*fvm::laplacian(M_a[1][1]*phi + M_l[1][1]*(1-phi), mu_2)
        -((dcdmu_a[1][0]*(mu_1 - B_a1) + dcdmu_a[1][1]*(mu_2 - B_a2)) - (dcdmu_l[1][0]*mu_1 + dcdmu_l[1][1]*mu_2))*dimt*fvc::ddt(phi)*30.0*phi*phi*(1.0-phi)*(1.0-phi) - Tdot*((dcdmu_a[1][0]*(-dB_a1dT) + dcdmu_adT[1][0]*(mu_1 - B_a1) + dcdmu_a[1][1]*(-dB_a2dT) + dcdmu_adT[1][1]*(mu_2 - B_a2))*hphi + (dcdmu_ldT[1][0]*mu_1 + dcdmu_ldT[1][1]*mu_2)*(1-hphi)) - anti_trap*epsilon*((dcdmu_a[1][0]*(mu_1 - B_a1) + dcdmu_a[1][1]*(mu_2 - B_a2)) - (dcdmu_l[1][0]*mu_1 + dcdmu_l[1][1]*mu_2))*dimx*fvc::div((n*dimt*fvc::ddt(phi)))
      );
      
      InitialResidual_2 = muEqn_2.solve().max().initialResidual();
      
    } while(InitialResidual_1 > 1e-5 || InitialResidual_2 > 1e-5);
    
    Info<< "Min/max mu_1:" << min(mu_1).value() << ' ' << max(mu_1).value() << endl;
    Info<< "Min/max mu_2:" << min(mu_2).value() << ' ' << max(mu_2).value() << endl;
   }
    
    //! Elastic stress, strain, displacement fields only for precipitate growth
    iCorr=0;
    if (swch == 2)
    {
        //nCorr = 0;
    
    
    do {
      fvVectorMatrix DEqn (
          dimt*dimt*fvm::d2dt2(D)
        ==
            dimx*dimx*fvm::laplacian(2*(mu1_elast*phi*phi*(3-2*phi) + mu2_elast*(1-phi)*(1-phi)*(1+2*phi)) 
          + lambda1*phi*phi*(3-2*phi) + lambda2*(1-phi)*(1-phi)*(1+2*phi), D, "laplacian(DD,D)")
          + dimx*dimx*divSigmaExp
          - dimx*dimx*fvc::div((2*mu1_elast*phi*phi*(3-2*phi) + 2*mu2_elast*(1-phi)*(1-phi)*(1+2*phi))*phi*phi*(3-2*phi)*cEigenStrain
          + (lambda1*phi*phi*(3-2*phi) + (1-phi)*(1-phi)*(1+2*phi)*lambda2)*I*tr(phi*phi*(3-2*phi)*cEigenStrain))
      );

      InitialResidual_3 = DEqn.solve().max().initialResidual();
    
      gradD  = fvc::grad(D);
      
      strain =((gradD-phi*phi*(3-2*phi)*cEigenStrain)&&symmTensor(1,0,0,0,0,0))*symmTensor(1,0,0,0,0,0)
             +((gradD-phi*phi*(3-2*phi)*cEigenStrain)&&symmTensor(0,0,0,1,0,0))*symmTensor(0,0,0,1,0,0)
             +((gradD-phi*phi*(3-2*phi)*cEigenStrain)&&symmTensor(0,0,0,0,0,1))*symmTensor(0,0,0,0,0,1);

      sigmaD = (mu1_elast*phi*phi*(3-2*phi)  + mu2_elast*(1-phi)*(1-phi)*(1+2*phi))*twoSymm(gradD) 
             + (lambda1*phi*phi*(3-2*phi)    + lambda2*(1-phi)*(1-phi)*(1+2*phi))*(I*tr(gradD))
             + (mu1_elast_*phi*phi*(3-2*phi) + mu2_elast_*(1-phi)*(1-phi)*(1+2*phi))*strain;
      
             
      divSigmaExp = fvc::div
                    (
                        sigmaD - (2*mu1_elast*phi*phi*(3-2*phi) + 2*mu2_elast*(1-phi)*(1-phi)*(1+2*phi)
			       + lambda1*phi*phi*(3-2*phi)      + (1-phi)*(1-phi)*(1+2*phi)*lambda2)*gradD,
                        "div(sigmaD)"
                    );
      
    }while(InitialResidual_3 > convergenceTolerance && ++iCorr < nCorr);
    
    Sigma = (2*(mu1_elast*phi*phi*(3-2*phi) + mu2_elast*(1-phi)*(1-phi)*(1+2*phi))*(symm(fvc::grad(D)) - phi*phi*(3-2*phi)*cEigenStrain) 
      + (lambda1*phi*phi*(3-2*phi)      + lambda2*(1-phi)*(1-phi)*(1+2*phi))*(I*tr(fvc::grad(D) - phi*phi*(3-2*phi)*cEigenStrain)))
      + (mu1_elast_*phi*phi*(3-2*phi)   + mu2_elast_*(1-phi)*(1-phi)*(1+2*phi))*strain;


    deltaSigmaD = ((mu1_elast-mu2_elast)*twoSymm(fvc::grad(D))           + (lambda1-lambda2)*(I*tr(fvc::grad(D))) 
            - 2*(mu1_elast-mu2_elast)*phi*phi*(3-2*phi)*cEigenStrain - (lambda1-lambda2)*(I*tr(phi*phi*(3-2*phi)*cEigenStrain)))
            + (mu1_elast_-mu2_elast_)*strain;
    
    deltaF = (0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi*phi*(3-2*phi)*cEigenStrain))-(Sigma && cEigenStrain));
    
    sigmaEq = sqrt((3.0/2.0)*magSqr(dev(Sigma)));
    }
    //#include "nucleateFields.H"
    
    avg_liq_conc[0][0] = fvc::domainIntegrate((dcdmu_l[0][0]*(mu_1) + dcdmu_l[0][1]*(mu_2))*(1-phi)).value()/fvc::domainIntegrate((1-phi)).value();

      avg_liq_conc[1][0] = fvc::domainIntegrate((dcdmu_l[1][0]*(mu_1) + dcdmu_l[1][1]*(mu_2))*(1-phi)).value()/fvc::domainIntegrate((1-phi)).value();
      
           
    //! Initial conditions for cooling simulation
if ((swcool == 1)&&(swch == 1))
{
      Tm.value() = calculate_melting_point(avg_liq_conc, Tm.value(), T1[np-1], T1[0]);
      
       Info << "Melting_point: "  << Tm.value() << " " << avg_liq_conc[0][0] << " " << avg_liq_conc[1][0] << endl;

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
            Info<< "xCenter, yCenter: " << xCenter << ", " << yCenter << endl;
            randTheta = randNumber.globalScalar01()*(pi.value()/2);
            Info<< "random theta: " << randTheta << endl;

            Info<< "Filling phi and theta fields in seeds" << endl;
            volScalarField gaussianSeed = (1-phi)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(seedRadius*seedRadius));

            //if (prob.value() > randNumber.globalScalar01()){
            theta = theta + randTheta*gaussianSeed;
            phi = phi + gaussianSeed;

              nucleation_event++;
          }
        }
      }
}
    
    //! Writing the results according to keywords in controlDict
    //if (runTime.writeTime()) {
      
      Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
      runTime.write();
    //}
        
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
    gsl_spline_free (spline11);
    gsl_spline_free (spline12);

    
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
    gsl_interp_accel_free (acc11);
    gsl_interp_accel_free (acc12);


    outpft.close();

  return 0;
}
