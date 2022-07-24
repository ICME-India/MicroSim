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
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    //#include "createDynamicFvMesh.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    #include "initializeFields.H"

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
    //volScalarField A_Sol = 0.0*phi;
    //volScalarField A_Liq = 0.0*phi;
    //volScalarField c_Sol = 0.0*phi;
    //volScalarField c_Liq = 0.0*phi;
    //volScalarField dcdmu_Liq = 0.0*phi;
    //volScalarField omega = 0.0*phi;

    
    //! Done reading thermodynamic database
    
    //! Allocation for GSL acceleration and interpolation
    
    gsl_spline *spline1 = gsl_spline_alloc (gsl_interp_cspline, np);
    gsl_interp_accel *acc1 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline1, T1, ASol, np);
    
    gsl_spline *spline2 = gsl_spline_alloc (gsl_interp_cspline, np);
    gsl_interp_accel *acc2 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline2, T1, ALiq, np);
    
    gsl_spline *spline3 = gsl_spline_alloc (gsl_interp_cspline, np);
    gsl_interp_accel *acc3 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline3, T1, cSol, np);
    
    gsl_spline *spline4 = gsl_spline_alloc (gsl_interp_cspline, np);
    gsl_interp_accel *acc4 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline4, T1, cLiq, np);

    gsl_spline *spline5 = gsl_spline_alloc (gsl_interp_cspline, np1);
    gsl_interp_accel *acc5 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline5, T2, Lf_arr, np1);

    gsl_spline *spline6 = gsl_spline_alloc (gsl_interp_cspline, np1);
    gsl_interp_accel *acc6 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline6, T2, Cp_arr, np1);
    
    
    dimensionedScalar T = initial;
    
    //! File to write temperature evolution
    const fileName pathToTemp = "tempvstime.dat";
    ofstream outpft(pathToTemp);
    if (!outpft.is_open())
    {
        Info << "No output file found" << endl;
    }
    

    dimensionedScalar c_a = gsl_spline_eval (spline3, T0.value(), acc3);
    dimensionedScalar c_l = gsl_spline_eval (spline4, T0.value(), acc4);
    
    //Info << T << "AS " << A_Sol << "AL " << A_Liq << "cS " << c_Sol << "cL " << c_Liq << endl;
    
    dimensionedScalar DeltaC = c_l-c_a;

    //! Total volume of solidus
   dimensionedScalar solidVolOld = fvc::domainIntegrate(phi).value();

   dimensionedScalar totalVol = fvc::domainIntegrate(1-0.0*phi).value();
    
    
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

            //! The imposed temperature field as a function of thermal gradient
            //! in the x direction, G and pulling velocity, v
            //T = G*( (1/dimx)* mesh.C().component(vector::X) - v*runTime.value()) + initial;
        
        
        dimensionedScalar Lf = gsl_spline_eval (spline5, T.value(), acc5);
        dimensionedScalar Cp = gsl_spline_eval (spline6, T.value(), acc6);
    
        dimensionedScalar solidVol = fvc::domainIntegrate(phi).value();
        dimensionedScalar T_old = T;
    
        if (swcool == 1)
        {
            T = T + (Lf/Cp)*(solidVol - solidVolOld)/totalVol - (qdot/Cp)*runTime.deltaTValue();
        
        }
    
        dimensionedScalar Tdot = (T-T_old)/runTime.deltaTValue();
    
        dimensionedScalar solidVolFrac = solidVol/totalVol;
    
        solidVolOld = solidVol;

        outpft << "time " << runTime.timeName() << " " << T.value() << " " << solidVolFrac.value() << nl << endl;
        Info << "Temperature: "  << T.value() << " K" << endl;

        
        dimensionedScalar A_Sol = gsl_spline_eval (spline1, T.value(), acc1);
        dimensionedScalar A_Liq = gsl_spline_eval (spline2, T.value(), acc2);
        
        dimensionedScalar A_SoldT = gsl_spline_eval_deriv (spline1, T.value(), acc1);
        dimensionedScalar A_LiqdT = gsl_spline_eval_deriv (spline2, T.value(), acc2);
        
        dimensionedScalar c_Sol = gsl_spline_eval (spline3, T.value(), acc3);
        dimensionedScalar c_Liq = gsl_spline_eval (spline4, T.value(), acc4);

        dimensionedScalar c_SoldT = gsl_spline_eval_deriv (spline3, T.value(), acc3);
        dimensionedScalar c_LiqdT = gsl_spline_eval_deriv (spline4, T.value(), acc4);
        
        
        dimensionedScalar dcdmu_Liq = 0.5/A_Liq;
            
        dimensionedScalar omega = epsilon*0.18*DeltaC*DeltaC/(dcdmu_Liq*diff_Liq*Vm);
    
        dimensionedScalar B_Sol = 2.0*(A_Liq*c_Liq - A_Sol*c_Sol);
        dimensionedScalar dBdT = 2.0*(A_LiqdT*c_Liq + A_Liq*c_LiqdT - A_SoldT*c_Sol - A_Sol*c_SoldT);
        
        dimensionedScalar B_Liq = 0.0; //("B_Liq", 0.*mu);

        dimensionedScalar D_Sol = -A_Liq*c_Liq*c_Liq + A_Sol*c_Sol*c_Sol;
        dimensionedScalar D_Liq = 0.0; //("D_Liq", 0.*mu);

          // For orthogonal correction of the finite volume mesh
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
    
    gsl_spline_free (spline1);
    gsl_spline_free (spline2);
    gsl_spline_free (spline3);
    gsl_spline_free (spline4);
    gsl_spline_free (spline5);
    gsl_spline_free (spline6);
    
    gsl_interp_accel_free (acc1);
    gsl_interp_accel_free (acc2);
    gsl_interp_accel_free (acc3);
    gsl_interp_accel_free (acc4);
    gsl_interp_accel_free (acc5);
    gsl_interp_accel_free (acc6);
    
    return 0;
}
