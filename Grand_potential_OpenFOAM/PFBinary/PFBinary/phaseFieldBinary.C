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

gsl_interp_accel *acc1;
gsl_interp_accel *acc2;
gsl_interp_accel *acc3;
gsl_interp_accel *acc4;
gsl_interp_accel *acc5;
gsl_interp_accel *acc6;

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
    //volScalarField A_Sol = 0.0*phi;
    //volScalarField A_Liq = 0.0*phi;
    //volScalarField c_Sol = 0.0*phi;
    //volScalarField c_Liq = 0.0*phi;
    //volScalarField dcdmu_Liq = 0.0*phi;
    //volScalarField omega = 0.0*phi;

    
    //! Done reading thermodynamic database
    
    //! Allocation for GSL acceleration and interpolation
    
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

    spline5 = gsl_spline_alloc (gsl_interp_cspline, np1);
    acc5 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline5, T2, Lf_arr, np1);

    spline6 = gsl_spline_alloc (gsl_interp_cspline, np1);
    acc6 = gsl_interp_accel_alloc ();
    gsl_spline_init (spline6, T2, Cp_arr, np1);
    
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
   scalar InitialResidual_1 = 0.0;
   scalar InitialResidual_2 = 0.0;
   scalar InitialResidual_3 = 0.0;
   
   scalar InitialDeltaT = runTime.deltaTValue();
   
   int iCorr = 0;
    
    dimensionedScalar c_a = gsl_spline_eval (spline3, Tm.value(), acc3);
    dimensionedScalar c_l = gsl_spline_eval (spline4, Tm.value(), acc4);
    
    //Info << T << "AS " << A_Sol << "AL " << A_Liq << "cS " << c_Sol << "cL " << c_Liq << endl;
    
    dimensionedScalar DeltaC = c_l-c_a;

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

	if (iter_num > 100 && max(max(InitialResidual_0,InitialResidual_1),max(InitialResidual_2,InitialResidual_3)) < Tol) //.value() && runTime.value() < 1000*InitialDeltaT
	{   
    	runTime.setDeltaT
        	(
            	runTime.deltaTValue() + dtf*InitialDeltaT
        	);
    	// set the iter_num = 0
    	iter_num = 0;
    	Info<< "deltaT increased to " <<  runTime.deltaTValue() << endl;
    	
    
	} else if (iter_num > 100 && max(max(InitialResidual_0,InitialResidual_1),max(InitialResidual_2,InitialResidual_3)) > Tol)
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
            dimensionedScalar Lf = gsl_spline_eval (spline5, T.value(), acc5);
            dimensionedScalar Cp = gsl_spline_eval (spline6, T.value(), acc6);
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
    
    dimensionedScalar avg_liq_conc = fvc::domainIntegrate((0.5*(mu-B_Liq)/A_Liq)*(1-phi)).value()/fvc::domainIntegrate((1-phi)).value();
    
    //Info << "avg_liq_conc: " << avg_liq_conc.value() << endl;
    
    //! Initial conditions for cooling simulation
if ((swcool == 1)&&(swch == 1))
{
    Tm.value() = calculate_melting_point(avg_liq_conc.value(), Tm.value(), T1[np-1], T1[0]);
      
       Info << "Melting_point: "  << Tm.value() << " " << avg_liq_conc.value() << endl;
       
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
    
    gsl_interp_accel_free (acc1);
    gsl_interp_accel_free (acc2);
    gsl_interp_accel_free (acc3);
    gsl_interp_accel_free (acc4);
    gsl_interp_accel_free (acc5);
    gsl_interp_accel_free (acc6);
    
    return 0;
}
