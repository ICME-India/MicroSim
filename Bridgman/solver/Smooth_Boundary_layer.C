/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
#include "mathematicalConstants.H"
#include "simpleControl.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

long num;
#define dt_0 0.1
#define dt_s 0.2

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
//     #include "initContinuityErrs.H"
//     #include "createControls.H"
//     #include "createTimeControls.H"
//     #include "CourantNo.H"
//     #include "createMoreFields.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating Phase-field profile\n" << endl;
    
//     scalar InitialResidual = 0.0; 
    
    double R               = seedRadius;
//     double r;
    scalar randTheta;
    scalar xCenter;
    scalar yCenter;
    
    volVectorField q          = dimx*fvc::grad(phi);
    volScalarField ac_01      = 0.0*phi;
    volVectorField dAdq01     = phi*vector(0,0,0);
    volVectorField dadgradPhi = q*0.0;
    int nucleation_event      = 0;
//volScalarField K = Ks*phi + Kl*(1-phi);
//volScalarField omega = epsilon*L*L*0.18/(Tm*K);
    
//     #include "initializeFields.H"
//     #include "nucleateSeeds.H"
    
//     double xMax = max(mesh.C().component(vector::X)).value();
//     double xMin = min(mesh.C().component(vector::X)).value();

//     double yMax = max(mesh.C().component(vector::Y)).value();
//     double yMin = min(mesh.C().component(vector::Y)).value();
    
    
    while (runTime.loop()) {
      
      label curTimeIndex = mesh.time().timeIndex();
      if (curTimeIndex <= 10) {
        #include "psiEqn_eq.H"
      } else {
      
        Info<< "Iteration: " << runTime.value() << nl << endl;
        Info<< "Time = " << runTime.timeName() << nl << endl;
        #include "TEqn.H"
        //#include "phiEqn.H"
//         #include "thetaEqn.H"
        //For nucleation
//         volScalarField deltaT    = mag(Tm-T+eps);
//         volScalarField Min_T_pos = pos(Tn - T);
//         volScalarField nprob = (0.5/((dt_s)*Foam::pow((2*M_PI),0.5)))*(1/(deltaT))*Foam::exp(-0.5*Foam::pow(((Foam::log(deltaT)-Foam::log(dt_0))/dt_s),2));
//         if (curTimeIndex%(int)(nucleation_interval)==0) {
//           
//           xCenter   = xbottom + Width*randNum.globalScalar01();
//           yCenter   = ybottom + Height*randNum.globalScalar01();
//           
//           randTheta = randNum.globalScalar01()*pi.value()/3;
//           
//           
//           volScalarField gaussianSeed = nprob*(1-phi)*(-theta+randTheta)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(R*R));
// 
//           volScalarField phiSeed = nprob*(1-phi)*exp(-((mesh.C().component(vector::X)/dimx-xCenter)*(mesh.C().component(vector::X)/dimx-xCenter) + (mesh.C().component(vector::Y)/dimx-yCenter)*(mesh.C().component(vector::Y)/dimx-yCenter))/(R*R));
// 
//           theta           = theta +  Min_T_pos*gaussianSeed;
//           phi             = phi   +  Min_T_pos*phiSeed;
//           
//         }
        //End nucleation
      }
      #include "write.H"
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
