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
//#include "fvOptions.H"
#include "simpleControl.H"
#include "Random.H"
//! Using dynamicFvMesh for dynamic interface refinement
#include "dynamicFvMesh.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating phase-field and chemical potential distributions\n" << endl;

   // The scalar and vector fields initialized below
   volVectorField q=dimx*fvc::grad(phi);
   volScalarField magq = 0.0*phi;
   volScalarField T = 0.0*phi + initial;
   volScalarField sumLHS = 0.0*phi;
   volVectorField q_4 = 0.0*q;
   volVectorField q_6 = 0.0*q;
   volScalarField ac_01 =0.0*phi;
   volVectorField dAdq01= phi*vector(0,0,0);
   volVectorField dadgradPhi=q*0.0;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

            //! The imposed temperature field as a function of thermal gradient
            //! in the x direction, G and pulling velocity, v
            T = G*( (1/dimx)* mesh.C().component(vector::X) - v*runTime.value()) + initial;

          // For orthogonal correction of the finite volume mesh
          while (simple.correctNonOrthogonal())
        {


          label curTimeIndex = mesh.time().timeIndex();

	for (int i=0; i<30; i++)
        {
		if(curTimeIndex == 0)
				{
                 //! The interpolation equation is used for smoothening of the
                 //! phase-field variable
                                    #include "preAlan.H"
                                }

	}

                mesh.update();
                //! Solving the phase-field and chemical potential equations after
                //! updating the mesh
                #include "alphaEqn.H"



	}


    //! Writing the results according to keywords in controlDict
    if (runTime.writeTime())
    {
      runTime.write();
    }
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}
