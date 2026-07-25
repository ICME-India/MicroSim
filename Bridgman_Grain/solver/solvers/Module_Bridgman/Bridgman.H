/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

Class
    Foam::solvers::solidDisplacement

Description
    Solver module for steady or transient segregated finite-volume solution of
    linear-elastic, small-strain deformation of a solid body, with optional
    thermal diffusion and thermal stresses.

    Solves for the displacement vector field D, also generating the stress
    tensor field sigma, including the thermal stress contribution if selected.

SourceFiles
    solidDisplacement.C

\*---------------------------------------------------------------------------*/

#ifndef Bridgman_H
#define Bridgman_H

#include "tmp.H"
#include "solid.H"
#include "OFstream.H"
#include "writeFile.H"
#include "OSspecific.H"
#include "bridgmanThermo.H"
#include "interfaceTopology.H"
#include "regionDataSharer.H"
#include "fvMeshMover.H"
#include "subCycle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace solvers
{

/*---------------------------------------------------------------------------*\
                              Class Bridgman Declaration
\*---------------------------------------------------------------------------*/

class Bridgman100
:
    public solid
{

protected:

    // Simulation controls

        //- Switch to simplify temperature equation, TEqn
        Switch simpleTEqn;

        //- Switch for domain parameter equation, psiEqn
        //  Used to smooth the ROI's boundary (crucible)
        Switch BPsi_diffusion;

        //- Flag to pick TEqn form
        //  value 0 is (least simplified), 1, 2, 3 (most simplified))
        int flag_1;

        //- No of subcycles of psiEqn() and STEqn

        label nPsiSubCycles,
              nStateSubCycles;

    //- Class containing bridgman thermophysical constants
    bridgmanThermo& thermo_;

    //- Settings manager for interface topology calculations
    interfaceTopologyUtils::interfaceTopologyManager topoManager_;

    //- Data sharer with neighbour region
    interfaceTopologyUtils::regionDataSharer dataSharer_;


    //- Current mesh point positions (updated each move step)
    pointField newPoints_;

    //- Simulation time at which the next mesh move is due
    scalar meshMoveTime_;

    //- Time interval between successive mesh move steps
    scalar meshMoveInterval_;

    // Constants

        //- Dimensions of space and time
        dimensioned<scalar> dimx, dimt;

        //- Smoothing parameters (used by psiEqn)
        dimensioned<scalar> omega, gamma, epsilon;

        //- Thermal properties of CMSX-4
        dimensioned<scalar> Ks, Kl, kc, Tm;

        //- Boundary condition parameters
        // Heat Transfer Coefficient and surrounding temperature
        dimensioned<scalar> flux, Tinf;

        //- Hot-zone temperature, pull velocity
        dimensioned<scalar> Bd,  Vpull, subCycleTime1;

        // Y coordinate of crucible bottom and associated references
            //- crucible bottom and masking region bounds
            scalar ybottom;

            //- Limit for solid-liquid interface tresholding
            scalar limT;

        //- Derived constants
        dimensioned<scalar> idimx, DCpp, Tls, pp;

    // Phase fields

        //- Domain parameter
        volScalarField BPsi_;

        //- Phase parameter 1
        volScalarField BPhi_;

        //- Masking field for dynamic meshing
        volScalarField BBPhi_;

        //- Grain temperature
        volScalarField grainT_;

        //- Field to display processor ownership of cells
        volScalarField procID_;

    // Cached temporary fields

        //- Gradient of BPsi
        autoPtr<volVectorField> GBPsi;

        //- Magnitude of BPsi gradient
        autoPtr<volScalarField> M_GBPsi;

        //- Multiplication of ((Bd*BPsi1)*GBPsi)
        autoPtr<volScalarField> BdBPGB;

    // Protected Member Functions

        //- Inline function used in temperature equation (TEqn)
        inline tmp<volScalarField> latent()
        {
          //Info << "TemperatureLatent :"
              // << max(grainT) << "/" <<  min(grainT) << endl;

          return DCpp - pp * exp(-pow((grainT - Tm) * Tls, 2));
        }

        //- Function containing equation to smoothen domain parameter
        void psiEqn();

        //- Function containing simplified forms of temperature equation (TEqn)
        void STEqn();

        //- Function making fvMatrix for domain matrix solving
        void solvePsiEqn();

        //- Return true if the solver's dependencies have been modified
        virtual bool dependenciesModified() const;

        //- Read controls
        virtual bool read();

public:

    // Public Data

        //- Reference to the Domain parameter
        volScalarField& BPsi;

        //- Reference to the Phase parameter 1
        volScalarField& BPhi;

        //- Reference to the Phase parameter 3
        volScalarField& BBPhi;

        //- Reference to modified MT
        volScalarField& grainT;

        //- Current simulation time
        scalar currentTime;

    //- Runtime type information
    TypeName("Bridgman100");

    // Constructors

        //- Construct from region mesh
        Bridgman100(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        Bridgman100(const Bridgman100&) = delete;


    //- Destructor
    virtual ~Bridgman100();

    // Member Functions

        //- Called at the beginning of the PIMPLE loop
        virtual void preSolve();

        //- Translate the mesh by one pull increment if the move
        //  deadline has elapsed; rebuild geometry cache and remap
        //  neighbour fields onto the updated domain
        virtual void moveMesh();

        //- Called at the beginning of the PIMPLE loop
        virtual void prePredictor();

        //- Update thermophysical and transport properties
        virtual void thermophysicalPredictor();

        //- Construct and solve the thermal equation to obtain the temperature
        virtual void pressureCorrector();

        //- Correct the thermophysical transport modelling
        virtual void postCorrector();

        //- Called after the PIMPLE loop at the end of the time-step
        virtual void postSolve();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Bridgman100&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
