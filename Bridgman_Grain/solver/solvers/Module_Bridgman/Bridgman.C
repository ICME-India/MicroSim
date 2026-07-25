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

\*---------------------------------------------------------------------------*/

#include "Bridgman.H"
#include "IOstream.H"
#include "fvcGrad.H"
#include "fvmDdt.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    defineTypeNameAndDebug(Bridgman100, 0);
    addToRunTimeSelectionTable(solver, Bridgman100, fvMesh);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::solvers::Bridgman100::dependenciesModified() const
{
    return solid::dependenciesModified() || mesh.solution().modified();
}


bool Foam::solvers::Bridgman100::read()
{
    solid::read();
     
    //- Solution control keys assignment
    simpleTEqn = pimple.dict().lookupOrDefault<Switch>("simpleTEqn", false);
    BPsi_diffusion = pimple.dict().lookupOrDefault<Switch>("BPsi_diffusion", true); 
    flag_1  = pimple.dict().lookupOrDefault<int>("flag_1", 0);


    nStateSubCycles =
        max
        (
            label(1),
            pimple.dict().lookupOrDefault<label>
            (
                "nTSubCycles",
                label(runTime.deltaTValue()/(4e-5))
            )
        );

    nPsiSubCycles =
        max
        (
            label(1),
            pimple.dict().lookupOrDefault<label>
            (
                "nPsiSubCycles",
                label(runTime.deltaTValue()/(4e-5))
            )
        );



    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::Bridgman100::Bridgman100(fvMesh& mesh)
:
    solid
    (
        mesh,
        autoPtr<solidThermo>(new bridgmanThermo(mesh))
    ),

    thermo_(refCast<bridgmanThermo>(solid::thermo_)),
    topoManager_(mesh),
    dataSharer_(mesh, topoManager_),

    newPoints_(mesh_.points()),
    meshMoveTime_(runTime.value()),
    meshMoveInterval_(0.05),

    dimx(thermo_.dimx()),
    dimt(thermo_.dimt()),
    omega(thermo_.omega()),
    gamma(thermo_.gamma()),
    epsilon(thermo_.epsilon()),
    Ks(thermo_.Ks()),
    Kl(thermo_.Kl()),
    kc(thermo_.kc()),
    Tm(thermo_.Tm()),
    flux(thermo_.flux()),
    Tinf(thermo_.Tinf()),
    Bd(thermo_.Bd()),
    Vpull(thermo_.Vpull()),
    subCycleTime1(runTime.value()),
    ybottom(thermo_.ybottom()),
    limT(2),


    idimx(1/dimx),
    DCpp( ((dimt/(dimx*dimx))*thermo_.Cpp()) ),
    Tls(4/(thermo_.Tl() - thermo_.Ts())),
    pp(-(thermo_.L() * (dimt/(dimx * dimx))) / 1.772 * Tls),


    BPsi_
    (
        IOobject
        (
            "BPsi",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    BPhi_
    (
        IOobject
        (
            "BPhi",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    BBPhi_
    (
        IOobject
        (
            "BBPhi",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    grainT_
    (
        IOobject
        (
            "grainT",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    procID_
    (
        IOobject
        (
            "procID",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Pstream::myProcNo())
    ),


    BPsi(BPsi_),

    BPhi(BPhi_),

    BBPhi(BBPhi_),

    grainT(grainT_),

    currentTime(runTime.value())
{
    //- Read the controls
    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::Bridgman100::~Bridgman100()
{}


// * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * //

void Foam::solvers::Bridgman100::preSolve()
{
    solid::preSolve();

    if ( mesh.changing() || (runTime.value() <= 48e-5) )
    {
        Info<< "Topo_change" << endl;

        //- Release storage of calculated terms
        GBPsi.clear();
        M_GBPsi.clear();
        BdBPGB.clear();

        const volVectorField gradBPsi = fvc::grad(BPsi);
        const volScalarField BPsiEps = BPsi + 1e-6;  


        // Terms calculated to simplify equation, 
        // TEqn, in thermophysicalPredictor function
 
        if (simpleTEqn) 
        {
            GBPsi = new volVectorField
            (
                "GBPsi", 
                (gradBPsi) / BPsiEps
            ); 
        }

        if (flag_1 >=1)
        {
            M_GBPsi = new volScalarField
            (
                "M_GBPsi", 
                mag(gradBPsi) * (flux * idimx)
            ); 
        }

        if (flag_1>2)
        {
            BdBPGB = new volScalarField
            (
                "BdBPGB",
                ( gradBPsi / BPsiEps ) & ( Bd * gradBPsi )
            );
        }

        // Store processor ID's
        procID_.primitiveFieldRef() = Pstream::myProcNo();
    }

    //- Calculate the distance traveled by crucible
    dimensioned<double> travel = runTime.value()*Vpull;

    //- Masking field used to highlight the hot-cold and solid-liquid interfaces
    BBPhi = (
                (pos(BPsi - 0.06) - pos(BPsi-0.94)) *
                ( pos(mesh.C().component(1)*idimx - (ybottom + -0.0007 + travel))
                 -pos(mesh.C().component(1)*idimx - (ybottom + 0.0007 + travel)))
            )
          + ( pos(BPsi-0.06) * ( pos(grainT-(Tm - limT))-pos(grainT-(Tm + limT)) ) );
/*
    const tmp<volScalarField> BPsiPos = pos(BPsi - 0.06);
    const auto y = mesh.C().component(vector::Y)*idimx;

    BBPhi =
    (
        (BPsiPos - pos(BPsi - 0.94))
       *
        (
            pos(y - (ybottom + yRef1 + travel))
          - pos(y - (ybottom + yRef2 + travel))
        )
    )
    +
    (
        BPsiPos
       *
        (
            pos(grainT - (Tm - limT))
          - pos(grainT - (Tm + limT))
        )
    );
*/
}



void Foam::solvers::Bridgman100::moveMesh()
{
    //solid::moveMesh();
}


//////////////////////////////////////////////////////////////////////////////////

void Foam::solvers::Bridgman100::prePredictor()
{
    //- Phase parameter evolution
    BPhi = ( 0.5 - 0.5 * Foam::erf( Tls * (grainT-Tm) ) );
}

//////////////////////////////////////////////////////////////////////////////////

void Foam::solvers::Bridgman100::thermophysicalPredictor()
{
    if(BPsi_diffusion) psiEqn();

    //- subCycle temperature equation, TEqn, form and solve

    List<volScalarField*> stateFieldPtrs({&grainT});

    label nStateSubCycles = (runTime.deltaTValue())/(4e-5);


    if (nStateSubCycles > 1)
    {
        List<volScalarField*> stateFieldPtrs({&grainT});

        for
        (
            subCycle<volScalarField, subCycleFields> stateFieldSubCycle
            (
                stateFieldPtrs,
                nStateSubCycles
            );
            !(++stateFieldSubCycle).end();
        )
        {
            STEqn();
        }
    }
    else
    {
        STEqn();
    }
}


void Foam::solvers::Bridgman100::pressureCorrector()
{
}

//////////////////////////////////////////////////////////////////////////////////

void Foam::solvers::Bridgman100::postCorrector()
{
}

//////////////////////////////////////////////////////////////////////////////////

void Foam::solvers::Bridgman100::postSolve()
{
    // Compute and export all configured entries
    // in shareData fvSolution/subDict
    dataSharer_.exportAll();

    // Access grain front coordinates of mesh/crucible region
    if ((runTime.value()) <= 0.01) return;

    const interfaceTopologyUtils::interfaceBoundingCoords grainFrontcoords =
        dataSharer_.importCoords("grainFrontCoords");

    // Access individual results
    Info<< "Grain front:" << nl
        << "    Top        : " << grainFrontcoords.maxBounds.y() << nl
        << "    Bottom     : " << grainFrontcoords.minBounds.y() << nl
        << "    Left       : " << grainFrontcoords.minBounds.x() << nl
        << "    Right      : " << grainFrontcoords.maxBounds.x() << nl
        << "    formation  : " << grainFrontcoords.formed << nl
        << endl;
}

// ************************************************************************* //
