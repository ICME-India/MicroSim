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

#include "grainFrontSim.H"
#include "IOstream.H"
#include "addToRunTimeSelectionTable.H"

#include "point.H"
#include "fvcSurfaceIntegrate.H"
#include "localEulerDdtScheme.H"
#include "surfaceFields.H"
#include "fvMatrix.H"
#include "ddtScheme.H"
#include "fvm.H"
#include "fvmLaplacian.H"
#include "fvc.H"
#include "quaternion.H"

#include "readPhysicalDict.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(grainFrontSim, 0);
    addToRunTimeSelectionTable(solver, grainFrontSim, fvMesh);

} // End namespace solvers
} // End namespace Foam


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //


void Foam::solvers::grainFrontSim::createFields()
{
    PhaseFields.setSize(NumPhases);
    forAll(PhaseFields, a)
    {
        if (a == NumPhases -1)
        {
            PhaseFields[a] = createVSF(mesh, "liq");
            break;
        }       
        PhaseFields[a] = createVSF(mesh, "g" + str(1+a));
    }

    TempGrain = createVSF(mesh, "temp");
    TempMod     = createVSF_NO(mesh, "temp_mod");
    aniso       = createVSF_NO(mesh, "aniso");
    Velocity    = createVSF_NO(mesh, "velocity");

    AllPhases   = createVSF(mesh, "allphases");
    SumPhases   = createVSF(mesh, "sumphases");
    // GradPhase   = createVVF(mesh, "gradphi");
    // PhiNCap     = createVVF(mesh, "phincap");
}


void Foam::solvers::grainFrontSim::createInitialConditions()
{  
    scalar s = 2e-6;
   
    scalar BaseX  = gMin
    (
        mesh.C().primitiveField().component(0)
    );
    
    scalar BaseY  = gMin
    (
        mesh.C().primitiveField().component(1)
    );
    
    scalar BaseZ  = gMin
    (
        mesh.C().primitiveField().component(2)
    );
   
    scalar BaseYinc = BaseY + s*6;
  
    scalar BaseZinc = BaseZ + s*2;
 

    fillCuboid(PHI(0), point(BaseX,   BaseY, BaseZ),  point(BaseX + s*69, BaseYinc, BaseZinc), 1);
    fillCuboid(PHI(1), point(BaseX + s*69,  BaseY, BaseZ),  point(BaseX + s*138, BaseYinc, BaseZinc), 1);
    fillCuboid(PHI(2), point(BaseX + s*138, BaseY, BaseZ),  point(BaseX + s*207, BaseYinc, BaseZinc), 1);


/*
    scalar BaseX = gMin(mesh.C().primitiveField().component(0));
    scalar BaseY = gMin(mesh.C().primitiveField().component(1));
    scalar BaseZ = gMin(mesh.C().primitiveField().component(2));

    scalar Lx = 160;
    scalar Lz = 160;

    label Nx = 4;
    label Nz = 4;

    scalar dx = Lx / Nx;
    scalar dz = Lz / Nz;

    label grainID = 0;

    for (label ix = 0; ix < Nx; ix++)
    {
        for (label iz = 0; iz < Nz; iz++)
        {
            point p1
            (
                BaseX + ix*dx,
                BaseY,
                BaseZ + iz*dz
            );

            point p2
            (
                BaseX + (ix+1)*dx,
                BaseY + s*6,
                BaseZ + (iz+1)*dz
            );

            fillCuboid(PHI(grainID), p1, p2, 1);

            grainID++;
        }
     }
*/
    fillLastPhase(PhaseFields);
    
    
    vector TreferencePoint( BaseX, BaseY, BaseZ);
    
    fillTemperatureGradient(TEMP(), gfc.Treference, gfc.TGradient, TreferencePoint);  
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::solvers::grainFrontSim::dependenciesModified() const
{
    return Bridgman1::dependenciesModified() || mesh.solution().modified();
}


bool Foam::solvers::grainFrontSim::read()
{
    Bridgman1::read();
           
    gfc.printConstants();  

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::grainFrontSim::grainFrontSim(fvMesh& mesh)
:
    Bridgman1(mesh),
    gfc
    (
        mpmsc::GrainFieldConstants_t
        (
            IOdictionary
            (
                IOobject
                (
                    "physicalProperties",
                    mesh.time().constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            ),
            mesh,
            mesh.time()
        )
    ),

    NumPhases(gfc.NumPhases),
    NumComps(4)  
{
    read();  
    createFields();
    createInitialConditions();  
    
    pointField staticPointsCopy(mesh_.points());
    mesh_.movePoints(staticPointsCopy);  
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::grainFrontSim::~grainFrontSim()
{}


// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::solvers::grainFrontSim::maxDeltaT() const
{
    scalar deltaT = min(fvModels().maxDeltaT(), vGreat);

    return deltaT;
}


void Foam::solvers::grainFrontSim::preSolve()
{
    Bridgman1::preSolve();
}

void Foam::solvers::grainFrontSim::moveMesh()
{

    shiftY_ = 0.0;

    // Translate mesh if levelSet is beyond meshMoveSwitch

    topoManager_.calc
    (
        "SumPhases1", 
        calcMode::interfaceBounds
    );

    interfaceBoundingCoords SumPhasesCoords_ =
        topoManager_.coordinates("SumPhases1");


    if (SumPhasesCoords_.formed)
    {
        const scalar maxY = gMax(mesh_.points().component(vector::Y));
        const scalar minY = gMin(mesh_.points().component(vector::Y));

        const scalar meshMoveSwitch =
            maxY - 0.2* (maxY - minY);

        const scalar meshMoveTarget =
            maxY - 0.5* (maxY - minY);

        if (SumPhasesCoords_.maxBounds.y() > meshMoveSwitch)
        {
            shiftY_ =
                SumPhasesCoords_.maxBounds.y()
              - meshMoveTarget;

            moveDue_ = true;

            Info << "levelSet tried to escape!" <<  meshMoveSwitch << "and"     
                 <<SumPhasesCoords_.maxBounds.y() << "and" << SumPhasesCoords_.minBounds.y() 
                 << "mesh" << maxY << "and" << minY << endl;
        }
    }


    Bridgman1::moveMesh();
}


void Foam::solvers::grainFrontSim::motionCorrector()
{}


void Foam::solvers::grainFrontSim::prePredictor()
{}


void Foam::solvers::grainFrontSim::momentumPredictor()
{}


void Foam::solvers::grainFrontSim::thermophysicalPredictor()
{


    const scalar BaseY  = gMin(mesh_.points().component(vector::Y));

    const scalar TopY  = gMax(mesh_.points().component(vector::Y));

    const scalar solveSwitch = BaseY + 0.01 * (TopY - BaseY);

    bool evolveGrains =
        liquidusTCoords_.formed
     && (liquidusTCoords_.minBounds.y() > solveSwitch);


    if (!evolveGrains)
    {
        SumPhases() = 0.0;

        forAll(PhaseFields, a)
        {
            if (a == (NumPhases-1)){continue;}
            SumPhases() += PHI(a);
        }
/*
        fvScalarMatrix
        (
            lapCoeff  * fvm::laplacian(PHI(0))
        );
*/
        Info << "Didn't solve for grains" << endl;
        
        return;
    }


    auto ddtCoeff = dt;
    auto lapCoeff = 2 * idx2 * gfc.Gamma * gfc.Epsilon;
    auto funCoeff = 9 * gfc.Gamma / gfc.Epsilon;
    auto sigCoeff = 2 * idx  * gfc.Gamma * gfc.Epsilon * gfc.Coeff_AntiCap;
    auto invCoeff = gfc.Epsilon / sqrt(2.0);

    auto&& PHI_LIQ = PHI(NumPhases-1);

    tmpVec3 grad_phil  = idx * fvc::grad(PHI_LIQ);

    //////////////////////////////////////////////////////////
    const word& meshName = nbrMeshNames_[0];
    const wordList& fieldNames = nbrFieldNames_[meshName];
    const word key1 = nbrMeshNames_[0] + "::" + fieldNames[0];

    Info << "meshSelected" << nbrMeshNames_[0] << endl;

    if (TempMod.valid() && mappedFieldPtrs_.found(key1))
    {
        volScalarField* ptr = mappedFieldPtrs_[key1];

        if (ptr)
        {
            TEMP() = *ptr;
        }
    }


    ///////////////////////////////////////////////////////////

    TempMod() = TEMP();

    if (gfc.UseTemperatureCorrector)
    {
        tmpVec3 ncap_phil  = (grad_phil() / (1e-9 + mag(grad_phil())));

        scalar tol = 1e-6;
        auto phi_mod = max(min(PHI_LIQ, 1.0-tol), tol);

        TempMod()   = TempMod()
                    - ((ncap_phil() & (idx * fvc::grad(TEMP())))
                    * invCoeff * atanh(2*phi_mod - 1));

    }

    tmpReal dTemp       =   mag(gfc.Tequlibrium - TempMod());

    tmpReal strength    =   gfc.Coeff_A1 * dTemp()
                        +   gfc.Coeff_A2 * pow(dTemp(), 2)
                        +   gfc.Coeff_A3 * pow(dTemp(), 3);

    tmpReal Mobility    =   (1.0-PHI_LIQ) * gfc.Mobility_sol
                        +   PHI_LIQ       * gfc.Mobility_liq;

    SumPhases() = 0.0;

        
    forAll(PhaseFields, a)
    {
        if (a == (NumPhases-1)){continue;}
        if (gSum(PHI(a)) < 0.01){continue;}


        tmpVec3 grad_phia = idx * fvc::grad(PHI(a));
        tmpReal magg_phia = mag(grad_phia());
        tmpVec3 ncap_phia = (grad_phia() / (1e-9 + magg_phia()));

        tmpVec3 qab_vec   = PHI(a) * grad_phil() - PHI_LIQ * grad_phia();
        tmpVec3 qab_cap   = (qab_vec() / (1e-9 + mag(qab_vec())));

        auto ax_xx  = gfc.Rotation[a] & dir_x;
        auto ax_yy  = gfc.Rotation[a] & dir_y;

        aniso()     = max(mag(ax_xx & qab_cap()), mag(ax_yy & qab_cap()));

        Velocity()  = strength() * aniso()
                    * fun_H_phi(PHI_LIQ) * gfc.Coeff_H_phi;


        scalar t0 = timer.elapsedTime();

         fvScalarMatrix
         (
             - ddtCoeff   * fvm::ddt(PHI(a))
             +   Mobility() *
             (
                  lapCoeff  * fvm::laplacian(PHI(a))
              -   funCoeff  * fun_g_prime(PHI(a))
              -   sigCoeff  * fvc::div(ncap_phia()) * magg_phia()
             )
             + Velocity() * magg_phia()
         ).solve();
         
         Info<< "grainC " << a << " solve = "
             << timer.elapsedTime()-t0
             << " s" << endl;

         SumPhases() += PHI(a);
    }

    PHI_LIQ = 1.0 - SumPhases();

    forAll(PhaseFields, a)
    {
        PHI(a) = max(PHI(a), 0.0);
        PHI(a) = min(PHI(a), 1.0);
    }


    AllPhases() = 0;
    forAll(PhaseFields, a)
    {
        if (a == (NumPhases-1)){continue;}
        AllPhases() += gfc.Angles_Z[a] * pos(PHI(a)-0.5);
    }

}

void Foam::solvers::grainFrontSim::pressureCorrector()
{    
    if (moveDue_)
    {
     	forAll(PhaseFields, a)
        {
            if (a == (NumPhases-1)){continue;}

            // SOURCE field = current alpha
            const volScalarField& alphaSrc = PHI(a);

            scalarField buffer(alphaSrc.primitiveField().size(),0.0);

            buffer =
                interpolatorPtrs_["bufferRegion"]().tgtToSrc
                (
                    alphaSrc,
                    buffer
                );

            // overwrite alpha
            volScalarField& alphaOut = PHI(a);

            alphaOut.primitiveFieldRef() = buffer;
            alphaOut.correctBoundaryConditions();
        }
    }
}


void Foam::solvers::grainFrontSim::postCorrector()
{}


void Foam::solvers::grainFrontSim::postSolve()
{    
    Bridgman1::postSolve();
    
    //remeshRegion
}

// ************************************************************************* //
