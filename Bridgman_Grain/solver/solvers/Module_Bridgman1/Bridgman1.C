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

#include "Bridgman1.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

    defineTypeNameAndDebug(Bridgman1, 0);
    addToRunTimeSelectionTable(solver, Bridgman1, fvMesh);

} // End namespace solvers
} // End namespace Foam


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::solvers::Bridgman1::dependenciesModified() const
{
    return runTime.controlDict().modified() || mesh.solution().modified();
}


bool Foam::solvers::Bridgman1::read()
{
    solver::read();

    meshMoveInterval_ =
        pimple.dict().lookupOrDefault<scalar>("meshMoveInterval", 0.05);

    nbrMapInterval_ =
        pimple.dict().lookupOrDefault<scalar>("nbrMapInterval", -1);

    // Read neighbour region and field configuration
    const dictionary& nbrRegionsDict =
        pimple.dict().subDict("neighbourRegions");


    nbrMeshNames_ = nbrRegionsDict.toc();
    

    forAll(nbrMeshNames_, i)
    {
        const word& meshName = nbrMeshNames_[i];
        const dictionary& meshDict = nbrRegionsDict.subDict(meshName);
        nbrFieldNames_.insert(meshName, wordList(meshDict.lookup("fields")));
    }
    
    // Read initial positioning config (optional block)
    if (pimple.dict().found("initialPositioning"))
    {
        const dictionary& posDict =
            pimple.dict().subDict("initialPositioning");

        refMeshName_  = posDict.lookupOrDefault<word>("region",  word::null);
        refFieldName_ = posDict.lookupOrDefault<word>("field",   word::null);

        refFieldThreshold_ =
            posDict.lookupOrDefault<scalar>("fieldThreshold", 0.96);
            
        meshInitialPosition_ =
            posDict.lookupOrDefault<vector>
            (
                "initialPosition",
                vector(GREAT, GREAT, GREAT)
            );

        meshInitialOffset_ =
            posDict.lookupOrDefault<vector>("initialOffset", vector::zero);
    }

    return true;
}


bool Foam::solvers::Bridgman1::checkNbrResources()
{
    // Check all requested meshes and fields are registered
    forAll(nbrMeshNames_, i)
    {
        const word& meshName = nbrMeshNames_[i];

        if (!runTime.foundObject<fvMesh>(meshName))
        {
            WarningInFunction
                << "Neighbour mesh '" << meshName
                << "' not found in object registry. "
                << "All mapping operations will be skipped."
                << endl;
            return false;
        }

        const fvMesh& nbr = runTime.lookupObject<fvMesh>(meshName);

        const wordList& fieldNames = nbrFieldNames_[meshName];
        forAll(fieldNames, j)
        {
            if (!nbr.foundObject<volScalarField>(fieldNames[j]))
            {
                WarningInFunction
                    << "Field '" << fieldNames[j]
                    << "' not found in neighbour mesh '" << meshName
                    << "'. All mapping operations will be skipped."
                    << endl;
                return false;
            }
        }
    }
    return true;
}


void Foam::solvers::Bridgman1::initialiseNbrResources()
{
    forAll(nbrMeshNames_, i)
    {
        const word& meshName = nbrMeshNames_[i];

        // Register neighbour mesh pointer
        const fvMesh& nbr = runTime.lookupObject<fvMesh>(meshName);
        nbrMeshPtrs_.insert(meshName, &nbr);

        // Create interpolator for this mesh pair
        // (no geometry computed yet — update() is called in mappingUpdate())
        interpolatorPtrs_.insert
        (
            meshName,
            cellsToCells::New("intersection")
        );

        // Register field pointers and create mapped storage fields
        const wordList& fieldNames = nbrFieldNames_[meshName];
        forAll(fieldNames, j)
        {
            const word key = meshName + "::" + fieldNames[j];

            // Pointer to source field on neighbour mesh
            nbrFieldPtrs_.insert
            (
                key,
                &nbr.lookupObject<volScalarField>(fieldNames[j])
            );

            // Mapped field on this mesh — read if restarting,
            // initialised to zero otherwise
            mappedFieldPtrs_.insert
            (
                key,
                new volScalarField
                (
                    IOobject
                    (
                        "mapped_" + meshName + "_" + fieldNames[j],
                        runTime.name(),
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedScalar(dimless, Zero)
                )
            );
        }
    }

    Info<< "Bridgman1: initialised " << nbrMeshNames_.size()
        << " neighbour mesh(es) and associated fields." << endl;
}


void Foam::solvers::Bridgman1::mappingUpdate()
{
    forAll(nbrMeshNames_, i)
    {
        const word& meshName = nbrMeshNames_[i];
        const fvMesh& nbr = *nbrMeshPtrs_[meshName];

        // Rebuild the intersection stencil only when mesh positions have
        // changed — avoids expensive geometric intersection every step
        if (stencilOutdated_)
        {
            interpolatorPtrs_[meshName]->update(mesh(), nbr);            
        }
        

        const wordList& fieldNames = nbrFieldNames_[meshName];
        
        
        forAll(fieldNames, j)
        {
            const word key = meshName + "::" + fieldNames[j];

            const volScalarField& nbrField = *nbrFieldPtrs_[key];
            volScalarField& mappedField    = *mappedFieldPtrs_[key];
            
            
            // Interpolate neighbour temperature into mappedField.
            // tgtToSrc maps from nbrMesh (target) onto this mesh (source).
            scalarField buffer(mappedField.primitiveField());
            buffer = interpolatorPtrs_[meshName]().tgtToSrc(nbrField, buffer);
            mappedField.primitiveFieldRef() = buffer;
                      

            // Manual linear interpolation
            
            if (meshName == "crucible")
            {
                        
                const vectorField& fineC = mesh.C();
            
                // Build the gradient field on coarse mesh (as you do)
                const volVectorField gradCoarse(fvc::grad(nbrField));

                // Map BOTH the coarse cell centre AND the coarse value
                // AND the coarse gradient onto the fine mesh 
                vectorField centreBuffer(mesh.C());      
                centreBuffer = 
                   interpolatorPtrs_[meshName]().tgtToSrc(nbr.C(), centreBuffer);
    
                vectorField gradBuffer(0*centreBuffer);
                gradBuffer = interpolatorPtrs_[meshName]().tgtToSrc(gradCoarse, gradBuffer);

                // LINEAR RECONSTRUCTION using only LOCAL,
                // already-correctly-mapped fine-mesh data 
                
                forAll(mappedField, celli)
                {
                    const vector& gradVal  = gradBuffer[celli];     // mapped coarse gradient

                    if (mag(gradVal) > 100)
                    {
                        mappedField[celli] = 
                            mappedField[celli] 
                          + (gradVal & (fineC[celli] - centreBuffer[celli]));
                    }
                }       
            }
            
            mappedField.correctBoundaryConditions();        
        }
    }

    stencilOutdated_ = false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::Bridgman1::Bridgman1(fvMesh& mesh)
:
    solver(mesh),
    
    topoManager_(mesh),
    dataSharer_(mesh, topoManager_),
    //meshPointsPredictor_(meshMotionControl::New(mesh)),
    
    newPoints_(mesh_.points()),
    predPoints_(newPoints_),
    meshMoveTime_(runTime.value()),
    shiftY_(0.0),
    fShift_(Foam::vector::zero),

    nbrResourcesAvailable_(false),
    
    stencilOutdated_(true),
    
    mapNbrFields_(false),
    nbrMapTime_(runTime.value()),  

    moveDue_(false),
    
    sProcID_
    (
        IOobject
        (
            "sProcID",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Pstream::myProcNo())
    ),
    
    sPhi_
    (
        IOobject
        (
            "sPhi",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    
    remeshRegion_
    (
        IOobject
        (
            "sPhi",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    
    liquidusTCoords_(),      
    solidusTCoords_(),  
      
    sPhi(sPhi_),
    remeshRegion(remeshRegion_)
{
    read();
    

    // Attempt to connect to all requested neighbour meshes & fields
    nbrResourcesAvailable_ = checkNbrResources();

    if (nbrResourcesAvailable_)
    {
        initialiseNbrResources();
        
        // Align mesh to domain parameter bottom before first mapping
        applyInitialMeshPosition();

        // Build initial stencils and map all fields
        mappingUpdate();
    }

    sPhi = 1;
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::Bridgman1::~Bridgman1()
{}


// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::solvers::Bridgman1::maxDeltaT() const
{
    scalar deltaT = min(fvModels().maxDeltaT(), vGreat);

    return deltaT;
}


void Foam::solvers::Bridgman1::preSolve()
{
    // Skip everything if neighbour resources are unavailable
    if (!nbrResourcesAvailable_) return;

    fvModels().preUpdateMesh();

    // Update the mesh for topology change, mesh to mesh mapping
    if (mesh().topoChanged())
    {
        mesh_.update();

        pointField staticPointsCopy(mesh_.points());
        mesh_.movePoints(staticPointsCopy);

        newPoints_ = mesh_.points();
        stencilOutdated_ = true;
    }

    // Neighbour mesh topology change also invalidates stencil
   // forAll(nbrMeshNames_, i)
   // {
        if (nbrMeshPtrs_[nbrMeshNames_[0]]->changing())
        {
            stencilOutdated_ = true;
            //break;
        }
   // }

    /*
    // Update mask for dynamic remeshing
    if (runTime.value() > 0.1 && runTime.value() < 0.13)
    {
       const dimensionedScalar idimx("idimx", dimless/dimLength,1.0);

       tmp<volScalarField> tCoordY = mesh.C().component(vector::Y)*idimx;

       const volScalarField& coordY = tCoordY();

       scalar MaxY= gMax(coordY);

       sPhi = - pos(coordY - (MaxY-0.0002))
            + pos(coordY - (MaxY-0.0004));

       Info<< "max small mesh point" << MaxY;
    }
    */

    sProcID_.primitiveFieldRef() = Pstream::myProcNo();

    fShift_ = 0*fShift_;

    moveDue_ = false;
    
    if (stencilOutdated_)
    {
        // Build initial stencils and map all fields
        mappingUpdate(); 
    }
}


void Foam::solvers::Bridgman1::motionCorrector()
{}


void Foam::solvers::Bridgman1::prePredictor()
{}


void Foam::solvers::Bridgman1::momentumPredictor()
{}


void Foam::solvers::Bridgman1::thermophysicalPredictor()
{}


void Foam::solvers::Bridgman1::pressureCorrector()
{}


void Foam::solvers::Bridgman1::postCorrector()
{}


void Foam::solvers::Bridgman1::postSolve()
{
    // Compute and export all configured entries
    // in shareData fvSolution/subDict
    dataSharer_.exportAll();


    //if (!dataSharer_.sharedEntryAvailable("liquidusTCoords")) return;

    // Access liquidus and solidus coordinates of nbrMesh/crucible region
    //if ((runTime.value()) <= 0.01) return;
    liquidusTCoords_ = dataSharer_.importCoords("liquidusTCoords");
    solidusTCoords_ = dataSharer_.importCoords("solidusTCoords");

    Info<< "liquidus interface:" << nl
        << "    Top    : " << liquidusTCoords_.maxBounds[1] << nl
        << "    Bottom : " << liquidusTCoords_.minBounds.y() << nl
        << "    Left   : " << liquidusTCoords_.minBounds.x() << nl
        << "    Right  : " << liquidusTCoords_.maxBounds.x() << nl
        << "    formation  : " << liquidusTCoords_.formed << nl
        << endl;

    Info<< "solidus interface:" << nl
        << "    Top        : " << solidusTCoords_.maxBounds.y() << nl
        << "    Bottom     : " << solidusTCoords_.minBounds.y() << nl
        << "    Left       : " << solidusTCoords_.minBounds.x() << nl
        << "    Right      : " << solidusTCoords_.maxBounds.x() << nl
        << "    formation  : " << solidusTCoords_.formed << nl
        << endl;
}

// ************************************************************************* //
