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

#include "grainBuffer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

    defineTypeNameAndDebug(grainBuffer, 0);
    addToRunTimeSelectionTable(solver, grainBuffer, fvMesh);

} // End namespace solvers
} // End namespace Foam


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::solvers::grainBuffer::dependenciesModified() const
{
    return runTime.controlDict().modified() || mesh.solution().modified();
}


bool Foam::solvers::grainBuffer::read()
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


bool Foam::solvers::grainBuffer::checkNbrResources()
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


void Foam::solvers::grainBuffer::initialiseNbrResources()
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

    Info<< "grainBuffer: initialised " << nbrMeshNames_.size()
        << " neighbour mesh(es) and associated fields." << endl;
}


void Foam::solvers::grainBuffer::mappingUpdate()
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
            
            mappedField.correctBoundaryConditions();
        }
    }

    stencilOutdated_ = false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::grainBuffer::grainBuffer(fvMesh& mesh)
:
    solver(mesh),
    
    
    topoManager_(mesh),
    dataSharer_(mesh, topoManager_),

    sharedData
    (
        IOobject
        (
            "sharedData",
            runTime.name(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    
    
    newPoints_(mesh_.points()),
    meshMoveTime_(runTime.value()),
    shiftY_(0.0),

    nbrResourcesAvailable_(false),
    
    stencilOutdated_(true),
    
    mapNbrFields_(false),
    nbrMapTime_(runTime.value()),  
    
    bufferField_
    (
        IOobject
        (
            "bufferField",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0)
    )
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
    
    vector F = vector(0, 0 , 0);
    
    sharedData.add("shiftVector", F, true); 
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::grainBuffer::~grainBuffer()
{}


// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::solvers::grainBuffer::maxDeltaT() const
{
    scalar deltaT = min(fvModels().maxDeltaT(), vGreat);

    return deltaT;
}


void Foam::solvers::grainBuffer::preSolve()
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
       
}


void Foam::solvers::grainBuffer::motionCorrector()
{}


void Foam::solvers::grainBuffer::prePredictor()
{}


void Foam::solvers::grainBuffer::momentumPredictor()
{}


void Foam::solvers::grainBuffer::thermophysicalPredictor()
{

}


void Foam::solvers::grainBuffer::pressureCorrector()
{}


void Foam::solvers::grainBuffer::postCorrector()
{}


void Foam::solvers::grainBuffer::postSolve()
{             
     
}

// ************************************************************************* //
