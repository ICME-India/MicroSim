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
#include "MULES.H"
#include "EulerDdtScheme.H"
#include "fvcDiv.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::solvers::Bridgman1::moveMesh()
{
    // Skip everything if neighbour resources are unavailable
    if (!nbrResourcesAvailable_) return;

     moveDue_ = runTime.value() > meshMoveTime_ || moveDue_;

    if (moveDue_)
    {

        // Before move — store old face centre
        // Translate mesh by one pull increment in the y-direction

        mesh_.movePoints(mesh_.points() + vector(0, shiftY_, 0));

        meshMoveTime_ += meshMoveInterval_;

        // Mesh moved — stencil is outdated
        stencilOutdated_ = true;


        fShift_ += vector(0, shiftY_, 0);

        dictionary& sharedData =
            runTime.lookupObjectRef<dictionary>("sharedData");

        sharedData.set("shiftVector", fShift_);

        const word& meshName = "bufferRegion";
        const fvMesh& nbr = *nbrMeshPtrs_[meshName];

        // Rebuild the intersection stencil only when mesh positions have
        // changed — avoids expensive geometric intersection every step
        interpolatorPtrs_[meshName]->update(mesh(), nbr);

	shiftY_ = 0;
    }

    // Determine whether to map neighbour fields this step.
    // Map if:
    //   - meshes moved or topology changed (stencil was just rebuilt), OR
    //   - the independent mapping deadline has elapsed, OR
    //   - nbrMapInterval_ < 0 (map every step)
    const bool mapDue = runTime.value() > nbrMapTime_;
    if (mapDue && nbrMapInterval_ >= 0)
    {
        nbrMapTime_ += nbrMapInterval_;
    }

    // Check overlap against all neighbour meshes
    bool anyOverlap = false;
    forAll(nbrMeshNames_, i)
    {
        if (mesh().bounds().overlaps(nbrMeshPtrs_[nbrMeshNames_[i]]->bounds()))
        {
            anyOverlap = true;
            break;
        }
    }

    mapNbrFields_ = anyOverlap
        && (stencilOutdated_ || mapDue || nbrMapInterval_ < 0);

    if (mapNbrFields_) mappingUpdate();
}

// ************************************************************************* //
