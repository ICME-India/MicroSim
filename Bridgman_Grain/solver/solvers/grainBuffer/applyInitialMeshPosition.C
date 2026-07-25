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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::solvers::grainBuffer::applyInitialMeshPosition()
{

    // Determine target position
    
    const bool positionPrescribed =
        mag(meshInitialPosition_) < mag(vector(GREAT, GREAT, GREAT));

    vector targetPosition = vector::zero;

    if (positionPrescribed)
    {
        targetPosition = meshInitialPosition_;

        Info<< "nbrRegionFieldMapper: using prescribed initial position "
            << targetPosition << endl;
    }
    else if
    (
        refMeshName_  != word::null
     && refFieldName_ != word::null
    )
    {
        if (!runTime.foundObject<fvMesh>(refMeshName_))
        {
            WarningInFunction
                << "Reference mesh '" << refMeshName_
                << "' not found. Initial positioning skipped."
                << endl;
            return;
        }

        const fvMesh& refMesh =
            runTime.lookupObject<fvMesh>(refMeshName_);

        if (!refMesh.foundObject<volScalarField>(refFieldName_))
        {
            WarningInFunction
                << "Reference field '" << refFieldName_
                << "' not found in mesh '" << refMeshName_
                << "'. Initial positioning skipped."
                << endl;
            return;
        }

        const volScalarField& refField =
            refMesh.lookupObject<volScalarField>(refFieldName_);

        const vectorField& cellCentres = refMesh.C();

        scalar targetY = GREAT;
        scalar sumX = 0.0;
        scalar sumZ = 0.0;
        scalar nInterface = 0.0;

        forAll(refField, celli)
        {
            if (refField[celli] > refFieldThreshold_)
            {
                const point& c = cellCentres[celli];

                targetY = min(targetY, c.y());
                sumX += c.x();
                sumZ += c.z();
                nInterface += 1.0;
            }
        }

        // Reduce FIRST — unconditionally on every processor
        targetY    = returnReduce(targetY, minOp<scalar>());
        sumX       = returnReduce(sumX, sumOp<scalar>());
        sumZ       = returnReduce(sumZ, sumOp<scalar>());
        nInterface = returnReduce(nInterface, sumOp<scalar>());

        // THEN check global value — safe across all processors
        if (nInterface < SMALL)
        {
            WarningInFunction
                << "No cells found where '" << refFieldName_
                << "' >= " << refFieldThreshold_
                << " in mesh '" << refMeshName_
                << "'. Initial positioning skipped."
                << endl;
            return;
        }

        const scalar targetX = sumX / nInterface;
        const scalar targetZ = sumZ / nInterface;

        targetPosition = vector(targetX, targetY, targetZ);

        Info<< "nbrRegionFieldMapper: auto-detected target position" << nl
            << "    interface bottom Y = " << targetY << nl
            << "    interface centre X = " << targetX << nl
            << "    interface centre Z = " << targetZ << endl;
    }
    else
    {
        return;
    }


    // Current mesh position 
    
    const scalar minX = gMin(newPoints_.component(vector::X));
    const scalar maxX = gMax(newPoints_.component(vector::X));

    const scalar minY = gMin(newPoints_.component(vector::Y));

    const scalar minZ = gMin(newPoints_.component(vector::Z));
    const scalar maxZ = gMax(newPoints_.component(vector::Z));

    const vector currentPosition
    (
        0.5*(minX + maxX),
        minY,
        0.5*(minZ + maxZ)
    );
        

    // Shifting mesh

    const vector shift = targetPosition - currentPosition
                        + meshInitialOffset_;

    if (mag(shift) < SMALL)
    {
        Info<< "nbrRegionFieldMapper: mesh already at target. "
             << "No shift applied." << endl;
            return;
    }

    newPoints_ += shift;
    mesh_.movePoints(newPoints_);
    
    //mesh_.clearOut();
    stencilOutdated_ = true;
    meshPositioned_  = true;

    Info<< "nbrRegionFieldMapper: mesh shifted by " << shift
        << " to reach position " << targetPosition << endl;
     
    // Update nbrMesh shift shift    
    fShift_ = shift;    
 Info << "I am here3" << endl; 
}

// ************************************************************************* //
