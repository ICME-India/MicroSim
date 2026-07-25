/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#include "interfaceTopology.H"
#include "globalMeshData.H"
#include "dimensionSets.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::interfaceTopologyUtils::interfaceBoundingCoords
Foam::interfaceTopologyUtils::loopCalcInterfaceBounds
(
    const volScalarField& BPsi,
    const volScalarField& Field,
    const interfaceTopologySettings& s
)
{
interfaceBoundingCoords result;

const fvMesh& mesh = Field.mesh();

const scalar Tlow  = s.Tm() - s.limT();
const scalar Thigh = s.Tm() + s.limT();
const scalar BPsiThreshold = s.interiorThreshold();

const vectorField& C = mesh.C();

vector minPos(GREAT, GREAT, GREAT);
vector maxPos(-GREAT, -GREAT, -GREAT);

scalar sumI = 0;

forAll(Field, celli)
{
    if
    (
        BPsi[celli] > BPsiThreshold
     && Field[celli] > Tlow
     && Field[celli] < Thigh
    )
    {
    
        const vector& p = C[celli];

        minPos.x() = min(minPos.x(), p.x());
        minPos.y() = min(minPos.y(), p.y());
        minPos.z() = min(minPos.z(), p.z());

        maxPos.x() = max(maxPos.x(), p.x());
        maxPos.y() = max(maxPos.y(), p.y());
        maxPos.z() = max(maxPos.z(), p.z());

        ++sumI;
    }
}

Info<< "iNTERFACE iNPUTS:" << Thigh << "  "<<BPsiThreshold << "  "<<Tlow << endl;

sumI = returnReduce(sumI, sumOp<scalar>());

if (sumI <= s.interfaceMinCells())
{
    result.formed = false;
    return result;
}

result.minBounds.x() = returnReduce(minPos.x(), minOp<scalar>());
result.minBounds.y() = returnReduce(minPos.y(), minOp<scalar>());
result.minBounds.z() = returnReduce(minPos.z(), minOp<scalar>());

result.maxBounds.x() = returnReduce(maxPos.x(), maxOp<scalar>());
result.maxBounds.y() = returnReduce(maxPos.y(), maxOp<scalar>());
result.maxBounds.z() = returnReduce(maxPos.z(), maxOp<scalar>());

result.formed = true;

return result;
}

// ************************************************************************* //
