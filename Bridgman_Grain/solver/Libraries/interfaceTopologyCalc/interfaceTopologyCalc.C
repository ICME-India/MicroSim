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

// * * * * * * * * * * * * * * Interface segmentation * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::interfaceTopologyUtils::calcInterfaceBand
(
    const volScalarField& Field,
    const interfaceTopologySettings& s
)
{
    return
        pos(Field - (s.Tm() - s.limT()))
      - pos(Field - (s.Tm() + s.limT()));
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceTopologyUtils::calcRegionI
(
    const volScalarField& BPsi,
    const volScalarField& band,
    const interfaceTopologySettings& s
)
{
    return pos(BPsi - s.interiorThreshold()) * band;
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceTopologyUtils::calcRegionI
(
    const volScalarField& BPsi,
    const volScalarField& Field,
    const interfaceTopologySettings& s,
    bool
)
{
    tmp<volScalarField> tBand = calcInterfaceBand(Field, s);
    return calcRegionI(BPsi, tBand(), s);
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceTopologyUtils::calcRegionII
(
    const volScalarField& BPsi,
    const volScalarField& band,
    const volScalarField& regionI,
    const interfaceTopologySettings& s
)
{
    const fvMesh& mesh = BPsi.mesh();

    return
        (
            pos(BPsi - s.boundaryThreshold()) * band - regionI
        )
       *pos
        (
            mesh.C().component(vector::Y) * s.idimx()
          - (s.ybottom() + s.yRef4())
        );
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceTopologyUtils::calcRegionII
(
    const volScalarField& BPsi,
    const volScalarField& Field,
    const interfaceTopologySettings& s,
    bool
)
{
    tmp<volScalarField> tBand    = calcInterfaceBand(Field, s);
    tmp<volScalarField> tRegionI = calcRegionI(BPsi, tBand(), s);
    return calcRegionII(BPsi, tBand(), tRegionI(), s);
}


// * * * * * * * * * * * * * * Interface geometry  * * * * * * * * * * * * * //

Foam::interfaceTopologyUtils::interfaceProperties
Foam::interfaceTopologyUtils::calcInterfaceGeometry
(
    const volScalarField& regionI,
    const volScalarField& regionII,
    const interfaceTopologySettings& s
)
{
    interfaceProperties result;

    const fvMesh& mesh = regionI.mesh();

    const scalar sumI  = gSum(regionI.primitiveField());
    const scalar sumII = gSum(regionII.primitiveField());

    if (sumI <= s.interfaceMinCells() || sumII <= s.interfaceMinCells())
        return result;

    const vector center_Ifv =
        mesh.C().weightedAverage(regionII).value();

    tmp<volScalarField> tDist
    (
        volScalarField::New
        (
            "interfaceDist",
            mag(mesh.C()*s.idimx() - center_Ifv)
          + (scalar(1) - pos(regionI - s.regionMaskEps()))*GREAT
        )
    );
    const volScalarField& Dist = tDist();

    result.height = gMin(Dist.primitiveField());

    const scalar cellSize =
        s.maskTolerance() * gAverage(pow(mesh.V().primitiveField(), 1.0/3.0));

    tmp<volScalarField> tMask(pos(cellSize - mag(Dist - result.height)));
    const volScalarField& Mask = tMask();
    const scalar sumMask = max(gSum(Mask.primitiveField()), SMALL);

    result.yCoord =
        gSum
        (
            Mask.primitiveField()
           *(mesh.C().primitiveField().component(vector::Y)
            *s.idimx().value())
        )
       /sumMask;

    const bool concaveUp =
        center_Ifv.component(vector::Y) > result.yCoord;

    result.concavity = concaveUp ? 3 : 1;
    result.height    = concaveUp ? -result.height : result.height;

    result.yTop = gMax
    (
        scalarField
        (
            regionI.primitiveField()
           *(mesh.C().primitiveField().component(vector::Y)
            *s.idimx().value())
          + (scalar(1) - regionI.primitiveField()) * (-GREAT)
        )
    );

    result.formed = true;
    return result;
}


Foam::interfaceTopologyUtils::interfaceProperties
Foam::interfaceTopologyUtils::calcInterfaceGeometry
(
    const volScalarField& BPsi,
    const volScalarField& Field,
    const interfaceTopologySettings& s,
    bool
)
{
    tmp<volScalarField> tBand    = calcInterfaceBand(Field, s);
    tmp<volScalarField> tRegionI = calcRegionI(BPsi, tBand(), s);
    tmp<volScalarField> tRegionII = calcRegionII(BPsi, tBand(), tRegionI(), s);
    return calcInterfaceGeometry(tRegionI(), tRegionII(), s);
}


Foam::interfaceTopologyUtils::interfaceProperties
Foam::interfaceTopologyUtils::calcInterface
(
    const volScalarField& BPsi,
    const volScalarField& Field,
    const interfaceTopologySettings& s
)
{
    tmp<volScalarField> tBand    = calcInterfaceBand(Field, s);
    tmp<volScalarField> tRegionI = calcRegionI(BPsi, tBand(), s);
    tmp<volScalarField> tRegionII = calcRegionII(BPsi, tBand(), tRegionI(), s);
    return calcInterfaceGeometry(tRegionI(), tRegionII(), s);
}


// * * * * * * * * * * * * * * Interface bounds   * * * * * * * * * * * * * //

Foam::interfaceTopologyUtils::interfaceBoundingCoords
Foam::interfaceTopologyUtils::findFieldBounds
(
    const volScalarField& regionI,     // plain ref — not tmp
    const interfaceTopologySettings& s
)
{
    interfaceBoundingCoords result;

    // Use plain ref — no () operator
    const scalarField& mask = regionI.primitiveField();
    const scalar sumI = gSum(mask);

    if (sumI <= s.interfaceMinCells())
    {
        result.formed = false;
        return result;
    }

    const scalarField invMask = (scalar(1) - mask) * GREAT;

    const volVectorField& C  = regionI.mesh().C();
    const vectorField&    Ci = C.primitiveField();

    result.minBounds.x() = gMin(Ci.component(vector::X) + invMask);
    result.minBounds.y()  = gMin(Ci.component(vector::Y) + invMask);
    result.minBounds.z()  = gMin(Ci.component(vector::Z) + invMask);
    result.maxBounds.x() = gMax(Ci.component(vector::X) - invMask);
    result.maxBounds.y() = gMax(Ci.component(vector::Y) - invMask);
    result.maxBounds.z() = gMax(Ci.component(vector::Z) - invMask);

    result.formed = true;
    return result;
}


Foam::interfaceTopologyUtils::interfaceBoundingCoords
Foam::interfaceTopologyUtils::calcInterfaceBounds
(
    const volScalarField& BPsi,
    const volScalarField& Field,
    const interfaceTopologySettings& s
)
{
    tmp<volScalarField> tBand    = calcInterfaceBand(Field, s);
    tmp<volScalarField> tRegionI = calcRegionI(BPsi, tBand(), s);
    tBand.clear();

    // Pass tRegionI() — const ref from tmp, valid until tRegionI cleared
    return findFieldBounds(tRegionI(), s);
}

// ************************************************************************* //
