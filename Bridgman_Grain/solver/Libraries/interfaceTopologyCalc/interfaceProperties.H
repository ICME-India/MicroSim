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
    
struct
    Foam::interfaceTopologyUtils::interfaceProperties

Description
    interfaceProperties holds results of interface topology calculation.
    Default constructed as unformed interface.

\*---------------------------------------------------------------------------*/

#ifndef interfaceTopologyData_H
#define interfaceTopologyData_H

#include "scalar.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace interfaceTopologyUtils
{

/*---------------------------------------------------------------------------*\
    interfaceProperties — result of interface topology calculation
\*---------------------------------------------------------------------------*/

struct interfaceProperties
{
    scalar height;      //- Interface height (negative if concave up)
    scalar yCoord;      //- Y-coordinate of interface centroid
    scalar yTop;        //- Y-coordinate of interface top
    label  concavity;   //- 1=concave down, 2=none, 3=concave up
    bool   formed;      //- True if interface was detected

    interfaceProperties()
    :
        height(0),
        yCoord(0),
        yTop(0),
        concavity(2),
        formed(false)
    {}
};


/*---------------------------------------------------------------------------*\
    interfaceBoundingCoords — axis-aligned bounding box of interface
\*---------------------------------------------------------------------------*/

struct interfaceBoundingCoords
{
    vector minBounds;
    vector maxBounds;
    bool   formed;

    interfaceBoundingCoords()
    :
        minBounds(vector::zero),
        maxBounds(vector::zero),
        formed(false)
    {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfaceTopologyUtils
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
