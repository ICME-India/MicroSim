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
    
Description
    Functions for identification and quantification of interface topology 
    in a field
    
    Two calling styles:

    Style A — Using pre-computed regions (regionI/regionII) 
              and dictionary (settings):

    Style B — Directly from domain parameter (BPsi), field (Field), and 
              dictionay (settings):

Usage
    Example usage:
    \verbatim    
        auto props   = calcInterfaceGeometry(regionI(), regionII(), settings);

        auto props = calcInterfaceGeometry(BPsi, Field, settings, true);
    \endverbatim
            
SourceFiles
    interfaceTopologyCalc.C
    
\*---------------------------------------------------------------------------*/

#ifndef interfaceTopology_H
#define interfaceTopology_H

#include "volFields.H"
#include "interfaceTopologyManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace interfaceTopologyUtils
{

    // Interface band
    tmp<volScalarField> calcInterfaceBand
    (const volScalarField& Field, const interfaceTopologySettings& s);

    // Region I
    tmp<volScalarField> calcRegionI
    (const volScalarField& BPsi, const volScalarField& band,
     const interfaceTopologySettings& s);

    tmp<volScalarField> calcRegionI
    (const volScalarField& BPsi, const volScalarField& Field,
     const interfaceTopologySettings& s, bool);

    // Region II
    tmp<volScalarField> calcRegionII
    (const volScalarField& BPsi, const volScalarField& band,
     const volScalarField& regionI, const interfaceTopologySettings& s);

    tmp<volScalarField> calcRegionII
    (const volScalarField& BPsi, const volScalarField& Field,
     const interfaceTopologySettings& s, bool);

    // Interface geometry
    interfaceProperties calcInterfaceGeometry
    (const volScalarField& regionI, const volScalarField& regionII,
     const interfaceTopologySettings& s);

    interfaceProperties calcInterfaceGeometry
    (const volScalarField& BPsi, const volScalarField& Field,
     const interfaceTopologySettings& s, bool);

    interfaceProperties calcInterface
    (const volScalarField& BPsi, const volScalarField& Field,
     const interfaceTopologySettings& s);

    // Interface bounds — field operations
    // Note: takes plain const ref not tmp
    interfaceBoundingCoords findFieldBounds
    (const volScalarField& regionI, const interfaceTopologySettings& s);

    interfaceBoundingCoords calcInterfaceBounds
    (const volScalarField& BPsi, const volScalarField& Field,
     const interfaceTopologySettings& s);

    // Interface bounds — forAll loop
    interfaceBoundingCoords loopCalcInterfaceBounds
    (const volScalarField& BPsi, const volScalarField& Field,
     const interfaceTopologySettings& s);
    
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfaceTopologyUtils
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
