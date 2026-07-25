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

#include "interfaceTopologyManager.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::interfaceTopologyUtils::interfaceTopologySettings&
Foam::interfaceTopologyUtils::interfaceTopologyManager::settings
(
    const word& subdictName
) const
{
    if (!settings_.found(subdictName))
    {
        FatalErrorInFunction
            << "No settings for '" << subdictName << "'" << nl
            << "Registered entries: " << entryNames_ << nl
            << abort(FatalError);
    }

    return *settings_[subdictName];
}


const Foam::interfaceTopologyUtils::interfaceProperties&
Foam::interfaceTopologyUtils::interfaceTopologyManager::properties
(
    const word& subdictName
) const
{
    if (!properties_.found(subdictName))
    {
        FatalErrorInFunction
            << "No properties for '" << subdictName << "'" << nl
            << "Registered entries: " << entryNames_ << nl
            << abort(FatalError);
    }

    return properties_[subdictName];
}


const Foam::interfaceTopologyUtils::interfaceBoundingCoords&
Foam::interfaceTopologyUtils::interfaceTopologyManager::coordinates
(
    const word& subdictName
) const
{
    if (!coordinates_.found(subdictName))
    {
        FatalErrorInFunction
            << "No coordinates for '" << subdictName << "'" << nl
            << "Registered entries    : " << entryNames_  << nl
            << "Registered combos     : " << comboNames_  << nl
            << "Computed coordinates  : " << coordinates_.toc() << nl
            << abort(FatalError);
    }

    return coordinates_[subdictName];
}

// ************************************************************************* //
