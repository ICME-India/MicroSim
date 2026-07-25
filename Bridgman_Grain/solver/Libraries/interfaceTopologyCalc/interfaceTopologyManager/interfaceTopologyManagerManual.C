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

void
Foam::interfaceTopologyUtils::interfaceTopologyManager::add
(
    const word& subdictName,
    const dictionary& dict
)
{
    Info<< "interfaceTopologyManager: adding manual entry '"
        << subdictName << "'" << endl;

    // Manual entries use no defaults — all values in dict
    addEntry(subdictName, dict, dictionary::null, true);
}


void
Foam::interfaceTopologyUtils::interfaceTopologyManager::add
(
    const word& subdictName,
    const word& fieldName,
    const word& BPsiName,
    const scalar Tm,
    const scalar limT,
    const scalar idimx,
    const scalar ybottom,
    const scalar yRef4,
    const scalar interiorThreshold,
    const scalar boundaryThreshold,
    const scalar interfaceMinCells,
    const scalar maskTolerance,
    const scalar regionMaskEps
)
{
    dictionary dict;
    dict.add("fieldName",          fieldName);
    dict.add("BPsiName",           BPsiName);
    dict.add("Tm",                 Tm);
    dict.add("limT",               limT);
    dict.add("idimx",              idimx);
    dict.add("ybottom",            ybottom);
    dict.add("yRef4",              yRef4);
    dict.add("interiorThreshold",  interiorThreshold);
    dict.add("boundaryThreshold",  boundaryThreshold);
    dict.add("interfaceMinCells",  interfaceMinCells);
    dict.add("maskTolerance",      maskTolerance);
    dict.add("regionMaskEps",      regionMaskEps);

    add(subdictName, dict);
}


void
Foam::interfaceTopologyUtils::interfaceTopologyManager::update
(
    const word& subdictName,
    const dictionary& dict
)
{
    if (!settings_.found(subdictName))
    {
        FatalErrorInFunction
            << "Entry '" << subdictName << "' not found." << nl
            << "Registered entries: " << entryNames_ << nl
            << abort(FatalError);
    }

    if (!isManual_[subdictName])
    {
        WarningInFunction
            << "Entry '" << subdictName << "' is file-based." << nl
            << "Manual update applied but will be overwritten"
            << " on next file reread." << endl;
    }

    // Manual updates use no defaults
    settings_[subdictName]->update(dict, dictionary::null);
}


void
Foam::interfaceTopologyUtils::interfaceTopologyManager::remove
(
    const word& subdictName
)
{
    if (!settings_.found(subdictName))
    {
        FatalErrorInFunction
            << "Entry '" << subdictName << "' not found." << nl
            << "Registered entries: " << entryNames_ << nl
            << abort(FatalError);
    }

    Info<< "interfaceTopologyManager: removing entry '"
        << subdictName << "'" << endl;

    settings_.erase(subdictName);
    properties_.erase(subdictName);
    coordinates_.erase(subdictName);
    isManual_.erase(subdictName);
    removeFromEntryNames(subdictName);
}


bool
Foam::interfaceTopologyUtils::interfaceTopologyManager::isManual
(
    const word& subdictName
) const
{
    if (!isManual_.found(subdictName))
    {
        FatalErrorInFunction
            << "Entry '" << subdictName << "' not found." << nl
            << "Registered entries: " << entryNames_ << nl
            << abort(FatalError);
    }

    return isManual_[subdictName];
}


// ************************************************************************* //
