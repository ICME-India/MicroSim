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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceTopologyUtils::interfaceTopologyManager::
interfaceTopologyManager
(
    const fvMesh& mesh,
    const bool readFile
)
:
    IOdictionary
    (
        IOobject
        (
            "interfaceSegmentation",
            mesh.time().constant(),
            mesh,
            readFile
                ? IOobject::MUST_READ_IF_MODIFIED
                : IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    settings_(),
    properties_(),
    entryNames_(),
    isManual_()
{
    if (readFile)
    {
        Info<< "\ninterfaceTopologyManager: reading "
            << "interfaceSegmentation" << endl;

        readAllSettings();
        readAllCombos(); 
    }
    else
    {
        Info<< "\ninterfaceTopologyManager: manual mode"
            << " — no file read" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::interfaceTopologyUtils::interfaceTopologyManager::readData
(
    Istream& is
)
{
    bool ok = IOdictionary::readData(is);

    if (ok)
    {
        Info<< "\ninterfaceTopologyManager: file changed"
            << " — updating file-based settings" << endl;

        // Re-read defaults from updated file
        const dictionary defaults =
            isDict("default")
                ? subDict("default")
                : dictionary();

        forAll(entryNames_, i)
        {
            const word& name = entryNames_[i];

            // Skip manual entries
            if (isManual_[name])
            {
                Info<< "    skipping manual entry: " << name << nl;
                continue;
            }

            // Get explicit subdict if exists — empty dict otherwise
            // No merging — pass directly to settings.update()
            const dictionary& entry =
                isDict(name)
                    ? subDict(name)
                    : dictionary::null;

            settings_[name]->update(entry, defaults);

            Info<< "    updated: " << name << nl;
        }
    }

    return ok;
}


// ************************************************************************* //
