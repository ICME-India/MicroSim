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
#include "interfaceTopology.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::interfaceTopologyUtils::interfaceTopologyManager::readAllCombos()
{
  
    if (!isDict("combos"))
    {
        Info<< "readAllCombos: no 'combos' subdict found" << endl;
        return;
    }
    
    const dictionary& combosDict = subDict("combos");

    Info<< "readAllCombos: combosDict.toc = " << combosDict.toc() << endl;

    const dictionary& defaults =
        found("default") ? subDict("default") : dictionary::null;

    const wordList comboKeys = combosDict.toc();

    forAll(comboKeys, i)
    {
        const word& name = comboKeys[i];

        Info<< "readAllCombos: processing key '" << name << "'" << endl;

        if (!combosDict.isDict(name))
        {
            WarningInFunction
                << "combos entry '" << name
                << "' is not a subdict — skipping"
                << endl;
            continue;
        }

        comboSettings_.insert
        (
            name,
            autoPtr<interfaceComboSettings>
            (
                new interfaceComboSettings
                (
                    name,
                    combosDict.subDict(name),
                    defaults
                )
            )
        );
        
        coordinates_.insert(name, interfaceBoundingCoords());

        comboNames_.append(name);

        Info<< "readAllCombos: registered combo '" << name
            << "' comboSettings_ now has: " << comboSettings_.toc()
            << endl;
    }

    Info<< "readAllCombos: done. comboNames_ = " << comboNames_ << endl;
}


bool
Foam::interfaceTopologyUtils::interfaceTopologyManager::validateCombo
(
    const word& comboName
) const
{
    const interfaceComboSettings& cs = *comboSettings_[comboName];

    bool valid = true;

    forAll(cs.factorEntries(), i)
    {
        const word& fname = settings_[cs.factorEntries()[i]]->fieldName();

        // Must be registered in settings_
        if (!settings_.found(cs.factorEntries()[i]))
        {
            WarningInFunction
                << "factorEntry '" << fname
                << "' in combo '" << comboName
                << "' has no settings entry. "
                << "Add it to the entries list in the dictionary."
                << endl;
            valid = false;
        }

        // Must exist in mesh registry
        if (!mesh_.foundObject<volScalarField>(fname))
        {
            WarningInFunction
                << "Field in factorEntry '" << fname
                << "' in combo '" << comboName
                << "' not found in mesh registry."
                << endl;
            valid = false;
        }
    }

    return valid;
}

// ************************************************************************* //
