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

using Foam::interfaceTopologyUtils::interfaceTopologyManager;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::interfaceTopologyUtils::interfaceTopologyManager::calcImpl
(
    const Foam::word& subdictName,
    interfaceTopologyManager::calcMode mode
)
{
    if (!dictionaryCheck(subdictName)) return;

    const interfaceTopologySettings& s = *settings_[subdictName];

    const volScalarField& field =
        mesh_.lookupObject<volScalarField>(s.fieldName());

    const volScalarField& BPsi =
        mesh_.lookupObject<volScalarField>(s.BPsiName());

    switch (mode)
    {
        case calcMode::interface:
            properties_[subdictName] = calcInterface(BPsi, field, s);
            break;

        case calcMode::interfaceBounds:
            coordinates_[subdictName] = calcInterfaceBounds(BPsi, field, s);
            break;
            
        case calcMode::interfaceBoundsViaLoop:
            coordinates_[subdictName] = loopCalcInterfaceBounds(BPsi, field, s);
            break;    
    }
}


bool
Foam::interfaceTopologyUtils::interfaceTopologyManager::dictionaryCheck
(
    const Foam::word& subdictName
)
{
    if (!settings_.found(subdictName))
    {
        WarningInFunction
            << "No settings for '" << subdictName << "'." << nl
            << "Registered entries: " << entryNames_ << nl
            << "Skipping."
            << endl;
        return false;
    }

    const interfaceTopologySettings& s = *settings_[subdictName];
    
    // Resolve field from mesh registry
    if (!mesh_.foundObject<volScalarField>(s.fieldName()))
    {
        WarningInFunction
            << "Field '" << s.fieldName()
            << "' not found in registry." << nl
            << "Skipping '" << subdictName << "'."
            << endl;
        return false;
    }

    // Resolve BPsi from mesh registry
    if (!mesh_.foundObject<volScalarField>(s.BPsiName()))
    {
        WarningInFunction
            << "BPsi field '" << s.BPsiName()
            << "' not found in registry." << nl
            << "Skipping '" << subdictName << "'."
            << endl;
        return false;
    }
    
    return true;
}


void
Foam::interfaceTopologyUtils::interfaceTopologyManager::calc
(
    const word& name,
    calcMode mode
)
{
    if (comboSettings_.found(name))
    {
        calcComboImpl(name, mode);
    }
    else
    {
        calcImpl(name, mode);
    }
}


void
Foam::interfaceTopologyUtils::interfaceTopologyManager::calcAll
(
    calcMode mode
)
{
    // Calculate all registered entries
    forAll(entryNames_, i)
    {
        calcImpl(entryNames_[i], mode);
    }

    // Calculate all registered combos
    forAll(comboNames_, i)
    {
        calcComboImpl(comboNames_[i], mode);
    }
}

// ************************************************************************* //
