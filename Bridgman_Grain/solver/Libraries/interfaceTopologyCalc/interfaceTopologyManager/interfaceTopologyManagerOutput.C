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
Foam::interfaceTopologyUtils::interfaceTopologyManager::printSettings() const
{
    Info<< " All registered settings:" << nl
        << "─────────────────────────────────────────" << nl;

    forAll(entryNames_, i)
    {
        const word& name = entryNames_[i];
        Info<< (isManual_[name] ? "[manual] " : "[file]   ");
        settings_[name]->print();
    }
}


void
Foam::interfaceTopologyUtils::interfaceTopologyManager::printProperties()
const
{
    Info<< "Interface properties summary:" << nl
        << "─────────────────────────────────────────" << nl;

    forAll(entryNames_, i)
    {
        const word&                      name = entryNames_[i];
        const interfaceTopologySettings& s    = *settings_[name];
        const interfaceProperties&       p    = properties_[name];

        Info<< "Entry     : " << name
            << (isManual_[name] ? " [manual]" : " [file]") << nl
            << "Field     : " << s.fieldName()              << nl
            << "BPsi      : " << s.BPsiName()               << nl;

        if (p.formed)
        {
            Info<< "    Height    : " << p.height    << nl
                << "    Y         : " << p.yCoord    << nl
                << "    Top       : " << p.yTop      << nl
                << "    Concavity : " << p.concavity << nl;
        }
        else
        {
            Info<< "    Interface did not form" << nl;
        }

        Info<< "─────────────────────────────────────────" << nl;
    }

    Info<< endl;
}


void
Foam::interfaceTopologyUtils::interfaceTopologyManager::writeProperties
(
    Ostream& os
) const
{
    forAll(entryNames_, i)
    {
        const word&                      name = entryNames_[i];
        const interfaceTopologySettings& s    = *settings_[name];
        const interfaceProperties&       p    = properties_[name];

        os  << name
            << (isManual_[name] ? "\tmanual" : "\tfile")
            << '\t' << s.fieldName()
            << '\t' << s.BPsiName()
            << '\t' << p.formed
            << '\t' << p.height
            << '\t' << p.yCoord
            << '\t' << p.yTop
            << '\t' << p.concavity
            << nl;
    }
}


// ************************************************************************* //
