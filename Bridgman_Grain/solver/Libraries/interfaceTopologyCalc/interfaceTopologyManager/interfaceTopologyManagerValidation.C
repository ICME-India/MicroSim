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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool
Foam::interfaceTopologyUtils::interfaceTopologyManager::validateRegistry()
const
{
    bool allValid = true;

    Info<< "interfaceTopologyManager: validating mesh registry" << nl
        << "─────────────────────────────────────────────────────" << nl;

    forAll(entryNames_, i)
    {
        const word&                      name = entryNames_[i];
        const interfaceTopologySettings& s    = *settings_[name];

        // Source tag
        const word source =
            isManual_[name] ? "[manual]" : "[file]  ";

        Info<< "  Entry: " << name << " " << source << nl;

        // ── Check field

        const bool fieldFound =
            mesh_.foundObject<volScalarField>(s.fieldName());

        Info<< "    [" << (fieldFound ? "OK  " : "FAIL") << "] field '"
            << s.fieldName() << "' "
            << (fieldFound ? "found" : "NOT found") << nl;

        allValid = allValid && fieldFound;

        // If field found — check dimensions and Tm range
        if (fieldFound)
        {
            const volScalarField& f =
                mesh_.lookupObject<volScalarField>(s.fieldName());

            Info<< "    [INFO] dimensions : " << f.dimensions() << nl;

            const scalar fMin = gMin(f.primitiveField());
            const scalar fMax = gMax(f.primitiveField());

            if (s.Tm() < fMin || s.Tm() > fMax)
            {
                Info<< "    [WARN] Tm (" << s.Tm()
                    << ") outside field range ["
                    << fMin << ", " << fMax << "]"
                    << " — interface band may be empty" << nl;
            }
            else
            {
                Info<< "    [OK  ] Tm (" << s.Tm()
                    << ") within field range ["
                    << fMin << ", " << fMax << "]" << nl;
            }
        }

        // ── Check BPsi 

        const bool BPsiFound =
            mesh_.foundObject<volScalarField>(s.BPsiName());

        Info<< "    [" << (BPsiFound ? "OK  " : "FAIL") << "] BPsi  '"
            << s.BPsiName() << "' "
            << (BPsiFound ? "found" : "NOT found") << nl;

        allValid = allValid && BPsiFound;

        // If BPsi found — check range and thresholds
        if (BPsiFound)
        {
            const volScalarField& BPsi =
                mesh_.lookupObject<volScalarField>(s.BPsiName());

            const scalar BPsiMin = gMin(BPsi.primitiveField());
            const scalar BPsiMax = gMax(BPsi.primitiveField());

            Info<< "    [INFO] BPsi range : ["
                << BPsiMin << ", " << BPsiMax << "]" << nl;

            // Warn if BPsi not in [0,1]
            if (BPsiMin < -SMALL || BPsiMax > 1.0 + SMALL)
            {
                Info<< "    [WARN] BPsi outside [0,1]"
                    << " — check field is a phase fraction" << nl;
            }

            // Warn if thresholds outside BPsi range
            if (s.interiorThreshold() > BPsiMax)
            {
                Info<< "    [WARN] interiorThreshold ("
                    << s.interiorThreshold()
                    << ") > BPsiMax (" << BPsiMax
                    << ") — Region I will be empty" << nl;
                allValid = false;
            }

            if (s.boundaryThreshold() > BPsiMax)
            {
                Info<< "    [WARN] boundaryThreshold ("
                    << s.boundaryThreshold()
                    << ") > BPsiMax (" << BPsiMax
                    << ") — Region II will be empty" << nl;
                allValid = false;
            }

            // Warn if thresholds wrong way around
            if (s.boundaryThreshold() >= s.interiorThreshold())
            {
                Info<< "    [WARN] boundaryThreshold ("
                    << s.boundaryThreshold()
                    << ") >= interiorThreshold ("
                    << s.interiorThreshold()
                    << ") — Region II will be empty" << nl;
                allValid = false;
            }
        }

        Info<< "─────────────────────────────────────────────────────" << nl;
    }

    Info<< "interfaceTopologyManager: validation "
        << (allValid ? "OK" : "FAILED") << nl << endl;

    return allValid;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::interfaceTopologyUtils::interfaceTopologyManager::validateOrFatal()
const
{
    if (!validateRegistry())
    {
        FatalErrorInFunction
            << "interfaceTopologyManager: registry validation failed." << nl
            << "Check fieldName and BPsiName entries in" << nl
            << "constant/interfaceSegmentation"
            << " or manually added entries." << nl
            << abort(FatalError);
    }
}


// ************************************************************************* //
