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


// * * * * * * * * * * * *  Compute field combination  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::interfaceTopologyUtils::interfaceTopologyManager::computeComboField
(
    const interfaceComboSettings& cs
) const
{
    const wordList& entries  = cs.factorEntries();
    const scalarList& w     = cs.weights();
    
    
    // Initialise with first factor field
    // Use tmp to avoid copy — ownership transferred through expressions
    tmp<volScalarField> tResult
    (
        new volScalarField
        (
            IOobject
            (
                cs.comboName(),
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_.lookupObject<volScalarField>(settings_[entries[0]]->fieldName())
        )
    );

    using op = interfaceComboSettings::operation;

    switch (cs.op())
    {
        case op::multiply:
        {
            // Pointwise product — weights not used for multiply
            for (label i = 1; i < entries.size(); i++)
            {
                tResult.ref() *=
                    mesh_.lookupObject<volScalarField>(settings_[entries[i]]->fieldName());
            }
            break;
        }

        case op::add:
        {
            // Weighted sum — apply first field weight then accumulate
            tResult.ref() *= w[0];
            for (label i = 1; i < entries.size(); i++)
            {
                tResult.ref() +=
                    w[i] * mesh_.lookupObject<volScalarField>
                    (settings_[entries[i]]->fieldName());
            }
            break;
        }

        case op::min:
        {
            // Pointwise minimum across all factor entries
            for (label i = 1; i < entries.size(); i++)
            {
                const volScalarField& fi =
                    mesh_.lookupObject<volScalarField>(settings_[entries[i]]->fieldName());

                tResult.ref().primitiveFieldRef() =
                    Foam::min
                    (
                        tResult().primitiveField(),
                        fi.primitiveField()
                    );
            }
            break;
        }

        case op::max:
        {
            // Pointwise maximum across all factor entries
            for (label i = 1; i < entries.size(); i++)
            {
                const volScalarField& fi =
                    mesh_.lookupObject<volScalarField>(settings_[entries[i]]->fieldName());

                tResult.ref().primitiveFieldRef() =
                    Foam::max
                    (
                        tResult().primitiveField(),
                        fi.primitiveField()
                    );
            }
            break;
        }
    }
    
    return tResult;
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceTopologyUtils::interfaceTopologyManager::computeBinarizedComboField
(
    const interfaceComboSettings& cs
) const
{
    const wordList& entries = cs.factorEntries();
    const scalarList& w    = cs.weights();

    // Lambda to get binarized field for a factorField
    // Uses that field's own settings (Tm, limT, interiorThreshold)
    // — each field binarized independently before combination
    auto getBinary = [&](const word& entryName) -> tmp<volScalarField>
    { 
        const volScalarField& field =
            mesh_.lookupObject<volScalarField>(settings_[entryName]->fieldName());

        const volScalarField& BPsi =
            mesh_.lookupObject<volScalarField>
            (
                settings_[entryName]->BPsiName()
            );

        const interfaceTopologySettings& s = *settings_[entryName];

        // calcInterfaceBand — 1 where field in [Tm-limT, Tm+limT], 0 elsewhere
        // calcRegionI       — masks band with BPsi interior region
        tmp<volScalarField> tBand   = calcInterfaceBand(field, s);
        tmp<volScalarField> tRegion = calcRegionI(BPsi, tBand(), s);
        tBand.clear();

        return tRegion;
    };

    // Initialise result with binarized first field
    tmp<volScalarField> tResult
    (
        new volScalarField
        (
            IOobject
            (
                cs.comboName() + "_binary",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            getBinary(entries[0])  // binarized first field
        )
    );

    using op = interfaceComboSettings::operation;

    switch (cs.op())
    {
        case op::multiply:
        {
            // AND-like operation — interface only where ALL fields
            // have interface simultaneously
            // result = binary1 * binary2 * ... (all must be 1)
            for (label i = 1; i < entries.size(); i++)
            {
                tmp<volScalarField> tBin = getBinary(entries[i]);
                tResult.ref() *= tBin();
                tBin.clear();
            }
            break;
        }

        case op::add:
        {
            // OR-like operation — interface where ANY field
            // has interface, weighted by contribution
            // result = w0*binary1 + w1*binary2 + ...
            tResult.ref() *= w[0];
            for (label i = 1; i < entries.size(); i++)
            {
                tmp<volScalarField> tBin = getBinary(entries[i]);
                tResult.ref() += w[i] * tBin();
                tBin.clear();
            }
            break;
        }

        case op::min:
        {
            // Conservative AND — takes minimum binary value
            // equivalent to multiply for pure 0/1 fields
            // but more robust for soft boundaries
            for (label i = 1; i < entries.size(); i++)
            {
                tmp<volScalarField> tBin = getBinary(entries[i]);
                tResult.ref().primitiveFieldRef() =
                    Foam::min
                    (
                        tResult().primitiveField(),
                        tBin().primitiveField()
                    );
                tBin.clear();
            }
            break;
        }

        case op::max:
        {
            // Liberal OR — takes maximum binary value
            // interface present if ANY field shows interface
            for (label i = 1; i < entries.size(); i++)
            {
                tmp<volScalarField> tBin = getBinary(entries[i]);
                tResult.ref().primitiveFieldRef() =
                    Foam::max
                    (
                        tResult().primitiveField(),
                        tBin().primitiveField()
                    );
                tBin.clear();
            }
            break;
        }
    }

    return tResult;
}


// * * * * * * * * * * *  Topology of field combination  * * * * * * * * * * //

void
Foam::interfaceTopologyUtils::interfaceTopologyManager::calcComboImpl
(
    const word& comboName,
    calcMode mode
)
{
    if (!comboSettings_.found(comboName))
    {
        Info<< "calcComboImpl: '" << comboName
            << "' NOT found in comboSettings_" << nl
            << "    comboSettings_ keys: " << comboSettings_.toc() << endl;

        FatalErrorInFunction
            << "No combo settings for '" << comboName << "'"
            << abort(FatalError);
    }

    const interfaceComboSettings& cs = *comboSettings_[comboName];

    if (!validateCombo(comboName))
    {
        Info<< "calcComboImpl: validateCombo FAILED for '"
            << comboName << "'" << endl;
        return;
    }

    if (!mesh_.foundObject<volScalarField>(cs.BPsiName()))
    {
        Info<< "calcComboImpl: BPsi '" << cs.BPsiName()
            << "' not found in registry" << endl;
        return;
    }

    const volScalarField& BPsi =
        mesh_.lookupObject<volScalarField>(cs.BPsiName());

    // In calcComboImpl
    const tmp<volScalarField> tCombo =
        cs.mode() == interfaceComboSettings::entryMode::binarized
      ? computeBinarizedComboField(cs)
      : computeComboField(cs);

    const interfaceTopologySettings& s = *settings_[comboName];
    
    switch (mode)
    {
        case calcMode::interface:
            properties_[comboName] = calcInterface(BPsi, tCombo(), s);
            break;

        case calcMode::interfaceBounds:
            coordinates_[comboName] = calcInterfaceBounds(BPsi, tCombo(), s);
            break;

        case calcMode::interfaceBoundsViaLoop:
            coordinates_[comboName] =
                loopCalcInterfaceBounds(BPsi, tCombo(), s);
            break;
    }
}

// ************************************************************************* //
