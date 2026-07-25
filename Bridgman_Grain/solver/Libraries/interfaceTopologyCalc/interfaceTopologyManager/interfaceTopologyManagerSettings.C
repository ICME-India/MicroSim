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

void
Foam::interfaceTopologyUtils::interfaceTopologyManager::removeFromEntryNames
(
    const word& subdictName
)
{
    // Rebuild preserving insertion order
    wordList newNames(entryNames_.size() - 1);
    label j = 0;

    forAll(entryNames_, i)
    {
        if (entryNames_[i] != subdictName)
            newNames[j++] = entryNames_[i];
    }

    entryNames_.transfer(newNames);
}


void
Foam::interfaceTopologyUtils::interfaceTopologyManager::addEntry
(
    const word& subdictName,
    const dictionary& entry,
    const dictionary& defaults,
    const bool manual
)
{
    if (settings_.found(subdictName))
    {
        WarningInFunction
            << "Entry '" << subdictName
            << "' already exists — overwriting." << endl;

        settings_.erase(subdictName);
        properties_.erase(subdictName);
        coordinates_.erase(subdictName);
        isManual_.erase(subdictName);
        removeFromEntryNames(subdictName);
    }

    // entry    = explicit subdict (T{}) or empty dict
    // defaults = default{} subdict
    // settings reads entry first then defaults for each key
    settings_.insert
    (
        subdictName,
        autoPtr<interfaceTopologySettings>
        (
            new interfaceTopologySettings(subdictName, entry, defaults)
        )
    );

    properties_.insert(subdictName, interfaceProperties());
    coordinates_.insert(subdictName, interfaceBoundingCoords());
    isManual_.insert(subdictName, manual);
    entryNames_.append(subdictName);
}


Foam::word
Foam::interfaceTopologyUtils::interfaceTopologyManager::resolveBPsiName
(
    const word& entryName,
    const dictionary& entry,
    const dictionary& defaults
) const
{
    // Priority 1 — explicit subdict BPsiName
    if (entry.found("BPsiName"))
    {
        const word name = entry.lookup<word>("BPsiName");

        Info<< "    [BPsiName] '" << entryName
            << "' → '" << name
            << "' (from explicit subdict)" << nl;

        return name;
    }

    // Priority 2 — BPsiNames{} mapping in defaults
    if (defaults.isDict("BPsiNames"))
    {
        const dictionary& BPsiNames = defaults.subDict("BPsiNames");

        if (BPsiNames.found(entryName))
        {
            const word name = BPsiNames.lookup<word>(entryName);

            Info<< "    [BPsiName] '" << entryName
                << "' → '" << name
                << "' (from default.BPsiNames)" << nl;

            return name;
        }
    }

    // Priority 3 — default BPsiName
    if (defaults.found("BPsiName"))
    {
        const word name = defaults.lookup<word>("BPsiName");

        Info<< "    [BPsiName] '" << entryName
            << "' → '" << name
            << "' (from default.BPsiName)" << nl;

        return name;
    }

    // Nothing found
    FatalErrorInFunction
        << "BPsiName not resolved for entry '" << entryName << "'." << nl
        << "Set one of:" << nl
        << "    1. BPsiName in field subdict T { BPsiName X; }" << nl
        << "    2. default { BPsiNames { " << entryName << " X; } }" << nl
        << "    3. default { BPsiName X; }" << nl
        << abort(FatalError);

    return word::null;
}


void
Foam::interfaceTopologyUtils::interfaceTopologyManager::readAllSettings()
{
    settings_.clear();
    properties_.clear();
    coordinates_.clear();
    entryNames_.clear();
    isManual_.clear();

    // ── Step 1: Read default subdict ─────────────────────────────────

    const dictionary defaults =
        isDict("default")
            ? subDict("default")
            : dictionary();

    if (!isDict("default"))
    {
        WarningInFunction
            << "No 'default' subdict found in interfaceSegmentation." << nl
            << "All values must be specified per entry." << endl;
    }
    else
    {
        Info<< "\ninterfaceTopologyManager: default subdict found" << nl;

        // Print BPsiName resolution info
        if (defaults.found("BPsiName"))
        {
            Info<< "    default BPsiName     : "
                << defaults.lookup<word>("BPsiName") << nl;
        }
        else
        {
            Info<< "    default BPsiName     : (not set)" << nl;
        }

        if (defaults.isDict("BPsiNames"))
        {
            Info<< "    per-field BPsiNames  : "
                << defaults.subDict("BPsiNames").toc() << nl;
        }
        else
        {
            Info<< "    per-field BPsiNames  : (not set)" << nl;
        }
    }

    // ── Step 2: Get field list from default ───────────────────────────

    wordList autoEntries;

    if (defaults.found("entries"))
    {
        autoEntries = defaults.lookup<wordList>("entries");
        Info<< "    entries               : " << autoEntries << nl << endl;
    }
    else
    {
        WarningInFunction
            << "No 'entries' list in default subdict." << nl
            << "No entries will be auto-registered." << nl
            << "Specify default { entries (T phi ...); }" << endl;
    }

    // ── Step 3: Register auto fields ─────────────────────────────────
    // For each field in default.fields:
    //   entry = explicit subdict if exists, otherwise empty dict
    //   defaults = default{} subdict
    //   settings reads entry first then defaults — no manual merging

    forAll(autoEntries, i)
    {
        const word& entryName = autoEntries[i];

        Info<< "\n    registering entry: " << entryName << nl;

        // Get explicit subdict if exists — empty dict otherwise
        // Do NOT manually merge — pass as-is to settings
        // settings.readFromDicts(entry, defaults) handles priority
        const dictionary& entry =
            isDict(entryName)
                ? subDict(entryName)
                : dictionary::null;

        if (isDict(entryName))
        {
            Info<< "    explicit subdict    : found — overrides apply" << nl;
        }
        else
        {
            Info<< "    explicit subdict    : not found — using defaults" << nl;
        }

        // Resolve BPsiName separately — needed for info output
        // settings also resolves it internally via readFromDicts
        const word BPsiName =
            resolveBPsiName(entryName, entry, defaults);

        Info<< "    resolved BPsiName   : " << BPsiName << nl;

        // false = from file
        addEntry(entryName, entry, defaults, false);
    }

    // ── Step 4: Register explicit subdicts not in auto list ───────────
    // Handles entries added explicitly without being in default.entries

    const wordList allEntries = toc();

    forAll(allEntries, i)
    {
        const word& name = allEntries[i];

        // Skip: non-subdicts, default subdict, already registered
        if
        (
            !isDict(name)
         || name == "default"
         || name == "combos"
         || name == "FoamFile"
         || settings_.found(name)
        )
        {
            Info<< "readAllSettings: skipping reserved key '"
            << name << "'" << endl;
            continue;
        }

        Info<< "\n    registering explicit entry: " << name << nl;

        // Pass explicit subdict directly — no merging
        const dictionary& entry = subDict(name);

        // Resolve and log BPsiName
        const word BPsiName =
            resolveBPsiName(name, entry, defaults);

        Info<< "    resolved BPsiName   : " << BPsiName << nl;

        // false = from file
        addEntry(name, entry, defaults, false);
    }

    Info<< "\ninterfaceTopologyManager: registered "
        << entryNames_.size() << " entries: "
        << entryNames_ << nl << endl;
}


// ************************************************************************* //
