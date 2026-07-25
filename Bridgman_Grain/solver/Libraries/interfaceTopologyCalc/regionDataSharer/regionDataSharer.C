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

#include "regionDataSharer.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceTopologyUtils::regionDataSharer::regionDataSharer
(
    const fvMesh& mesh,
    interfaceTopologyManager& topoMananger
)
:        
    mesh_(mesh),
    runTime_(mesh.time()),
    sharedDictPtr_(nullptr),
    topoManagerPtr_(&topoMananger)    
{
    readConfiguration();
    initialiseSharedEntries();
} 


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void 
Foam::interfaceTopologyUtils::regionDataSharer::readConfiguration()
{
    const dictionary& fvSolution =
        mesh_.lookupObject<IOdictionary>("fvSolution");
        
    if (!fvSolution.found("shareData"))
    {    
        Info<< "regionDataSharer [" << mesh_.name()
            << "]: no shareData in fvSolution — nothing shared"
            << endl;
            
        return;
    }
    

    const dictionary& shareDict =
         fvSolution.subDict("shareData");
         

    // ── Read exports ──────────────────────────────────────
    if (shareDict.isDict("export"))
    {
        const dictionary& exportDict =
            shareDict.subDict("export");

        const wordList exportKeys = exportDict.toc();

        forAll(exportKeys, i)
        {
            const word& name = exportKeys[i];
            const dictionary& entry = exportDict.subDict(name);

            exportEntry e;
            e.sharedName = name;
            e.type       = parseDataType
            (
                entry.lookup<word>("type")
            );
            e.source     = entry.lookup<word>("source");
            e.component  = entry.lookupOrDefault<word>
            (
                "component", word::null
            );
            e.calcmode   = parseCalcMode
            (
                entry.lookupOrDefault<word>
                (
                    "calcMode", "interfaceBounds"
                )
            );

                exports_.append(e);

                Info<< "regionDataSharer [" << mesh_.name()
                    << "]: export '" << name
                    << "' type=" << entry.lookup<word>("type")
                    << " source=" << e.source
                    << endl;
        }
    }

    // ── Read imports ──────────────────────────────────────
    if (shareDict.isDict("import"))
    {
        const dictionary& importDict =
            shareDict.subDict("import");

        const wordList importKeys = importDict.toc();

        forAll(importKeys, i)
        {
            const word& name = importKeys[i];
            const dictionary& entry = importDict.subDict(name);

            importEntry imp;
            imp.sharedName = name;
            imp.type       = parseDataType
            (
                entry.lookup<word>("type")
            );

            imports_.append(imp);

            Info<< "regionDataSharer [" << mesh_.name()
                << "]: import '" << name
                << "' type=" << entry.lookup<word>("type")
                << endl;
        }
    }
}
        
        
void
Foam::interfaceTopologyUtils::regionDataSharer::exportAll()
{
    if (exports_.empty()) return;

    IOdictionary& dict = sharedDict();
                
    forAll(exports_, i)
    {
        const exportEntry& e = exports_[i];

        // Ensure topoManager has computed this entry
        topoManagerPtr_->calc(e.source, e.calcmode);

        switch (e.type)
        {
            case dataType::boundingCoords:
            {
                const interfaceBoundingCoords& coords =
                    topoManagerPtr_->coordinates(e.source);

                // Update subdict
                dictionary& sub =
                    dict.subDict(e.sharedName);

                sub.set("formed", coords.formed);
                sub.set("minBounds",   coords.minBounds); 
                sub.set("maxBounds",   coords.maxBounds);


                break;
            }

            case dataType::scalar_:
            {
                // Extract named component from boundingCoords
                const interfaceBoundingCoords& coords =
                    topoManagerPtr_->coordinates(e.source);

                scalar val = extractComponent(coords, e.component);
                dict.set(e.sharedName, val);

                break;
            }

            case dataType::vector_:
            {
                // Export min vector 
                const interfaceBoundingCoords& coords =
                      topoManagerPtr_->coordinates(e.source);

                dict.set
                (
                    e.sharedName,
                    coords.minBounds
                );
                break;
            }
        }
    }
}  


Foam::HashTable<Foam::scalar>  
Foam::interfaceTopologyUtils::regionDataSharer::importedScalars() const
{
    HashTable<scalar> result;

    const IOdictionary& dict =
        runTime_.lookupObject<IOdictionary>("regionSharedData");

    forAll(imports_, i)
    {
        const importEntry& imp = imports_[i];

        if (imp.type == dataType::scalar_)
        {
            if (dict.found(imp.sharedName))
            {
                result.insert
                (
                    imp.sharedName,
                    dict.lookup<scalar>(imp.sharedName)
                );
            }
        }
    }

    return result;
}


Foam::interfaceTopologyUtils::interfaceBoundingCoords 
Foam::interfaceTopologyUtils::regionDataSharer::importCoords
(
    const word& sharedName,
    const bool retry
) const
{
    interfaceBoundingCoords result;

    if (!runTime_.foundObject<IOdictionary>("regionSharedData"))
    {
   
        // Not available yet — return default and retry next step
        WarningInFunction
            << "Shared registry not available yet for '"
            << sharedName << "' — will retry next timestep"
            << endl;

        const_cast<HashTable<bool>&>(retryPending_)
            .set(sharedName, true);     
            
            
        return interfaceBoundingCoords();
    }

    const IOdictionary& dict =
        runTime_.lookupObject<IOdictionary>("regionSharedData");

    if (!dict.isDict(sharedName))
    {
        if (retry)
        {
            WarningInFunction
                << "'" << sharedName
                << "' not yet in shared registry — "
                << "will retry next timestep"
                << endl;

            const_cast<HashTable<bool>&>(retryPending_)
                .set(sharedName, true);

            return interfaceBoundingCoords();
        }
    }
    
    // Clear retry flag — successfully imported
    const_cast<HashTable<bool>&>(retryPending_)
        .set(sharedName, false);

    const dictionary& sub = dict.subDict(sharedName);

    result.formed = sub.lookup<bool>  ("formed");
    result.minBounds   = sub.lookup<vector>("minBounds");
    result.maxBounds   = sub.lookup<vector>("maxBounds");


    return result;
}

// ************************************************************************* //
