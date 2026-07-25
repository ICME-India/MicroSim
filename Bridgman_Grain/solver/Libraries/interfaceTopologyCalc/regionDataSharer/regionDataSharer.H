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
    
class
    Foam::interfaceTopologyUtils::regionDataSharer

Description
    Manages export and import of computed interface data between
    solver regions via the shared time object registry.

    Export — region computes data and writes to time registry
    Import — region reads data from time registry

    Both sides configured via fvSolution shareData subdict.
    Registry object created once at construction and updated
    in-place each step — no allocation per timestep.

\*---------------------------------------------------------------------------*/

#ifndef regionDataSharer_H
#define regionDataSharer_H

#include "IOdictionary.H"
#include "fvMesh.H"
#include "HashTable.H"
#include "autoPtr.H"            //////////
#include "interfaceTopologyManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace interfaceTopologyUtils
{

/*---------------------------------------------------------------------------*\
                   class regionDataSharer Declaration
\*---------------------------------------------------------------------------*/

class regionDataSharer
{

public:

    //- Supported shared data types
    enum class dataType
    {
        boundingCoords,     // Bounding Coords of interface
        scalar_,          
        vector_
    };
    
    //- Export entry - what is computed by current region 
    //  and shared with neighbouring region
    struct exportEntry
    {
        word sharedName;    // key in time registry
        dataType type;         
        word source;        // topoManager entry name     
        word component;     // for scalar_ type — which component
        interfaceTopologyManager::calcMode calcmode;    // how to compute
    };
    
    //- Import entry - what current region reads from
    //  neighbouring region
    struct importEntry
    {
        word     sharedName;   // key in time registry
        dataType type;
    };


private:    
    
    const fvMesh& mesh_;
    const Time&   runTime_;
    
    //- Export entries read from fvSolution
    List<exportEntry> exports_;
    
    //- Import entries read from fvSolution
    List<importEntry> imports_;
    
    //- Shared IOdictionary in time registry
    IOdictionary* sharedDictPtr_;
    
    //- Reference to topoManager - for export computation
    interfaceTopologyManager* topoManagerPtr_;
    
    //- Track which imports were unavailable last step
    HashTable<bool> retryPending_;
    
      
    // Private member functions
    
        //- Parse dataType from word
        static dataType parseDataType(const word& w)
        {
            if (w == "boundingCoords") return dataType::boundingCoords;
            if (w == "scalar")         return dataType::scalar_;
            if (w == "vector")         return dataType::vector_;
        
        
            FatalErrorInFunction
                << "Unknown shared data type '" << w << "'" << nl
                << "Valid types: boundingCoords scalar vector"
                << abort(FatalError);
                
            return dataType::scalar_;    
        }
        
        //- Parse calcMode from word
        static interfaceTopologyManager::calcMode
        parseCalcMode(const word& w)
        {
            using cm = interfaceTopologyManager::calcMode;

            if (w == "interface")             return cm::interface;
            if (w == "interfaceBounds")       return cm::interfaceBounds;
            if (w == "interfaceBoundsViaLoop")return cm::interfaceBoundsViaLoop;

            FatalErrorInFunction
                << "Unknown calcMode '" << w << "'" << nl
                << "Valid: interface interfaceBounds interfaceBoundsViaLoop"
                << abort(FatalError);

            return cm::interfaceBounds;
        }  
        
        //- Get or create IOdictionary in time registry
        IOdictionary& sharedDict()
        {
            const word dictName = "regionSharedData";

            //if (role_ == role::importer)
            //{
                // Importer never creates — only accesses
                //if (!runTime_.foundObject<IOdictionary>(dictName))
                /*
                {
                    FatalErrorInFunction
                        << "Shared registry '" << dictName
                        << "' not found in time registry." << nl
                        << "Exporter region must be constructed first." << nl
                        << "Check multiRegion solver construction order."
                        << abort(FatalError);
                }
                
                return runTime_.lookupObjectRef<IOdictionary>(dictName);
                */
            //}

            // Exporter — create if not exists
            if (!runTime_.foundObject<IOdictionary>(dictName))
            {
                IOdictionary* ptr = new IOdictionary
                (
                    IOobject
                    (
                        dictName,
                        runTime_.name(),
                        runTime_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    )
                );
                ptr->checkIn();
                sharedDictPtr_ = ptr;
            }
            else
            {
                sharedDictPtr_ =
                    &runTime_.lookupObjectRef<IOdictionary>(dictName);
            }

            return *sharedDictPtr_;
        }
        
        //- Initialise shared Dict entries with defaults
        void initialiseSharedEntries()
        {
            // Only exporter creates the registry
            // Importer waits until first import call
            //if (role_ != role::exporter) return;
            IOdictionary& dict = sharedDict();
        
            forAll(exports_, i)
            {
                const exportEntry& e = exports_[i];
                
                if (dict.found(e.sharedName)) continue;
                
                switch (e.type)
                {
                    case dataType::boundingCoords:
                    {
                        // Store as subdict — one entry per component
                        dictionary subDict;
                        subDict.add("formed", false);
                        subDict.add
                        (
                            "minBounds",   
                            vector( -Foam::great, -Foam::great, -Foam::great)
                        );
                        subDict.add
                        (
                            "maxBounds",   
                            vector( Foam::great, Foam::great, Foam::great)
                        );
                        dict.add(e.sharedName, subDict);
                        break;
                    }

                    case dataType::scalar_:
                        dict.add(e.sharedName, scalar(-Foam::great));
                        break;

                    case dataType::vector_:
                        dict.add(e.sharedName, vector::zero);
                        break;
                }
            }
            
            Info<< "regionDataSharer [" << mesh_.name()
                << "]: shared registry created with "
                << exports_.size() << " export entries" << endl;            
        }
        
public:

    // Constructor
    
        regionDataSharer
        (
            const fvMesh& mesh,
            interfaceTopologyManager& topoMananger
        );
         
         
        // Member functions — configuration

        //- Read export/import config from fvSolution shareData subdict
        void readConfiguration();
            
            
    // Member functions - runtime
    
        //- Compute and export all registered entries to shared registry
        //  Call after topoManager_.calcAll() in postSolve
        void exportAll();

        //- Read all imported entries from shared registry
        //  Returns dict of imported values — call in postSolve
        //  after exporting region has run
        HashTable<scalar>  importedScalars() const;

        //- Import a specific boundingCoords entry
        interfaceBoundingCoords importCoords
        (
            const word& sharedName,
            const bool retry = true
        ) const;

        //- Import a specific scalar entry
        scalar importScalar(const word& sharedName) const
        {
            const IOdictionary& dict =
                runTime_.lookupObject<IOdictionary>("regionSharedData");

            return dict.lookup<scalar>(sharedName);
        }

        //- Check if shared registry exists and has entry
        bool sharedEntryAvailable(const word& sharedName) const
        {
            if (!runTime_.foundObject<IOdictionary>("regionSharedData"))
                return false;

            return runTime_
                .lookupObject<IOdictionary>("regionSharedData")
                .found(sharedName);
        }


    // Private helpers

        static scalar extractComponent
        (
            const interfaceBoundingCoords& coords,
            const word& component
        )
        {
            if      (component == "minX") return coords.minBounds.x();
            else if (component == "minY") return coords.minBounds.y();
            else if (component == "minZ") return coords.minBounds.z();
            else if (component == "maxX") return coords.maxBounds.x();
            else if (component == "maxY") return coords.maxBounds.y();
            else if (component == "maxZ") return coords.maxBounds.z();
            else
            {
                FatalErrorInFunction
                    << "Unknown component '" << component << "'" << nl
                    << "Valid: minX minY minZ maxX maxY maxZ"
                    << abort(FatalError);
                return 0;
            }
        }  
        
        //- True if any imports are still pending from last step
        bool hasRetryPending() const
        {
            forAllConstIter(HashTable<bool>, retryPending_, iter)
            {
                if (iter()) return true;
            }
            return false;
        }               
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfaceTopologyUtils
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
