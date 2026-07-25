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

Class
    Foam::interfaceTopologyUtils::interfaceTopologyManager

Description   
    A helper class for creating dictionary for each field  
     
    Supports three construction modes:

    Mode 1 — File based:
        interfaceTopologyManager manager(mesh);
        Reads constant/interfaceSegmentationDict automatically.
        MUST_READ_IF_MODIFIED — auto-reread on file change.

    Mode 2 — Manual only:
        interfaceTopologyManager manager(mesh, false);
        No file read. Add entries manually via add().

    Mode 3 — Mixed:
        interfaceTopologyManager manager(mesh);
        manager.add("extra", dict);
        Combines file entries with manual entries.

    Manual entries are NOT auto-reread — update() must be
    called explicitly if needed.
    
    Supports following member function calls:
    
    // Verbose way 
    manager.calcInterfaceExtremaAll();

    // Simple way
    manager.calcAll(calcMode::interfaceExtrema);

    // Single subdict
    manager.calc("myField", calcMode::interfaceTop);
    
Usage
    Example usage:
    \verbatim     
    default
    {
        // Required entries
        fields      (T Phi);
        
        BPsiName            BPsi;    // Default domain parameter
    
        BPsiNames   { T BPsi1; }     // Selecting domain parameter for a field

        Tm                  1681;
        limT                5.41;    // Interface band [Tm - limT, Tm + limT]

        idimx               2;       //- Inverse length scale [1/m]

        ybottom             0.01231; //- Bottom reference Y coordinate [m]

        yRef4               0.00415; //- Y reference offset [m]

        // Optional entries 
        interiorThreshold   0.91;    //- BPsi interior region threshold

        boundaryThreshold   0.055;   //- BPsi boundary region threshold
  
        interfaceMinCells   2;       //- Minimum cells required to 
                                     //  declare interface formed  
     }

     T                               //- Default values over-written
     {
         Tm                  1683; 
         yRef4               0.00412;      
     }   
    \endverbatim    

SourceFiles
    interfaceTopologyManagerIO.C
    interfaceTopologyManagerSettings.C
    interfaceTopologyManagerManual.C
    interfaceTopologyManagerValidation.C
    interfaceTopologyManagerAccess.C
    interfaceTopologyManagerCalc.C
    interfaceTopologyManagerOutput.C
    
See also
    Foam::interfaceTopologyUtils::interfaceTopologySettings    

\*---------------------------------------------------------------------------*/

#ifndef interfaceTopologyManager_H
#define interfaceTopologyManager_H

#include "IOdictionary.H"
#include "HashTable.H"
#include "wordList.H"
#include "autoPtr.H"
#include "volFields.H"
#include "interfaceProperties.H"
#include "interfaceTopologySettings.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace interfaceTopologyUtils
{

/*---------------------------------------------------------------------------*\
              Class interfaceTopologyManager Declaration
\*---------------------------------------------------------------------------*/

class interfaceTopologyManager
:
    public IOdictionary
{

public:

    //- Calculation mode
    enum class calcMode
    {
        interface,
        interfaceBounds,
        interfaceBoundsViaLoop
    };
    
    
private: 

        const fvMesh& mesh_;

        //- Settings per subdict name
        HashTable<autoPtr<interfaceTopologySettings>> settings_;

        //- Computed properties per subdict name
        HashTable<interfaceProperties> properties_;
        
        //- Computed interface bounding box coordinates
        //  per subdict name
        HashTable<interfaceBoundingCoords> coordinates_;

        //- Subdict names in insertion order
        wordList entryNames_;

        //- Track which entries are manual (not from file)
        //  Manual entries not updated on file reread
        HashTable<bool> isManual_;
        
        // Field combo support
        
            //- Combo settings keyed by combo name
            HashTable<autoPtr<interfaceComboSettings>> comboSettings_;

            //- Combo names in insertion order
            wordList comboNames_;


    // Private member functions
    
        //- Functions for a field

            // Remove entry from entryNames_ preserving order
            void removeFromEntryNames(const word& subdictName);
        
            // Internal add — used by both file and manual paths
            void addEntry
            (
                const word& subdictName,
                const dictionary& entry,
                const dictionary& defaults,   
                const bool manual
            );
        
            // Resolve BPsiName for a field using three-level priority:
            //  1. Explicit subdict BPsiName
            //  2. BPsiNames{} mapping in defaults
            //  3. default BPsiName
            //  FatalError if none found
            word resolveBPsiName
            (
                const word& entryName,
                const dictionary& fieldEntry,
                const dictionary& defaults
            ) const;
        
            // Read all subdicts from file
            void readAllSettings();

            // Validate all fields in mesh registry
            bool validateRegistry() const;
        
            // Validate all entries in subdict 
            bool dictionaryCheck(const word& subdictName);
        
            // Common lookup and dispatch
            void calcImpl(const word& subdictName, calcMode calcmode); 
            
        //- Functions for field combinations 
           
            // Read all combos from 'combos' subdict
            void readAllCombos();

            // Validate all factorFields of a combo exist in 
            //  settings_ and mesh
            bool validateCombo(const word& comboName) const;

            // Compute combined field from factorFields at runtime
            //  Returns tmp — no permanent storage
            tmp<volScalarField> computeComboField
            (
                const interfaceComboSettings& cs
            ) const;
            
            // Compute combined field from binarized factorFields 
            // at runtime - Returns tmp 
            tmp<volScalarField> computeBinarizedComboField
            (
                const interfaceComboSettings& cs
            ) const;

            //- Dispatch calc for a combo entry
            void calcComboImpl
            (
                const word& comboName,
                calcMode mode
            );
    

public:

    // Constructors

        //- Construct from mesh — reads interfaceSegmentationDict
        //  readFile=true  → read constant/interfaceSegmentationDict
        //  readFile=false → no file read, add entries manually
        explicit interfaceTopologyManager
        (
            const fvMesh& mesh,
            const bool readFile = true
        );


    // Member functions — IO

        //- Override readData — called automatically on file change
        virtual bool readData(Istream& is);


    // Member functions — manual entry management

        //- Add entry from dictionary
        //  subdictName — key to identify this entry
        //  dict        — dictionary containing all settings
        //                must contain fieldName, BPsiName, Tm etc
        void add( const word& subdictName, const dictionary& dict);

        //- Add entry by building dict from individual parameters
        void add
        (
            const word& subdictName,
            const word& fieldName,
            const word& BPsiName,
            const scalar Tm,
            const scalar limT,
            const scalar idimx,
            const scalar ybottom,
            const scalar yRef4,
            const scalar interiorThreshold = 0.96,
            const scalar boundaryThreshold = 0.05,
            const scalar interfaceMinCells = 5,
            const scalar maskTolerance     = 2.0,
            const scalar regionMaskEps     = 0.001
        );

        //- Update existing manual entry from new dictionary
        //  FatalError if entry not found or not manual
        void update
        (
            const word& subdictName,
            const dictionary& dict
        );

        //- Remove entry
        //  FatalError if entry not found
        void remove(const word& subdictName);

        //- Return true if entry is manual
        bool isManual(const word& subdictName) const;
        
        
    // combo management
    
        void addCombo
        (
            const word& comboName,
            const dictionary& dict
        );

        void removeCombo(const word& comboName);

 
    // Member functions — validation

        //- Validate and abort if any field missing
        void validateOrFatal() const;


    // Member functions — access
    
        //- Names of entries and settings access

           const wordList& entryNames() const 
            { 
                return entryNames_;
            }

            bool found(const word& subdictName) const
            {
                return settings_.found(subdictName);
            }

            const interfaceTopologySettings& settings
            (
                const word& subdictName
            ) const;

            const interfaceProperties& properties
            (
                const word& subdictName
            ) const;
        
            const interfaceBoundingCoords& coordinates
            (
                const word& subdictName
            ) const;

            const HashTable<interfaceProperties>& allProperties() const
            {
                return properties_;
            }
        
            const HashTable<interfaceBoundingCoords>& allCoordinates() const
            {
                return coordinates_;
            }
        
        //- Field combos name and setting access
        
            bool isCombo(const word& name) const
            {
                return comboSettings_.found(name);
            }

            const interfaceComboSettings& comboSettings
            (
                const word& comboName
            ) const;


    // Member functions — computation
    

        //- Loop over all entryNames_ with given mode
        void calcAll(calcMode mode);
        
        //- Calculate interface property for a given 
        //  fiels _with given mode
        void calc(const word& subdictName, calcMode mode);
        
    // Member functions — output

        void printSettings()              const;
        void printProperties()            const;
        void writeProperties(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfaceTopologyUtils
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
