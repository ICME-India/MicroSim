/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::solvers::solidDisplacement

Description
    Solver module for steady or transient segregated finite-volume solution of
    linear-elastic, small-strain deformation of a solid body, with optional
    thermal diffusion and thermal stresses.

    Solves for the displacement vector field D, also generating the stress
    tensor field sigma, including the thermal stress contribution if selected.

SourceFiles
    solidDisplacement.C

\*---------------------------------------------------------------------------*/

#ifndef Bridgman1_H
#define Bridgman1_H

#include "solver.H"
#include "fvcSurfaceIntegrate.H"
#include "fvMeshMover.H"
#include "volFieldValue.H"
#include "cellsToCells.H"
#include "solidDisplacementThermoo.H"
#include "interfaceTopology.H"
#include "regionDataSharer.H"
//#include "meshMotionControl.H"
#include "meshSearch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                         Class Bridgman1 Declaration

    Solid solver implementing the Bridgman crystal growth process.
    Moves the solid mesh at regular intervals to simulate crystal pull,
    and maps temperature from a neighbouring mesh (e.g. crucible/furnace)
    onto the solid domain via mesh-to-mesh intersection interpolation.
\*---------------------------------------------------------------------------*/

class Bridgman1
:
    public solver
{

public:

    using interfaceBoundingCoords =
        interfaceTopologyUtils::interfaceBoundingCoords;

protected:

    //- Settings manager for interface topology calculations
    interfaceTopologyUtils::interfaceTopologyManager topoManager_; 
    
    //- Data sharer with neighbour region
    interfaceTopologyUtils::regionDataSharer dataSharer_; 
 
    //- Mesh points predictor
    //autoPtr<meshMotionControl> meshPointsPredictor_; 

    // Mesh motion

        //- Current mesh point positions (updated each move step)
        pointField newPoints_;
        
        //- Predicted mesh point positions (updated each move step)
        pointField predPoints_;

        //- Simulation time at which the next mesh move is due
        scalar meshMoveTime_;

        //- Time interval between successive mesh move steps
        scalar meshMoveInterval_;
        
        //- Shift in Y direction to be applied to mesh   
        scalar shiftY_;

        //- Storing shift for nbrMesh 
        vector fShift_;

        
    // Neighbour mesh and field availability
    
        //- True if all required neighbour resources are available
        bool nbrResourcesAvailable_;


    // Neighbour mesh/field selection (read from fvSolution)

        //- Names of neighbour meshes to map from
        wordList nbrMeshNames_;

        //- Per-mesh list of field names to map
        //  e.g. nbrFieldNames_["crucible"] = {"T", "alpha1"}
        HashTable<wordList> nbrFieldNames_;

        //- Pointers to neighbour meshes (populated in constructor)
        HashTable<const fvMesh*> nbrMeshPtrs_;

        //- Pointers to neighbour fields, keyed by "meshName::fieldName"
        HashTable<const volScalarField*> nbrFieldPtrs_;

        //- Mapped field storage, keyed by "meshName::fieldName"
        HashTable<volScalarField*> mappedFieldPtrs_;

        //- Per-mesh interpolation objects
        mutable HashTable<autoPtr<cellsToCells>> interpolatorPtrs_;  
        
    
    // Initial mesh positioning

        //- Name of the neighbour mesh used as position reference
        word refMeshName_;

        //- Name of the field used to detect the reference interface
        word refFieldName_;

        //- Field threshold to detect the interface 
        scalar refFieldThreshold_;

        //- Absolute 3D position to place mesh centroid at (GREAT = not set)
        vector meshInitialPosition_;

        //- Additional 3D offset applied after positioning
        vector meshInitialOffset_;

        //- True once the initial mesh position has been applied
        bool meshPositioned_;      
     
        
    // Neighbour fields mapping control

        //- True if meshes moved or changed topology this step —
        //  triggers stencil rebuild in mappingUpdate()
        bool stencilOutdated_;

        //- True if neighbour fields should be mapped this step
        bool mapNbrFields_;

        //- Simulation time at which the next field mapping is due
        scalar nbrMapTime_;

        //- Time interval between neighbour field mapping operations.
        //  Negative value means map every step.
        scalar nbrMapInterval_;     

        //- True if mesh has to move
        bool moveDue_; 
        
        //- List of cells coupled from mesh to crucible
        labelList coarseCellAddressing_;

        //- Field to display processor ownership of cells 
        volScalarField sProcID_;
        
        //- Sample domain parameter  
        volScalarField sPhi_;
        
        //- Indicator field for remeshing
        volScalarField remeshRegion_;
        
    // Neighbour fields level set data  
      
        interfaceBoundingCoords liquidusTCoords_;  
        
        interfaceBoundingCoords solidusTCoords_;  
        
        
    // Protected member functions

        //- Return true if solver dependencies have been modified
        virtual bool dependenciesModified() const;
        
        //- Read solver controls from fvSolution
        virtual bool read();
        
        //- Check if all resources exist
        virtual bool checkNbrResources();
        
        //- Populate neighbour mesh/field pointers, interpolators and
        //  mapped field storage.
        void initialiseNbrResources();
        
        //- Align mesh bottom to the interface position detected in the
        //  reference neighbour field, or to a prescribed Y coordinate.
        //  sets meshPositioned_ = true on success.
        void applyInitialMeshPosition();
        
        //- Rebuild intersection stencil (if needed) and remap
        //  neighbour temperature onto mappedT
        virtual void mappingUpdate();
        

public:

    //- Runtime type information
    TypeName("Bridgman1");


    // Constructors

        //- Construct from region mesh
        Bridgman1(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        Bridgman1(const Bridgman1&) = delete;


    //- Destructor
    virtual ~Bridgman1();


    // Public Data

        //- Reference to Sample domain parameter 
        volScalarField& sPhi;
    
        //- Reference indicator field for remeshing
        volScalarField& remeshRegion;

    // Member functions

        //- Return the current maximum time-step for stable solution
        virtual scalar maxDeltaT() const;

        //- Synchronise newPoints_ if mesh topology changed;
        //  called before moveMesh() each time step
        virtual void preSolve();

        //- Translate the mesh by one pull increment if the move
        //  deadline has elapsed; rebuild geometry cache and remap
        //  neighbour fields onto the updated domain
        virtual void moveMesh();
        
        //- Corrections that follow mesh motion
        virtual void motionCorrector();

        //- Called at the beginning of each PIMPLE iteration
        virtual void prePredictor();
        
        //- Construct and optionally solve the momentum equation
        virtual void momentumPredictor();

        //- Solve the energy equation to obtain temperature
        virtual void thermophysicalPredictor();

        //- Pressure corrector (unused for solid; present for interface)
        virtual void pressureCorrector();

        //- Post-correction hook (unused for solid; present for interface)
        virtual void postCorrector();

        //- Called after the PIMPLE loop at the end of the time step
        virtual void postSolve();


    // Member operators

        //- Disallow default bitwise assignment
        void operator=(const Bridgman1&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*
struct LevelSetAxisState
{
    State  stateMax;   // state for max side
    State  stateMin;   // state for min side
    label  waitCount;
    scalar offsetMax;  // trigger offset from mesh max
    scalar offsetMin;  // trigger offset from mesh min
    scalar travel;     // target = trigger ± travel
    label  priority;   // per axis priority (lower = higher priority)
};

meshTracking
{
    // Default offsets — overridden per level set
    default
    {
        offsetMax   (GREAT 0.002 GREAT);  // only Y tracked by default
        offsetMin   (GREAT 0.002 GREAT);
        travel      (0     0.001 0    );
        minCells    2;
        waitSteps   5;
        maxStep     (GREAT 0.01 GREAT);
    }

    levelSets
    {
        liquidus
        {
            // inherits defaults — override selectively
            offsetMax   (GREAT 0.003 GREAT);
            offsetMin   (GREAT 0.003 GREAT);
            travel      (0     0.001 0    );
            priority    (0     0     0   );  // per axis — 0 = highest
            minCells    2;
        }

        solidus
        {
            offsetMax   (GREAT 0.005 GREAT);
            offsetMin   (GREAT 0.005 GREAT);
            travel      (0     0.001 0    );
            priority    (0     1     0   );  // lower priority than liquidus in Y
            minCells    2;
        }
    }
}

meshTracking
{
    primaryLevelSet     liquidus;   // ← only this triggers reposition

    levelSets
    {
        liquidus { ... }
        solidus  { ... }
    }
}
*/

