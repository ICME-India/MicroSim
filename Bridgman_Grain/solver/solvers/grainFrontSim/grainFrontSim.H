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
    Foam::solvers::grainFrontSim

Description


SourceFiles
    grainFrontSim.C

\*---------------------------------------------------------------------------*/

#ifndef grainFrontSim_H
#define grainFrontSim_H

#include "Bridgman1.H"

#include "readPhysicalDict.H"
#include "clockTime.H"

using namespace Foam::mpmsc;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                         Class grainFrontSim Declaration

\*---------------------------------------------------------------------------*/

class grainFrontSim
:
    public Bridgman1
{
    private:

    using calcMode =
       Foam::interfaceTopologyUtils::interfaceTopologyManager::calcMode;

        mpmsc::GrainFieldConstants_t gfc;

        const label 
            NumPhases, 
            NumComps;            

       FieldArr 
            PhaseFields;
      

        FieldPtr 
            TempGrain,
            AllPhases,
            aniso,
            Velocity,
            SumPhases, 
            TempMod;


        clockTime timer;

        // Field3Ptr 
        //     GradPhase,
        //     PhiNCap;
                
        void createFields();
        void createInitialConditions();


protected:
  
        
    // Protected member functions

        //- Return true if solver dependencies have been modified
        virtual bool dependenciesModified() const;
        
        //- Read solver controls from fvSolution
        virtual bool read();
        

public:



    //- Runtime type information
    TypeName("grainFrontSim");


    // Constructors

        //- Construct from region mesh
        grainFrontSim(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        grainFrontSim(const grainFrontSim&) = delete;


    //- Destructor
    virtual ~grainFrontSim();
    
    
    // Member functions
    inline auto&& PHI(label a) noexcept
    {
        return *(PhaseFields[a]);
    }

    inline auto&& TEMP() noexcept
    {
        return *(TempGrain);
    }


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
        void operator=(const grainFrontSim&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
