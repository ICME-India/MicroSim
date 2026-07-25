/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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
    Foam::bridgmanThermo

Description
    Fundamental bridgman thermodynamic properties

SourceFiles
    bridgmanThermo.C

\*---------------------------------------------------------------------------*/

#ifndef bridgmanThermo_H
#define bridgmanThermo_H

#include "constSolidThermoo.H"
#include "PhysicalPropertiesThermo.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class bridgmanThermo Declaration
\*---------------------------------------------------------------------------*/

class bridgmanThermo
:public constSolidThermoo
{
    // Private data
    
        //- Distance and time scales      
        dimensionedScalar dimx_;        
        dimensionedScalar dimt_;
       
        //- Solid-liquid interface smoothing parameters
        dimensionedScalar omega_; 
        dimensionedScalar gamma_;
        dimensionedScalar epsilon_; 
       
        //- Crucible heat transfer and motion parameters
        dimensionedScalar flux_; 
        dimensionedScalar Tinf_;         
        dimensionedScalar Bd_; 
        dimensionedScalar Vpull_;       
        scalar ybottom_;
       
        //- CMSX4 thermophysical properties
        dimensionedScalar Cpp_;  
        dimensionedScalar L_;        
        dimensionedScalar Ks_;
        dimensionedScalar Kl_; 
        dimensionedScalar kc_;       
        dimensionedScalar Tm_;        
        dimensionedScalar Tl_;       
        dimensionedScalar Ts_;         
    
        //- List of sine series coefficients 
        //  for furnace temperature interpolation        
        List<FixedList<scalar, 3>> coeffsTinf2_;        
   
public:

    //- Runtime type information
    TypeName("bridgmanThermo");


    // Constructors

        //- Construct from mesh and phase name
        bridgmanThermo
        (
            const fvMesh&,
            const word& phaseName = word::null
        );


    //- Destructor
    virtual ~bridgmanThermo();


    // Member Functions
   
       //dimensioned scalars      
       dimensionedScalar dimx() const
       {
           return dimx_;
       }
               
       dimensionedScalar dimt() const
       {
           return dimt_;
       }
       
       dimensionedScalar omega() const
       {
           return omega_;
       }
               
       dimensionedScalar gamma() const
       {
           return gamma_;
       }
       
       dimensionedScalar epsilon() const
       {
           return epsilon_;
       }
       
       dimensionedScalar flux() const
       {
           return flux_;
       }
             
       dimensionedScalar Cpp() const
       {
           return Cpp_; 
       } 
        
       dimensionedScalar L() const
       {
           return L_;
        }
                
       dimensionedScalar Ks() const
       {
           return Ks_;
       }
       
       dimensionedScalar Kl() const
       {
           return Kl_;
       }
        
       dimensionedScalar kc() const
       {
           return kc_;
       }
              
       dimensionedScalar Tm() const
       {
           return Tm_;
       }
               
       dimensionedScalar Tl() const
       {
           return Tl_;
       }
              
       dimensionedScalar Ts() const
       {
           return Ts_;
       }         
       
       dimensionedScalar Tinf() const
       {
           return Tinf_;
       }
                
       dimensionedScalar Bd() const
       { 
           return Bd_;
       }    
       
       List<FixedList<scalar, 3>> coeffsTinf2() const
       {
           return coeffsTinf2_;
       }
        
       dimensionedScalar Vpull() const
       {
            return Vpull_;
       }       
       
       scalar ybottom() const
       {
            return ybottom_;
       }
            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
