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
    Foam::solidDisplacementThermoo

Description
    Fundamental solid thermodynamic properties

SourceFiles
    solidDisplacementThermoo.C

\*---------------------------------------------------------------------------*/

#ifndef solidDisplacementThermoo_H
#define solidDisplacementThermoo_H

#include "constSolidThermoo.H"
#include "PhysicalPropertiesThermo.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class solidDisplacementThermoo Declaration
\*---------------------------------------------------------------------------*/

class solidDisplacementThermoo
:public constSolidThermoo
{
    // Private data

        //- Switch to enable plane stress
        Switch planeStress_;

        //- Youngs modulus [Pa]
        volScalarField REf_;
   
 ///////////////////////////////////////////////////////////////////    
    //Bridgman properties
       //dimensioned scalars      
       dimensioned<scalar> dimx_;        
       dimensioned<scalar> dimt_;
       dimensioned<scalar> omega_; 
       dimensioned<scalar> gamma_;
       dimensioned<scalar> epsilon_; 
       dimensioned<scalar> flux_; 
       dimensioned<scalar> emissivity_;
       dimensioned<scalar> Cpp_;  
       dimensioned<scalar> L_;        
       dimensioned<scalar> Ks_;
       dimensioned<scalar> Kl_; 
       dimensioned<scalar> kc_;       
       dimensioned<scalar> Tm_;        
       dimensioned<scalar> Tl_;       
       dimensioned<scalar> Ts_;         
       dimensioned<scalar> Tinf_;         
       dimensioned<scalar> Bd_;         
       dimensioned<scalar> Vpull_;

       //DNM
       dimensioned<scalar> liquidD_;
       
       //scalars
       scalar ybottom_;
       scalar yRef1_;
       scalar yRef2_;
       scalar yRef3_;
       scalar yRef4_;
       scalar limT_;
        
 //////////////////////////////////////////////////////////////////          
public:

    //- Runtime type information
    TypeName("solidDisplacementThermoo");


    // Constructors

        //- Construct from mesh and phase name
        solidDisplacementThermoo
        (
            const fvMesh&,
            const word& phaseName = word::null
        );


    //- Destructor
    virtual ~solidDisplacementThermoo();


    // Member Functions

        //- Returns true to enable plane stress
        bool planeStress() const
        {
            return planeStress_;
        }

        // Mechanical state

            //- Youngs modulus [Pa]
            virtual const volScalarField& REf() const;

            //- Youngs modulus for a patch [Pa]
            virtual const scalarField& REf(const label patchi) const;
            
    // Bridgman
             
 //////////////////////////////////////////////////////////////////
 
       //dimensioned scalars      
       dimensioned<scalar> dimx() const
       {
           return dimx_;
       }
               
       dimensioned<scalar> dimt() const
       {
           return dimt_;
       }
       
       dimensioned<scalar> omega() const
       {
           return omega_;
       }
               
       dimensioned<scalar> gamma() const
       {
           return gamma_;
       }
       
       dimensioned<scalar> epsilon() const
       {
           return epsilon_;
       }
       
       dimensioned<scalar> flux() const
       {
           return flux_;
       }
       
       dimensioned<scalar> emissivity() const
       {
           return emissivity_;
       }
         
       dimensioned<scalar> Cpp() const
       {
           return Cpp_; 
        } 
        
       dimensioned<scalar> L() const
       {
           return L_;
        }
                
       dimensioned<scalar> Ks() const
       {
           return Ks_;
       }
       
       dimensioned<scalar> Kl() const
       {
           return Kl_;
       }
        
       dimensioned<scalar> kc() const
       {
           return kc_;
       }
              
       dimensioned<scalar> Tm() const
       {
           return Tm_;
       }
               
       dimensioned<scalar> Tl() const
       {
           return Tl_;
       }
              
       dimensioned<scalar> Ts() const
       {
           return Ts_;
       }         
       
       dimensioned<scalar> Tinf() const
       {
           return Tinf_;
       }
                
       dimensioned<scalar> Bd() const
       { 
           return Bd_;
       }         
        
       dimensioned<scalar> Vpull() const
       {
            return Vpull_;
       }    

       dimensioned<scalar> liquidD() const
       {
            return liquidD_;
       }   
       
       //scalars
        scalar ybottom() const
        {
            return ybottom_;
        }
        
        scalar yRef1() const
        {
            return yRef1_;
        }
        
        scalar yRef2() const
        {
            return yRef2_;
        }
        
        scalar yRef3() const
        {
            return yRef3_;
        }
        
        scalar yRef4() const
        {
            return yRef4_;
        }
        
          scalar limT() const
        {
            return limT_;
        }
            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
