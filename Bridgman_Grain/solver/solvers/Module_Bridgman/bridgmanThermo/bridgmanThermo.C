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

\*---------------------------------------------------------------------------*/

#include "bridgmanThermo.H"

/* * * * * * * * * * * * * * * Private Static Data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(bridgmanThermo, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bridgmanThermo::bridgmanThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    constSolidThermoo(mesh, phaseName),
     
    dimx_("dimx",dimLength,lookupOrDefault<scalar>("dimx",1)), 
    dimt_("dimt",dimTime,lookupOrDefault<scalar>("dimt",1)),   
     
    omega_
    (
    	"omega",
    	dimless,
    	subDict("interfaceSL").lookupOrDefault<scalar>("omega",1e5)
    ),       
    gamma_
    (
    	"gamma",
    	dimless,
    	subDict("interfaceSL").lookupOrDefault<scalar>("gamma",2)
    ),    
    epsilon_
    (
    	"epsilon",
    	dimless,
    	subDict("interfaceSL").lookupOrDefault<scalar>("epsilon",200e-6)
    ),  
    
    flux_
    (
    	"flux",
    	dimless,
    	subDict("crucible").lookup<scalar>("flux")
    ), 
    Tinf_
    (
    	"Tinf",
    	dimless,
    	subDict("crucible").lookupOrDefault<scalar>("Tinf",300)
    ),  
    Bd_
    (
    	"Bd",
    	dimless,
    	subDict("crucible").lookupOrDefault<scalar>("Bd",1700)
    ),   
    Vpull_
    (
    	"Vpull",
    	dimless,
    	subDict("crucible").lookupOrDefault<scalar>("Vpull",5e-5)
    ), 
    ybottom_(subDict("crucible").lookupOrDefault<scalar>("ybottom",0)),  
      
    Cpp_
    (
    	"Cpp",
     	dimless,
     	subDict("CMSX4").lookupOrDefault<scalar>("Cpp",5471550)
    ),     
    L_
    (
    	"L",
    	dimless,
    	subDict("CMSX4").lookupOrDefault<scalar>("L",1945440000)
    ),    
    Ks_
    (
    	"Ks",
    	dimless,
    	subDict("CMSX4").lookupOrDefault<scalar>("Ks",90.7)
    ),     
    Kl_
    (
    	"Kl",
    	dimless,
    	subDict("CMSX4").lookupOrDefault<scalar>("Kl",18.14)
    ),     
    kc_
    (
    	"kc",
    	dimless,
    	subDict("CMSX4").lookupOrDefault<scalar>("kc",1.31)
    ),     
    Tm_
    (
    	"Tm",
    	dimless,
    	subDict("CMSX4").lookupOrDefault<scalar>("Tm",1623)
    ),     
    Tl_
    (
    	"Tl",
    	dimless,
    	subDict("CMSX4").lookupOrDefault<scalar>("Tl",1653)
    ),    
    Ts_
    (
    	"Ts",
    	dimless,
    	subDict("CMSX4").lookupOrDefault<scalar>("Ts",1593)
    ),    

    
    coeffsTinf2_
    (    
	lookupOrDefault<List<FixedList<scalar, 3>>>
        (
		"coeffsT",          
		List<FixedList<scalar, 3>>
		{
			{2643.90000, 9.18290, -0.28100},
			{1128.50000, 18.2673, 1.213929},
			{302.534000, 34.7846, 1.293900},
			{153.779100, 51.4367, 1.292700},
			{142.741700, 69.1533, 1.058700},
			{112.999100, 78.0118, 2.525000},
			{1.785430e4, 97.9142, 1.793700},
			{1.782715e4, 97.9285, 4.932500}
		}
	)   
    )  
   
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bridgmanThermo::~bridgmanThermo()
{}

// ************************************************************************* //
