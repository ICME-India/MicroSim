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
//constSolidThermoTemplatess.C
#include "solidDisplacementThermoo.H"

/* * * * * * * * * * * * * * * Private Static Data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(solidDisplacementThermoo, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidDisplacementThermoo::solidDisplacementThermoo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    constSolidThermoo(mesh, phaseName),
    planeStress_(lookup("planeStress")),
    REf_(readProperty<scalar>("REf", dimPressure)),
    dimx_(dimensioned<scalar>(dimensionSet(0, 1, 0, 0, 0, 0, 0),lookup<scalar>("dimx"))), 
    dimt_(dimensioned<scalar>(dimensionSet(0, 0, 1, 0, 0, 0, 0),lookup<scalar>("dimt"))), 
    omega_(dimensioned<scalar>(dimensionSet(0, 0, 0, 0, 0, 0, 0),lookup<scalar>("omega"))), 
    gamma_(dimensioned<scalar>(dimensionSet(0, 0, 0, 0, 0, 0, 0),lookup<scalar>("gamma"))), 
    epsilon_(dimensioned<scalar>(dimensionSet(0, 0, 0, 0, 0, 0, 0),lookup<scalar>("epsilon"))), 
    flux_(dimensioned<scalar>(dimensionSet(0, 0, 0, 0, 0, 0, 0),lookup<scalar>("flux"))),
    emissivity_
        (
            dimensioned<scalar>(dimensionSet(0, 0, 0, 0, 0, 0, 0),
            lookup<scalar>("emissivity"))
        ),
    Cpp_(dimensioned<scalar>(dimensionSet(0, 0, 0, 0, 0, 0, 0),lookup<scalar>("Cpp"))), 
    L_(dimensioned<scalar>(dimensionSet(0, 0, 0, 0, 0, 0, 0),lookup<scalar>("L"))),
    Ks_(dimensioned<scalar>(dimensionSet(0, 0, 0, 0, 0, 0, 0),lookup<scalar>("Ks"))), 
    Kl_(dimensioned<scalar>(dimensionSet(0, 0, 0, 0, 0, 0, 0),lookup<scalar>("Kl"))), 
    kc_(dimensioned<scalar>(dimensionSet(0, 0, 0, 0, 0, 0, 0),lookup<scalar>("kc"))), 
    Tm_(dimensioned<scalar>(dimensionSet(0, 0, 0, 0, 0, 0, 0),lookup<scalar>("Tm"))), 
    Tl_(dimensioned<scalar>(dimensionSet(0, 0, 0, 0, 0, 0, 0),lookup<scalar>("Tl"))),
    Ts_(dimensioned<scalar>(dimensionSet(0, 0, 0, 0, 0, 0, 0),lookup<scalar>("Ts"))), 
    Tinf_(dimensioned<scalar>(dimensionSet(0, 0, 0, 0, 0, 0, 0),lookup<scalar>("Tinf"))), 
    Bd_(dimensioned<scalar>(dimensionSet(0, 0, 0, 0, 0, 0, 0),lookup<scalar>("Bd"))), 
    Vpull_(dimensioned<scalar>(dimensionSet(0, 0, 0, 0, 0, 0, 0),lookup<scalar>("Vpull"))), 
    liquidD_(dimensioned<scalar>(dimensionSet(0, 0, 0, 0, 0, 0, 0),lookup<scalar>("liquidD"))),
    ybottom_(lookup<scalar>("ybottom")),
    yRef1_(lookup<scalar>("yRef1")),
    yRef2_(lookup<scalar>("yRef2")),
    yRef3_(lookup<scalar>("yRef3")),
    yRef4_(lookup<scalar>("yRef4")),
    limT_(lookup<scalar>("limT"))     
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidDisplacementThermoo::~solidDisplacementThermoo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volScalarField& Foam::solidDisplacementThermoo::REf() const
{
    return REf_;
}


const Foam::scalarField& Foam::solidDisplacementThermoo::REf
(
    const label patchi
) const
{
    return REf_.boundaryField()[patchi];
}

//Bridgman

















// ************************************************************************* //
