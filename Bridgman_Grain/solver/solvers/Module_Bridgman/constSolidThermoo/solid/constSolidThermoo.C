/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "constSolidThermoo.H"
#include "addToRunTimeSelectionTable.H"

/* * * * * * * * * * * * * * * Private Static Data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(constSolidThermoo, 0);
    addToRunTimeSelectionTable(basicThermo, constSolidThermoo, fvMesh);
    addToRunTimeSelectionTable(solidThermo, constSolidThermoo, fvMesh);
}


// * * * * * * * * * * * * * * Protected Constructors  * * * * * * * * * * * //

Foam::constSolidThermoo::constSolidThermoo
(
    const fvMesh& mesh,
    const bool readKappa,
    const word& phaseName
)
:
    PhysicalPropertiesThermo<solidThermo::composite>(mesh, phaseName),
    Cv_(readProperty<scalar>("Cv", dimEnergy/dimMass/dimTemperature)),
    e_
    (
        IOobject
        (
            phasePropertyName("e"),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        T_,
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes(),
        this->heSourcesTypes()
    )
{

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constSolidThermoo::constSolidThermoo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    constSolidThermoo(mesh, true, phaseName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constSolidThermoo::~constSolidThermoo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::constSolidThermoo::W() const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermoo::W(const label patchi) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


const Foam::volScalarField& Foam::constSolidThermoo::he() const
{
    return e_;
}


Foam::volScalarField& Foam::constSolidThermoo::he()
{
    return e_;
}


const Foam::volScalarField& Foam::constSolidThermoo::Cp() const
{
    return T_; 
}


const Foam::volScalarField& Foam::constSolidThermoo::Cv() const
{
    return T_; 
}


const Foam::volScalarField& Foam::constSolidThermoo::Cpv() const
{
    return T_; 
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermoo::he
(
    const scalarField& T,
    const labelList& cells
) const
{
    return scalarField(T_, cells)*T;
}


Foam::tmp<Foam::volScalarField> Foam::constSolidThermoo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return T;
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermoo::he
(
    const scalarField& T,
    const label patchi
) const
{
    return T_.boundaryField()[patchi]*T; 
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermoo::he
(
    const scalarField& T,
    const fvSource& source
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::constSolidThermoo::hs() const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::constSolidThermoo::hs
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermoo::hs
(
    const scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermoo::hs
(
    const scalarField& T,
    const labelList& cells
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::constSolidThermoo::ha() const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::constSolidThermoo::ha
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermoo::ha
(
    const scalarField& T,
    const label patchi
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermoo::ha
(
    const scalarField& T,
    const labelList& cells
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermoo::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    return T_.boundaryField()[patchi]; 
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermoo::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    return T_.boundaryField()[patchi]; 
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermoo::Cpv
(
    const scalarField& T,
    const label patchi
) const
{
    return T_.boundaryField()[patchi]; 
}


 Foam::tmp<Foam::volScalarField> Foam::constSolidThermoo::The
(
    const volScalarField& h,
    const volScalarField& p,
    const volScalarField& T0
) const
{
    NotImplemented;
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermoo::The
(
    const scalarField& he,
    const scalarField& T0,
    const labelList& cells
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermoo::The
(
    const scalarField& he,
    const scalarField& T0,
    const label patchi
) const
{
    NotImplemented;
    return tmp<scalarField>(nullptr);
} 


const Foam::volVectorField& Foam::constSolidThermoo::Kappa() const
{
    NotImplemented;
    return volVectorField::null();
}


void Foam::constSolidThermoo::correct()
{
    T_ = e_;
}


// ************************************************************************* //
