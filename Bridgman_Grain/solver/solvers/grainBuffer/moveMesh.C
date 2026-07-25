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

\*---------------------------------------------------------------------------*/

#include "grainBuffer.H"
#include "fvcDiv.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::solvers::grainBuffer::moveMesh()
{
    // Skip everything if neighbour resources are unavailable
    if (!nbrResourcesAvailable_) return;

    bool moveDue = false;

    vector F =
        sharedData.lookup<vector>("shiftVector");

    fShift_.value() = F;

    Info << "shift value iddddd" << fShift_.value() << endl;


    if (mag(fShift_.value())>0) moveDue = true;


    // Translate mesh if liquidus is beyond mesh max offset


    if (moveDue)
    {  
        // Translate mesh by one pull increment in the y-direction        
        mesh_.movePoints(mesh_.points() + fShift_.value());  
        
        sharedData.set("shiftVector", vector(0, 0, 0));
    }

}

// ************************************************************************* //
