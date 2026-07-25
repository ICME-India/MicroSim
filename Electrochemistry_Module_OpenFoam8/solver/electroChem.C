/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Application
    laplacianFoam

Description
    Solves a simple Laplace equation, e.g. for thermal diffusion in a solid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//template <typename T>
//inline T call_g(T& phi) {
//	return phi * (1-phi);
//}
#include "Functions.H"

//inline auto call_g(auto& phi) {
//	return phi * (1-phi);


int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating eta distribution\n" << endl;
    
        
    while (simple.loop(runTime))
    {
        
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            volScalarField h_e = 30*eta*eta*(1-eta)*(1-eta);
            volScalarField c_Lip = ((exp((mu-(e_l*R*T))/(R*T)))/(1+exp((mu-(e_l*R*T))/(R*T))))*(1-(eta*eta*eta*(6*eta*eta-15*eta+10)));
            volScalarField c_l= (exp((mu-(e_l*R*T))/(R*T))/(1+exp((mu-(e_l*R*T))/(R*T))));
            volScalarField c_s= (exp((mu-(e_s*R*T))/(R*T))/(1+exp((mu-(e_s*R*T))/(R*T))));
            volScalarField chi= (1/(R*T))*(((exp((mu-(e_l*R*T))/(R*T)))/((1+exp((mu-(e_l*R*T))/(R*T)))*(1+exp((mu-(e_l*R*T))/(R*T)))))*(1-(eta*eta*eta*(6*eta*eta-15*eta+10)))+((C_m_s/C_m_l)*((exp((mu-(e_s*R*T))/(R*T)))/((1+exp((mu-(e_s*R*T))/(R*T)))*(1+exp((mu-(e_s*R*T))/(R*T)))))*(eta*eta*eta*(6*eta*eta-15*eta+10))));
            volScalarField sigma= S_s*(eta*eta*eta*(6*eta*eta-15*eta+10))+S_l*(1-(eta*eta*eta*(6*eta*eta-15*eta+10)));
            volScalarField Deff = D_s*(eta*eta*eta*(6*eta*eta-15*eta+10))+D_l*(1-(eta*eta*eta*(6*eta*eta-15*eta+10)));
            
            fvScalarMatrix etaEqn
            (
              dimt*fvm::ddt(eta) == -L_s*(-dimx*dimx*(1)*fvm::laplacian(k,eta) + 2*W*eta*(1-eta)*(1-2*eta))-L_e*h_e*(exp((1-alpha)*n*(F*phi)/(R*T))-((C_m_l*c_Lip)/C_0)*exp(-(alpha)*n*(F*phi)/(R*T)))
	     );
            etaEqn.solve();
            
            fvScalarMatrix muEqn
             (
               (chi)*dimt*fvm::ddt(mu) == dimx*dimx*fvm::laplacian(((D*c_Lip)/(R*T)),mu)+dimx*dimx*fvc::laplacian(((D*c_Lip*n*F)/(R*T)),phi) - dimt*(h_e)*(c_s*(C_m_s/C_m_l)-c_l)*fvc::ddt(eta)             
	     );            
             muEqn.solve();
                       
            solve
            (
               (1/(dx*dx*1e-12))*fvm::laplacian(sigma,phi) == (n*F*C_m_s_r/(4000))*fvc::ddt(eta)                
            );
        }
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
