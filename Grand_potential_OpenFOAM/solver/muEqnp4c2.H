

//! Chemical potential equation with approximate slope of c-mu curve according to the parabolic approximation for free energy
if (phases == 4)
{
n2=dimx*fvc::grad(phi_2)/(1E-20+mag(dimx*fvc::grad(phi_2)));
n3=dimx*fvc::grad(phi_3)/(1E-20+mag(dimx*fvc::grad(phi_3)));
n4=dimx*fvc::grad(phi_4)/(1E-20+mag(dimx*fvc::grad(phi_4)));

hphi = phi_1*phi_1*(3-2*phi_1) + 2*phi_1*(phi_2*phi_3 + phi_2*phi_4 + phi_3*phi_4);
hphi2 = phi_2*phi_2*(3-2*phi_2) + 2*phi_2*(phi_1*phi_3 + phi_1*phi_4 + phi_3*phi_4);
hphi3 = phi_3*phi_3*(3-2*phi_3) + 2*phi_3*(phi_1*phi_2 + phi_1*phi_4 + phi_2*phi_4);
hphi4 = phi_4*phi_4*(3-2*phi_4) + 2*phi_4*(phi_1*phi_2 + phi_1*phi_3 + phi_2*phi_3);
if (components == 2)
{
fvScalarMatrix muEqn
(
    (0.5)*((hphi/A[0]) + (hphi2/A[1]) + (hphi3/A[2]) + (hphi4/A[int(phases)-1]))*dimt*fvm::ddt(mu_1) == dimx*dimx*fvm::laplacian(diff_1*0.5*(1/A[0])*phi_1 + diff_2*0.5*(1/A[1])*phi_2 + diff_3*0.5*(1/A[2])*phi_3 + diff_Liq*0.5*(1/A[int(phases)-1])*phi_4, mu_1) 
    - 0.5*((mu_1 - B_bin[0])*dimt*fvc::ddt(phi_1)*(6.0*phi_1*(1-phi_1) + 2*(phi_2*phi_3 + phi_2*phi_4 + phi_3*phi_4))/A[0]
     + (mu_1 - B_bin[1])*dimt*fvc::ddt(phi_2)*(6.0*phi_2*(1-phi_2) + 2*(phi_1*phi_3 + phi_1*phi_4 + phi_3*phi_4))/A[1]
      + (mu_1 - B_bin[2])*dimt*fvc::ddt(phi_3)*(6.0*phi_3*(1-phi_3) + 2*(phi_1*phi_2 + phi_1*phi_4 + phi_2*phi_4))/A[2]
       + (mu_1 - B_bin[int(phases)-1])*dimt*fvc::ddt(phi_4)*(6.0*phi_4*(1-phi_4) + 2*(phi_1*phi_2 + phi_1*phi_3 + phi_2*phi_3))/A[int(phases)-1]) 
    - Tdot*(((-B_bin_dT[0])/(2*A[0]) + (mu_1 - B_bin[0])*A_dT[0]/(-2*A[0]*A[0]))*hphi + ((-B_bin_dT[1])/(2*A[1]) + (mu_1 - B_bin[1])*A_dT[1]/(-2*A[1]*A[1]))*hphi2 + ((-B_bin_dT[2])/(2*A[2]) + (mu_1 - B_bin[2])*A_dT[2]/(-2*A[2]*A[2]))*hphi3 + (mu_1 - B_bin[int(phases)-1])*A_dT[int(phases)-1]*hphi4/(-2*A[int(phases)-1]*A[int(phases)-1]))
     //- anti_trap*epsilon*0.5*((mu_1 - B_Sol)/A_Sol - (mu_1 - B_Liq)/A_Liq)*dimx*fvc::div((n*dimt*fvc::ddt(phi_1)))
     + 0.5*anti_trap*epsilon *((mu_1 - B_bin[int(phases)-1])/A[int(phases)-1] - (mu_1 - B_bin[0])/A[0])*dimx*fvc::div(n&n4*dimt*fvc::ddt(phi_1)*-1*n) // see this last term
     + 0.5*anti_trap*epsilon *((mu_1 - B_bin[int(phases)-1])/A[int(phases)-1] - (mu_1 - B_bin[1])/A[1])* dimx*fvc::div(n2&n4* dimt*fvc::ddt(phi_2)* -1*n2)
     + 0.5*anti_trap*epsilon *((mu_1 - B_bin[int(phases)-1])/A[int(phases)-1] - (mu_1 - B_bin[2])/A[2])* dimx*fvc::div(n3&n4* dimt*fvc::ddt(phi_3)* -1*n3)
);

//muEqn.solve();

InitialResidual_10 = muEqn.solve().max().initialResidual();

    Info<< "Min/max mu:" << min(mu_1).value() << ' ' << max(mu_1).value() << endl;
}

}

