

//! Chemical potential equation with approximate slope of c-mu curve according to the parabolic approximation for free energy
if (phases == 3)
{
n2=dimx*fvc::grad(phi_2)/(1E-20+mag(dimx*fvc::grad(phi_2)));
n3=dimx*fvc::grad(phi_3)/(1E-20+mag(dimx*fvc::grad(phi_3)));

hphi = 0.25*phi_1*phi_1*(15.0*(1-phi_1)*(1.0 + phi_1 - (phi_2 - phi_3)*(phi_2 - phi_3))
+ phi_1*(9.0*phi_1*phi_1 - 5.0));

hphi2 = 0.25*phi_2*phi_2*(15.0*(1-phi_2)*(1.0 + phi_2 - (phi_1 - phi_3)*(phi_1 - phi_3))
+ phi_2*(9.0*phi_2*phi_2 - 5.0));

hphi3 = 0.25*phi_3*phi_3*(15.0*(1-phi_3)*(1.0 + phi_3 - (phi_1 - phi_2)*(phi_1 - phi_2))
+ phi_3*(9.0*phi_3*phi_3 - 5.0));
                                                                                                        
if (components == 2)
{
fvScalarMatrix muEqn
(
    (0.5)*((hphi/A[0]) + (hphi2/A[1]) + (hphi3/A[int(phases)-1]))*dimt*fvm::ddt(mu_1) == dimx*dimx*fvm::laplacian(diff_1*0.5*(1/A[0])*phi_1 + diff_2*0.5*(1/A[1])*phi_2 + diff_Liq*0.5*(1/A[int(phases)-1])*phi_3, mu_1) 
    - 0.5*((mu_1 - B_bin[0])*dimt*fvc::ddt(phi_1)*0.25*15*(3.0*phi_1*phi_1*phi_1*phi_1 - 4.0*phi_1*phi_1*phi_1 - phi_1*phi_1 + 2.0*phi_1 - (phi_2 - phi_3)*(phi_2 - phi_3)*(2.0*phi_1 - 3.0*phi_1*phi_1))/A[0]
     + (mu_1 - B_bin[1])*dimt*fvc::ddt(phi_2)*0.25*15*(3.0*phi_2*phi_2*phi_2*phi_2 - 4.0*phi_2*phi_2*phi_2 - phi_2*phi_2 + 2.0*phi_2 - (phi_3 - phi_1)*(phi_3 - phi_1)*(2.0*phi_2 - 3.0*phi_2*phi_2))/A[1]
      + (mu_1 - B_bin[int(phases)-1])*dimt*fvc::ddt(phi_3)*0.25*15*(3.0*phi_3*phi_3*phi_3*phi_3 - 4.0*phi_3*phi_3*phi_3 - phi_3*phi_3 + 2.0*phi_3 - (phi_1 - phi_2)*(phi_1 - phi_2)*(2.0*phi_3 - 3.0*phi_3*phi_3))/A[int(phases)-1]) 
    - Tdot*(((-B_bin_dT[0])/(2*A[0]) + (mu_1 - B_bin[0])*A_dT[0]/(-2*A[0]*A[0]))*hphi + ((-B_bin_dT[1])/(2*A[1]) + (mu_1 - B_bin[1])*A_dT[1]/(-2*A[1]*A[1]))*hphi2 + (mu_1 - B_bin[int(phases)-1])*A_dT[int(phases)-1]*hphi3/(-2*A[int(phases)-1]*A[int(phases)-1]))
     //- anti_trap*epsilon*0.5*((mu_1 - B_Sol)/A_Sol - (mu_1 - B_Liq)/A_Liq)*dimx*fvc::div((n*dimt*fvc::ddt(phi_1)))
     + 0.5*anti_trap*epsilon *((mu_1 - B_bin[int(phases)-1])/A[int(phases)-1] - (mu_1 - B_bin[0])/A[0])*dimx*fvc::div(n&n3*dimt*fvc::ddt(phi_1)*-1*n) // see this last term
     + 0.5*anti_trap*epsilon *((mu_1 - B_bin[int(phases)-1])/A[int(phases)-1] - (mu_1 - B_bin[1])/A[1])* dimx*fvc::div(n2&n3* dimt*fvc::ddt(phi_2)* -1*n2)
);

//muEqn.solve();

InitialResidual_10 = muEqn.solve().max().initialResidual();

    Info<< "Min/max mu:" << min(mu_1).value() << ' ' << max(mu_1).value() << endl;
}

}

