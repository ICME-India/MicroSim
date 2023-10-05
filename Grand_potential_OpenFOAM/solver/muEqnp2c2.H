
//! The unit normal vector to the interface with a small number in denominator to prevent solution from diverging
n=dimx*fvc::grad(phi_1)/(1E-20+mag(dimx*fvc::grad(phi_1)));

//! Chemical potential equation with approximate slope of c-mu curve according to the parabolic approximation for free energy
if (phases == 2)
{
hphi = phi_1*phi_1*phi_1*(6*phi_1*phi_1 - 15*phi_1 +10);
if (components == 2)
{
fvScalarMatrix muEqn
(
    (0.5)*((hphi/A[0]) + ((1-hphi)/A[int(phases)-1]))*dimt*fvm::ddt(mu_1) == dimx*dimx*fvm::laplacian(diff_1*0.5*(1/A[0])*phi_1 + diff_Liq*0.5*(1/A[int(phases)-1])*(1-phi_1), mu_1) - 0.5*((mu_1 - B_bin[0])/A[0] - (mu_1 - B_bin[int(phases)-1])/A[int(phases)-1])*dimt*fvc::ddt(phi_1)*30.0*phi_1*phi_1*(1.0-phi_1)*(1.0-phi_1) - Tdot*(((-B_bin_dT[0])/(2*A[0]) + (mu_1 - B_bin[0])*A_dT[0]/(-2*A[0]*A[0]))*hphi + (mu_1 - B_bin[int(phases)-1])*A_dT[int(phases)-1]*(1-hphi)/(-2*A[int(phases)-1]*A[int(phases)-1])) - anti_trap*epsilon*0.5*((mu_1 - B_bin[0])/A[0] - (mu_1 - B_bin[int(phases)-1])/A[int(phases)-1])*dimx*fvc::div((n*dimt*fvc::ddt(phi_1)))
);

//muEqn.solve();

InitialResidual_10 = muEqn.solve().max().initialResidual();

    Info<< "Min/max mu:" << min(mu_1).value() << ' ' << max(mu_1).value() << endl;
}

}

