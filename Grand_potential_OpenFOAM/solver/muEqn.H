
//! The unit normal vector to the interface with a small number in denominator to prevent solution from diverging
volVectorField n=dimx*fvc::grad(phi)/(1E-20+mag(dimx*fvc::grad(phi)));

volScalarField hphi = phi*phi*phi*(6*phi*phi - 15*phi +10); //-phi*phi*phi*(10*phi*phi*phi - 36*phi*phi + 45*phi - 20);

//! Chemical potential equation with approximate slope of c-mu curve according to the parabolic approximation for free energy
if (components == 2)
{
fvScalarMatrix muEqn
(
    (0.5)*((hphi/A_Sol) + ((1-hphi)/A_Liq))*dimt*fvm::ddt(mu_1) == dimx*dimx*fvm::laplacian(diff_Sol*0.5*(1/A_Sol)*phi + diff_Liq*0.5*(1/A_Liq)*(1-phi), mu_1) - 0.5*((mu_1 - B_Sol)/A_Sol - (mu_1 - B_Liq)/A_Liq)*dimt*fvc::ddt(phi)*30.0*phi*phi*(1.0-phi)*(1.0-phi) - Tdot*(((-dBdT)/(2*A_Sol) + (mu_1 - B_Sol)*A_SoldT/(-2*A_Sol*A_Sol))*hphi + (mu_1 - B_Liq)*A_LiqdT*(1-hphi)/(-2*A_Liq*A_Liq)) - anti_trap*epsilon*0.5*((mu_1 - B_Sol)/A_Sol - (mu_1 - B_Liq)/A_Liq)*dimx*fvc::div((n*dimt*fvc::ddt(phi)))
);

//muEqn.solve();

InitialResidual_10 = muEqn.solve().max().initialResidual();

    Info<< "Min/max mu:" << min(mu_1).value() << ' ' << max(mu_1).value() << endl;
}

if (components == 3)
{
    do {
      
      fvScalarMatrix muEqn_1 (
        (dcdmu_a[0][0]*hphi + dcdmu_l[0][0]*(1-hphi))*dimt*fvm::ddt(mu_1) +  (dcdmu_a[0][1]*hphi + dcdmu_l[0][1]*(1-hphi))*dimt*fvc::ddt(mu_2) 
        == dimx*dimx*fvm::laplacian(M_a[0][0]*phi + M_l[0][0]*(1-phi), mu_1) + dimx*dimx*fvc::laplacian(M_a[0][1]*phi + M_l[0][1]*(1-phi), mu_2)
        -((dcdmu_a[0][0]*(mu_1 - B_a1) + dcdmu_a[0][1]*(mu_2 - B_a2)) - (dcdmu_l[0][0]*mu_1 + dcdmu_l[0][1]*mu_2))* dimt*fvc::ddt(phi)*30.0*phi*phi*(1.0-phi)*(1.0-phi) - Tdot*((dcdmu_a[0][0]*(-dB_a1dT) + dcdmu_adT[0][0]*(mu_1 - B_a1) + dcdmu_a[0][1]*(-dB_a2dT) + dcdmu_adT[0][1]*(mu_2 - B_a2))*hphi + (dcdmu_ldT[0][0]*mu_1 + dcdmu_ldT[0][1]*mu_2)*(1-hphi)) - anti_trap*epsilon*((dcdmu_a[0][0]*(mu_1 - B_a1) + dcdmu_a[0][1]*(mu_2 - B_a2)) - (dcdmu_l[0][0]*mu_1 + dcdmu_l[0][1]*mu_2))*dimx*fvc::div((n*dimt*fvc::ddt(phi)))
      );
      
      InitialResidual_10 = muEqn_1.solve().max().initialResidual();
      
      fvScalarMatrix muEqn_2 (
        (dcdmu_a[1][0]*hphi + dcdmu_l[1][0]*(1-hphi))*dimt*fvc::ddt(mu_1) +  (dcdmu_a[1][1]*hphi + dcdmu_l[1][1]*(1-hphi))*dimt*fvm::ddt(mu_2) == 
        dimx*dimx*fvc::laplacian(M_a[1][0]*phi + M_l[1][0]*(1-phi), mu_1) + dimx*dimx*fvm::laplacian(M_a[1][1]*phi + M_l[1][1]*(1-phi), mu_2)
        -((dcdmu_a[1][0]*(mu_1 - B_a1) + dcdmu_a[1][1]*(mu_2 - B_a2)) - (dcdmu_l[1][0]*mu_1 + dcdmu_l[1][1]*mu_2))*dimt*fvc::ddt(phi)*30.0*phi*phi*(1.0-phi)*(1.0-phi) - Tdot*((dcdmu_a[1][0]*(-dB_a1dT) + dcdmu_adT[1][0]*(mu_1 - B_a1) + dcdmu_a[1][1]*(-dB_a2dT) + dcdmu_adT[1][1]*(mu_2 - B_a2))*hphi + (dcdmu_ldT[1][0]*mu_1 + dcdmu_ldT[1][1]*mu_2)*(1-hphi)) - anti_trap*epsilon*((dcdmu_a[1][0]*(mu_1 - B_a1) + dcdmu_a[1][1]*(mu_2 - B_a2)) - (dcdmu_l[1][0]*mu_1 + dcdmu_l[1][1]*mu_2))*dimx*fvc::div((n*dimt*fvc::ddt(phi)))
      );
      
      InitialResidual_11 = muEqn_2.solve().max().initialResidual();
      
    } while(InitialResidual_10 > 1e-5 || InitialResidual_11 > 1e-5);
    
    Info<< "Min/max mu_1:" << min(mu_1).value() << ' ' << max(mu_1).value() << endl;
    Info<< "Min/max mu_2:" << min(mu_2).value() << ' ' << max(mu_2).value() << endl;
}


