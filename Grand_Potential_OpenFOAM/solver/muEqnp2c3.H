
//! Chemical potential equation with approximate slope of c-mu curve according to the parabolic approximation for free energy
if (phases == 2)
{

if (components == 3)
{
    do {
      
      fvScalarMatrix muEqn_1 (
        (dcdmu[0][0]*hphi + dcdmu[2*(int(phases)-1)][0]*(1-hphi))*dimt*fvm::ddt(mu_1) +  (dcdmu[0][1]*hphi + dcdmu[2*(int(phases)-1)][1]*(1-hphi))*dimt*fvc::ddt(mu_2) 
        == dimx*dimx*fvm::laplacian(Mob[0][0]*phi_1 + Mob[2*(int(phases)-1)][0]*(1-phi_1), mu_1) + dimx*dimx*fvc::laplacian(Mob[0][1]*phi_1 + Mob[2*(int(phases)-1)][1]*(1-phi_1), mu_2)
        -((dcdmu[0][0]*(mu_1 - B_ter1[0]) + dcdmu[0][1]*(mu_2 - B_ter2[0])) - (dcdmu[2*(int(phases)-1)][0]*mu_1 + dcdmu[2*(int(phases)-1)][1]*mu_2))* dimt*fvc::ddt(phi_1)*30.0*phi_1*phi_1*(1.0-phi_1)*(1.0-phi_1) - Tdot*((dcdmu[0][0]*(-B_ter1_dT[0]) + dcdmu_dT[0][0]*(mu_1 - B_ter1[0]) + dcdmu[0][1]*(-B_ter2_dT[0]) + dcdmu_dT[0][1]*(mu_2 - B_ter2[0]))*hphi + (dcdmu_dT[2*(int(phases)-1)][0]*mu_1 + dcdmu_dT[2*(int(phases)-1)][1]*mu_2)*(1-hphi)) - anti_trap*epsilon*((dcdmu[0][0]*(mu_1 - B_ter1[0]) + dcdmu[0][1]*(mu_2 - B_ter2[0])) - (dcdmu[2*(int(phases)-1)][0]*mu_1 + dcdmu[2*(int(phases)-1)][1]*mu_2))*dimx*fvc::div((n*dimt*fvc::ddt(phi_1)))
      );
      
      InitialResidual_10 = muEqn_1.solve().max().initialResidual();
      
      fvScalarMatrix muEqn_2 (
        (dcdmu[1][0]*hphi + dcdmu[2*(int(phases)-1)+1][0]*(1-hphi))*dimt*fvc::ddt(mu_1) +  (dcdmu[1][1]*hphi + dcdmu[2*(int(phases)-1)+1][1]*(1-hphi))*dimt*fvm::ddt(mu_2) == 
        dimx*dimx*fvc::laplacian(Mob[1][0]*phi_1 + Mob[2*(int(phases)-1)+1][0]*(1-phi_1), mu_1) + dimx*dimx*fvm::laplacian(Mob[1][1]*phi_1 + Mob[2*(int(phases)-1)+1][1]*(1-phi_1), mu_2)
        -((dcdmu[1][0]*(mu_1 - B_ter1[0]) + dcdmu[1][1]*(mu_2 - B_ter2[0])) - (dcdmu[2*(int(phases)-1)+1][0]*mu_1 + dcdmu[2*(int(phases)-1)+1][1]*mu_2))*dimt*fvc::ddt(phi_1)*30.0*phi_1*phi_1*(1.0-phi_1)*(1.0-phi_1) - Tdot*((dcdmu[1][0]*(-B_ter1_dT[0]) + dcdmu_dT[1][0]*(mu_1 - B_ter1[0]) + dcdmu[1][1]*(-B_ter2_dT[0]) + dcdmu_dT[1][1]*(mu_2 - B_ter2[0]))*hphi + (dcdmu_dT[2*(int(phases)-1)+1][0]*mu_1 + dcdmu_dT[2*(int(phases)-1)+1][1]*mu_2)*(1-hphi)) - anti_trap*epsilon*((dcdmu[1][0]*(mu_1 - B_ter1[0]) + dcdmu[1][1]*(mu_2 - B_ter2[0])) - (dcdmu[2*(int(phases)-1)+1][0]*mu_1 + dcdmu[2*(int(phases)-1)+1][1]*mu_2))*dimx*fvc::div((n*dimt*fvc::ddt(phi_1)))
      );
      
      InitialResidual_11 = muEqn_2.solve().max().initialResidual();
      
    } while(InitialResidual_10 > 1e-5 || InitialResidual_11 > 1e-5);
    
    Info<< "Min/max mu_1:" << min(mu_1).value() << ' ' << max(mu_1).value() << endl;
    Info<< "Min/max mu_2:" << min(mu_2).value() << ' ' << max(mu_2).value() << endl;
}
}

