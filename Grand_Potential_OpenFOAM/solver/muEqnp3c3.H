
//! Chemical potential equation with approximate slope of c-mu curve according to the parabolic approximation for free energy
if (phases == 3)
{

if (components == 3)
{
    do {
      
      fvScalarMatrix muEqn_1 (
        (dcdmu[0][0]*hphi + dcdmu[2][0]*hphi2 + dcdmu[2*(int(phases)-1)][0]*hphi3)*dimt*fvm::ddt(mu_1) +  (dcdmu[0][1]*hphi + dcdmu[2][1]*hphi2 + dcdmu[2*(int(phases)-1)][1]*hphi3)*dimt*fvc::ddt(mu_2) 
        == dimx*dimx*fvm::laplacian(Mob[0][0]*phi_1 + Mob[2][0]*phi_2 + Mob[2*(int(phases)-1)][0]*phi_3, mu_1) + dimx*dimx*fvc::laplacian(Mob[0][1]*phi_1 + Mob[2][1]*phi_2 + Mob[2*(int(phases)-1)][1]*phi_3, mu_2)
        -((dcdmu[0][0]*(mu_1 - B_ter1[0]) + dcdmu[0][1]*(mu_2 - B_ter2[0]))* dimt*fvc::ddt(phi_1)*0.25*15*(3.0*phi_1*phi_1*phi_1*phi_1 - 4.0*phi_1*phi_1*phi_1 - phi_1*phi_1 + 2.0*phi_1 - (phi_2 - phi_3)*(phi_2 - phi_3)*(2.0*phi_1 - 3.0*phi_1*phi_1))
        + (dcdmu[2][0]*(mu_1 - B_ter1[1]) + dcdmu[2][1]*(mu_2 - B_ter2[1]))* dimt*fvc::ddt(phi_2)*0.25*15*(3.0*phi_2*phi_2*phi_2*phi_2 - 4.0*phi_2*phi_2*phi_2 - phi_2*phi_2 + 2.0*phi_2 - (phi_3 - phi_1)*(phi_3 - phi_1)*(2.0*phi_2 - 3.0*phi_2*phi_2))
         + (dcdmu[2*(int(phases)-1)][0]*mu_1 + dcdmu[2*(int(phases)-1)][1]*mu_2)* dimt*fvc::ddt(phi_3)*0.25*15*(3.0*phi_3*phi_3*phi_3*phi_3 - 4.0*phi_3*phi_3*phi_3 - phi_3*phi_3 + 2.0*phi_3 - (phi_1 - phi_2)*(phi_1 - phi_2)*(2.0*phi_3 - 3.0*phi_3*phi_3)))
         - Tdot*((dcdmu[0][0]*(-B_ter1_dT[0]) + dcdmu_dT[0][0]*(mu_1 - B_ter1[0]) + dcdmu[0][1]*(-B_ter2_dT[0]) + dcdmu_dT[0][1]*(mu_2 - B_ter2[0]))*hphi + (dcdmu[2][0]*(-B_ter1_dT[1]) + dcdmu_dT[2][0]*(mu_1 - B_ter1[1]) + dcdmu[2][1]*(-B_ter2_dT[1]) + dcdmu_dT[2][1]*(mu_2 - B_ter2[1]))*hphi2 + (dcdmu_dT[2*(int(phases)-1)][0]*mu_1 + dcdmu_dT[2*(int(phases)-1)][1]*mu_2)*hphi3)
          + anti_trap*epsilon *((dcdmu[0][0]*(mu_1 - B_ter1[0]) + dcdmu[0][1]*(mu_2 - B_ter2[0])) - (dcdmu[2*(int(phases)-1)][0]*mu_1 + dcdmu[2*(int(phases)-1)][1]*mu_2))*dimx*fvc::div(n&n3*dimt*fvc::ddt(phi_1)*-1*n) // see this last term
     + anti_trap*epsilon *((dcdmu[2][0]*(mu_1 - B_ter1[1]) + dcdmu[2][1]*(mu_2 - B_ter2[1])) - (dcdmu[2*(int(phases)-1)][0]*mu_1 + dcdmu[2*(int(phases)-1)][1]*mu_2))* dimx*fvc::div(n2&n3* dimt*fvc::ddt(phi_2)* -1*n2)
      );
      
      InitialResidual_10 = muEqn_1.solve().max().initialResidual();
      
      fvScalarMatrix muEqn_2 (
        (dcdmu[1][0]*hphi + dcdmu[3][0]*hphi2 + dcdmu[2*(int(phases)-1)+1][0]*hphi3)*dimt*fvc::ddt(mu_1) +  (dcdmu[1][1]*hphi + dcdmu[3][1]*hphi2 + dcdmu[2*(int(phases)-1)+1][1]*hphi3)*dimt*fvm::ddt(mu_2) == 
        dimx*dimx*fvc::laplacian(Mob[1][0]*phi_1 + Mob[3][0]*phi_2 + Mob[2*(int(phases)-1)+1][0]*phi_3, mu_1) + dimx*dimx*fvm::laplacian(Mob[1][1]*phi_1 + Mob[3][1]*phi_2 + Mob[2*(int(phases)-1)+1][1]*phi_3, mu_2)
        -((dcdmu[1][0]*(mu_1 - B_ter1[0]) + dcdmu[1][1]*(mu_2 - B_ter2[0]))*dimt*fvc::ddt(phi_1)*0.25*15*(3.0*phi_1*phi_1*phi_1*phi_1 - 4.0*phi_1*phi_1*phi_1 - phi_1*phi_1 + 2.0*phi_1 - (phi_2 - phi_3)*(phi_2 - phi_3)*(2.0*phi_1 - 3.0*phi_1*phi_1))
        + (dcdmu[3][0]*(mu_1 - B_ter1[1]) + dcdmu[3][1]*(mu_2 - B_ter2[1]))*dimt*fvc::ddt(phi_2)*0.25*15*(3.0*phi_2*phi_2*phi_2*phi_2 - 4.0*phi_2*phi_2*phi_2 - phi_2*phi_2 + 2.0*phi_2 - (phi_3 - phi_1)*(phi_3 - phi_1)*(2.0*phi_2 - 3.0*phi_2*phi_2))
         + (dcdmu[2*(int(phases)-1)+1][0]*mu_1 + dcdmu[2*(int(phases)-1)+1][1]*mu_2)*dimt*fvc::ddt(phi_3)*0.25*15*(3.0*phi_3*phi_3*phi_3*phi_3 - 4.0*phi_3*phi_3*phi_3 - phi_3*phi_3 + 2.0*phi_3 - (phi_1 - phi_2)*(phi_1 - phi_2)*(2.0*phi_3 - 3.0*phi_3*phi_3)))
         - Tdot*((dcdmu[1][0]*(-B_ter1_dT[0]) + dcdmu_dT[1][0]*(mu_1 - B_ter1[0]) + dcdmu[1][1]*(-B_ter2_dT[0]) + dcdmu_dT[1][1]*(mu_2 - B_ter2[0]))*hphi + (dcdmu[3][0]*(-B_ter1_dT[1]) + dcdmu_dT[3][0]*(mu_1 - B_ter1[1]) + dcdmu[3][1]*(-B_ter2_dT[1]) + dcdmu_dT[3][1]*(mu_2 - B_ter2[1]))*hphi2 + (dcdmu_dT[2*(int(phases)-1)+1][0]*mu_1 + dcdmu_dT[2*(int(phases)-1)+1][1]*mu_2)*hphi3)
          + anti_trap*epsilon *((dcdmu[1][0]*(mu_1 - B_ter1[0]) + dcdmu[1][1]*(mu_2 - B_ter2[0])) - (dcdmu[2*(int(phases)-1)+1][0]*mu_1 + dcdmu[2*(int(phases)-1)+1][1]*mu_2))*dimx*fvc::div(n&n3*dimt*fvc::ddt(phi_1)*-1*n) // see this last term
     + anti_trap*epsilon *((dcdmu[3][0]*(mu_1 - B_ter1[1]) + dcdmu[3][1]*(mu_2 - B_ter2[1])) - (dcdmu[2*(int(phases)-1)+1][0]*mu_1 + dcdmu[2*(int(phases)-1)+1][1]*mu_2))* dimx*fvc::div(n2&n3* dimt*fvc::ddt(phi_2)* -1*n2)
      );
      
      InitialResidual_11 = muEqn_2.solve().max().initialResidual();
      
    } while(InitialResidual_10 > 1e-5 || InitialResidual_11 > 1e-5);
    
    Info<< "Min/max mu_1:" << min(mu_1).value() << ' ' << max(mu_1).value() << endl;
    Info<< "Min/max mu_2:" << min(mu_2).value() << ' ' << max(mu_2).value() << endl;
}
}

