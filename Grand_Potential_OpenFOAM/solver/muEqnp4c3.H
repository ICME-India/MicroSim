
//! Chemical potential equation with approximate slope of c-mu curve according to the parabolic approximation for free energy
if (phases == 4)
{

if (components == 3)
{
    do {
      
      fvScalarMatrix muEqn_1 (
        (dcdmu[0][0]*hphi + dcdmu[2][0]*hphi2 + dcdmu[4][0]*hphi3 + dcdmu[2*(int(phases)-1)][0]*hphi4)*dimt*fvm::ddt(mu_1) +  (dcdmu[0][1]*hphi + dcdmu[2][1]*hphi2 + dcdmu[4][1]*hphi3 + dcdmu[2*(int(phases)-1)][1]*hphi4)*dimt*fvc::ddt(mu_2) 
        == dimx*dimx*fvm::laplacian(Mob[0][0]*phi_1 + Mob[2][0]*phi_2 + Mob[4][0]*phi_3 + Mob[2*(int(phases)-1)][0]*phi_4, mu_1) + dimx*dimx*fvc::laplacian(Mob[0][1]*phi_1 + Mob[2][1]*phi_2 + Mob[4][1]*phi_3 + Mob[2*(int(phases)-1)][1]*phi_4, mu_2)
        -((dcdmu[0][0]*(mu_1 - B_ter1[0]) + dcdmu[0][1]*(mu_2 - B_ter2[0]))* dimt*fvc::ddt(phi_1)*(6.0*phi_1*(1-phi_1) + 2*(phi_2*phi_3 + phi_2*phi_4 + phi_3*phi_4))
        + (dcdmu[2][0]*(mu_1 - B_ter1[1]) + dcdmu[2][1]*(mu_2 - B_ter2[1]))* dimt*fvc::ddt(phi_2)*(6.0*phi_2*(1-phi_2) + 2*(phi_1*phi_3 + phi_1*phi_4 + phi_3*phi_4))
        + (dcdmu[4][0]*(mu_1 - B_ter1[2]) + dcdmu[4][1]*(mu_2 - B_ter2[2]))* dimt*fvc::ddt(phi_3)*(6.0*phi_3*(1-phi_3) + 2*(phi_1*phi_2 + phi_1*phi_4 + phi_2*phi_4))
         + (dcdmu[2*(int(phases)-1)][0]*mu_1 + dcdmu[2*(int(phases)-1)][1]*mu_2)* dimt*fvc::ddt(phi_4)*(6.0*phi_4*(1-phi_4) + 2*(phi_1*phi_2 + phi_1*phi_3 + phi_2*phi_3)))
         - Tdot*((dcdmu[0][0]*(-B_ter1_dT[0]) + dcdmu_dT[0][0]*(mu_1 - B_ter1[0]) + dcdmu[0][1]*(-B_ter2_dT[0]) + dcdmu_dT[0][1]*(mu_2 - B_ter2[0]))*hphi + (dcdmu[2][0]*(-B_ter1_dT[1]) + dcdmu_dT[2][0]*(mu_1 - B_ter1[1]) + dcdmu[2][1]*(-B_ter2_dT[1]) + dcdmu_dT[2][1]*(mu_2 - B_ter2[1]))*hphi2 + (dcdmu[4][0]*(-B_ter1_dT[2]) + dcdmu_dT[4][0]*(mu_1 - B_ter1[2]) + dcdmu[4][1]*(-B_ter2_dT[2]) + dcdmu_dT[4][1]*(mu_2 - B_ter2[2]))*hphi3 + (dcdmu_dT[2*(int(phases)-1)][0]*mu_1 + dcdmu_dT[2*(int(phases)-1)][1]*mu_2)*hphi4)
          + anti_trap*epsilon *((dcdmu[0][0]*(mu_1 - B_ter1[0]) + dcdmu[0][1]*(mu_2 - B_ter2[0])) - (dcdmu[2*(int(phases)-1)][0]*mu_1 + dcdmu[2*(int(phases)-1)][1]*mu_2))*dimx*fvc::div(n&n4*dimt*fvc::ddt(phi_1)*-1*n) // see this last term
     + anti_trap*epsilon *((dcdmu[2][0]*(mu_1 - B_ter1[1]) + dcdmu[2][1]*(mu_2 - B_ter2[1])) - (dcdmu[2*(int(phases)-1)][0]*mu_1 + dcdmu[2*(int(phases)-1)][1]*mu_2))* dimx*fvc::div(n2&n4* dimt*fvc::ddt(phi_2)* -1*n2)
     + anti_trap*epsilon *((dcdmu[4][0]*(mu_1 - B_ter1[2]) + dcdmu[4][1]*(mu_2 - B_ter2[2])) - (dcdmu[2*(int(phases)-1)][0]*mu_1 + dcdmu[2*(int(phases)-1)][1]*mu_2))* dimx*fvc::div(n3&n4* dimt*fvc::ddt(phi_3)* -1*n3)
      );
      
      InitialResidual_10 = muEqn_1.solve().max().initialResidual();
      
      fvScalarMatrix muEqn_2 (
        (dcdmu[1][0]*hphi + dcdmu[3][0]*hphi2 + dcdmu[5][0]*hphi3 + dcdmu[2*(int(phases)-1)+1][0]*hphi4)*dimt*fvc::ddt(mu_1) +  (dcdmu[1][1]*hphi + dcdmu[3][1]*hphi2 + dcdmu[5][1]*hphi3 + dcdmu[2*(int(phases)-1)+1][1]*hphi4)*dimt*fvm::ddt(mu_2) == 
        dimx*dimx*fvc::laplacian(Mob[1][0]*phi_1 + Mob[3][0]*phi_2 + Mob[5][0]*phi_3 + Mob[2*(int(phases)-1)+1][0]*phi_4, mu_1) + dimx*dimx*fvm::laplacian(Mob[1][1]*phi_1 + Mob[3][1]*phi_2 + Mob[5][1]*phi_3 + Mob[2*(int(phases)-1)+1][1]*phi_4, mu_2)
        -((dcdmu[1][0]*(mu_1 - B_ter1[0]) + dcdmu[1][1]*(mu_2 - B_ter2[0]))*dimt*fvc::ddt(phi_1)*(6.0*phi_1*(1-phi_1) + 2*(phi_2*phi_3 + phi_2*phi_4 + phi_3*phi_4))
        + (dcdmu[3][0]*(mu_1 - B_ter1[1]) + dcdmu[3][1]*(mu_2 - B_ter2[1]))*dimt*fvc::ddt(phi_2)*(6.0*phi_2*(1-phi_2) + 2*(phi_1*phi_3 + phi_1*phi_4 + phi_3*phi_4))
        + (dcdmu[5][0]*(mu_1 - B_ter1[2]) + dcdmu[5][1]*(mu_2 - B_ter2[2]))*dimt*fvc::ddt(phi_3)*(6.0*phi_3*(1-phi_3) + 2*(phi_1*phi_2 + phi_1*phi_4 + phi_2*phi_4))
         + (dcdmu[2*(int(phases)-1)+1][0]*mu_1 + dcdmu[2*(int(phases)-1)+1][1]*mu_2)*dimt*fvc::ddt(phi_4)*(6.0*phi_4*(1-phi_4) + 2*(phi_1*phi_2 + phi_1*phi_3 + phi_2*phi_3)))
         - Tdot*((dcdmu[1][0]*(-B_ter1_dT[0]) + dcdmu_dT[1][0]*(mu_1 - B_ter1[0]) + dcdmu[1][1]*(-B_ter2_dT[0]) + dcdmu_dT[1][1]*(mu_2 - B_ter2[0]))*hphi + (dcdmu[3][0]*(-B_ter1_dT[1]) + dcdmu_dT[3][0]*(mu_1 - B_ter1[1]) + dcdmu[3][1]*(-B_ter2_dT[1]) + dcdmu_dT[3][1]*(mu_2 - B_ter2[1]))*hphi2 + (dcdmu[5][0]*(-B_ter1_dT[2]) + dcdmu_dT[5][0]*(mu_1 - B_ter1[2]) + dcdmu[5][1]*(-B_ter2_dT[2]) + dcdmu_dT[5][1]*(mu_2 - B_ter2[2]))*hphi3 + (dcdmu_dT[2*(int(phases)-1)+1][0]*mu_1 + dcdmu_dT[2*(int(phases)-1)+1][1]*mu_2)*hphi4)
          + anti_trap*epsilon *((dcdmu[1][0]*(mu_1 - B_ter1[0]) + dcdmu[1][1]*(mu_2 - B_ter2[0])) - (dcdmu[2*(int(phases)-1)+1][0]*mu_1 + dcdmu[2*(int(phases)-1)+1][1]*mu_2))*dimx*fvc::div(n&n4*dimt*fvc::ddt(phi_1)*-1*n) // see this last term
     + anti_trap*epsilon *((dcdmu[3][0]*(mu_1 - B_ter1[1]) + dcdmu[3][1]*(mu_2 - B_ter2[1])) - (dcdmu[2*(int(phases)-1)+1][0]*mu_1 + dcdmu[2*(int(phases)-1)+1][1]*mu_2))* dimx*fvc::div(n2&n4* dimt*fvc::ddt(phi_2)* -1*n2)
     + anti_trap*epsilon *((dcdmu[5][0]*(mu_1 - B_ter1[2]) + dcdmu[5][1]*(mu_2 - B_ter2[2])) - (dcdmu[2*(int(phases)-1)+1][0]*mu_1 + dcdmu[2*(int(phases)-1)+1][1]*mu_2))* dimx*fvc::div(n3&n4* dimt*fvc::ddt(phi_3)* -1*n3)
      );
      
      InitialResidual_11 = muEqn_2.solve().max().initialResidual();
      
    } while(InitialResidual_10 > 1e-5 || InitialResidual_11 > 1e-5);
    
    Info<< "Min/max mu_1:" << min(mu_1).value() << ' ' << max(mu_1).value() << endl;
    Info<< "Min/max mu_2:" << min(mu_2).value() << ' ' << max(mu_2).value() << endl;
}
}

