
if (phases == 4)
{
if (components == 3)
{
//For phi_1
    
    fvScalarMatrix phi_1Eqn
    (
        dimt*omega*epsilon*fvm::ddt(phi_1) -2.0*gamma_1*epsilon*dimx*dimx*fvm::laplacian(phi_1) == -(18.0*gamma_1*phi_1*(1.0 -phi_1)*(1.0 -2.0*phi_1))/(epsilon)
        - ((-0.5*H11[0]*(dcdmu[0][0]*(mu_1 - B_ter1[0]) + dcdmu[0][1]*(mu_2 - B_ter2[0]))*(dcdmu[0][0]*(mu_1 - B_ter1[0]) + dcdmu[0][1]*(mu_2 - B_ter2[0])) - H12[0]*(dcdmu[0][0]*(mu_1 - B_ter1[0]) + dcdmu[0][1]*(mu_2 - B_ter2[0]))*(dcdmu[1][0]*(mu_1 - B_ter1[0]) + dcdmu[1][1]*(mu_2 - B_ter2[0])) - 0.5*H22[0]*(dcdmu[1][0]*(mu_1 - B_ter1[0]) + dcdmu[1][1]*(mu_2 - B_ter2[0]))*(dcdmu[1][0]*(mu_1 - B_ter1[0]) + dcdmu[1][1]*(mu_2 - B_ter2[0])) + D_ter[0])*(6.0*phi_1*(1-phi_1) + 2*(phi_2*phi_3 + phi_2*phi_4 + phi_3*phi_4))
         + (-0.5*H11[1]*(dcdmu[2][0]*(mu_1 - B_ter1[1]) + dcdmu[2][1]*(mu_2 - B_ter2[1]))*(dcdmu[2][0]*(mu_1 - B_ter1[1]) + dcdmu[2][1]*(mu_2 - B_ter2[1])) - H12[1]*(dcdmu[2][0]*(mu_1 - B_ter1[1]) + dcdmu[2][1]*(mu_2 - B_ter2[1]))*(dcdmu[3][0]*(mu_1 - B_ter1[1]) + dcdmu[3][1]*(mu_2 - B_ter2[1])) - 0.5*H22[1]*(dcdmu[3][0]*(mu_1 - B_ter1[1]) + dcdmu[3][1]*(mu_2 - B_ter2[1]))*(dcdmu[3][0]*(mu_1 - B_ter1[1]) + dcdmu[3][1]*(mu_2 - B_ter2[1])) + D_ter[1])*2*phi_2*(phi_3 + phi_4)
         + (-0.5*H11[2]*(dcdmu[4][0]*(mu_1 - B_ter1[2]) + dcdmu[4][1]*(mu_2 - B_ter2[2]))*(dcdmu[4][0]*(mu_1 - B_ter1[2]) + dcdmu[4][1]*(mu_2 - B_ter2[2])) - H12[2]*(dcdmu[4][0]*(mu_1 - B_ter1[2]) + dcdmu[4][1]*(mu_2 - B_ter2[2]))*(dcdmu[5][0]*(mu_1 - B_ter1[2]) + dcdmu[5][1]*(mu_2 - B_ter2[2])) - 0.5*H22[2]*(dcdmu[5][0]*(mu_1 - B_ter1[2]) + dcdmu[5][1]*(mu_2 - B_ter2[2]))*(dcdmu[5][0]*(mu_1 - B_ter1[2]) + dcdmu[5][1]*(mu_2 - B_ter2[2])) + D_ter[2])*2*phi_3*(phi_2 + phi_4)
           + (-0.5*H11[int(phases)-1]*(dcdmu[2*(int(phases)-1)][0]*mu_1 + dcdmu[2*(int(phases)-1)][1]*mu_2)*(dcdmu[2*(int(phases)-1)][0]*mu_1 + dcdmu[2*(int(phases)-1)][1]*mu_2) - H12[int(phases)-1]*(dcdmu[2*(int(phases)-1)][0]*mu_1 + dcdmu[2*(int(phases)-1)][1]*mu_2)*(dcdmu[2*(int(phases)-1)+1][0]*mu_1 + dcdmu[2*(int(phases)-1)+1][1]*mu_2) - 0.5*H22[int(phases)-1]*(dcdmu[2*(int(phases)-1)+1][0]*mu_1 + dcdmu[2*(int(phases)-1)+1][1]*mu_2)*(dcdmu[2*(int(phases)-1)+1][0]*mu_1 + dcdmu[2*(int(phases)-1)+1][1]*mu_2))*2*phi_4*(phi_2 + phi_3))*(1 + noise_mag*randNum.globalScalar01())/Vm
             - (6.0*phi_1*(1-phi_1)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 - phi_2*phi_2*(3-2*phi_2)*cEigenStrain2
                                                                    - phi_3*phi_3*(3-2*phi_3)*cEigenStrain3))-(Sigma && cEigenStrain1)))
    );
//    phi_1Eqn -= 6.0*phi_1*(1-phi_1)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 - phi_2*phi_2*(3-2*phi_2)*cEigenStrain2 
//								    - phi_3*phi_3*(3-2*phi_3)*cEigenStrain3))-(Sigma && cEigenStrain1));

    
//For phi_2        
    
    fvScalarMatrix phi_2Eqn
    (
        dimt*omega*epsilon*fvm::ddt(phi_2) -2.0*gamma_2*epsilon*dimx*dimx*fvm::laplacian(phi_2) == -(18.0*gamma_2*phi_2*(1.0 -phi_2)*(1.0 -2.0*phi_2))/(epsilon)
        - ((-0.5*H11[0]*(dcdmu[0][0]*(mu_1 - B_ter1[0]) + dcdmu[0][1]*(mu_2 - B_ter2[0]))*(dcdmu[0][0]*(mu_1 - B_ter1[0]) + dcdmu[0][1]*(mu_2 - B_ter2[0])) - H12[0]*(dcdmu[0][0]*(mu_1 - B_ter1[0]) + dcdmu[0][1]*(mu_2 - B_ter2[0]))*(dcdmu[1][0]*(mu_1 - B_ter1[0]) + dcdmu[1][1]*(mu_2 - B_ter2[0])) - 0.5*H22[0]*(dcdmu[1][0]*(mu_1 - B_ter1[0]) + dcdmu[1][1]*(mu_2 - B_ter2[0]))*(dcdmu[1][0]*(mu_1 - B_ter1[0]) + dcdmu[1][1]*(mu_2 - B_ter2[0])) + D_ter[0])*2*phi_1*(phi_3 + phi_4)
         + (-0.5*H11[1]*(dcdmu[2][0]*(mu_1 - B_ter1[1]) + dcdmu[2][1]*(mu_2 - B_ter2[1]))*(dcdmu[2][0]*(mu_1 - B_ter1[1]) + dcdmu[2][1]*(mu_2 - B_ter2[1])) - H12[1]*(dcdmu[2][0]*(mu_1 - B_ter1[1]) + dcdmu[2][1]*(mu_2 - B_ter2[1]))*(dcdmu[3][0]*(mu_1 - B_ter1[1]) + dcdmu[3][1]*(mu_2 - B_ter2[1])) - 0.5*H22[1]*(dcdmu[3][0]*(mu_1 - B_ter1[1]) + dcdmu[3][1]*(mu_2 - B_ter2[1]))*(dcdmu[3][0]*(mu_1 - B_ter1[1]) + dcdmu[3][1]*(mu_2 - B_ter2[1])) + D_ter[1])*(6.0*phi_2*(1-phi_2) + 2*(phi_1*phi_3 + phi_1*phi_4 + phi_3*phi_4))
         + (-0.5*H11[2]*(dcdmu[4][0]*(mu_1 - B_ter1[2]) + dcdmu[4][1]*(mu_2 - B_ter2[2]))*(dcdmu[4][0]*(mu_1 - B_ter1[2]) + dcdmu[4][1]*(mu_2 - B_ter2[2])) - H12[2]*(dcdmu[4][0]*(mu_1 - B_ter1[2]) + dcdmu[4][1]*(mu_2 - B_ter2[2]))*(dcdmu[5][0]*(mu_1 - B_ter1[2]) + dcdmu[5][1]*(mu_2 - B_ter2[2])) - 0.5*H22[2]*(dcdmu[5][0]*(mu_1 - B_ter1[2]) + dcdmu[5][1]*(mu_2 - B_ter2[2]))*(dcdmu[5][0]*(mu_1 - B_ter1[2]) + dcdmu[5][1]*(mu_2 - B_ter2[2])) + D_ter[2])*2*phi_3*(phi_1 + phi_4)
          + (-0.5*H11[int(phases)-1]*(dcdmu[2*(int(phases)-1)][0]*mu_1 + dcdmu[2*(int(phases)-1)][1]*mu_2)*(dcdmu[2*(int(phases)-1)][0]*mu_1 + dcdmu[2*(int(phases)-1)][1]*mu_2) - H12[int(phases)-1]*(dcdmu[2*(int(phases)-1)][0]*mu_1 + dcdmu[2*(int(phases)-1)][1]*mu_2)*(dcdmu[2*(int(phases)-1)+1][0]*mu_1 + dcdmu[2*(int(phases)-1)+1][1]*mu_2) - 0.5*H22[int(phases)-1]*(dcdmu[2*(int(phases)-1)+1][0]*mu_1 + dcdmu[2*(int(phases)-1)+1][1]*mu_2)*(dcdmu[2*(int(phases)-1)+1][0]*mu_1 + dcdmu[2*(int(phases)-1)+1][1]*mu_2))*2*phi_4*(phi_1 + phi_3))*(1 + noise_mag*randNum.globalScalar01())/Vm
		- (6.0*phi_2*(1-phi_2)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 - phi_2*phi_2*(3-2*phi_2)*cEigenStrain2
                                                                    - phi_3*phi_3*(3-2*phi_3)*cEigenStrain3))-(Sigma && cEigenStrain2)))
    );
//    phi_2Eqn -= 6.0*phi_2*(1-phi_2)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 - phi_2*phi_2*(3-2*phi_2)*cEigenStrain2 
//								    - phi_3*phi_3*(3-2*phi_3)*cEigenStrain3))-(Sigma && cEigenStrain2));

    
//For phi_3  
    
    fvScalarMatrix phi_3Eqn
    (
        dimt*omega*epsilon*fvm::ddt(phi_3) -2.0*gamma_3*epsilon*dimx*dimx*fvm::laplacian(phi_3) == -(18.0*gamma_3*phi_3*(1.0 -phi_3)*(1.0 -2.0*phi_3))/(epsilon)
        - ((-0.5*H11[0]*(dcdmu[0][0]*(mu_1 - B_ter1[0]) + dcdmu[0][1]*(mu_2 - B_ter2[0]))*(dcdmu[0][0]*(mu_1 - B_ter1[0]) + dcdmu[0][1]*(mu_2 - B_ter2[0])) - H12[0]*(dcdmu[0][0]*(mu_1 - B_ter1[0]) + dcdmu[0][1]*(mu_2 - B_ter2[0]))*(dcdmu[1][0]*(mu_1 - B_ter1[0]) + dcdmu[1][1]*(mu_2 - B_ter2[0])) - 0.5*H22[0]*(dcdmu[1][0]*(mu_1 - B_ter1[0]) + dcdmu[1][1]*(mu_2 - B_ter2[0]))*(dcdmu[1][0]*(mu_1 - B_ter1[0]) + dcdmu[1][1]*(mu_2 - B_ter2[0])) + D_ter[0])*2*phi_1*(phi_2 + phi_4)
         + (-0.5*H11[1]*(dcdmu[2][0]*(mu_1 - B_ter1[1]) + dcdmu[2][1]*(mu_2 - B_ter2[1]))*(dcdmu[2][0]*(mu_1 - B_ter1[1]) + dcdmu[2][1]*(mu_2 - B_ter2[1])) - H12[1]*(dcdmu[2][0]*(mu_1 - B_ter1[1]) + dcdmu[2][1]*(mu_2 - B_ter2[1]))*(dcdmu[3][0]*(mu_1 - B_ter1[1]) + dcdmu[3][1]*(mu_2 - B_ter2[1])) - 0.5*H22[1]*(dcdmu[3][0]*(mu_1 - B_ter1[1]) + dcdmu[3][1]*(mu_2 - B_ter2[1]))*(dcdmu[3][0]*(mu_1 - B_ter1[1]) + dcdmu[3][1]*(mu_2 - B_ter2[1])) + D_ter[1])*2*phi_2*(phi_1 + phi_4)
         + (-0.5*H11[2]*(dcdmu[4][0]*(mu_1 - B_ter1[2]) + dcdmu[4][1]*(mu_2 - B_ter2[2]))*(dcdmu[4][0]*(mu_1 - B_ter1[2]) + dcdmu[4][1]*(mu_2 - B_ter2[2])) - H12[2]*(dcdmu[4][0]*(mu_1 - B_ter1[2]) + dcdmu[4][1]*(mu_2 - B_ter2[2]))*(dcdmu[5][0]*(mu_1 - B_ter1[2]) + dcdmu[5][1]*(mu_2 - B_ter2[2])) - 0.5*H22[2]*(dcdmu[5][0]*(mu_1 - B_ter1[2]) + dcdmu[5][1]*(mu_2 - B_ter2[2]))*(dcdmu[5][0]*(mu_1 - B_ter1[2]) + dcdmu[5][1]*(mu_2 - B_ter2[2])) + D_ter[2])*(6.0*phi_3*(1-phi_3) + 2*(phi_1*phi_2 + phi_1*phi_4 + phi_2*phi_4))
         + (-0.5*H11[int(phases)-1]*(dcdmu[2*(int(phases)-1)][0]*mu_1 + dcdmu[2*(int(phases)-1)][1]*mu_2)*(dcdmu[2*(int(phases)-1)][0]*mu_1 + dcdmu[2*(int(phases)-1)][1]*mu_2) - H12[int(phases)-1]*(dcdmu[2*(int(phases)-1)][0]*mu_1 + dcdmu[2*(int(phases)-1)][1]*mu_2)*(dcdmu[2*(int(phases)-1)+1][0]*mu_1 + dcdmu[2*(int(phases)-1)+1][1]*mu_2) - 0.5*H22[int(phases)-1]*(dcdmu[2*(int(phases)-1)+1][0]*mu_1 + dcdmu[2*(int(phases)-1)+1][1]*mu_2)*(dcdmu[2*(int(phases)-1)+1][0]*mu_1 + dcdmu[2*(int(phases)-1)+1][1]*mu_2))*2*phi_4*(phi_1 + phi_2))*(1 + noise_mag*randNum.globalScalar01())/Vm
			- (6.0*phi_3*(1-phi_3)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 - phi_2*phi_2*(3-2*phi_2)*cEigenStrain2
                                                                    - phi_3*phi_3*(3-2*phi_3)*cEigenStrain3))-(Sigma && cEigenStrain3)))
    );
//    phi_3Eqn -= 6.0*phi_3*(1-phi_3)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 - phi_2*phi_2*(3-2*phi_2)*cEigenStrain2 
//								    - phi_3*phi_3*(3-2*phi_3)*cEigenStrain3))-(Sigma && cEigenStrain3));    

    
//For phi_4     
    
     deltaSigmaD = ((mu2_elast)*twoSymm(fvc::grad(D)) + (lambda2)*(I*tr(fvc::grad(D))) 
		   - 2*(mu2_elast)*(phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 + phi_2*phi_2*(3-2*phi_2)*cEigenStrain2 +phi_3*phi_3*(3-2*phi_3)*cEigenStrain3) 
		   - (lambda2)*(I*tr((phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 + phi_2*phi_2*(3-2*phi_2)*cEigenStrain2 +phi_3*phi_3*(3-2*phi_3)*cEigenStrain3))))
		   + (mu2_elast_)*strain;
    
    fvScalarMatrix phi_4Eqn
    (
        dimt*omega*epsilon*fvm::ddt(phi_4) -2.0*gamma_4*epsilon*dimx*dimx*fvm::laplacian(phi_4) == -(18.0*gamma_4*phi_4*(1.0 -phi_4)*(1.0 -2.0*phi_4))/(epsilon)
        - ((-0.5*H11[0]*(dcdmu[0][0]*(mu_1 - B_ter1[0]) + dcdmu[0][1]*(mu_2 - B_ter2[0]))*(dcdmu[0][0]*(mu_1 - B_ter1[0]) + dcdmu[0][1]*(mu_2 - B_ter2[0])) - H12[0]*(dcdmu[0][0]*(mu_1 - B_ter1[0]) + dcdmu[0][1]*(mu_2 - B_ter2[0]))*(dcdmu[1][0]*(mu_1 - B_ter1[0]) + dcdmu[1][1]*(mu_2 - B_ter2[0])) - 0.5*H22[0]*(dcdmu[1][0]*(mu_1 - B_ter1[0]) + dcdmu[1][1]*(mu_2 - B_ter2[0]))*(dcdmu[1][0]*(mu_1 - B_ter1[0]) + dcdmu[1][1]*(mu_2 - B_ter2[0])) + D_ter[0])*2*phi_1*(phi_2 + phi_3)
         + (-0.5*H11[1]*(dcdmu[2][0]*(mu_1 - B_ter1[1]) + dcdmu[2][1]*(mu_2 - B_ter2[1]))*(dcdmu[2][0]*(mu_1 - B_ter1[1]) + dcdmu[2][1]*(mu_2 - B_ter2[1])) - H12[1]*(dcdmu[2][0]*(mu_1 - B_ter1[1]) + dcdmu[2][1]*(mu_2 - B_ter2[1]))*(dcdmu[3][0]*(mu_1 - B_ter1[1]) + dcdmu[3][1]*(mu_2 - B_ter2[1])) - 0.5*H22[1]*(dcdmu[3][0]*(mu_1 - B_ter1[1]) + dcdmu[3][1]*(mu_2 - B_ter2[1]))*(dcdmu[3][0]*(mu_1 - B_ter1[1]) + dcdmu[3][1]*(mu_2 - B_ter2[1])) + D_ter[1])*2*phi_2*(phi_1 + phi_3)
         + (-0.5*H11[2]*(dcdmu[4][0]*(mu_1 - B_ter1[2]) + dcdmu[4][1]*(mu_2 - B_ter2[2]))*(dcdmu[4][0]*(mu_1 - B_ter1[2]) + dcdmu[4][1]*(mu_2 - B_ter2[2])) - H12[2]*(dcdmu[4][0]*(mu_1 - B_ter1[2]) + dcdmu[4][1]*(mu_2 - B_ter2[2]))*(dcdmu[5][0]*(mu_1 - B_ter1[2]) + dcdmu[5][1]*(mu_2 - B_ter2[2])) - 0.5*H22[2]*(dcdmu[5][0]*(mu_1 - B_ter1[2]) + dcdmu[5][1]*(mu_2 - B_ter2[2]))*(dcdmu[5][0]*(mu_1 - B_ter1[2]) + dcdmu[5][1]*(mu_2 - B_ter2[2])) + D_ter[2])*2*phi_3*(phi_1 + phi_2)
         + (-0.5*H11[int(phases)-1]*(dcdmu[2*(int(phases)-1)][0]*mu_1 + dcdmu[2*(int(phases)-1)][1]*mu_2)*(dcdmu[2*(int(phases)-1)][0]*mu_1 + dcdmu[2*(int(phases)-1)][1]*mu_2) - H12[int(phases)-1]*(dcdmu[2*(int(phases)-1)][0]*mu_1 + dcdmu[2*(int(phases)-1)][1]*mu_2)*(dcdmu[2*(int(phases)-1)+1][0]*mu_1 + dcdmu[2*(int(phases)-1)+1][1]*mu_2) - 0.5*H22[int(phases)-1]*(dcdmu[2*(int(phases)-1)+1][0]*mu_1 + dcdmu[2*(int(phases)-1)+1][1]*mu_2)*(dcdmu[2*(int(phases)-1)+1][0]*mu_1 + dcdmu[2*(int(phases)-1)+1][1]*mu_2))*(6.0*phi_4*(1-phi_4) + 2*(phi_1*phi_2 + phi_1*phi_3 + phi_2*phi_3)))*(1 + noise_mag*randNum.globalScalar01())/Vm
			- 6.0*phi_4*(1-phi_4)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 - phi_2*phi_2*(3-2*phi_2)*cEigenStrain2
                                                                    - phi_3*phi_3*(3-2*phi_3)*cEigenStrain3)))
    );
    
    
    phi_1Eqn.solve();
    phi_2Eqn.solve();
    phi_3Eqn.solve();
    phi_4Eqn.solve();

  Info<< "Min/max sigmaproduct:" << min((Sigma && cEigenStrain2)()).value() << '\t'<< max((Sigma && cEigenStrain2)()).value() << endl;
 
//    scalar numphases = 4.0;   
    
    volScalarField lambda_phase = dimt*((fvc::ddt(phi_1) + fvc::ddt(phi_2) + fvc::ddt(phi_3) + fvc::ddt(phi_4))/phases);
    Info<< "Min/max lambda_phase:" << min(lambda_phase()).value() << '\t'<< max(lambda_phase()).value() << endl;


//=============== Substracting lambda_phase from evolution equations ==============================//  
    
  //Update the phase-values such that sum upto 1
    phi_1 -=  lambda_phase*runTime.deltaTValue();
    phi_2 -=  lambda_phase*runTime.deltaTValue();
    phi_3 -=  lambda_phase*runTime.deltaTValue();
    phi_4 -=  lambda_phase*runTime.deltaTValue();
}
phase_index = phi_1 + 2*phi_2 + 3*phi_3 + 4*phi_4;
}

    // Writing into the log file
    Info<< "Min/max phi_1:" << min(phi_1).value() << ' ' << max(phi_1).value() << endl;

if (phases == 3 || phases == 4)
{
Info<< "Min/max phi_2:" << min(phi_2).value() << ' ' << max(phi_2).value() << endl;
Info<< "Min/max phi_3:" << min(phi_3).value() << ' ' << max(phi_3).value() << endl;
}

if (phases == 4)
{
Info<< "Min/max phi_4:" << min(phi_4).value() << ' ' << max(phi_4).value() << endl;
}

