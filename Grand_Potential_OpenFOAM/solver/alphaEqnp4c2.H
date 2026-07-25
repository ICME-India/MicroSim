if (swch == 2)
    {
    if (phases == 4)
    {
   strain=((gradD-((phi_1*phi_1*(3-2*phi_1)*cEigenStrain1)+(phi_2*phi_2*(3-2*phi_2)*cEigenStrain2)+(phi_3*phi_3*(3-2*phi_3)*cEigenStrain3)))
		&&symmTensor(1,0,0,0,0,0))*symmTensor(1,0,0,0,0,0)
	  +((gradD-((phi_1*phi_1*(3-2*phi_1)*cEigenStrain1)+(phi_2*phi_2*(3-2*phi_2)*cEigenStrain2)+(phi_3*phi_3*(3-2*phi_3)*cEigenStrain3)))
                &&symmTensor(0,0,0,1,0,0))*symmTensor(0,0,0,1,0,0)
	  +((gradD-((phi_1*phi_1*(3-2*phi_1)*cEigenStrain1)+(phi_2*phi_2*(3-2*phi_2)*cEigenStrain2)+(phi_3*phi_3*(3-2*phi_3)*cEigenStrain3)))
                &&symmTensor(0,0,0,0,0,1))*symmTensor(0,0,0,0,0,1);
//    vector eStrain = (T0*T0*(3-2*T0)*cEigenStrain0 + T1*T1*(3-2*T1)*cEigenStrain1 +T2*T2*(3-2*T2)*cEigenStrain2);
           
   Sigma = (2*(mu2_elast*phi_4*phi_4*(3-2*phi_4) + mu1_elast*(1-phi_4)*(1-phi_4)*(1+2*phi_4))*(symm(fvc::grad(D)) 
	      - (phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 + phi_2*phi_2*(3-2*phi_2)*cEigenStrain2 +phi_3*phi_3*(3-2*phi_3)*cEigenStrain3)) 
	  + (lambda2*phi_4*phi_4*(3-2*phi_4) + lambda1*(1-phi_4)*(1-phi_4)*(1+2*phi_4))*(I*tr(fvc::grad(D) 
	      - (phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 + phi_2*phi_2*(3-2*phi_2)*cEigenStrain2 +phi_3*phi_3*(3-2*phi_3)*cEigenStrain3))))
	  + (mu2_elast_*phi_4*phi_4*(3-2*phi_4) + mu1_elast_*(1-phi_4)*(1-phi_4)*(1+2*phi_4))*strain;		
		
   deltaSigmaD = ((mu1_elast-mu2_elast)*twoSymm(fvc::grad(D)) + (lambda1-lambda2)*(I*tr(fvc::grad(D))) 
		- 2*(mu1_elast-mu2_elast)*(phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 + phi_2*phi_2*(3-2*phi_2)*cEigenStrain2 +phi_3*phi_3*(3-2*phi_3)*cEigenStrain3) 
		- (lambda1-lambda2)*(I*tr((phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 + phi_2*phi_2*(3-2*phi_2)*cEigenStrain2 +phi_3*phi_3*(3-2*phi_3)*cEigenStrain3))))
		+ (mu1_elast_-mu2_elast_)*strain;
}
}


if (phases == 4)
{
if (components == 2)
{
//For phi_1
    
    fvScalarMatrix phi_1Eqn
    (
        dimt*omega*epsilon*fvm::ddt(phi_1) -2.0*gamma_1*epsilon*dimx*dimx*fvm::laplacian(phi_1) == -(18.0*gamma_1*phi_1*(1.0 -phi_1)*(1.0 -2.0*phi_1))/(epsilon)
        - ((-1.0/4.0*(mu_1 - B_bin[0])*(mu_1 - B_bin[0])/A[0] + D_bin[0])*(6.0*phi_1*(1-phi_1) + 2*(phi_2*phi_3 + phi_2*phi_4 + phi_3*phi_4)) + (-1.0/4.0*(mu_1 - B_bin[1])*(mu_1 - B_bin[1])/A[1] + D_bin[1])*2*phi_2*(phi_3 + phi_4) + (-1.0/4.0*(mu_1 - B_bin[2])*(mu_1 - B_bin[2])/A[2] + D_bin[2])*2*phi_3*(phi_2 + phi_4) + (-1.0/4.0*(mu_1 - B_bin[int(phases)-1])*(mu_1 - B_bin[int(phases)-1])/A[int(phases)-1] + D_bin[int(phases)-1])*2*phi_4*(phi_2 + phi_3))/Vm
             - (6.0*phi_1*(1-phi_1)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 - phi_2*phi_2*(3-2*phi_2)*cEigenStrain2
                                                                    - phi_3*phi_3*(3-2*phi_3)*cEigenStrain3))-(Sigma && cEigenStrain1)))
    );
//    phi_1Eqn -= 6.0*phi_1*(1-phi_1)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 - phi_2*phi_2*(3-2*phi_2)*cEigenStrain2 
//								    - phi_3*phi_3*(3-2*phi_3)*cEigenStrain3))-(Sigma && cEigenStrain1));

    
//For phi_2        
    
    fvScalarMatrix phi_2Eqn
    (
        dimt*omega*epsilon*fvm::ddt(phi_2) -2.0*gamma_2*epsilon*dimx*dimx*fvm::laplacian(phi_2) == -(18.0*gamma_2*phi_2*(1.0 -phi_2)*(1.0 -2.0*phi_2))/(epsilon)
        - ((-1.0/4.0*(mu_1 - B_bin[0])*(mu_1 - B_bin[0])/A[0] + D_bin[0])*2*phi_1*(phi_3 + phi_4) + (-1.0/4.0*(mu_1 - B_bin[1])*(mu_1 - B_bin[1])/A[1] + D_bin[1])*(6.0*phi_2*(1-phi_2) + 2*(phi_1*phi_3 + phi_1*phi_4 + phi_3*phi_4)) + (-1.0/4.0*(mu_1 - B_bin[2])*(mu_1 - B_bin[2])/A[2] + D_bin[2])*2*phi_3*(phi_1 + phi_4) + (-1.0/4.0*(mu_1 - B_bin[int(phases)-1])*(mu_1 - B_bin[int(phases)-1])/A[int(phases)-1] + D_bin[int(phases)-1])*2*phi_4*(phi_1 + phi_3))/Vm
		- (6.0*phi_2*(1-phi_2)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 - phi_2*phi_2*(3-2*phi_2)*cEigenStrain2
                                                                    - phi_3*phi_3*(3-2*phi_3)*cEigenStrain3))-(Sigma && cEigenStrain2)))
    );
//    phi_2Eqn -= 6.0*phi_2*(1-phi_2)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 - phi_2*phi_2*(3-2*phi_2)*cEigenStrain2 
//								    - phi_3*phi_3*(3-2*phi_3)*cEigenStrain3))-(Sigma && cEigenStrain2));

    
//For phi_3  
    
    fvScalarMatrix phi_3Eqn
    (
        dimt*omega*epsilon*fvm::ddt(phi_3) -2.0*gamma_3*epsilon*dimx*dimx*fvm::laplacian(phi_3) == -(18.0*gamma_3*phi_3*(1.0 -phi_3)*(1.0 -2.0*phi_3))/(epsilon)
        - ((-1.0/4.0*(mu_1 - B_bin[0])*(mu_1 - B_bin[0])/A[0] + D_bin[0])*2*phi_1*(phi_2 + phi_4) + (-1.0/4.0*(mu_1 - B_bin[1])*(mu_1 - B_bin[1])/A[1] + D_bin[1])*2*phi_2*(phi_1 + phi_4) + (-1.0/4.0*(mu_1 - B_bin[2])*(mu_1 - B_bin[2])/A[2] + D_bin[2])*(6.0*phi_3*(1-phi_3) + 2*(phi_1*phi_2 + phi_1*phi_4 + phi_2*phi_4)) + (-1.0/4.0*(mu_1 - B_bin[int(phases)-1])*(mu_1 - B_bin[int(phases)-1])/A[int(phases)-1] + D_bin[int(phases)-1])*2*phi_4*(phi_1 + phi_2))/Vm
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
        - ((-1.0/4.0*(mu_1 - B_bin[0])*(mu_1 - B_bin[0])/A[0] + D_bin[0])*2*phi_1*(phi_2 + phi_3) + (-1.0/4.0*(mu_1 - B_bin[1])*(mu_1 - B_bin[1])/A[1] + D_bin[1])*2*phi_2*(phi_1 + phi_3) + (-1.0/4.0*(mu_1 - B_bin[2])*(mu_1 - B_bin[2])/A[2] + D_bin[2])*2*phi_3*(phi_1 + phi_2) + (-1.0/4.0*(mu_1 - B_bin[int(phases)-1])*(mu_1 - B_bin[int(phases)-1])/A[int(phases)-1] + D_bin[int(phases)-1])*(6.0*phi_4*(1-phi_4) + 2*(phi_1*phi_2 + phi_1*phi_3 + phi_2*phi_3)))/Vm
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
}


