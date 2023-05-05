    //! Elastic stress, strain, displacement fields only for precipitate growth
    iCorr=0;
    
    if (swch == 2)
    {
        //nCorr = 0;
    
    do {
      fvVectorMatrix DEqn (
          dimt*dimt*fvm::d2dt2(D)
        ==
            dimx*dimx*fvm::laplacian(2*(mu1_elast*phi*phi*(3-2*phi) + mu2_elast*(1-phi)*(1-phi)*(1+2*phi)) 
          + lambda1*phi*phi*(3-2*phi) + lambda2*(1-phi)*(1-phi)*(1+2*phi), D, "laplacian(DD,D)")
          + dimx*dimx*divSigmaExp
          - dimx*dimx*fvc::div((2*mu1_elast*phi*phi*(3-2*phi) + 2*mu2_elast*(1-phi)*(1-phi)*(1+2*phi))*phi*phi*(3-2*phi)*cEigenStrain
          + (lambda1*phi*phi*(3-2*phi) + (1-phi)*(1-phi)*(1+2*phi)*lambda2)*I*tr(phi*phi*(3-2*phi)*cEigenStrain))
      );

      InitialResidual_3 = DEqn.solve().max().initialResidual();
    
      gradD  = fvc::grad(D);
      
      strain =((gradD-phi*phi*(3-2*phi)*cEigenStrain)&&symmTensor(1,0,0,0,0,0))*symmTensor(1,0,0,0,0,0)
             +((gradD-phi*phi*(3-2*phi)*cEigenStrain)&&symmTensor(0,0,0,1,0,0))*symmTensor(0,0,0,1,0,0)
             +((gradD-phi*phi*(3-2*phi)*cEigenStrain)&&symmTensor(0,0,0,0,0,1))*symmTensor(0,0,0,0,0,1);

      sigmaD = (mu1_elast*phi*phi*(3-2*phi)  + mu2_elast*(1-phi)*(1-phi)*(1+2*phi))*twoSymm(gradD) 
             + (lambda1*phi*phi*(3-2*phi)    + lambda2*(1-phi)*(1-phi)*(1+2*phi))*(I*tr(gradD))
             + (mu1_elast_*phi*phi*(3-2*phi) + mu2_elast_*(1-phi)*(1-phi)*(1+2*phi))*strain;
      
             
      divSigmaExp = fvc::div
                    (
                        sigmaD - (2*mu1_elast*phi*phi*(3-2*phi) + 2*mu2_elast*(1-phi)*(1-phi)*(1+2*phi)
			       + lambda1*phi*phi*(3-2*phi)      + (1-phi)*(1-phi)*(1+2*phi)*lambda2)*gradD,
                        "div(sigmaD)"
                    );
      
    }while(InitialResidual_3 > convergenceTolerance && ++iCorr < nCorr);
    
    Sigma = (2*(mu1_elast*phi*phi*(3-2*phi) + mu2_elast*(1-phi)*(1-phi)*(1+2*phi))*(symm(fvc::grad(D)) - phi*phi*(3-2*phi)*cEigenStrain) 
      + (lambda1*phi*phi*(3-2*phi)      + lambda2*(1-phi)*(1-phi)*(1+2*phi))*(I*tr(fvc::grad(D) - phi*phi*(3-2*phi)*cEigenStrain)))
      + (mu1_elast_*phi*phi*(3-2*phi)   + mu2_elast_*(1-phi)*(1-phi)*(1+2*phi))*strain;


    deltaSigmaD = ((mu1_elast-mu2_elast)*twoSymm(fvc::grad(D))           + (lambda1-lambda2)*(I*tr(fvc::grad(D))) 
            - 2*(mu1_elast-mu2_elast)*phi*phi*(3-2*phi)*cEigenStrain - (lambda1-lambda2)*(I*tr(phi*phi*(3-2*phi)*cEigenStrain)))
            + (mu1_elast_-mu2_elast_)*strain;
    
    deltaF = (0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi*phi*(3-2*phi)*cEigenStrain))-(Sigma && cEigenStrain));
    
    sigmaEq = sqrt((3.0/2.0)*magSqr(dev(Sigma)));
    
    	Info<< "Min/max Dx:" << min(D.component(0)).value() << ' ' << max(D.component(0)).value() << endl;
	Info<< "Min/max Dy:" << min(D.component(1)).value() << ' ' << max(D.component(1)).value() << endl;
	if (dimensions == 3)
	{
	Info<< "Min/max Dz:" << min(D.component(2)).value() << ' ' << max(D.component(2)).value() << endl;
	}
    }
    
