fvScalarMatrix psiEqn_eq
(
omega*(epsilon)*dimt*fvm::ddt(psi)
	-2.0*(epsilon)*gamma*dimx*dimx*fvm::laplacian(psi)
	+18.0*gamma*psi*(psi-1.0)*(2.0*psi-1.0)/(epsilon)		
);

psiEqn_eq.solve();

Info<< "Min/max psi:" << min(psi()).value() << ' '
    << max(psi()).value() << endl;

// volScalarField deltapsi = dimt*fvc::ddt(psi)*runTime.deltaTValue();
//     
// scalar sum;
// scalar sum2;
// scalar sum3;
// 
// sum  = gSum(psi());
// sum2 = gSum((psi*psi)());
// sum3 = gSum(deltapsi());
//     
// 
// scalar lambdas = sum3/(6.0*(sum-sum2));  
// 
// psi -= 6.0*lambdas*psi*(1.0-psi);
