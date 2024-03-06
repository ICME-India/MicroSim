// volScalarField de = 1.0/(psi*(1-psi)+1e-6);

// volScalarField e = psi*(1-psi);

volVectorField normal = fvc::grad(psi)/(mag(fvc::grad(psi))+1.0e-9/dimx);

fvScalarMatrix psiEqn
(
//  	dimt*(fvm::ddt(gamma)
//         - (massDepRate/rhoMetal)*(mag(fvc::grad(gamma)))
// 	)
  omega*dimt*fvm::ddt(psi)
 	==
	gamma*
	(
		2.0*dimx*dimx*fvm::laplacian(psi)
    - 2.0*dimx*dimx*mag(fvc::grad(psi))*fvc::div(normal)
		- (18.0/(epsilon*epsilon))*(1.0-2*psi)*psi*(1.0-psi)
    
//     - 2.0*dimx*(3*psi*(1-psi)/epsilon)*fvc::div(fvc::grad(psi)/(mag(fvc::grad(psi))+1.0e-9/dimx))
//     - 2.0*dimx*dimx*e*fvm::laplacian(de,psi)
	)
);
InitialResidual = psiEqn.solve().max().initialResidual();

// deltapsi = dimt*fvc::ddt(psi)*runTime.deltaTValue();
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
// lambdas = sum3/(6.0*(sum-sum2));  
// 
// psi -= 6.0*lambdas*psi*(1.0-psi);
