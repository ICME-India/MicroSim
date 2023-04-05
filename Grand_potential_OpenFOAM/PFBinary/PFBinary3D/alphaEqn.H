// The phase-field and chemical potential equations to be used in main
//! Random number generation for adding noise to induce side branching
#include "clock.H"
Random randNum(clock::getTime());
//Random obj(1);
//const scalar randNumber(obj.scalar01());
// Header file to obtain some specific terms in the phase-field equation
//
//

qt0 = (Foam::cos(0.5*thetax)*Foam::cos(0.5*thetay)*Foam::cos(0.5*theta) + Foam::sin(0.5*thetax)*Foam::sin(0.5*thetay)*Foam::sin(0.5*theta));
qt1 = (Foam::sin(0.5*thetax)*Foam::cos(0.5*thetay)*Foam::cos(0.5*theta) - Foam::cos(0.5*thetax)*Foam::sin(0.5*thetay)*Foam::sin(0.5*theta));
qt2 = (Foam::cos(0.5*thetax)*Foam::sin(0.5*thetay)*Foam::cos(0.5*theta) + Foam::sin(0.5*thetax)*Foam::cos(0.5*thetay)*Foam::sin(0.5*theta));
qt3 = (Foam::cos(0.5*thetax)*Foam::cos(0.5*thetay)*Foam::sin(0.5*theta) - Foam::sin(0.5*thetax)*Foam::sin(0.5*thetay)*Foam::cos(0.5*theta));

volVectorField grad_qt0 = dimx*fvc::grad(qt0);
volVectorField grad_qt1 = dimx*fvc::grad(qt1);
volVectorField grad_qt2 = dimx*fvc::grad(qt2);
volVectorField grad_qt3 = dimx*fvc::grad(qt3);

volScalarField grad_qt_sqr = magSqr(grad_qt0) + magSqr(grad_qt1) + magSqr(grad_qt2) + magSqr(grad_qt3);

//
//#include "dAdgradPhiMod.H"
#include "dAdgradPhiMod3D.H"

//volVectorField grad_theta = dimx*fvc::grad(theta);
//volScalarField grad_theta_sqr = magSqr(grad_theta);

if (swch == 2)
    {
strain = ((gradD-phi*phi*(3-2*phi)*cEigenStrain)&&symmTensor(1,0,0,0,0,0))*symmTensor(1,0,0,0,0,0)
       + ((gradD-phi*phi*(3-2*phi)*cEigenStrain)&&symmTensor(0,0,0,1,0,0))*symmTensor(0,0,0,1,0,0)
       + ((gradD-phi*phi*(3-2*phi)*cEigenStrain)&&symmTensor(0,0,0,0,0,1))*symmTensor(0,0,0,0,0,1);

Sigma = (2*(mu1_elast*phi*phi*(3-2*phi) + mu2_elast*(1-phi)*(1-phi)*(1+2*phi))*(symm(fvc::grad(D)) - phi*phi*(3-2*phi)*cEigenStrain) 
      + (lambda1*phi*phi*(3-2*phi)      + lambda2*(1-phi)*(1-phi)*(1+2*phi))*(I*tr(fvc::grad(D) - phi*phi*(3-2*phi)*cEigenStrain)))
      + (mu1_elast_*phi*phi*(3-2*phi)   + mu2_elast_*(1-phi)*(1-phi)*(1+2*phi))*strain;

//Sigma=Sigma-T*cEigenStrain;

deltaSigmaD = ((mu1_elast-mu2_elast)*twoSymm(fvc::grad(D))           + (lambda1-lambda2)*(I*tr(fvc::grad(D))) 
            - 2*(mu1_elast-mu2_elast)*phi*phi*(3-2*phi)*cEigenStrain - (lambda1-lambda2)*(I*tr(phi*phi*(3-2*phi)*cEigenStrain)))
            + (mu1_elast_-mu2_elast_)*strain;
}            

//! Implicit discretization using fvm class for time derivative and laplacian.
//! Explicit discretization using fvc class for divergence.
//! Phase-field equation with approximate relation between c, mu and T according to the parabolic approximation for free energy
fvScalarMatrix alphaEqn
(
/*omega*epsilon*dimt*fvm::ddt(alpha)
	-epsilon*gamma*dimx*dimx*fvm::laplacian(alpha)
	+2.0*alpha*(alpha-1.0)*(2.0*alpha-1.0)/epsilon		*/

/*omega*epsilon*dimt*fvm::ddt(phi)     //alpha=1 liquid phase
                ==
		2.0*epsilon*gamma*dimx*dimx*fvm::laplacian(ac_01*ac_01,phi) + 2*gamma*epsilon*dimx*fvc::div(dadgradPhi) - 18.0*(gamma/epsilon)*(phi)*(1-phi)*(1-2.0*phi) + ((-1.0/4.0*(mu - B_Liq)*(mu - B_Liq)/A_Liq + D_Liq) - (-1.0/4.0*(mu - B_Sol)*(mu - B_Sol)/A_Sol + D_Sol))*(-60)*phi*phi*(-1+3*phi-3*phi*phi+phi*phi*phi)/Vm - gamma*epsilon*mu_th*mu_th*grad_theta_sqr*((3*(7*phi*phi+8*phi*phi*phi)*(1-phi) + 2*phi*phi*phi*(7-6*phi))/(1+eps-phi)) + 6*noise_mag*phi*(1.0-phi)*phi*(1-phi)*randNum.globalScalar01()

*/

omega*epsilon*dimt*fvm::ddt(phi)     //alpha=1 liquid phase
                ==
		2.0*epsilon*gamma*dimx*dimx*fvm::laplacian(ac_01*ac_01,phi) + 2*gamma*epsilon*dimx*fvc::div(dadgradPhi) - 18.0*(gamma/epsilon)*(phi)*(1-phi)*(1-2.0*phi) + ((-1.0/4.0*(mu - B_Liq)*(mu - B_Liq)/A_Liq + D_Liq) - (-1.0/4.0*(mu - B_Sol)*(mu - B_Sol)/A_Sol + D_Sol))*(30.0*phi*phi*(1.0-phi)*(1.0-phi))*(1 + noise_mag*randNum.globalScalar01())/Vm - epsilon*mu_th*mu_th*grad_qt_sqr*((3*(7*phi*phi+8*phi*phi*phi)*(1-phi) + 2*phi*phi*phi*(7-6*phi))/(1+eps-phi)) - 6.0*phi*(1-phi)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi*phi*(3-2*phi)*cEigenStrain))-(Sigma && cEigenStrain))
//+ 30.0*noise_mag*phi*(1.0-phi)*phi*(1-phi)*randNum.globalScalar01()


);

//alphaEqn.solve();

InitialResidual_0 = alphaEqn.solve().max().initialResidual();

    // Writing into the log file
    Info<< "Min/max phi:" << min(phi).value() << ' ' << max(phi).value() << endl;

//! The unit normal vector to the interface with a small number in denominator to prevent solution from diverging
volVectorField n=dimx*fvc::grad(phi)/(1E-20+mag(dimx*fvc::grad(phi)));

volScalarField hphi = phi*phi*phi*(6*phi*phi - 15*phi +10); //-phi*phi*phi*(10*phi*phi*phi - 36*phi*phi + 45*phi - 20);

//! Chemical potential equation with approximate slope of c-mu curve according to the parabolic approximation for free energy
fvScalarMatrix muEqn
(
    (0.5)*((hphi/A_Sol) + ((1-hphi)/A_Liq))*dimt*fvm::ddt(mu) == dimx*dimx*fvm::laplacian(diff_Sol*0.5*(1/A_Sol)*phi + diff_Liq*0.5*(1/A_Liq)*(1-phi), mu) - 0.5*((mu - B_Sol)/A_Sol - (mu - B_Liq)/A_Liq)*dimt*fvc::ddt(phi)*30.0*phi*phi*(1.0-phi)*(1.0-phi) - Tdot*(((-dBdT)/(2*A_Sol) + (mu - B_Sol)*A_SoldT/(-2*A_Sol*A_Sol))*hphi + (mu - B_Liq)*A_LiqdT*(1-hphi)/(-2*A_Liq*A_Liq)) - anti_trap*epsilon*0.5*((mu - B_Sol)/A_Sol - (mu - B_Liq)/A_Liq)*dimx*fvc::div((n*dimt*fvc::ddt(phi)))
);

//muEqn.solve();

InitialResidual_1 = muEqn.solve().max().initialResidual();

    Info<< "Min/max mu:" << min(mu).value() << ' ' << max(mu).value() << endl;

//! Orientation field equation is solved only for multigrain problems
if (swch == 1)
{
    
volScalarField g_phi = phi*phi*phi*(7-6*phi)/((1+eps-phi)*(1+eps-phi));

volScalarField M_theta_func = M_Scale*epsilon*(1-phi*phi*phi*(10-15*phi+6*phi*phi))*(1-phi)*(1-phi);

/*forAll(grad_theta_sqr, i)
{
    Sout << Pstream::myProcNo << " " << g_phi[i] << " " << M_theta_func[i] << endl;
}*/

//! Orientation field evolution equation
fvScalarMatrix qt0Eqn
(
    dimt*fvm::ddt(qt0)  == 2*M_theta_func*mu_th*mu_th*(dimx*dimx*fvm::laplacian(g_phi,qt0))
);

//thetaEqn.solve();

InitialResidual_20 = qt0Eqn.solve().max().initialResidual();

fvScalarMatrix qt1Eqn
(
    dimt*fvm::ddt(qt1)  == 2*M_theta_func*mu_th*mu_th*(dimx*dimx*fvm::laplacian(g_phi,qt1))
);

//thetaEqn.solve();

InitialResidual_21 = qt1Eqn.solve().max().initialResidual();

fvScalarMatrix qt2Eqn
(
    dimt*fvm::ddt(qt2)  == 2*M_theta_func*mu_th*mu_th*(dimx*dimx*fvm::laplacian(g_phi,qt2))
);

//thetaEqn.solve();

InitialResidual_22 = qt2Eqn.solve().max().initialResidual();

fvScalarMatrix qt3Eqn
(
    dimt*fvm::ddt(qt3)  == 2*M_theta_func*mu_th*mu_th*(dimx*dimx*fvm::laplacian(g_phi,qt3))
);

//thetaEqn.solve();

InitialResidual_23 = qt3Eqn.solve().max().initialResidual();


volScalarField lambda = qt0*dimt*fvc::ddt(qt0) + qt1*dimt*fvc::ddt(qt1) + qt2*dimt*fvc::ddt(qt2) + qt3*dimt*fvc::ddt(qt3);

    qt0 -= lambda*qt0*runTime.deltaTValue();
    qt1 -= lambda*qt1*runTime.deltaTValue();
    qt2 -= lambda*qt2*runTime.deltaTValue();
    qt3 -= lambda*qt3*runTime.deltaTValue();

    // Writing into the log file
    Info<< "Min/max qt0:" << min(qt0).value() << ' ' << max(qt0).value() << endl;
    Info<< "Min/max qt1:" << min(qt1).value() << ' ' << max(qt1).value() << endl;
    Info<< "Min/max qt2:" << min(qt2).value() << ' ' << max(qt2).value() << endl;
    Info<< "Min/max qt3:" << min(qt3).value() << ' ' << max(qt3).value() << endl;

thetax = Foam::atan2(2*(qt2*qt3 + qt0*qt1), (qt0*qt0 - qt1*qt1 - qt2*qt2 + qt3*qt3));
thetay = Foam::asin(2*(qt0*qt2 - qt1*qt3));
theta = Foam::atan2(2*(qt1*qt2 + qt0*qt3), (qt0*qt0 + qt1*qt1 - qt2*qt2 - qt3*qt3));

}

Info<< "Min/max thetax:" << min(thetax).value() << ' ' << max(thetax).value() << endl;
Info<< "Min/max thetay:" << min(thetay).value() << ' ' << max(thetay).value() << endl;
Info<< "Min/max theta:" << min(theta).value() << ' ' << max(theta).value() << endl;
