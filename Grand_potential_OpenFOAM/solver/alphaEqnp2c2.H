// The phase-field and chemical potential equations to be used in main
//! Random number generation for adding noise to induce side branching
#include "clock.H"
Random randNum(clock::getTime());
//Random obj(1);
//const scalar randNumber(obj.scalar01());
// Header file to obtain some specific terms in the phase-field equation
//
//
if (phases == 2)
{
//! For 2D simulations
if (dimensions == 2)
{
qt0 = theta.component(2);

volVectorField grad_qt0 = dimx*fvc::grad(qt0);

volScalarField grad_qt_sqr = magSqr(grad_qt0);
}

//! For 3D simulations
if (dimensions == 3)
{
qt0 = (Foam::cos(0.5*theta.component(0))*Foam::cos(0.5*theta.component(1))*Foam::cos(0.5*theta.component(2)) + Foam::sin(0.5*theta.component(0))*Foam::sin(0.5*theta.component(1))*Foam::sin(0.5*theta.component(2)));
qt1 = (Foam::sin(0.5*theta.component(0))*Foam::cos(0.5*theta.component(1))*Foam::cos(0.5*theta.component(2)) - Foam::cos(0.5*theta.component(0))*Foam::sin(0.5*theta.component(1))*Foam::sin(0.5*theta.component(2)));
qt2 = (Foam::cos(0.5*theta.component(0))*Foam::sin(0.5*theta.component(1))*Foam::cos(0.5*theta.component(2)) + Foam::sin(0.5*theta.component(0))*Foam::cos(0.5*theta.component(1))*Foam::sin(0.5*theta.component(2)));
qt3 = (Foam::cos(0.5*theta.component(0))*Foam::cos(0.5*theta.component(1))*Foam::sin(0.5*theta.component(2)) - Foam::sin(0.5*theta.component(0))*Foam::sin(0.5*theta.component(1))*Foam::cos(0.5*theta.component(2)));

volVectorField grad_qt0 = dimx*fvc::grad(qt0);
volVectorField grad_qt1 = dimx*fvc::grad(qt1);
volVectorField grad_qt2 = dimx*fvc::grad(qt2);
volVectorField grad_qt3 = dimx*fvc::grad(qt3);

volScalarField grad_qt_sqr = magSqr(grad_qt0) + magSqr(grad_qt1) + magSqr(grad_qt2) + magSqr(grad_qt3);
}
}

#include "dAdgradPhiMod0.H"

//volVectorField grad_theta = dimx*fvc::grad(theta);
//volScalarField grad_theta_sqr = magSqr(grad_theta);

if (swch == 2)
    {
    if (phases == 2)
    {
strain = ((gradD-phi_1*phi_1*(3-2*phi_1)*cEigenStrain1)&&symmTensor(1,0,0,0,0,0))*symmTensor(1,0,0,0,0,0)
       + ((gradD-phi_1*phi_1*(3-2*phi_1)*cEigenStrain1)&&symmTensor(0,0,0,1,0,0))*symmTensor(0,0,0,1,0,0)
       + ((gradD-phi_1*phi_1*(3-2*phi_1)*cEigenStrain1)&&symmTensor(0,0,0,0,0,1))*symmTensor(0,0,0,0,0,1);

Sigma = (2*(mu1_elast*phi_1*phi_1*(3-2*phi_1) + mu2_elast*(1-phi_1)*(1-phi_1)*(1+2*phi_1))*(symm(fvc::grad(D)) - phi_1*phi_1*(3-2*phi_1)*cEigenStrain1) 
      + (lambda1*phi_1*phi_1*(3-2*phi_1)      + lambda2*(1-phi_1)*(1-phi_1)*(1+2*phi_1))*(I*tr(fvc::grad(D) - phi_1*phi_1*(3-2*phi_1)*cEigenStrain1)))
      + (mu1_elast_*phi_1*phi_1*(3-2*phi_1)   + mu2_elast_*(1-phi_1)*(1-phi_1)*(1+2*phi_1))*strain;

//Sigma=Sigma-T*cEigenStrain;

deltaSigmaD = ((mu1_elast-mu2_elast)*twoSymm(fvc::grad(D))           + (lambda1-lambda2)*(I*tr(fvc::grad(D))) 
            - 2*(mu1_elast-mu2_elast)*phi_1*phi_1*(3-2*phi_1)*cEigenStrain1 - (lambda1-lambda2)*(I*tr(phi_1*phi_1*(3-2*phi_1)*cEigenStrain1)))
            + (mu1_elast_-mu2_elast_)*strain;
}
}

//! Implicit discretization using fvm class for time derivative and laplacian.
//! Explicit discretization using fvc class for divergence.
//! Phase-field equation with approximate relation between c, mu and T according to the parabolic approximation for free energy
if (phases == 2)
{
if (components == 2)
{
fvScalarMatrix alphaEqn
(
/*omega*epsilon*dimt*fvm::ddt(alpha)
	-epsilon*gamma*dimx*dimx*fvm::laplacian(alpha)
	+2.0*alpha*(alpha-1.0)*(2.0*alpha-1.0)/epsilon		*/

/*omega*epsilon*dimt*fvm::ddt(phi)     //alpha=1 liquid phase
                ==
		2.0*epsilon*gamma*dimx*dimx*fvm::laplacian(ac_01*ac_01,phi) + 2*gamma*epsilon*dimx*fvc::div(dadgradPhi) - 18.0*(gamma/epsilon)*(phi)*(1-phi)*(1-2.0*phi) + ((-1.0/4.0*(mu - B_Liq)*(mu - B_Liq)/A_Liq + D_Liq) - (-1.0/4.0*(mu - B_Sol)*(mu - B_Sol)/A_Sol + D_Sol))*(-60)*phi*phi*(-1+3*phi-3*phi*phi+phi*phi*phi)/Vm - gamma*epsilon*mu_th*mu_th*grad_theta_sqr*((3*(7*phi*phi+8*phi*phi*phi)*(1-phi) + 2*phi*phi*phi*(7-6*phi))/(1+eps-phi)) + 6*noise_mag*phi*(1.0-phi)*phi*(1-phi)*randNum.globalScalar01()

*/

omega*epsilon*dimt*fvm::ddt(phi_1)     //alpha=1 liquid phase
                ==
		2.0*epsilon*gamma_1*dimx*dimx*fvm::laplacian(ac_01*ac_01,phi_1) + 2*gamma_1*epsilon*dimx*fvc::div(dadgradPhi) - 18.0*(gamma_1/epsilon)*(phi_1)*(1-phi_1)*(1-2.0*phi_1) + ((-1.0/4.0*(mu_1 - B_bin[int(phases)-1])*(mu_1 - B_bin[int(phases)-1])/A[int(phases)-1] + D_bin[int(phases)-1]) - (-1.0/4.0*(mu_1 - B_bin[0])*(mu_1 - B_bin[0])/A[0] + D_bin[0]))*(30.0*phi_1*phi_1*(1.0-phi_1)*(1.0-phi_1))*(1 + noise_mag*randNum.globalScalar01())/Vm - epsilon*mu_th*mu_th*grad_qt_sqr*((3*(7*phi_1*phi_1+8*phi_1*phi_1*phi_1)*(1-phi_1) + 2*phi_1*phi_1*phi_1*(7-6*phi_1))/(1+eps-phi_1)) - 6.0*phi_1*(1-phi_1)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi_1*phi_1*(3-2*phi_1)*cEigenStrain1))-(Sigma && cEigenStrain1))
//+ 30.0*noise_mag*phi*(1.0-phi)*phi*(1-phi)*randNum.globalScalar01()


);

//alphaEqn.solve();

InitialResidual_0 = alphaEqn.solve().max().initialResidual();
}

}


