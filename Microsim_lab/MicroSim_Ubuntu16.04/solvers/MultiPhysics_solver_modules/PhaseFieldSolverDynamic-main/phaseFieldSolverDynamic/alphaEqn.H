// The phase-field and chemical potential equations to be used in main
//! Random number generation for adding noise to induce side branching
Random obj(1);
const scalar randNumber(obj.scalar01());
// Header file to obtain a specific terms in the phase-field equation
#include "dAdgradPhiMod.H"
//! Implicit discretization using fvm class for time derivative and laplacian.
//! Explicit discretization using fvc class for divergence.
//! Phase-field equation with approximate relation between c, mu and T according to the parabolic approximation for free energy
fvScalarMatrix alphaEqn
(
/*omega*epsilon*dimt*fvm::ddt(alpha)
	-epsilon*gamma*dimx*dimx*fvm::laplacian(alpha)
	+2.0*alpha*(alpha-1.0)*(2.0*alpha-1.0)/epsilon		*/

omega*epsilon*dimt*fvm::ddt(phi)     //alpha=1 liquid phase
                ==
		2.0*epsilon*gamma*dimx*dimx*fvm::laplacian(ac_01*ac_01,phi) + 2*gamma*epsilon*dimx*fvc::div(dadgradPhi) - 18.0*(gamma/epsilon)*(phi)*(1-phi)*(1-2.0*phi) + (c_Sol-c_Liq)*(mu - (2*A*c_eq - 2*A*(T0-T)/m_1))*30.0*phi*phi*(1.0-phi)*(1.0-phi) + 6*noise_mag*phi*(1.0-phi)*phi*(1-phi)*randNumber


);

alphaEqn.solve();

//! The unit normal vector to the interface with a small number in denominator to prevent solution from diverging
volVectorField n=dimx*fvc::grad(phi)/(1E-20+mag(dimx*fvc::grad(phi)));

//! Chemical potential equation with approximate slope of c-mu curve according to the parabolic approximation for free energy
fvScalarMatrix muEqn
(
    (0.5)*(1/A)*dimt*fvm::ddt(mu) == diff_Liq*0.5*(1/A)*dimx*dimx*fvm::laplacian((1-phi),mu) - (c_Sol-c_Liq)* dimt*fvc::ddt(phi)*30.0*phi*phi*(1.0-phi)*(1.0-phi) - anti_trap*epsilon*(c_Sol-c_Liq)*dimx*fvc::div((n*dimt*fvc::ddt(phi)))
);

muEqn.solve();


    // Writing into the log file
    Info<< "Min/max phi:" << min(phi).value() << ' '
    << max(phi).value() << endl;

