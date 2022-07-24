// The phase-field and chemical potential equations to be used in main
//! Random number generation for adding noise to induce side branching
#include "clock.H"
Random randNum(clock::getTime());
//Random obj(1);
//const scalar randNumber(obj.scalar01());
// Header file to obtain some specific terms in the phase-field equation
//
//
//
#include "dAdgradPhiMod.H"

volVectorField grad_theta = dimx*fvc::grad(theta);
volScalarField grad_theta_sqr = magSqr(grad_theta);


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
		2.0*epsilon*gamma*dimx*dimx*fvm::laplacian(ac_01*ac_01,phi) + 2*gamma*epsilon*dimx*fvc::div(dadgradPhi) - 18.0*(gamma/epsilon)*(phi)*(1-phi)*(1-2.0*phi) + ((-1.0/4.0*(mu - B_Liq)*(mu - B_Liq)/A_Liq + D_Liq) - (-1.0/4.0*(mu - B_Sol)*(mu - B_Sol)/A_Sol + D_Sol))*(30.0*phi*phi*(1.0-phi)*(1.0-phi))*(1 + noise_mag*randNum.globalScalar01())/Vm - epsilon*mu_th*mu_th*grad_theta_sqr*((3*(7*phi*phi+8*phi*phi*phi)*(1-phi) + 2*phi*phi*phi*(7-6*phi))/(1+eps-phi))
//+ 30.0*noise_mag*phi*(1.0-phi)*phi*(1-phi)*randNum.globalScalar01()


);

alphaEqn.solve();

    // Writing into the log file
    Info<< "Min/max phi:" << min(phi).value() << ' ' << max(phi).value() << endl;

//! The unit normal vector to the interface with a small number in denominator to prevent solution from diverging
volVectorField n=dimx*fvc::grad(phi)/(1E-20+mag(dimx*fvc::grad(phi)));

volScalarField hphi = phi*phi*phi*(6*phi*phi - 15*phi +10); //-phi*phi*phi*(10*phi*phi*phi - 36*phi*phi + 45*phi - 20);

//! Chemical potential equation with approximate slope of c-mu curve according to the parabolic approximation for free energy
fvScalarMatrix muEqn
(
    (0.5)*((hphi/A_Sol) + ((1-hphi)/A_Liq))*dimt*fvm::ddt(mu) == diff_Liq*0.5*(1/A_Liq)*dimx*dimx*fvm::laplacian((1-phi),mu) - 0.5*((mu - B_Sol)/A_Sol - (mu - B_Liq)/A_Liq)*dimt*fvc::ddt(phi)*30.0*phi*phi*(1.0-phi)*(1.0-phi) - Tdot*(((-dBdT)/(2*A_Sol) + (mu - B_Sol)*A_SoldT/(-2*A_Sol*A_Sol))*hphi + (mu - B_Liq)*A_LiqdT*(1-hphi)/(-2*A_Liq*A_Liq)) - anti_trap*epsilon*0.5*((mu - B_Sol)/A_Sol - (mu - B_Liq)/A_Liq)*dimx*fvc::div((n*dimt*fvc::ddt(phi)))
);

muEqn.solve();

    Info<< "Min/max mu:" << min(mu).value() << ' ' << max(mu).value() << endl;

if (swch != 0)
{
    
volScalarField g_phi = phi*phi*phi*(7-6*phi)/((1+eps-phi)*(1+eps-phi));

volScalarField M_theta_func = M_Scale*epsilon*(1-phi*phi*phi*(10-15*phi+6*phi*phi))*(1-phi)*(1-phi);

/*forAll(grad_theta_sqr, i)
{
    Sout << Pstream::myProcNo << " " << g_phi[i] << " " << M_theta_func[i] << endl;
}*/

//! Orientation field evolution equation
fvScalarMatrix thetaEqn
(
    dimt*fvm::ddt(theta)  == 2*M_theta_func*mu_th*mu_th*(dimx*dimx*fvm::laplacian(g_phi,theta))
);

thetaEqn.solve();
}

    // Writing into the log file
    Info<< "Min/max theta:" << min(theta).value() << ' ' << max(theta).value() << endl;

