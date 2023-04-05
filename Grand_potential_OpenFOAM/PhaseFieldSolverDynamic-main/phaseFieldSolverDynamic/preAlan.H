//! The interpolation equation for smoothening of the phase-field variable.
//! Implicit discretization using fvm class for time derivative and laplacian.
fvScalarMatrix alpha1Eqn
(
omega*epsilon*dimt*fvm::ddt(phi)
	-epsilon*gamma*dimx*dimx*fvm::laplacian(phi)
	+2.0*phi*(phi-1.0)*(2.0*phi-1.0)/(epsilon)		

/*omega*epsilon*dimt*fvm::ddt(alpha)     //alpha=1 liquid phase
                ==
		2.0*epsilon*gamma*dimx*dimx*fvm::laplacian(alpha) - 18.0*(gamma/epsilon)*(alpha)*(1-alpha)*(1-2.0*alpha)
               // -(c_Sol-c_Liq)*(mu - (2*A*c_eq - 2*A*(T0-T)/m_1))*30.0*alpha*alpha*(1.0-alpha)*(1.0-alpha)*/

		
);

alpha1Eqn.solve();


/*fvScalarMatrix muEqn
(
    (0.5)*dimt*fvm::ddt(mu) == diff_Liq*0.5*dimx*dimx*fvm::laplacian(alpha,mu) + (c_Sol-c_Liq)* dimt*fvc::ddt(alpha)*30.0*alpha*alpha*(1.0-alpha)*(1.0-alpha)
);

muEqn.solve();
*/

    // Writing into the log file
    Info<< "Min/max phi:" << min(phi).value() << ' '
    << max(phi).value() << endl;

