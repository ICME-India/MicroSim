// scalar sumphi  = 0.0;
// scalar sum2phi = 0.0;
// scalar sum3phi = 0.0;
// Random obj(1);
// const scalar randNumber(obj.scalar01());//For initiating dendritic growth in slightly larger systems
// Random obj(1);
// const scalar randNumber(obj.scalar01());//For initiating dendritic growth in slightly larger systems
#include "clock.H"
Random randNum(clock::getTime());

// #include "dAdgradPhi_mod.H"

// volVectorField grad_theta     = dimx*fvc::grad(theta);
// volScalarField grad_theta_sqr = magSqr(grad_theta);

fvScalarMatrix phiEqn
    (
      omega*epsilon*dimt*fvm::ddt(phi)     //alpha=1 liquid phase
                == (1/(psi + 1e-6))*2.0*epsilon*gamma*dimx*dimx*fvm::laplacian(psi, phi) 
//                 + 2*gamma*epsilon*dimx*fvc::div(psi*dadgradPhi)*(1/(psi + 1e-6)) 
                - 18.0*(gamma/epsilon)*(phi)*(1-phi)*(1-2.0*phi) + L*((Tm-T)/Tm)*30.0*phi*phi*(1.0-phi)*(1.0-phi)
//                 - epsilon*mu_t*mu_t*grad_theta_sqr*((3*(7*phi*phi+8*phi*phi*phi)*(1-phi) + 2*phi*phi*phi*(7-6*phi))/(1+eps-phi))
//                 - epsilon*mu_t*mu_t*grad_theta_sqr*((3*(7*phi*phi+8*phi*phi*phi)*(1-phi) + 2*phi*phi*phi*(7-6*phi))/(1+eps-phi))
//                 + dimx*costheta*(mag(fvc::grad(psi))/(psi + 1e-6))*(6.0*gamma*phi*(1.0-phi)) 
//                 + 6*noise_mag*phi*(1.0-phi)*phi*(1-phi)*randNumber
    );    
    phiEqn.solve();

Info<< "Min/max phi:" << min(phi()).value() << ' '
    << max(phi()).value() << endl;


