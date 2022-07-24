// Random obj(1);
// const scalar randNumber(obj.scalar01());
#include "clock.H"
Random randNum(clock::getTime());

//! Implicit discretization using fvm class for time derivative and laplacian.
//! Explicit discretization using fvc class for divergence.
/*
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
*/
            
//! Grand potential of solidus (alpha) and liquidus for each components
/*
psi_a = (-0.5*dmudc_1[0][0]*(dcdmu_1[0][0]*(mu_1 - B_a1) + dcdmu_1[0][1]*(mu_2 - B_a2))*(dcdmu_1[0][0]*(mu_1 - B_a1) + dcdmu_1[0][1]*(mu_2 - B_a2)) - dmudc_1[0][1]*(dcdmu_1[0][0]*(mu_1 - B_a1) + dcdmu_1[0][1]*(mu_2 - B_a2))*(dcdmu_1[1][0]*(mu_1 - B_a1) + dcdmu_1[1][1]*(mu_2 - B_a2)) - 0.5*dmudc_1[1][1]*(dcdmu_1[1][0]*(mu_1 - B_a1) + dcdmu_1[1][1]*(mu_2 - B_a2))*(dcdmu_1[1][0]*(mu_1 - B_a1) + dcdmu_1[1][1]*(mu_2 - B_a2)) + D_a);

psi_l = (-0.5*dmudc_2[0][0]*(dcdmu_2[0][0]*mu_1 + dcdmu_2[0][1]*mu_2)*(dcdmu_2[0][0]*mu_1 + dcdmu_2[0][1]*mu_2) - dmudc_2[0][1]*(dcdmu_2[0][0]*mu_1 + dcdmu_2[0][1]*mu_2)*(dcdmu_2[1][0]*mu_1 + dcdmu_2[1][1]*mu_2) - 0.5*dmudc_2[1][1]*(dcdmu_2[1][0]*mu_1 + dcdmu_2[1][1]*mu_2)*(dcdmu_2[1][0]*mu_1 + dcdmu_2[1][1]*mu_2));

*/
#include "dAdgradPhiMod.H"

volVectorField grad_theta = dimx*fvc::grad(theta);
volScalarField grad_theta_sqr = magSqr(grad_theta);

//! Implicit discretization using fvm class for time derivative and laplacian.
//! Explicit discretization using fvc class for divergence.
//! Phase-field equation with approximate relation between c, mu and T according to the parabolic approximation for free energy
fvScalarMatrix alphaEqn (
/*omega*epsilon*dimt*fvm::ddt(alpha)
	-epsilon*gamma*dimx*dimx*fvm::laplacian(alpha)
	+2.0*alpha*(alpha-1.0)*(2.0*alpha-1.0)/epsilon		*/

omega*epsilon*dimt*fvm::ddt(phi)     //alpha=1 liquid phase
                ==
// 		2.0*epsilon*gamma*dimx*dimx*fvm::laplacian(ac_01*ac_01,phi) + 2*gamma*epsilon*dimx*fvc::div(dadgradPhi) - 18.0*(gamma/epsilon)*(phi)*(1-phi)*(1-2.0*phi) 
//                 + ((ceq_1[0][1] - ceq_2[0][1])*(mu_1 - mu_eq[0][1]) + (ceq_1[1][1] - ceq_2[1][1])*(mu_2 - mu_eq[1][1]))*30.0*phi*phi*(1.0-phi)*(1.0-phi)
    2.0*epsilon*gamma*dimx*dimx*fvm::laplacian(ac_01*ac_01,phi) + 2*gamma*epsilon*dimx*fvc::div(dadgradPhi) - 18.0*(gamma/epsilon)*(phi)*(1-phi)*(1-2.0*phi) 
    + ((-0.5*dmudc_l[0][0]*(dcdmu_l[0][0]*mu_1 + dcdmu_l[0][1]*mu_2)*(dcdmu_l[0][0]*mu_1 + dcdmu_l[0][1]*mu_2) - dmudc_l[0][1]*(dcdmu_l[0][0]*mu_1 + dcdmu_l[0][1]*mu_2)*(dcdmu_l[1][0]*mu_1 + dcdmu_l[1][1]*mu_2) - 0.5*dmudc_l[1][1]*(dcdmu_l[1][0]*mu_1 + dcdmu_l[1][1]*mu_2)*(dcdmu_l[1][0]*mu_1 + dcdmu_l[1][1]*mu_2)) - (-0.5*dmudc_a[0][0]*(dcdmu_a[0][0]*(mu_1 - B_a1) + dcdmu_a[0][1]*(mu_2 - B_a2))*(dcdmu_a[0][0]*(mu_1 - B_a1) + dcdmu_a[0][1]*(mu_2 - B_a2)) - dmudc_a[0][1]*(dcdmu_a[0][0]*(mu_1 - B_a1) + dcdmu_a[0][1]*(mu_2 - B_a2))*(dcdmu_a[1][0]*(mu_1 - B_a1) + dcdmu_a[1][1]*(mu_2 - B_a2)) - 0.5*dmudc_a[1][1]*(dcdmu_a[1][0]*(mu_1 - B_a1) + dcdmu_a[1][1]*(mu_2 - B_a2))*(dcdmu_a[1][0]*(mu_1 - B_a1) + dcdmu_a[1][1]*(mu_2 - B_a2)) + DD_a))*30.0*phi*phi*(1.0-phi)*(1.0-phi)*(1 + noise_mag*randNum.globalScalar01())/Vm - epsilon*mu_th*mu_th*grad_theta_sqr*((3*(7*phi*phi+8*phi*phi*phi)*(1-phi) + 2*phi*phi*phi*(7-6*phi))/(1+eps-phi))
    //- 6.0*phi*(1-phi)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi*phi*(3-2*phi)*cEigenStrain))-(Sigma && cEigenStrain))
);

alphaEqn.solve();

//deltaF = (6.0*phi*(1-phi)*0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi*phi*(3-2*phi)*cEigenStrain))-(Sigma && cEigenStrain));


Info<< "Min/max phi:" << min(phi).value() << ' '
<< max(phi).value() << endl;

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
