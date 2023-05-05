
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

//! For 2D simulations
if (dimensions == 2)
{
theta = qt0*vector(0,0,1);
}

//! For 3D simulations
if (dimensions == 3)
{
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

theta = Foam::atan2(2*(qt2*qt3 + qt0*qt1), (qt0*qt0 - qt1*qt1 - qt2*qt2 + qt3*qt3))*vector(1,0,0)
 + Foam::asin(2*(qt0*qt2 - qt1*qt3))*vector(0,1,0)
  + Foam::atan2(2*(qt1*qt2 + qt0*qt3), (qt0*qt0 + qt1*qt1 - qt2*qt2 - qt3*qt3))*vector(0,0,1);
}

}

    // Writing into the log file
    //! For 3D simulations
if (dimensions == 3)
{
Info<< "Min/max thetax:" << min(theta.component(0)).value() << ' ' << max(theta.component(0)).value() << endl;
Info<< "Min/max thetay:" << min(theta.component(1)).value() << ' ' << max(theta.component(1)).value() << endl;
}
Info<< "Min/max thetaz:" << min(theta.component(2)).value() << ' ' << max(theta.component(2)).value() << endl;

