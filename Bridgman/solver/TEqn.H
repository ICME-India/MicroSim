//  fvScalarMatrix TEqn
//   (
//     Cp*dimt*fvm::ddt(T)- (dimx*dimx/(psi + 1e-6))*fvm::laplacian(psi*K, T) == dimx*(mag(fvc::grad(psi))/(psi + 1e-6))*flux*(T)
//   );
// 
//   TEqn.solve();

volScalarField height = pos(mesh.C().component(1)/dimx-(ybottom + 30*50e-6 + runTime.value()*Vpull));

// volScalarField gradpsi2 = (fvc::grad(psi))&(fvc::grad(psi));
phi = (1 - (0.5*Foam::erf(4.0*(T-Tm)/(Tl-Ts))+scalar(0.5)));

volScalarField K = Ks*phi + Kl*(1-phi);
//volScalarField omega = epsilon*L*L*0.18/(Tm*K);


 fvScalarMatrix TEqn
  (
    Cp*dimt*fvm::ddt(T)- (dimx*dimx/(psi + 1e-6))*fvm::laplacian(psi*K, T) == dimx*(mag(fvc::grad(psi))/(psi + 1e-6))*flux*(T-Tinf)*(1-height) - dimx*dimx*((fvc::grad(psi))&(fvc::grad(psi*T) - Bd*fvc::grad(psi)))*height*(K/((psi+1e-6)*(psi+1e-6))) -L*4.0*exp(-pow(4.0*(T-Tm)/(Tl-Ts),2))/Foam::sqrt(pi)/(Tl-Ts)*(dimt*fvm::ddt(T))
  );

  TEqn.solve();
  
Info<< "Min/max T:" << min(T()).value() << ' '
    << max(T()).value() << endl;
