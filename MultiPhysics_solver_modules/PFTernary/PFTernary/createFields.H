  //! Creating phase-field with the option to write
   volScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    //! Creating chemical potential fields with the option to write
    volScalarField mu_1
    (
        IOobject
        (
            "mu_1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    volScalarField mu_2
    (
        IOobject
        (
            "mu_2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    //! Creating orientation field with the option to write
    volScalarField theta
    (
        IOobject
        (
            "theta",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
     #include "readTransportProperties.H"
     #include "readTDB.H"
//Elasticity
//! Creating elastic strain fields with the option to write

volVectorField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volTensorField gradD
(
    IOobject
    (
        "gradD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
//    mesh
     fvc::grad(D)
);

volSymmTensorField strain
(
    IOobject
    (
        "strain",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
//    mesh
   ((gradD-phi*phi*(3-2*phi)*cEigenStrain)&&symmTensor(1,0,0,0,0,0))*symmTensor(1,0,0,0,0,0)
  +((gradD-phi*phi*(3-2*phi)*cEigenStrain)&&symmTensor(0,0,0,1,0,0))*symmTensor(0,0,0,1,0,0)
  +((gradD-phi*phi*(3-2*phi)*cEigenStrain)&&symmTensor(0,0,0,0,0,1))*symmTensor(0,0,0,0,0,1)
    //volSymmTensorField("strain", dimensionSet(0,0,0,0,0,0,0), Foam::symmTensor(0,0,0,0,0,0))
   // Foam::symmTensor(0,0,0,0,0,0)
);

volSymmTensorField sigmaD
(
    IOobject
    (
        "sigmaD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
//     (mu1*T+mu2*(1-T))*twoSymm(fvc::grad(D)) + (lambda1*T+lambda2*(1-T))*(I*tr(fvc::grad(D)))
    (mu1_elast*phi*phi*(3-2*phi) + mu2_elast*(1-phi)*(1-phi)*(1+2*phi))*twoSymm(gradD) + (lambda1*phi*phi*(3-2*phi)    + lambda2*(1-phi)*(1-phi)*(1+2*phi))*(I*tr(gradD))
);

volSymmTensorField deltaSigmaD
(
    IOobject
    (
        "deltaSigmaD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
//      ((mu1-mu2)*twoSymm(fvc::grad(D)) + (lambda1-lambda2)*(I*tr(fvc::grad(D))) - 2*(mu1-mu2)*T*cEigenStrain - (lambda1-lambda2)*(I*tr(T*cEigenStrain)))
     ((mu1_elast-mu2_elast)*twoSymm(fvc::grad(D)) + (lambda1-lambda2)*(I*tr(fvc::grad(D))) - 2*(mu1_elast-mu2_elast)*phi*phi*(3-2*phi)*cEigenStrain 
     - (lambda1-lambda2)*(I*tr(phi*phi*(3-2*phi)*cEigenStrain)))
     + (mu1_elast_-mu2_elast_)*strain
);


volSymmTensorField Sigma
(
    IOobject
    (
        "Sigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
//     (2*(mu1*T+mu2*(1-T))*(0.5*twoSymm(fvc::grad(D)) - T*cEigenStrain) + (lambda1*T+lambda2*(1-T))*(I*tr(fvc::grad(D) - T*cEigenStrain)))
    (2*(mu1_elast*phi*phi*(3-2*phi) + mu2_elast*(1-phi)*(1-phi)*(1+2*phi))*(symm(fvc::grad(D)) - phi*phi*(3-2*phi)*cEigenStrain) + (lambda1*phi*phi*(3-2*phi) + lambda2*(1-phi)*(1-phi)*(1+2*phi))*(I*tr(fvc::grad(D) - phi*phi*(3-2*phi)*cEigenStrain))) + (mu1_elast_*phi*phi*(3-2*phi) + mu2_elast_*(1-phi)*(1-phi)*(1+2*phi))*strain
);

 volScalarField deltaF
    (
        IOobject
        (
            "deltaF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
//         mesh
        6.0*phi*(1-phi)*(0.5*(deltaSigmaD && (symm(fvc::grad(D))-phi*phi*(3-2*phi)*cEigenStrain))-(Sigma && cEigenStrain))
    );


Info<< "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
volVectorField divSigmaExp
(
    IOobject
    (
        "divSigmaExp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(sigmaD)
);

volScalarField sigmaEq
(
    IOobject
    (
        "sigmaEq",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sqrt((3.0/2.0)*magSqr(dev(Sigma)))
);


mesh.setFluxRequired(D.name());
    
    
    

    
