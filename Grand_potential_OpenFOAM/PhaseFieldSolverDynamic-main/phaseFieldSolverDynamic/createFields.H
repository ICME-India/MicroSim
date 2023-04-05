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
      
  //! Creating chemical potential field with the option to write
    volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
 // The input properties are imported
 #include "readTransportProperties.H"
    
