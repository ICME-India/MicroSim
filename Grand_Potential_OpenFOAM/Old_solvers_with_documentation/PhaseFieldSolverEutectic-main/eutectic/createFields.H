//! Creating phase-fields with the option to write
 volScalarField phi_alpha
  (
      IOobject
      (
          "phi_alpha",
          runTime.timeName(),
          mesh,
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
      ),
      mesh
  );
volScalarField phi_beta
  (
      IOobject
      (
          "phi_beta",
          runTime.timeName(),
          mesh,
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
      ),
      mesh
  );
volScalarField phi_liq
  (
      IOobject
      (
          "phi_liq",
          runTime.timeName(),
          mesh,
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
      ),
      mesh
  );
//! Creating temperature field with the option to write
volScalarField T
    (
        IOobject
        (
            "T",
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

#include "readTransportProperties.H"
