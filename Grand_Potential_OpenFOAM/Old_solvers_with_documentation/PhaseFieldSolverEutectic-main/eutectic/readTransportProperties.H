//! The input properties to be exported to createFields
   IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

//! The input properties are read from constant/transportProperties dictionary
	//dimensionedScalar pi = constant::mathematical::pi;
        //!Dimension of time
        dimensionedScalar dimt(transportProperties.lookup("dimt"));
        //!Dimension of position
        dimensionedScalar dimx(transportProperties.lookup("dimx"));
        //!Slope liq-alpha and alpha-liq
        dimensionedScalar ms_alpha(transportProperties.lookup("ms_alpha"));
        //!Slope liq-beta and beta-liq
        dimensionedScalar ms_beta(transportProperties.lookup("ms_beta"));
        //!Eutectic composition of liquid phase
        dimensionedScalar c_eq_liq(transportProperties.lookup("c_eq_liq"));
        //!Eutectic composition of alpha phase
        dimensionedScalar c_eq_alpha(transportProperties.lookup("c_eq_alpha"));
        //!Eutectic composition of beta phase
        dimensionedScalar c_eq_beta(transportProperties.lookup("c_eq_beta"));
        //!Thermal gradient
        dimensionedScalar G(transportProperties.lookup("G"));
        //!Velocity
        dimensionedScalar v(transportProperties.lookup("v"));
        dimensionedScalar A(transportProperties.lookup("A"));
        //!Diffusivity in liquid
        dimensionedScalar D(transportProperties.lookup("D"));
        //!Eutectic temperature
        dimensionedScalar T_eut(transportProperties.lookup("T_eut"));
        //!Constant value from temperature profile
        dimensionedScalar initial(transportProperties.lookup("initial"));
        //!Relaxation coefficient 
        dimensionedScalar tau(transportProperties.lookup("tau"));
        //!Interface energy parameter
        dimensionedScalar gamma(transportProperties.lookup("gamma"));
        //!Interface width parameter
        dimensionedScalar epsilon(transportProperties.lookup("epsilon"));
