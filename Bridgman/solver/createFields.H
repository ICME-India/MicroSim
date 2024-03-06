//     Info<< "Reading field phi\n" << endl;

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
//     volScalarField phiOld
//     (
//         IOobject
//         (
//             "phiOld",
//             runTime.timeName(),
//             mesh,
//             IOobject::NO_READ,
//             IOobject::NO_WRITE
//         ),
//         mesh
//     );
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
    volScalarField psi
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
//     volScalarField theta
//     (
//         IOobject
//         (
//             "theta",
//             runTime.timeName(),
//             mesh,
//             IOobject::MUST_READ,
//             IOobject::AUTO_WRITE
//         ),
//         mesh
//     );


    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );


    Info<< "Reading Conductivity\n" << endl;

//     dimensionedScalar Dphi
//     (
//         transportProperties.lookup("K")
//     );

//     dimensionedScalar omega
//     (
//         transportProperties.lookup("omega")
//     );
    dimensionedScalar pi = constant::mathematical::pi;
    dimensionedScalar dimt(transportProperties.lookup("dimt"));             //dimension of time
    dimensionedScalar dimx(transportProperties.lookup("dimx"));             //diemsion of position
    dimensionedScalar omega(transportProperties.lookup("omega"));           //relaxation coefficient for alpha (order parameter)
    dimensionedScalar gamma(transportProperties.lookup("gamma"));           //Surface Energy
    dimensionedScalar epsilon(transportProperties.lookup("epsilon"));       //Interface Width
    dimensionedScalar flux(transportProperties.lookup("flux"));             //Flux
    dimensionedScalar Cp(transportProperties.lookup("Cp"));                 //Specific heat per unit volume
    dimensionedScalar L(transportProperties.lookup("L"));                   //Latent heat per unit volume
    dimensionedScalar Ks(transportProperties.lookup("Ks"));                   //Thermal Conductivity(solid)
    dimensionedScalar Kl(transportProperties.lookup("Kl"));                   //Thermal Conductivity(liquid)
    dimensionedScalar Vpull(transportProperties.lookup("Vpull"));                   //Thermal Conductivity
    dimensionedScalar Tm(transportProperties.lookup("Tm"));                 //Melting point
    dimensionedScalar Tl(transportProperties.lookup("Tl"));                 //Liquidus 
    dimensionedScalar Ts(transportProperties.lookup("Ts"));                 //Solidus
    dimensionedScalar Tinf(transportProperties.lookup("Tinf"));                 //Ambient temperature
    dimensionedScalar Bd(transportProperties.lookup("Bd"));                 //Dirichlet boundary temperature
//     dimensionedScalar delta_01(transportProperties.lookup("delta_01"));     //strength of anisotropy
    //dimensionedScalar A(transportProperties.lookup("A"));                   //G=Ac^2;d2G/dc^2=A
//     dimensionedScalar mu_t(transportProperties.lookup("mu_t")); 
    //dimensionedScalar M_theta(transportProperties.lookup("M_theta")); 
//     dimensionedScalar eps(transportProperties.lookup("eps"));
    //!M_theta function scaling with time step
//     scalar M_Scale(readScalar(transportProperties.lookup("M_Scale")));
//     dimensionedScalar noise_mag(transportProperties.lookup("noise_mag"));   //noise magnitude
    
    scalar costheta;
//     scalar center_x;
//     scalar center_y;
// //     scalar delta;
//     scalar radius_psi;
    scalar seedRadius;
    scalar xbottom;
    scalar ybottom;
    scalar Height;
    scalar Width;
    
    costheta = readScalar(transportProperties.lookup("costheta"));          //cos of the contact angle
//     center_x = readScalar(transportProperties.lookup("center_x"));          //Center of the domain
//     center_y = readScalar(transportProperties.lookup("center_y"));          //Center of the domain
//     delta    = readScalar(transportProperties.lookup("delta"));             //Width of interface
//     radius_psi = readScalar(transportProperties.lookup("delta"));           //radius of the powder particle
    scalar Tn(readScalar(transportProperties.lookup("Nucleation_undercooling")));
    scalar nucleation_interval(readScalar(transportProperties.lookup("Nucleation_Interval")));
    seedRadius           = readScalar(transportProperties.lookup("seedRadius"));
    xbottom              = readScalar(transportProperties.lookup("xbottom"));
    ybottom              = readScalar(transportProperties.lookup("ybottom"));
    Height               = readScalar(transportProperties.lookup("Height"));
    Width                = readScalar(transportProperties.lookup("Width"));
//     int numSeeds(readScalar(transportProperties.lookup("numSeeds")));
