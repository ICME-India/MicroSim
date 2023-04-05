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
    
   IOdictionary MechProperties
    (
        IOobject
        (
            "MechProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
//! The input properties are read from constant/transportProperties dictionary
	dimensionedScalar pi = constant::mathematical::pi;
	    //!Dimension of time
        dimensionedScalar dimt(transportProperties.lookup("dimt"));
        //!Dimension of position
        dimensionedScalar dimx(transportProperties.lookup("dimx"));
	    //!Number of dimensions
        dimensionedScalar dimensions(transportProperties.lookup("dimensions"));
        //!Number of phases
        dimensionedScalar phases(transportProperties.lookup("phases"));
        //!Number of components
        dimensionedScalar components(transportProperties.lookup("components"));
        //!Slope liquidus
        //dimensionedScalar m_1(transportProperties.lookup("m_1"));
        //!Slope solidus
        //dimensionedScalar m_0(transportProperties.lookup("m_0"));
        //Relaxation coefficient for phi (order parameter)
        //dimensionedScalar omega(transportProperties.lookup("omega"));
        //!Surface Energy
        dimensionedScalar gamma(transportProperties.lookup("gamma"));
        //!Interface Width
        dimensionedScalar epsilon(transportProperties.lookup("epsilon"));
        //!Composition of solid in equilibrium with liquid
        //dimensionedScalar c_Sol(transportProperties.lookup("c_Sol"));
        //!Composition of liquid in equilibrium with solid
        //dimensionedScalar c_Liq(transportProperties.lookup("c_Liq"));
        //Equilibrium composition or average composition of alloy
        //dimensionedScalar c_eq(transportProperties.lookup("c_eq"));
        //!Anti-trapping coefficient
        dimensionedScalar anti_trap(transportProperties.lookup("anti_trap"));
        //!Diffusivity in solid
        dimensionedScalar diff_Sol(transportProperties.lookup("diff_Sol"));
        //!Diffusivity in liquid
        dimensionedScalar diff_Liq(transportProperties.lookup("diff_Liq"));
        //!Thermal gradient
        dimensionedScalar G(transportProperties.lookup("G"));
        //!Velocity
        dimensionedScalar v(transportProperties.lookup("v"));
        //!Strength of anisotropy
        dimensionedScalar delta_01(transportProperties.lookup("delta_01"));
        //!Constant orientation in radian about x
        dimensionedScalar theta_x(transportProperties.lookup("theta_x"));
        //!Constant orientation in radian about y
        dimensionedScalar theta_y(transportProperties.lookup("theta_y"));
        //!Constant orientation in radian about z
        dimensionedScalar theta_val(transportProperties.lookup("theta_val"));
        //!mu^2 is the strength of the coupling in f_ori equation
        dimensionedScalar mu_th(transportProperties.lookup("mu_th"));             
        //!Small value to avoid singularity
        dimensionedScalar eps(transportProperties.lookup("eps"));
	//!Molar volume
        dimensionedScalar Vm(transportProperties.lookup("Vm"));
        //!Properties from TDB
        //dimensionedScalar A_Liq(transportProperties.lookup("A_Liq"));
        //dimensionedScalar B(transportProperties.lookup("B"));//G=Ac^2;d2G/dc^2=A
        //dimensionedScalar A_Sol(transportProperties.lookup("A_Sol"));
        //dimensionedScalar B_Sol(transportProperties.lookup("B_Sol"));
        //dimensionedScalar D_Sol(transportProperties.lookup("D_Sol"));
        //!Melting temperature
        dimensionedScalar T0(transportProperties.lookup("T0"));
        //!Noise magnitude
        dimensionedScalar noise_mag(transportProperties.lookup("noise_mag"));
        //!Constant value from temperature profile
        dimensionedScalar initial(transportProperties.lookup("initial"));
       
        //!Heat dissipation in J/mol-s
        dimensionedScalar qdot(transportProperties.lookup("qdot"));
        
        
        //!Switch = 1 to enable cooling, else 0
        scalar swcool(readScalar(transportProperties.lookup("swcool")));
        
        //!Switch = 1 to enable orientaton field, 2 to model precipitate growth,
        //!else 0 for single dendrite solidification
        scalar swch(readScalar(transportProperties.lookup("swch")));
        
        //!M_theta function scaling with time step
        scalar M_Scale(readScalar(transportProperties.lookup("M_Scale")));
 
        scalar seedRadius(readScalar(transportProperties.lookup("seedRadius")));
        
        //! Maximum number depending on probability of nucleation
        scalar numSeeds(readScalar(transportProperties.lookup("numSeeds")));

        //! Temperature parameters in probability of nucleation
        scalar dt_0(readScalar(transportProperties.lookup("dt_0")));
        scalar dt_s(readScalar(transportProperties.lookup("dt_s")));
        
        scalar nucleation_interval(readScalar(transportProperties.lookup("nucleation_interval")));
        
        //!Properties for incorporating elasticity
        
        dimensionedSymmTensor cEigenStrain(MechProperties.lookup("cEigenStrain"));
        dimensionedScalar     mu1_elast(MechProperties.lookup("mu1_elast"));
        dimensionedScalar     mu2_elast(MechProperties.lookup("mu2_elast"));
        dimensionedScalar     mu1_elast_(MechProperties.lookup("mu1_elast_")); 
        dimensionedScalar     mu2_elast_(MechProperties.lookup("mu2_elast_"));
        dimensionedScalar     lambda1(MechProperties.lookup("lambda1"));
        dimensionedScalar     lambda2(MechProperties.lookup("lambda2"));
        
        
