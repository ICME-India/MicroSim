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
        scalar dimensions(readScalar(transportProperties.lookup("dimensions")));
        //!Number of phases
        scalar phases(readScalar(transportProperties.lookup("phases")));
        //!Number of components
        scalar components(readScalar(transportProperties.lookup("components")));
        //!Slope liquidus
        //dimensionedScalar m_1(transportProperties.lookup("m_1"));
        //!Slope solidus
        //dimensionedScalar m_0(transportProperties.lookup("m_0"));
        //Relaxation coefficient for phi (order parameter)
        //dimensionedScalar omega(transportProperties.lookup("omega"));
        //!Surface Energy
        dimensionedScalar gamma_1(transportProperties.lookup("gamma_1"));
        dimensionedScalar gamma_2(transportProperties.lookup("gamma_2"));
        dimensionedScalar gamma_3(transportProperties.lookup("gamma_3"));
        dimensionedScalar gamma_4(transportProperties.lookup("gamma_4"));
        
        //dimensionedScalar DTs(transportProperties.lookup("DTs"));
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
        dimensionedScalar diff_1(transportProperties.lookup("diff_Sol"));
        dimensionedScalar diff_2(transportProperties.lookup("diff_Sol"));
        dimensionedScalar diff_3(transportProperties.lookup("diff_Sol"));
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
        dimensionedScalar theta_z(transportProperties.lookup("theta_z"));
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
 
        //!Reading Filling.in
        //scalar seedRadius(readScalar(transportProperties.lookup("seedRadius")));

        //scalar volume_fraction(readScalar(transportProperties.lookup("volume_fraction")));

        //scalar spread(readScalar(transportProperties.lookup("spread")));

        //scalar centerX(readScalar(transportProperties.lookup("centerX")));
        //scalar centerY(readScalar(transportProperties.lookup("centerY")));
        //scalar centerZ(readScalar(transportProperties.lookup("centerZ")));
        
        scalarRectangularMatrix centerX((phases-1),1,0);
        scalarRectangularMatrix centerY((phases-1),1,0);
        scalarRectangularMatrix centerZ((phases-1),1,0);
        scalarRectangularMatrix seedRadius((phases-1),1,0);
        scalarRectangularMatrix volumeFraction((phases-1),1,0);
        scalarRectangularMatrix shieldDist((phases-1),1,0);
        scalarRectangularMatrix spread((phases-1),1,0);
        
        centerX[0][0]   = readScalar(transportProperties.lookup("centerX[0]"));      
        centerY[0][0]   = readScalar(transportProperties.lookup("centerY[0]"));      
        centerZ[0][0]   = readScalar(transportProperties.lookup("centerZ[0]"));
        seedRadius[0][0] = readScalar(transportProperties.lookup("seedRadius[0]")); 
        volumeFraction[0][0] = readScalar(transportProperties.lookup("volumeFraction[0]"));  
        shieldDist[0][0] = readScalar(transportProperties.lookup("shieldDist[0]")); 
        spread[0][0] = readScalar(transportProperties.lookup("spread[0]")); 
	if (phases != 2)
    	{
        centerX[1][0]   = readScalar(transportProperties.lookup("centerX[1]"));
        centerY[1][0]   = readScalar(transportProperties.lookup("centerY[1]"));
        centerZ[1][0]   = readScalar(transportProperties.lookup("centerZ[1]"));
        seedRadius[1][0] = readScalar(transportProperties.lookup("seedRadius[1]")); 
        volumeFraction[1][0] = readScalar(transportProperties.lookup("volumeFraction[1]")); 
        shieldDist[1][0] = readScalar(transportProperties.lookup("shieldDist[1]")); 
        spread[1][0] = readScalar(transportProperties.lookup("spread[1]")); 
        }
	if (phases == 4)
    	{
        centerX[2][0]   = readScalar(transportProperties.lookup("centerX[2]"));      
        centerY[2][0]   = readScalar(transportProperties.lookup("centerY[2]"));      
        centerZ[2][0]   = readScalar(transportProperties.lookup("centerZ[2]"));
        seedRadius[2][0] = readScalar(transportProperties.lookup("seedRadius[2]")); 
        volumeFraction[2][0] = readScalar(transportProperties.lookup("volumeFraction[2]")); 
        shieldDist[2][0] = readScalar(transportProperties.lookup("shieldDist[2]")); 
        spread[2][0] = readScalar(transportProperties.lookup("spread[2]")); 
        }
        
        //! Maximum number depending on probability of nucleation
        scalar numSeeds(readScalar(transportProperties.lookup("numSeeds")));

        //! Temperature parameters in probability of nucleation
        scalar dt_0(readScalar(transportProperties.lookup("dt_0")));
        scalar dt_s(readScalar(transportProperties.lookup("dt_s")));
        
        scalar nucleation_interval(readScalar(transportProperties.lookup("nucleation_interval")));
       

        scalarRectangularMatrix diff_ter(2*phases,2,0);
        
        
        if (components == 3)
    	{
        diff_ter[0][0]   = readScalar(transportProperties.lookup("D_a1[0][0]"));      
        diff_ter[0][1]   = readScalar(transportProperties.lookup("D_a1[0][1]"));
        diff_ter[1][0]   = readScalar(transportProperties.lookup("D_a1[1][0]"));
        diff_ter[1][1]   = readScalar(transportProperties.lookup("D_a1[1][1]"));
        
        if (phases == 2)
    	{
    	diff_ter[2][0]   = readScalar(transportProperties.lookup("D_l[0][0]"));      
        diff_ter[2][1]   = readScalar(transportProperties.lookup("D_l[0][1]"));
        diff_ter[3][0]   = readScalar(transportProperties.lookup("D_l[1][0]"));
        diff_ter[3][1]   = readScalar(transportProperties.lookup("D_l[1][1]"));
    	}
        
    	if (phases == 3)
    	{
        diff_ter[2][0]   = readScalar(transportProperties.lookup("D_a2[0][0]"));      
        diff_ter[2][1]   = readScalar(transportProperties.lookup("D_a2[0][1]"));
        diff_ter[3][0]   = readScalar(transportProperties.lookup("D_a2[1][0]"));
        diff_ter[3][1]   = readScalar(transportProperties.lookup("D_a2[1][1]"));
        
        diff_ter[4][0]   = readScalar(transportProperties.lookup("D_l[0][0]"));      
        diff_ter[4][1]   = readScalar(transportProperties.lookup("D_l[0][1]"));
        diff_ter[5][0]   = readScalar(transportProperties.lookup("D_l[1][0]"));
        diff_ter[5][1]   = readScalar(transportProperties.lookup("D_l[1][1]"));
        }
        
        if (phases == 4)
    	{        
        diff_ter[2][0]   = readScalar(transportProperties.lookup("D_a2[0][0]"));      
        diff_ter[2][1]   = readScalar(transportProperties.lookup("D_a2[0][1]"));
        diff_ter[3][0]   = readScalar(transportProperties.lookup("D_a2[1][0]"));
        diff_ter[3][1]   = readScalar(transportProperties.lookup("D_a2[1][1]"));
        
        diff_ter[4][0]   = readScalar(transportProperties.lookup("D_a3[0][0]"));      
        diff_ter[4][1]   = readScalar(transportProperties.lookup("D_a3[0][1]"));
        diff_ter[5][0]   = readScalar(transportProperties.lookup("D_a3[1][0]"));
        diff_ter[5][1]   = readScalar(transportProperties.lookup("D_a3[1][1]"));
        
        diff_ter[6][0]   = readScalar(transportProperties.lookup("D_l[0][0]"));      
        diff_ter[6][1]   = readScalar(transportProperties.lookup("D_l[0][1]"));
        diff_ter[7][0]   = readScalar(transportProperties.lookup("D_l[1][0]"));
        diff_ter[7][1]   = readScalar(transportProperties.lookup("D_l[1][1]"));
        }
    	}
        
        //!Properties for incorporating elasticity
        
        dimensionedSymmTensor cEigenStrain1(MechProperties.lookup("cEigenStrain1"));
        dimensionedSymmTensor cEigenStrain2(MechProperties.lookup("cEigenStrain2"));
        dimensionedSymmTensor cEigenStrain3(MechProperties.lookup("cEigenStrain3"));

        dimensionedScalar     mu1_elast(MechProperties.lookup("mu1_elast"));
        dimensionedScalar     mu2_elast(MechProperties.lookup("mu2_elast"));
        dimensionedScalar     mu1_elast_(MechProperties.lookup("mu1_elast_")); 
        dimensionedScalar     mu2_elast_(MechProperties.lookup("mu2_elast_"));
        dimensionedScalar     lambda1(MechProperties.lookup("lambda1"));
        dimensionedScalar     lambda2(MechProperties.lookup("lambda2"));
        
        
