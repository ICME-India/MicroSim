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

    dimensionedScalar dimt(transportProperties.lookup("dimt"));    //diemsion of time
    dimensionedScalar dimx(transportProperties.lookup("dimx"));    //diemsion of position
	dimensionedScalar L_s(transportProperties.lookup("L_s"));
    dimensionedScalar L_s_r(transportProperties.lookup("L_s_r"));                    
	dimensionedScalar L_e(transportProperties.lookup("L_e"));       
    dimensionedScalar L_e_r(transportProperties.lookup("L_e_r"));
    dimensionedScalar k(transportProperties.lookup("k"));    
	dimensionedScalar k_r(transportProperties.lookup("k_r"));
    dimensionedScalar W(transportProperties.lookup("W"));
    dimensionedScalar W_r(transportProperties.lookup("W_r"));
	dimensionedScalar alpha(transportProperties.lookup("alpha"));
	dimensionedScalar C_m_s(transportProperties.lookup("C_m_s"));
	dimensionedScalar C_m_l(transportProperties.lookup("C_m_l"));
    dimensionedScalar C_m_l_r(transportProperties.lookup("C_m_l_r"));
    dimensionedScalar C_m_s_r(transportProperties.lookup("C_m_s_r"));
    dimensionedScalar C_0(transportProperties.lookup("C_0"));
    dimensionedScalar C_0_r(transportProperties.lookup("C_0_r"));
	dimensionedScalar S_s(transportProperties.lookup("S_s"));
	dimensionedScalar S_l(transportProperties.lookup("S_l"));
	dimensionedScalar D(transportProperties.lookup("D"));
    dimensionedScalar D_s(transportProperties.lookup("D_s"));
    dimensionedScalar D_l(transportProperties.lookup("D_l"));
	dimensionedScalar n(transportProperties.lookup("n"));
	dimensionedScalar F(transportProperties.lookup("F"));
	dimensionedScalar R(transportProperties.lookup("R"));
	dimensionedScalar T(transportProperties.lookup("T"));
	dimensionedScalar e_s(transportProperties.lookup("e_s"));
	dimensionedScalar e_l(transportProperties.lookup("e_l"));
    dimensionedScalar dt(transportProperties.lookup("dt"));    
    dimensionedScalar dx(transportProperties.lookup("dx"));    
    //dimensionedScalar dt_r(transportProperties.lookup("dt_r"));    //diemsion of time
    //dimensionedScalar dx_r(transportProperties.lookup("dx_r"));    //diemsion of position
	//dimensionedScalar W(transportProperties.lookup("W"));
	//dimensionedScalar kappa(transportProperties.lookup("kappa"));


        
//         dimensionedScalar m_1(transportProperties.lookup("m_1"));      //slope liquidus
//         dimensionedScalar m_0(transportProperties.lookup("m_0"));      //slope solidus
//         dimensionedScalar omega(transportProperties.lookup("omega"));  //relaxation coefficient for alpha (order parameter)
//         dimensionedScalar gamma(transportProperties.lookup("gamma"));           //Surface Energy
//         dimensionedScalar epsilon(transportProperties.lookup("epsilon"));       //Interface Width 
//         dimensionedScalar c_Sol(transportProperties.lookup("c_Sol"));           //composition of solid in equilibrium with liquid
//         dimensionedScalar c_Liq(transportProperties.lookup("c_Liq"));          //composition of liquid in quilibrium with solid
//         dimensionedScalar c_eq(transportProperties.lookup("c_eq"));            //equilibrium composition or avg composition of alloy
//         dimensionedScalar anti_trap(transportProperties.lookup("anti_trap"));           //antitrapping coefficient
//         dimensionedScalar diff_Sol(transportProperties.lookup("diff_Sol"));             //diffusivity in solid
//         dimensionedScalar diff_Liq(transportProperties.lookup("diff_Liq"));             //diffusivity in liquid
//         dimensionedScalar G(transportProperties.lookup("G"));                           //thermal gradient
//         dimensionedScalar v(transportProperties.lookup("v"));                           //velocity
//         dimensionedScalar delta_01(transportProperties.lookup("delta_01"));             //strength of anisotropy
//         dimensionedScalar A(transportProperties.lookup("A")); 
//         dimensionedScalar B(transportProperties.lookup("B"));//G=Ac^2;d2G/dc^2=A
//         dimensionedScalar T0(transportProperties.lookup("T0"));                         //Melting Temp
// 	//dimensionedScalar T(transportProperties.lookup("T"));                         //operating Temp
//         dimensionedScalar noise_mag(transportProperties.lookup("noise_mag"));           //noise magnitude
//         dimensionedScalar initial(transportProperties.lookup("initial"));               //constant value from temperature profile

