        dimensionedScalar T_old = T;
        
    solidVolFrac = fvc::domainIntegrate(phi_1)/totalVol;
    
    dsolidVolFracOld = dsolidVolFrac;
    dsolidVolFrac = (solidVolFrac - solidVolFracOld);
    
        if (swcool == 1)
        {
            dimensionedScalar Lf = gsl_spline_eval (spline_L, T.value(), acc_L);
            dimensionedScalar Cp = gsl_spline_eval (spline_C, T.value(), acc_C);
            T = T + (Lf/Cp)*dsolidVolFrac - (qdot/Cp)*runTime.deltaTValue();
        
        }

        if (T.value() > Tm.value()) {
        T.value() = Tm.value();
    	}
    	else if (T.value() < T1[0]) {
        T.value() = T1[0];
        break;
    	}
    
        solidVolFracOld = solidVolFrac;

        outpft << "time " << runTime.timeName() << " " << T.value() << " " << solidVolFrac.value() << nl << endl;
        Info << "Temperature: "  << T.value() << " K" << endl;
        
        dimensionedScalar Tdot = (T-T_old)/runTime.deltaTValue();
