        //! Timestep increment for accelerating simulation.
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//! Checking monotonicity for stability of phi equation
	//! Change of average phi showing no oscillation
	if (sign(dsolidVolFrac.value()) == sign(dsolidVolFracOld.value()))
	//if (dsolidVolFrac.value() >= 0)
	{
	iter_num += 1;

	if (iter_num > 100 && max(max(InitialResidual_0,max(InitialResidual_10,InitialResidual_11)),max(max(max(InitialResidual_20,InitialResidual_21),max(InitialResidual_22,InitialResidual_23)),InitialResidual_3)) < Tol) //.value() && runTime.value() < 1000*InitialDeltaT
	{   
    	runTime.setDeltaT
        	(
            	runTime.deltaTValue() + dtf*InitialDeltaT
        	);
    	// set the iter_num = 0
    	iter_num = 0;
    	Info<< "deltaT increased to " <<  runTime.deltaTValue() << endl;
    	
    
	} else if (iter_num > 100 && max(max(InitialResidual_0,max(InitialResidual_10,InitialResidual_11)),max(max(max(InitialResidual_20,InitialResidual_21),max(InitialResidual_22,InitialResidual_23)),InitialResidual_3)) > Tol)
	{
    
	runTime.setDeltaT
        	(
            	runTime.deltaTValue() - dtf*InitialDeltaT 
        	);
	if (runTime.deltaTValue() < InitialDeltaT) //|| mag(Tdot) > )
	  {
          	runTime.setDeltaT(InitialDeltaT);  
          }
    	iter_num =0;
    	Info<< "deltaT decreased to " <<  runTime.deltaTValue() << endl;

	}

	Info<< "iter_num = " <<  iter_num << endl;
	}
	//! Change of average phi showing oscillation
	if (sign(dsolidVolFrac.value()) != sign(dsolidVolFracOld.value()))
	//if (dsolidVolFrac.value() < 0)
	{
	runTime.setDeltaT(InitialDeltaT);
    	Info<< "deltaT decreased to " <<  runTime.deltaTValue() << endl;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
