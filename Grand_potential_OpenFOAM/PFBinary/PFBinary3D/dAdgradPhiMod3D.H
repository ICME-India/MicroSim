// Header file to obtain a specific terms in the phase-field equation (check userGuide for the equations)
    //! Normal vector to the interface, q
     volVectorField q=dimx*fvc::grad(phi);
     
     // Writing into the log file
     Info<< "Min/max magq:" << min(mag(q)).value() << ' '
     << max(mag(q)).value() << endl;
     
     //scalar theta1 = theta*pi/180.0;

     //! Rotated field   
     q = (q.component(0)*(Foam::cos(thetay)*Foam::cos(theta)) + q.component(1)*(Foam::sin(thetax)*Foam::sin(thetay)*Foam::cos(theta) - Foam::cos(thetax)*Foam::sin(theta)) + q.component(2)*(Foam::cos(thetax)*Foam::sin(thetay)*Foam::cos(theta) + Foam::sin(thetax)*Foam::sin(theta)))*vector(1,0,0) 
     + (q.component(0)*(Foam::cos(thetay)*Foam::sin(theta)) + q.component(1)*(Foam::sin(thetax)*Foam::sin(thetay)*Foam::sin(theta) + Foam::cos(thetax)*Foam::cos(theta)) + q.component(2)*(Foam::cos(thetax)*Foam::sin(thetay)*Foam::sin(theta) - Foam::sin(thetax)*Foam::cos(theta)))*vector(0,1,0) 
     + (q.component(0)*(-Foam::sin(thetay)) + q.component(1)*(Foam::sin(thetax)*Foam::cos(thetay)) + q.component(2)*(Foam::cos(thetax)*Foam::cos(thetay)))*vector(0,0,1);

        //Info<< "Min/max mag qDash:" << min(mag(qDash)).value() << ", " << max(mag(qDash)).value() << endl;

     //! Anisotropy, ac(q)
        ac_01 = (1.0-delta_01*(3.0-4.0*(pow(q.component(0),4) + pow(q.component(1), 4) + pow(q.component(2), 4))
                                                        / ( 1E-20 + pow(mag(q),4) )));

     //! dAdq01 = Derivative of ac(q)
         dAdq01
            = 16.0*delta_01*(( pow(q.component(0),3)/(1E-20 +pow(mag(q),4))
                 - q.component(0)*(pow(q.component(0),4) + pow(q.component(1), 4) + pow(q.component(2), 4))/(1E-20 +pow(mag(q),6))) *vector(1,0,0)
                                     + (pow(q.component(1),3)/(1E-20 +pow(mag(q),4))
                                     - q.component(1)*(pow(q.component(0),4) + pow(q.component(1), 4) + pow(q.component(2), 4) )/(1E-20 +pow(mag(q),6))) *vector(0,1,0)
                                     + (pow(q.component(2),3)/(1E-20 +pow(mag(q),4))
                                     - q.component(2)*(pow(q.component(0),4) + pow(q.component(1), 4) + pow(q.component(2), 4) )/(1E-20 +pow(mag(q),6))) *vector(0,0,1)
                          );
	dAdq01
	= (dAdq01.component(0)*(Foam::cos(thetay)*Foam::cos(theta)) + dAdq01.component(1)*(Foam::cos(thetay)*Foam::sin(theta)) + dAdq01.component(2)*(-Foam::sin(thetay))) *vector(1,0,0)
	+ (dAdq01.component(0)*(Foam::sin(thetax)*Foam::sin(thetay)*Foam::cos(theta) - Foam::cos(thetax)*Foam::sin(theta)) + dAdq01.component(1)*(Foam::sin(thetax)*Foam::sin(thetay)*Foam::sin(theta) + Foam::cos(thetax)*Foam::cos(theta)) + dAdq01.component(2)*(Foam::sin(thetax)*Foam::cos(thetay))) *vector(0,1,0) 
	+ (dAdq01.component(0)*(Foam::cos(thetax)*Foam::sin(thetay)*Foam::cos(theta) + Foam::sin(thetax)*Foam::sin(theta)) + dAdq01.component(1)*(Foam::cos(thetax)*Foam::sin(thetay)*Foam::sin(theta) - Foam::sin(thetax)*Foam::cos(theta)) + dAdq01.component(2)*(Foam::cos(thetax)*Foam::cos(thetay)))*vector(0,0,1);

 //! A term in the derivative of gradient energy density, a(q) w.r.t q in the phase-field equation
 dadgradPhi = (ac_01*dAdq01*magSqr(q));
 
 //! Writing into the log file
 Info<< "Min/max ac_01:" << min(ac_01).value() << ' '
     << max(ac_01).value() << endl;
     
//! Writing into the log file
Info<< "Min/max dadgradPhi:" << min(dadgradPhi).value() << ' '
     << max(dadgradPhi).value() << endl;
