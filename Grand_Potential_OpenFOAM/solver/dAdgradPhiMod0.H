// Header file to obtain a specific terms in the phase-field equation (check userGuide for the equations)
    //! Normal vector to the interface, q
     volVectorField q=dimx*fvc::grad(phi_1);
     
     // Writing into the log file
     Info<< "Min/max magq:" << min(mag(q)).value() << ' '
     << max(mag(q)).value() << endl;
     
     //scalar theta1 = theta*pi/180.0;
if (phases == 2)
{
    //! For 2D simulations
    if (dimensions== 2)
    {
     //! Rotated field   
     q = (q.component(0)*Foam::cos(theta.component(2))-q.component(1)*Foam::sin(theta.component(2)))*vector(1,0,0) + (q.component(0)*Foam::sin(theta.component(2))+q.component(1)*Foam::cos(theta.component(2)))*vector(0,1,0);

        //Info<< "Min/max mag qDash:" << min(mag(qDash)).value() << ", " << max(mag(qDash)).value() << endl;

     //! Anisotropy, ac(q)
        ac_01 = (1.0-delta_01*(3.0-4.0*(pow(q.component(0),4) + pow(q.component(1), 4))
                                                        / ( 1E-20 + pow(mag(q),4) )));

     //! dAdq01 = Derivative of ac(q)
         dAdq01
            = 16.0*delta_01*(( pow(q.component(0),3)/(1E-20 +pow(mag(q),4))
                 - q.component(0)*(pow(q.component(0),4) + pow(q.component(1), 4))/(1E-20 +pow(mag(q),6))) *vector(1,0,0)
                                     + (pow(q.component(1),3)/(1E-20 +pow(mag(q),4))
                                     - q.component(1)*(pow(q.component(0),4) + pow(q.component(1), 4) )/(1E-20 +pow(mag(q),6))) *vector(0,1,0)
                          );
	dAdq01
	= (dAdq01.component(0)*Foam::cos(theta.component(2)) + dAdq01.component(1)*Foam::sin(theta.component(2))) *vector(1,0,0)
	+ (-dAdq01.component(0)*Foam::sin(theta.component(2)) + dAdq01.component(1)*Foam::cos(theta.component(2))) *vector(0,1,0);

}

//! For 3D simulations
if (dimensions == 3)
{
     //! Rotated field
     q = (q.component(0)*(Foam::cos(theta.component(1))*Foam::cos(theta.component(2))) + q.component(1)*(Foam::sin(theta.component(0))*Foam::sin(theta.component(1))*Foam::cos(theta.component(2)) - Foam::cos(theta.component(0))*Foam::sin(theta.component(2))) + q.component(2)*(Foam::cos(theta.component(0))*Foam::sin(theta.component(1))*Foam::cos(theta.component(2)) + Foam::sin(theta.component(0))*Foam::sin(theta.component(2))))*vector(1,0,0)
     + (q.component(0)*(Foam::cos(theta.component(1))*Foam::sin(theta.component(2))) + q.component(1)*(Foam::sin(theta.component(0))*Foam::sin(theta.component(1))*Foam::sin(theta.component(2)) + Foam::cos(theta.component(0))*Foam::cos(theta.component(2))) + q.component(2)*(Foam::cos(theta.component(0))*Foam::sin(theta.component(1))*Foam::sin(theta.component(2)) - Foam::sin(theta.component(0))*Foam::cos(theta.component(2))))*vector(0,1,0)
     + (q.component(0)*(-Foam::sin(theta.component(1))) + q.component(1)*(Foam::sin(theta.component(0))*Foam::cos(theta.component(1))) + q.component(2)*(Foam::cos(theta.component(0))*Foam::cos(theta.component(1))))*vector(0,0,1);

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
	= (dAdq01.component(0)*(Foam::cos(theta.component(1))*Foam::cos(theta.component(2))) + dAdq01.component(1)*(Foam::cos(theta.component(1))*Foam::sin(theta.component(2))) + dAdq01.component(2)*(-Foam::sin(theta.component(1)))) *vector(1,0,0)
	+ (dAdq01.component(0)*(Foam::sin(theta.component(0))*Foam::sin(theta.component(1))*Foam::cos(theta.component(2)) - Foam::cos(theta.component(0))*Foam::sin(theta.component(2))) + dAdq01.component(1)*(Foam::sin(theta.component(0))*Foam::sin(theta.component(1))*Foam::sin(theta.component(2)) + Foam::cos(theta.component(0))*Foam::cos(theta.component(2))) + dAdq01.component(2)*(Foam::sin(theta.component(0))*Foam::cos(theta.component(1)))) *vector(0,1,0)
	+ (dAdq01.component(0)*(Foam::cos(theta.component(0))*Foam::sin(theta.component(1))*Foam::cos(theta.component(2)) + Foam::sin(theta.component(0))*Foam::sin(theta.component(2))) + dAdq01.component(1)*(Foam::cos(theta.component(0))*Foam::sin(theta.component(1))*Foam::sin(theta.component(2)) - Foam::sin(theta.component(0))*Foam::cos(theta.component(2))) + dAdq01.component(2)*(Foam::cos(theta.component(0))*Foam::cos(theta.component(1))))*vector(0,0,1);
}
 //! A term in the derivative of gradient energy density, a(q) w.r.t q in the phase-field equation
 dadgradPhi = (ac_01*dAdq01*magSqr(q));

 //! Writing into the log file
 Info<< "Min/max ac_01:" << min(ac_01).value() << ' '
     << max(ac_01).value() << endl;

//! Writing into the log file
Info<< "Min/max dadgradPhi:" << min(dadgradPhi).value() << ' '
     << max(dadgradPhi).value() << endl;
     }
