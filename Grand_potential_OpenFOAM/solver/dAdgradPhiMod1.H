// Header file to obtain a specific terms in the phase-field equation (check userGuide for the equations)
    //! Normal vector to the interface, q
     //volVectorField q=dimx*fvc::grad(phi_1);
     
     // Writing into the log file
     //Info<< "Min/max magq:" << min(mag(q)).value() << ' '
     //<< max(mag(q)).value() << endl;
     
     //scalar theta1 = theta*pi/180.0;
if (phases == 5)
{
q2=dimx*fvc::grad(phi_2);
q3=dimx*fvc::grad(phi_3);
q4=dimx*fvc::grad(phi_4);
     
     // Writing into the log file
     Info<< "Min/max magq2:" << min(mag(q2)).value() << ' '
     << max(mag(q2)).value() << endl;
     Info<< "Min/max magq3:" << min(mag(q3)).value() << ' '
     << max(mag(q3)).value() << endl;
     Info<< "Min/max magq4:" << min(mag(q4)).value() << ' '
     << max(mag(q4)).value() << endl;

    //! For 2D simulations
    if (dimensions== 2)
    {
     //! Rotated field   
     q2 = (q2.component(0)*Foam::cos(theta.component(2))-q2.component(1)*Foam::sin(theta.component(2)))*vector(1,0,0) + (q2.component(0)*Foam::sin(theta.component(2))+q2.component(1)*Foam::cos(theta.component(2)))*vector(0,1,0);

     q3 = (q3.component(0)*Foam::cos(theta.component(2))-q3.component(1)*Foam::sin(theta.component(2)))*vector(1,0,0) + (q3.component(0)*Foam::sin(theta.component(2))+q3.component(1)*Foam::cos(theta.component(2)))*vector(0,1,0);

     q4 = (q4.component(0)*Foam::cos(theta.component(2))-q4.component(1)*Foam::sin(theta.component(2)))*vector(1,0,0) + (q4.component(0)*Foam::sin(theta.component(2))+q4.component(1)*Foam::cos(theta.component(2)))*vector(0,1,0);

        //Info<< "Min/max mag qDash:" << min(mag(qDash)).value() << ", " << max(mag(qDash)).value() << endl;

     //! Anisotropy, ac(q)
        ac_02 = (1.0-delta_01*(3.0-4.0*(pow(q2.component(0),4) + pow(q2.component(1), 4))
                                                        / ( 1E-20 + pow(mag(q2),4) )));

        ac_03 = (1.0-delta_01*(3.0-4.0*(pow(q3.component(0),4) + pow(q3.component(1), 4))
                                                        / ( 1E-20 + pow(mag(q3),4) )));

        ac_04 = (1.0-delta_01*(3.0-4.0*(pow(q4.component(0),4) + pow(q4.component(1), 4))
                                                        / ( 1E-20 + pow(mag(q4),4) )));

     //! dAdq01 = Derivative of ac(q)
        /* dAdq01
            = 16.0*delta_01*(( pow(q.component(0),3)/(1E-20 +pow(mag(q),4))
                 - q.component(0)*(pow(q.component(0),4) + pow(q.component(1), 4))/(1E-20 +pow(mag(q),6))) *vector(1,0,0)
                                     + (pow(q.component(1),3)/(1E-20 +pow(mag(q),4))
                                     - q.component(1)*(pow(q.component(0),4) + pow(q.component(1), 4) )/(1E-20 +pow(mag(q),6))) *vector(0,1,0)
                          );
	dAdq01
	= (dAdq01.component(0)*Foam::cos(theta.component(2)) + dAdq01.component(1)*Foam::sin(theta.component(2))) *vector(1,0,0)
	+ (-dAdq01.component(0)*Foam::sin(theta.component(2)) + dAdq01.component(1)*Foam::cos(theta.component(2))) *vector(0,1,0);*/

}

//! For 3D simulations
if (dimensions == 3)
{
     //! Rotated field
     q2 = (q2.component(0)*(Foam::cos(theta.component(1))*Foam::cos(theta.component(2))) + q2.component(1)*(Foam::sin(theta.component(0))*Foam::sin(theta.component(1))*Foam::cos(theta.component(2)) - Foam::cos(theta.component(0))*Foam::sin(theta.component(2))) + q2.component(2)*(Foam::cos(theta.component(0))*Foam::sin(theta.component(1))*Foam::cos(theta.component(2)) + Foam::sin(theta.component(0))*Foam::sin(theta.component(2))))*vector(1,0,0)
     + (q2.component(0)*(Foam::cos(theta.component(1))*Foam::sin(theta.component(2))) + q2.component(1)*(Foam::sin(theta.component(0))*Foam::sin(theta.component(1))*Foam::sin(theta.component(2)) + Foam::cos(theta.component(0))*Foam::cos(theta.component(2))) + q2.component(2)*(Foam::cos(theta.component(0))*Foam::sin(theta.component(1))*Foam::sin(theta.component(2)) - Foam::sin(theta.component(0))*Foam::cos(theta.component(2))))*vector(0,1,0)
     + (q2.component(0)*(-Foam::sin(theta.component(1))) + q2.component(1)*(Foam::sin(theta.component(0))*Foam::cos(theta.component(1))) + q2.component(2)*(Foam::cos(theta.component(0))*Foam::cos(theta.component(1))))*vector(0,0,1);

     q3 = (q3.component(0)*(Foam::cos(theta.component(1))*Foam::cos(theta.component(2))) + q3.component(1)*(Foam::sin(theta.component(0))*Foam::sin(theta.component(1))*Foam::cos(theta.component(2)) - Foam::cos(theta.component(0))*Foam::sin(theta.component(2))) + q3.component(2)*(Foam::cos(theta.component(0))*Foam::sin(theta.component(1))*Foam::cos(theta.component(2)) + Foam::sin(theta.component(0))*Foam::sin(theta.component(2))))*vector(1,0,0)
     + (q3.component(0)*(Foam::cos(theta.component(1))*Foam::sin(theta.component(2))) + q3.component(1)*(Foam::sin(theta.component(0))*Foam::sin(theta.component(1))*Foam::sin(theta.component(2)) + Foam::cos(theta.component(0))*Foam::cos(theta.component(2))) + q3.component(2)*(Foam::cos(theta.component(0))*Foam::sin(theta.component(1))*Foam::sin(theta.component(2)) - Foam::sin(theta.component(0))*Foam::cos(theta.component(2))))*vector(0,1,0)
     + (q3.component(0)*(-Foam::sin(theta.component(1))) + q3.component(1)*(Foam::sin(theta.component(0))*Foam::cos(theta.component(1))) + q3.component(2)*(Foam::cos(theta.component(0))*Foam::cos(theta.component(1))))*vector(0,0,1);

     q4 = (q4.component(0)*(Foam::cos(theta.component(1))*Foam::cos(theta.component(2))) + q4.component(1)*(Foam::sin(theta.component(0))*Foam::sin(theta.component(1))*Foam::cos(theta.component(2)) - Foam::cos(theta.component(0))*Foam::sin(theta.component(2))) + q4.component(2)*(Foam::cos(theta.component(0))*Foam::sin(theta.component(1))*Foam::cos(theta.component(2)) + Foam::sin(theta.component(0))*Foam::sin(theta.component(2))))*vector(1,0,0)
     + (q4.component(0)*(Foam::cos(theta.component(1))*Foam::sin(theta.component(2))) + q4.component(1)*(Foam::sin(theta.component(0))*Foam::sin(theta.component(1))*Foam::sin(theta.component(2)) + Foam::cos(theta.component(0))*Foam::cos(theta.component(2))) + q4.component(2)*(Foam::cos(theta.component(0))*Foam::sin(theta.component(1))*Foam::sin(theta.component(2)) - Foam::sin(theta.component(0))*Foam::cos(theta.component(2))))*vector(0,1,0)
     + (q4.component(0)*(-Foam::sin(theta.component(1))) + q4.component(1)*(Foam::sin(theta.component(0))*Foam::cos(theta.component(1))) + q4.component(2)*(Foam::cos(theta.component(0))*Foam::cos(theta.component(1))))*vector(0,0,1);

        //Info<< "Min/max mag qDash:" << min(mag(qDash)).value() << ", " << max(mag(qDash)).value() << endl;

     //! Anisotropy, ac(q)
        ac_02 = (1.0-delta_01*(3.0-4.0*(pow(q2.component(0),4) + pow(q2.component(1), 4) + pow(q2.component(2), 4))
                                                        / ( 1E-20 + pow(mag(q2),4) )));

        ac_03 = (1.0-delta_01*(3.0-4.0*(pow(q3.component(0),4) + pow(q3.component(1), 4) + pow(q3.component(2), 4))
                                                        / ( 1E-20 + pow(mag(q3),4) )));

        ac_04 = (1.0-delta_01*(3.0-4.0*(pow(q4.component(0),4) + pow(q4.component(1), 4) + pow(q4.component(2), 4))
                                                        / ( 1E-20 + pow(mag(q4),4) )));

     //! dAdq01 = Derivative of ac(q)
        /* dAdq01
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
	+ (dAdq01.component(0)*(Foam::cos(theta.component(0))*Foam::sin(theta.component(1))*Foam::cos(theta.component(2)) + Foam::sin(theta.component(0))*Foam::sin(theta.component(2))) + dAdq01.component(1)*(Foam::cos(theta.component(0))*Foam::sin(theta.component(1))*Foam::sin(theta.component(2)) - Foam::sin(theta.component(0))*Foam::cos(theta.component(2))) + dAdq01.component(2)*(Foam::cos(theta.component(0))*Foam::cos(theta.component(1))))*vector(0,0,1);*/
}
 //! A term in the derivative of gradient energy density, a(q) w.r.t q in the phase-field equation
 //dadgradPhi = (ac_01*dAdq01*magSqr(q));

 //! Writing into the log file
 Info<< "Min/max ac_02:" << min(ac_02).value() << ' '
     << max(ac_02).value() << endl;

 Info<< "Min/max ac_03:" << min(ac_03).value() << ' '
     << max(ac_03).value() << endl;

 Info<< "Min/max ac_04:" << min(ac_04).value() << ' '
     << max(ac_04).value() << endl;

//! Writing into the log file
//Info<< "Min/max dadgradPhi:" << min(dadgradPhi).value() << ' '
     //<< max(dadgradPhi).value() << endl;
     }
