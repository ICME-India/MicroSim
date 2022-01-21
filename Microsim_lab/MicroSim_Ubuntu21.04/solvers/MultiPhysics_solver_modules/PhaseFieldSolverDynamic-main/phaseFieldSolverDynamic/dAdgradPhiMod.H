// Header file to obtain a specific terms in the phase-field equation (check userGuide for the equations)
    //! Normal vector to the interface, q
     volVectorField q=dimx*fvc::grad(phi);
     
     // Writing into the log file
     Info<< "Min/max magq:" << min(mag(q)).value() << ' '
     << max(mag(q)).value() << endl;
     
     //! Anisotropy, ac(q)
      ac_01 = 1.0 - delta_01*(3.0 - 4.0*(pow(q.component(0),4) + pow(q.component(1), 4) + pow(q.component(2), 4))/(1E-20 +pow(mag(q),4)));
     
     //! dAdq01 = Derivative of ac(q) with respect to q
      dAdq01 = 16.0*delta_01*((pow(q.component(0),3)/(1E-20 +pow(mag(q),4)) - q.component(0)*(pow(q.component(0),4) + pow(q.component(1), 4)+ pow(q.component(2), 4))/(1E-20 +pow(mag(q),6)))*vector(1,0,0)
      
      +(pow(q.component(1),3)/(1E-20 +pow(mag(q),4)) - q.component(1)*(pow(q.component(0),4) + pow(q.component(1), 4) + pow(q.component(2), 4))/(1E-20 +pow(mag(q),6)))*vector(0,1,0)
      
      +(pow(q.component(2),3)/(1E-20 +pow(mag(q),4)) - q.component(2)*(pow(q.component(0),4) + pow(q.component(1), 4) + pow(q.component(2), 4))/(1E-20 +pow(mag(q),6)))*vector(0,0,1));

 //! A term in the derivative of gradient energy density, a(q) w.r.t q in the phase-field equation
 dadgradPhi = (ac_01*dAdq01*mag(q)*mag(q));
 
 //! Writing into the log file
 Info<< "Min/max ac_01:" << min(ac_01).value() << ' '
     << max(ac_01).value() << endl;
     
//! Writing into the log file
Info<< "Min/max dadgradPhi:" << min(dadgradPhi).value() << ' '
     << max(dadgradPhi).value() << endl;
