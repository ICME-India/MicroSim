//#include "hdh_phi.H" 
dimensionedScalar initial_residual_alpha ;
dimensionedScalar initial_residual_beta ;
dimensionedScalar initial_residual_liq;
scalar counter= 0;
//! Implicit discretization using fvm class for time derivative and laplacian.
//! Explicit discretization using fvc class for gradient and divergence.
//! Phase-field equations with approximate relation between c, mu and T according to the parabolic approximation for free energy
do{
    volVectorField grad_alpha =dimx*fvc::grad(phi_alpha);
    volVectorField grad_beta  =dimx*fvc::grad(phi_beta);
    volVectorField grad_liq   =dimx*fvc::grad(phi_liq);

    fvScalarMatrix phi_alphaEqn
    (
        tau*epsilon*dimt*fvm::ddt(phi_alpha) - v*grad_alpha.component(vector::X)

                ==
                (-2.0/3.0)*(            -2.0*gamma*epsilon * dimx*dimx*fvm::laplacian(phi_alpha)
                                + 	9.0*gamma/epsilon*2.0*(phi_alpha)*(1.0-phi_alpha)*(1.0-2.0*phi_alpha)
                                +		(-1.0/4.0*(mu - B_alpha)*(mu - B_alpha)/A + D_alpha) * 0.25*15*(

                                                3.0*phi_alpha*phi_alpha*phi_alpha*phi_alpha - 4.0*phi_alpha*phi_alpha*phi_alpha
                                            -   phi_alpha*phi_alpha + 2.0*phi_alpha - (phi_beta - phi_liq)*(phi_beta - phi_liq)*(2.0*phi_alpha - 3.0*phi_alpha*phi_alpha)

                                                                                                                )

                                +		( -1.0/4.0*(mu - B_beta)*(mu - B_beta)/A + D_beta ) * -0.5*15*(

                                                       (phi_beta*phi_beta - phi_beta*phi_beta*phi_beta)*(phi_alpha - phi_liq)

                                                                                                         )

                                +		( -1.0/4.0*(mu - B_liq)*(mu - B_liq)/A + D_liq ) * -0.5*15*(

                                                      (phi_liq*phi_liq - phi_liq*phi_liq*phi_liq)*(phi_alpha - phi_beta)

                                                                                                          )
                           )

             +  (1.0/3.0)*(            -2.0*gamma*epsilon * dimx*dimx*fvc::laplacian(phi_liq)
                                + 	9.0*gamma/epsilon*2.0*(phi_liq)*(1.0-phi_liq)*(1.0-2.0*phi_liq)
                                +		( -1.0/4.0*(mu - B_alpha)*(mu - B_alpha)/A + D_alpha) * -0.5*15*(

                                                        (phi_alpha*phi_alpha - phi_alpha*phi_alpha*phi_alpha)*(phi_liq - phi_beta)

                                                                                                         )

                                +		( -1.0/4.0*(mu - B_beta)*(mu - B_beta)/A + D_beta ) * -0.5*15*(

                                                        (phi_beta*phi_beta - phi_beta*phi_beta*phi_beta)*(phi_liq - phi_alpha)

                                                                                                          )

                                +		( -1.0/4.0*(mu - B_liq)*(mu - B_liq)/A + D_liq ) *0.25*15* (

                                                        3.0*phi_liq*phi_liq*phi_liq*phi_liq - 4.0*phi_liq*phi_liq*phi_liq - phi_liq*phi_liq
                                                    +   2.0*phi_liq - (phi_alpha - phi_beta)*(phi_alpha - phi_beta)*(2.0*phi_liq - 3.0*phi_liq*phi_liq)


                                                                                                        )
                           )

             +  (1.0/3.0)*(            -2.0*gamma*epsilon * dimx*dimx*fvc::laplacian(phi_beta)
                                + 	9.0*gamma/epsilon*2.0*(phi_beta)*(1.0-phi_beta)*(1.0-2.0*phi_beta)
                                +		( -1.0/4.0*(mu - B_alpha)*(mu - B_alpha)/A + D_alpha) * -0.5*15*(

                                                                (phi_alpha*phi_alpha - phi_alpha*phi_alpha*phi_alpha)*(phi_beta - phi_liq)

                                                                                                           )

                                +		( -1.0/4.0*(mu - B_beta)*(mu - B_beta)/A + D_beta ) * 0.25*15*(

                                                        3.0*phi_beta*phi_beta*phi_beta*phi_beta - 4.0*phi_beta*phi_beta*phi_beta  - phi_beta*phi_beta
                                                    +   2.0*phi_beta - (phi_liq - phi_alpha)*(phi_liq - phi_alpha)*(2.0*phi_beta - 3.0*phi_beta*phi_beta)
                                                                                                        )

                                +		( -1.0/4.0*(mu - B_liq)*(mu - B_liq)/A + D_liq ) * -0.5*15*(

                                                        (phi_liq*phi_liq - phi_liq*phi_liq*phi_liq)*(phi_beta - phi_alpha)

                                                                                                       )
                           )

    );

    initial_residual_alpha = phi_alphaEqn.solve().max().initialResidual();

    fvScalarMatrix phi_betaEqn
    (tau
        *epsilon*dimt*fvm::ddt(phi_beta) - v*grad_beta.component(vector::X)

                ==
                   (-2.0/3.0)*(         -2.0*gamma*epsilon * dimx*dimx*fvm::laplacian(phi_beta)
                                + 	9.0*gamma/epsilon*2.0*(phi_beta)*(1.0-phi_beta)*(1.0-2.0*phi_beta)
                                +		( -1.0/4.0*(mu - B_alpha)*(mu - B_alpha)/A + D_alpha)* -0.5*15*(

                                                     (phi_alpha*phi_alpha - phi_alpha*phi_alpha*phi_alpha)*(phi_beta - phi_liq)
                                                                                                           )

                                +		( -1.0/4.0*(mu - B_beta)*(mu - B_beta)/A + D_beta ) *  0.25*15*(

                                                       3.0*phi_beta*phi_beta*phi_beta*phi_beta - 4.0*phi_beta*phi_beta*phi_beta  - phi_beta*phi_beta
                                                    +   2.0*phi_beta - (phi_liq - phi_alpha)*(phi_liq - phi_alpha)*(2.0*phi_beta - 3.0*phi_beta*phi_beta)

                                                                                                        )

                                +		( -1.0/4.0*(mu - B_liq)*(mu - B_liq)/A + D_liq ) * -0.5*15*(

                                                        (phi_liq*phi_liq - phi_liq*phi_liq*phi_liq)*(phi_beta - phi_alpha)

                                                                                                       )

                           )

                + (1.0/3.0)*(          -2.0*gamma*epsilon * dimx*dimx*fvc::laplacian(phi_liq)
                                + 	9.0*gamma/epsilon*2.0*(phi_liq)*(1.0-phi_liq)*(1.0-2.0*phi_liq)
                                +		( -1.0/4.0*(mu - B_alpha)*(mu - B_alpha)/A + D_alpha) *-0.5*15* (

                                                        (phi_alpha*phi_alpha - phi_alpha*phi_alpha*phi_alpha)*(phi_liq - phi_beta)
                                                                                                         )

                                +		( -1.0/4.0*(mu - B_beta)*(mu - B_beta)/A + D_beta ) * -0.5*15* (


                                                        (phi_beta*phi_beta - phi_beta*phi_beta*phi_beta)*(phi_liq - phi_alpha)
                                                                                                          )

                                +		( -1.0/4.0*(mu - B_liq)*(mu - B_liq)/A + D_liq ) *  0.25*15*(


                                                        3.0*phi_liq*phi_liq*phi_liq*phi_liq - 4.0*phi_liq*phi_liq*phi_liq - phi_liq*phi_liq
                                                    +   2.0*phi_liq - (phi_alpha - phi_beta)*(phi_alpha - phi_beta)*(2.0*phi_liq - 3.0*phi_liq*phi_liq)

                                                                                                        )

                             )


                + (1.0/3.0)*(          -2.0*gamma*epsilon * dimx*dimx*fvc::laplacian(phi_alpha)
                                + 	9.0*gamma/epsilon*2.0*(phi_alpha)*(1.0-phi_alpha)*(1.0-2.0*phi_alpha)
                                +		( -1.0/4.0*(mu - B_alpha)*(mu - B_alpha)/A + D_alpha) *  0.25*15*(


                                                3.0*phi_alpha*phi_alpha*phi_alpha*phi_alpha - 4.0*phi_alpha*phi_alpha*phi_alpha
                                            -   phi_alpha*phi_alpha + 2.0*phi_alpha - (phi_beta - phi_liq)*(phi_beta - phi_liq)*(2.0*phi_alpha - 3.0*phi_alpha*phi_alpha)
                                                                                                        )

                                +		( -1.0/4.0*(mu - B_beta)*(mu - B_beta)/A + D_beta ) * -0.5*15*(
                                                              (phi_beta*phi_beta - phi_beta*phi_beta*phi_beta)*(phi_alpha - phi_liq)


                                                                                                         )

                                +		( -1.0/4.0*(mu - B_liq)*(mu - B_liq)/A + D_liq ) *-0.5*15* (


                                                                             (phi_liq*phi_liq - phi_liq*phi_liq*phi_liq)*(phi_alpha - phi_beta)
                                                                                                          )

                            )
    );

    initial_residual_beta = phi_betaEqn.solve().max().initialResidual();

    fvScalarMatrix phi_liqEqn
    (tau
        *epsilon*dimt*fvm::ddt(phi_liq) - v*grad_liq.component(vector::X)

                ==
                (-2.0/3.0)*(           -2.0*gamma*epsilon * dimx*dimx*fvm::laplacian(phi_liq)
                                + 	9.0*gamma/epsilon*2.0*(phi_liq)*(1.0-phi_liq)*(1.0-2.0*phi_liq)
                                +		( -1.0/4.0*(mu - B_alpha)*(mu - B_alpha)/A + D_alpha) * -0.5*15*(

                                                        (phi_alpha*phi_alpha - phi_alpha*phi_alpha*phi_alpha)*(phi_liq - phi_beta)
                                                                                                         )

                                +		( -1.0/4.0*(mu - B_beta)*(mu - B_beta)/A + D_beta ) * -0.5*15*(

                                                       (phi_beta*phi_beta - phi_beta*phi_beta*phi_beta)*(phi_liq - phi_alpha)

                                                                                                          )

                                +		( -1.0/4.0*(mu - B_liq)*(mu - B_liq)/A + D_liq ) * 0.25*15* (

                                                    3.0*phi_liq*phi_liq*phi_liq*phi_liq - 4.0*phi_liq*phi_liq*phi_liq - phi_liq*phi_liq
                                                    +   2.0*phi_liq - (phi_alpha - phi_beta)*(phi_alpha - phi_beta)*(2.0*phi_liq - 3.0*phi_liq*phi_liq)

                                                                                                        )

                )

                + (1.0/3.0)*(         -2.0*gamma*epsilon * dimx*dimx*fvc::laplacian(phi_alpha)
                                + 	9.0*gamma/epsilon*2.0*(phi_alpha)*(1.0-phi_alpha)*(1.0-2.0*phi_alpha)
                                +		( -1.0/4.0*(mu - B_alpha)*(mu - B_alpha)/A + D_alpha) *  0.25*15*(

                                                    3.0*phi_alpha*phi_alpha*phi_alpha*phi_alpha - 4.0*phi_alpha*phi_alpha*phi_alpha
                                            -   phi_alpha*phi_alpha + 2.0*phi_alpha - (phi_beta - phi_liq)*(phi_beta - phi_liq)*(2.0*phi_alpha - 3.0*phi_alpha*phi_alpha)
                                                                                                        )

                                +		( -1.0/4.0*(mu - B_beta)*(mu - B_beta)/A + D_beta ) * -0.5*15*(


                                                              (phi_beta*phi_beta - phi_beta*phi_beta*phi_beta)*(phi_alpha - phi_liq)

                                                                                                         )

                                +		( -1.0/4.0*(mu - B_liq)*(mu - B_liq)/A + D_liq ) * -0.5*15*(

                                                          (phi_liq*phi_liq - phi_liq*phi_liq*phi_liq)*(phi_alpha - phi_beta)

                                                                                                          )
                )

                + (1.0/3.0)*(         -2.0*gamma*epsilon * dimx*dimx*fvc::laplacian(phi_beta)
                                + 	9.0*gamma/epsilon*2.0*(phi_beta)*(1.0-phi_beta)*(1.0-2.0*phi_beta)
                                +		( -1.0/4.0*(mu - B_alpha)*(mu - B_alpha)/A + D_alpha) *-0.5*15* (

                                                    (phi_alpha*phi_alpha - phi_alpha*phi_alpha*phi_alpha)*(phi_beta - phi_liq)
                                                                                                           )

                                +		( -1.0/4.0*(mu - B_beta)*(mu - B_beta)/A + D_beta ) *  0.25*15*(

                                                       3.0*phi_beta*phi_beta*phi_beta*phi_beta - 4.0*phi_beta*phi_beta*phi_beta  - phi_beta*phi_beta
                                                    +   2.0*phi_beta - (phi_liq - phi_alpha)*(phi_liq - phi_alpha)*(2.0*phi_beta - 3.0*phi_beta*phi_beta)

                                                                                                        )

                                +		( -1.0/4.0*(mu - B_liq)*(mu - B_liq)/A + D_liq ) * -0.5*15*(

                                                         (phi_liq*phi_liq - phi_liq*phi_liq*phi_liq)*(phi_beta - phi_alpha)
                                                                                                       )

                )
    );

    initial_residual_liq = phi_liqEqn.solve().max().initialResidual();
    counter = counter +1; 

}   while(( initial_residual_alpha.value() > Tol || initial_residual_beta.value() > Tol || initial_residual_liq.value() > Tol ) && counter < 100);
//Info<< "Tol value = " << Tol <<  endl;
Info<< "counter =" << counter <<  endl;

//! Timestep reduction for improving convergence.
//////////////////////////////////////////////////////////////////////////////////////////////////////////
iter_num += 1;

if (initial_residual_alpha.value() > Tol)
{   
    runTime.setDeltaT
        (
            runTime.deltaTValue() - 0.1
        );
	if (runTime.deltaTValue() < 0)
	  {
          	runTime.setDeltaT(0.01);  
          }
    // set the iter_num = 0
    iter_num = 0;
    Info<< "deltaT decreased and = " <<  runTime.deltaTValue() << endl;
    
} else if (iter_num > 100 && initial_residual_alpha.value() < Tol && runTime.deltaTValue() < 1 )
{
    
runTime.setDeltaT
        (
            runTime.deltaTValue() + 0.01
	   
        );
    iter_num =0;
    Info<< "deltaT increased and = " <<  runTime.deltaTValue() << endl;

}

Info<< "iter_num = " <<  iter_num << endl;
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//! Gradient term associated with advection in chemical potential equation
volVectorField grad_c =dimx*fvc::grad( ((mu - B_alpha)/(2*A))*
                                                    (
                                                      0.25*phi_alpha*phi_alpha*
                                                          (
                                                                 15.0*(1-phi_alpha)*( 1.0 + phi_alpha - (phi_beta - phi_liq)*(phi_beta - phi_liq)   )
                                                             +   phi_alpha*( 9.0 * phi_alpha * phi_alpha  - 5.0     )
                                                          )
                                                    )
                                      +
                                      ((mu - B_beta )/(2*A)) *
                                                   (
                                                     0.25*phi_beta*phi_beta*
                                                        (
                                                                15.0*(1-phi_beta)*( 1.0 + phi_beta - (phi_alpha - phi_liq)*(phi_alpha - phi_liq)   )
                                                            +     phi_beta*( 9.0 * phi_beta * phi_beta  - 5.0     )
                                                        )

                                                   )
                                      +
                                        ((mu - B_liq )/(2*A)) *
                                                  (
                                                    0.25*phi_liq*phi_liq*
                                                       (
                                                               15.0*(1-phi_liq)*( 1.0 + phi_liq - (phi_alpha - phi_beta)*(phi_alpha - phi_beta)   )
                                                           +     phi_liq*( 9.0 * phi_liq * phi_liq  - 5.0     )

                                                  )


                                                ));
//! The unit normal vectors to the interfaces with a small number in denominator to prevent solution from diverging
volVectorField n_alpha = -1*dimx*fvc::grad(phi_alpha)/(1E-20+mag(dimx*fvc::grad(phi_alpha)));
volVectorField n_beta  = -1*dimx*fvc::grad(phi_beta) /(1E-20+mag(dimx*fvc::grad(phi_beta)));
volVectorField n_liq   = -1*dimx*fvc::grad(phi_liq)  /(1E-20+mag(dimx*fvc::grad(phi_liq)));

//! Chemical potential equation with approximation according to the parabolic form of free energy
fvScalarMatrix muEqn
(
         (0.5* 1.0/A* (
                                  0.25*phi_alpha*phi_alpha*
                                  (
                                         15.0*(1-phi_alpha)*( 1.0 + phi_alpha - (phi_beta - phi_liq)*(phi_beta - phi_liq)   )
                                     +   phi_alpha*( 9.0 * phi_alpha * phi_alpha  - 5.0     )
                                  )


                              +   0.25*phi_liq*phi_liq*
                                 (
                                         15.0*(1-phi_liq)*( 1.0 + phi_liq - (phi_alpha - phi_beta)*(phi_alpha - phi_beta)   )
                                     +     phi_liq*( 9.0 * phi_liq * phi_liq  - 5.0     )
                                 )


                              +   0.25*phi_beta*phi_beta*
                                 (
                                         15.0*(1-phi_beta)*( 1.0 + phi_beta - (phi_alpha - phi_liq)*(phi_alpha - phi_liq)   )
                                     +     phi_beta*( 9.0 * phi_beta * phi_beta  - 5.0     )
                                 )


                               ))  * dimt * fvm::ddt(mu) - v*grad_c.component(vector::X)


			+  ((mu - B_alpha)/(2*A)) * (
                                            dimt*fvc::ddt(phi_alpha)*0.25*15*(

                                                      3.0*phi_alpha*phi_alpha*phi_alpha*phi_alpha - 4.0*phi_alpha*phi_alpha*phi_alpha
                                            -   phi_alpha*phi_alpha + 2.0*phi_alpha - (phi_beta - phi_liq)*(phi_beta - phi_liq)*(2.0*phi_alpha - 3.0*phi_alpha*phi_alpha)

                                                                                                        )

                                        +   dimt*fvc::ddt(phi_beta)*-0.5*15*(

                                                       (phi_alpha*phi_alpha - phi_alpha*phi_alpha*phi_alpha)*(phi_beta - phi_liq)

                                                                                                           )

                                        +   dimt*fvc::ddt(phi_liq)*-0.5*15*(

                                                                         (phi_alpha*phi_alpha - phi_alpha*phi_alpha*phi_alpha)*(phi_liq - phi_beta)


                                                                                                         )

                                     )

			+  ((mu - B_beta )/(2*A)) *  (
                                            dimt*fvc::ddt(phi_alpha)*-0.5*15*(

                                                        (phi_beta*phi_beta - phi_beta*phi_beta*phi_beta)*(phi_alpha - phi_liq)

                                                                                                         )

                                        +   dimt*fvc::ddt(phi_beta)*0.25*15* (

                                                     3.0*phi_beta*phi_beta*phi_beta*phi_beta - 4.0*phi_beta*phi_beta*phi_beta  - phi_beta*phi_beta
                                                    +   2.0*phi_beta - (phi_liq - phi_alpha)*(phi_liq - phi_alpha)*(2.0*phi_beta - 3.0*phi_beta*phi_beta)

                                                                                                        )

                                        +   dimt*fvc::ddt(phi_liq)*-0.5*15*(

                                                         (phi_beta*phi_beta - phi_beta*phi_beta*phi_beta)*(phi_liq - phi_alpha)


                                                                                                          )

                                     )

			+  ( (mu - B_liq  )/(2*A))  *  (
                                            dimt*fvc::ddt(phi_alpha)*-0.5*15*(
                                                                                     (phi_liq*phi_liq - phi_liq*phi_liq*phi_liq)*(phi_alpha - phi_beta)
                                                                                                          )

			                +   dimt*fvc::ddt(phi_beta)*-0.5*15*(

                                                      (phi_liq*phi_liq - phi_liq*phi_liq*phi_liq)*(phi_beta - phi_alpha)
                                                                                                       )

			                +   dimt*fvc::ddt(phi_liq)* 0.25*15*(

                                                 3.0*phi_liq*phi_liq*phi_liq*phi_liq - 4.0*phi_liq*phi_liq*phi_liq - phi_liq*phi_liq
                                                    +   2.0*phi_liq - (phi_alpha - phi_beta)*(phi_alpha - phi_beta)*(2.0*phi_liq - 3.0*phi_liq*phi_liq)
                                                                                                        )

                                                         )

				 == 0.5*dimx*dimx*fvm::laplacian(D*phi_liq, mu) 
                                 
                            + 0.5*0.707*epsilon *(c_eq_liq - c_eq_alpha)*dimx*fvc::div(n_alpha&n_liq*dimt*fvc::ddt(phi_alpha)*-1*n_alpha) // see this last term
                                 
                            + 0.5*0.707*epsilon *(c_eq_liq - c_eq_beta)* dimx*fvc::div(n_beta&n_liq* dimt*fvc::ddt(phi_beta)* -1*n_beta)    // see this last term


);

muEqn.solve();

 Info<< "Min/max phi_alpha:" << min(phi_alpha).value() << ' '
    << max(phi_alpha).value() << endl;

 Info<< "Min/max phi_beta:" << min(phi_beta).value() << ' '
    << max(phi_beta).value() << endl;

 Info<< "Min/max phi_liq:" << min(phi_liq).value() << ' '
    << max(phi_liq).value() << endl;

 Info<< "Min/max mu:" << min(mu).value() << ' '
    << max(mu).value() << endl;
