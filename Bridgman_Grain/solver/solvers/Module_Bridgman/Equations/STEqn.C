#include "Bridgman.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvmD2dt2.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::Bridgman100::STEqn()
{

    // Computing temperory fields

        BPhi = ( 0.5 - 0.5 * Foam::erf( Tls * (grainT-Tm) ) );

        //- Field height has value 1 in hot zone and 0 in cold zone
        tmp<volScalarField> height
        (
            pos
            (
                mesh.C().component(vector::Y)*idimx
              - (ybottom + runTime.value()*Vpull)
            )
         );

    //- Conductivity field K assigned using distribution
        //  of phase parameters (BPhi)
        tmp<volScalarField> K
        (
            ( (Ks - Kl) * BPhi ) + Kl
        );

        tmp<volScalarField> BPsiEps = BPsi + scalar(1e-6);

    const volVectorField& GBPsi_R = GBPsi();

    // LHS coefficient
    tmp<volScalarField> lhsCoeff =
        latent() * BPsiEps;


    // Neumann Boundary source
    tmp<volScalarField> boundaryNeumannTerm =
        (flag_1 == 0)
      ? mag(GBPsi_R) * (flux/dimx) * (grainT - Tinf) * (1 - height())   // case 0
      : M_GBPsi()   *              (grainT - Tinf) * (1 - height());    // case 1,2,3

    // Dirichlet Boundary condition terms
    tmp<volScalarField> boundaryDirictletTerm =
        (flag_1 > 2)
     ?
        // case 3 — precomputed
        (
            (GBPsi_R & fvc::grad(BPsi*grainT))
          - BdBPGB()
        ) * height() * K()
     :
        // case 0,1,2 — expanded form
        (
            GBPsi_R
          & (
                fvc::grad(BPsi*grainT)
              - (Bd * BPsiEps) * GBPsi_R
            )
       ) * height() * K();

    // Clear memory-1
    BPsiEps.clear();
    height.clear();

    // Assemble equation
    fvScalarMatrix grainTEqn
    (
        lhsCoeff() * fvm::ddt(grainT)
      - fvm::laplacian(BPsi*K(), grainT)
     ==
        boundaryNeumannTerm()
      - boundaryDirictletTerm()
    );

    // Clear memory-2
    K.clear();
    lhsCoeff.clear();
    boundaryNeumannTerm.clear();
    boundaryDirictletTerm.clear();

    //- Solving for T
    if (runTime.value() >= 56e-5) grainTEqn.solve();

    Info << "Temperature2 :" << max(grainT) << "/" <<  min(grainT) << endl;
}


/*
    // Boundary coefficient — coefficient of T separated from Tinf source
    // boundarySource = coeff*(T - Tinf)*(1-height)
    //                = coeff*(1-height)*T  ← implicit part
    //                - coeff*(1-height)*Tinf ← explicit source
    const tmp<volScalarField> boundaryCoeff =
        (flag_1 == 0)
      ? mag(GBPsi_R) * (flux/dimx) * (1 - height)
      : M_GBPsi()                  * (1 - height);

    // Assemble equation — boundary T term is implicit (fvm::Sp)
    // improves diagonal dominance and convergence
    fvScalarMatrix TEqn
    (
        lhsCoeff        * fvm::ddt(T)
      - fvm::laplacian(BPsi*K, T)
      + fvm::Sp(boundaryCoeff, T)          // implicit — strengthens diagonal
     ==
        boundaryCoeff * Tinf               // explicit source part of boundary term
      - interfaceSource
    );
*/









