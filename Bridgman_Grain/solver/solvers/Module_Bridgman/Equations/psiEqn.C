#include "Bridgman.H"
#include "fvcLaplacian.H"
#include "fvmD2dt2.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::Bridgman100::psiEqn()
{

    if (runTime.value() < 5e-4)
    {
        //- Domain parameter (BPsi) equation

        if(nPsiSubCycles > 1)
        {

            List<volScalarField*> state1FieldPtrs({&BPsi});


            for
	    (
                subCycle<volScalarField, subCycleFields> stateFieldSubCycle1
                (
             	    state1FieldPtrs,
                    nPsiSubCycles
                );
                !(++stateFieldSubCycle1).end();
            )
	    {
                solvePsiEqn();
            }

        }
        else
        {
            solvePsiEqn();
        }

        //- BPsi extrema display
        Info << "Min/max psi:" << min(BPsi()).value()
             << ' ' << max(BPsi()).value()  << endl;

        Info << "BPsisum :" << sum(BPsi)
             << "smoothtime :" << runTime.deltaTValue() << endl;

    }
    else
    {
     	//- Dummy matrix (for proper dynamic mesh functioning)
        fvScalarMatrix psiEqn_eq
        (
             fvm::laplacian(BPsi)
        );
    }
}

void Foam::solvers::Bridgman100::solvePsiEqn()
{
    fvScalarMatrix psiEqn_eq
    (
        omega * epsilon * dimt * fvm::ddt(BPsi)
      - 2.0 * epsilon * gamma * dimx * dimx * fvm::laplacian(BPsi)
      + 18.0 * gamma * BPsi * (BPsi - 1.0) * (2.0 * BPsi - 1.0) / epsilon
    );

    psiEqn_eq.solve();
}
