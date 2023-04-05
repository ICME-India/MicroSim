const dictionary& get_Tol_from_this = mesh.solutionDict().subDict("Tol_is_defined_here");
//! Initial tolerance "Tol" within subDict "Tol_is_defined_here" in fvSolution
scalar Tol(readScalar(get_Tol_from_this.lookup("Tol")));

const dictionary& get_dtf_from_this = mesh.solutionDict().subDict("deltaTFactor");
//! Time step changing factor "dtf" within subDict "deltaTFactor" in fvSolution
scalar dtf(readScalar(get_dtf_from_this.lookup("dtf")));
