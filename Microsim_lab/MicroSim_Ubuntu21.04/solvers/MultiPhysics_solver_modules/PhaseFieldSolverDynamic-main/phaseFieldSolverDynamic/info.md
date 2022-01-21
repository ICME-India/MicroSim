@mainpage Code-guide for solver phaseFieldDynamic

The solver requires *dynamicInterfaceRefineFvMesh* library, and has been successfully tested using OpenFOAM v6. Before compiling the solver, the library must be obtained in the following procedure:

* You can compile the lib where ever you want. This is just an example:
> mkdir -p $FOAM_RUN/../OpenFOAM_extensions

* Switch to this directory
> cd $FOAM_RUN/../OpenFOAM_extensions

* Clone the repository
> git clone https://bitbucket.org/shor-ty/dynamicinterfacerefinefvmesh.git

* Move to the new library folder
> cd dynamicinterfacerefinefvmesh

* Checkout the openfoam version you need (e.g. using 5.x)
> git checkout OpenFOAM-5.x

* Go to the library source
> cd src/dynamicFvMesh

* Compile the lib
> wmake libso

* Finally you have to include the libary into your solver set-up. Therefore add the following line to your dynamicMeshDict
> dynamicFvMeshLibs ( "libdynamicInterfaceRefineFvMesh.so" );

* The best way is to copy the dynamicMeshDict to your case and modify it.

@section s1 Compiling the solver

* Following commands should create the executable of the solver
> cd $FOAM_RUN/phaseFieldSolverDynamic/phaseFieldSolverDynamic

> wclean

> wmake

* The solver can be run by following the instructions in *userGuide*.

@section s2 Further details

The implementation, client and header files of the solver have been written following OpenFOAM conventions. These are explained next with flow charts generated from the source code using Doxygen. It must be noted that the solver is based on [laplacianFoam](https://www.openfoam.com/documentation/guides/latest/doc/guide-applications-solvers-basic-laplacianFoam.html "laplacianFoam") solver within OpenFOAM. Hence, it may be helpful for the user to become familiar with [OpenFOAM Programmer’s Guide](http://foam.sourceforge.net/docs/Guides-a4/ProgrammersGuide.pdf "OpenFOAM Programmer’s Guide") and [laplacianFoam](https://www.openfoam.com/documentation/guides/latest/doc/guide-applications-solvers-basic-laplacianFoam.html "laplacianFoam") beforehand.
