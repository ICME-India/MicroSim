#include <AMReX_Gpu.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_BCRec.H>
#include <AMReX_BCUtil.H>
#include <AMReX_Utility.H>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

#include "Variables.h"
#include "Tau.h"
#include "Advance_sol.h"
#include "Read_input.h"
#include "Boundary_conditions.h"
#include "Filling_fv.h"
#include "Initialise_functions.h"
#include "Function_F4.h"
#include "Chkpnt.h"
#include "Calc_part.h"


using namespace std;

using namespace amrex; 

void GPotential();

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

	{	
		BL_PROFILE("main()");	
		
        //Start time of the simulation-------------------------------------------------------------------------
		strt_time = ParallelDescriptor::second();
		
		//Read the input parameters from the infile-------------------------------------------------------------------------
    	readinput(); 
		
		//Calling the Grand Potential function-------------------------------------------------------------------------  
		GPotential();

        //Calculating end time-------------------------------------------------------------------------
		stop_time = ParallelDescriptor::second()-strt_time+rst_time;
	 	const int IOProc = ParallelDescriptor::IOProcessorNumber();
	 	ParallelDescriptor::ReduceRealMax(stop_time, IOProc);

        //Print stop time------------------------------------------------------------------------------
	 	amrex::Print()<<"Run time = "<<stop_time<<"\n";
	}

    amrex::Finalize();
}


void GPotential()
{	
	//Declare maxgrid variables
	int maxy_x{0};
	int maxy_y{0};
	int maxy_z{0};

	//Calculating number of partions of the Box depending on the number of procs-------------------------------------------------------------------------
	Calc_partition(maxy_x,maxy_y,maxy_z);
	
	//Defining maxgrid
    IntVect maxgrid(AMREX_D_DECL(maxy_x,maxy_y,maxy_z));

    //Restart from a pre-existing checkpoint file if advised so-------------------------------------------------------------------------
	if(restart==1){
    restart_chkfile = chk_file + to_string(startt);
    }

	//bc_lo is the boundary condition at the lower edge(eg. X-, Y-, Z-)-------------------------------------------------------------------------
	Vector<int> bc_lo_phi(AMREX_SPACEDIM,0);			//Boundary condition for phi
	Vector<int> bc_lo_mu(AMREX_SPACEDIM,0);				//Boundary condition for chemical potential
	//Vector<int> bc_lo_comp(AMREX_SPACEDIM,0);			//Boundary condition for composition

	//bc_hi is the boundary condition at the upper edge(eg. X+, Y+, Z+)-------------------------------------------------------------------------
    Vector<int> bc_hi_phi(AMREX_SPACEDIM,0);			//Boundary condition for phi
	Vector<int> bc_hi_mu(AMREX_SPACEDIM,0);				//Boundary condition for chemical potential
	//Vector<int> bc_hi_comp(AMREX_SPACEDIM,0);			//Boundary condition for composition

	//Reading the boundary condition from the infile into bc_lo and bc_hi for phi, mu and comp-------------------------------------------------------------------------
	bc_hi_phi[X] = stod(bound[0][1]);
	bc_hi_phi[Y] = stod(bound[0][3]);
	bc_lo_phi[X] = stod(bound[0][2]);
	bc_lo_phi[Y] = stod(bound[0][4]);

	#if(AMREX_SPACEDIM>2)
	bc_hi_phi[Z] = stod(bound[0][5]);
	bc_lo_phi[Z] = stod(bound[0][6]);
	#endif

	bc_hi_mu[X] = stod(bound[1][1]);
	bc_hi_mu[Y] = stod(bound[1][3]);
	bc_lo_mu[X] = stod(bound[1][2]);
	bc_lo_mu[Y] = stod(bound[1][4]);

	#if(AMREX_SPACEDIM>2)
	bc_hi_mu[Z] = stod(bound[1][5]);
	bc_lo_mu[Z] = stod(bound[1][6]);
	#endif

	// bc_hi_comp[X] = stod(bound[2][1]);
	// bc_hi_comp[Y] = stod(bound[2][3]);
	// bc_lo_comp[X] = stod(bound[2][2]);
	// bc_lo_comp[Y] = stod(bound[2][4]);

	// #if(AMREX_SPACEDIM>2)
	// bc_hi_comp[Z] = stod(bound[2][5]);
	// bc_lo_comp[Z] = stod(bound[2][6]);
	// #endif
	
	//Parameter for periodic boundary condition-------------------------------------------------------------------------
	Vector<int> is_periodic(AMREX_SPACEDIM,0);
	for (int idim=0; idim < AMREX_SPACEDIM; ++idim)
	{
		if (bc_lo_phi[idim] == 3 && bc_hi_phi[idim] == 3)
		{
		    is_periodic[idim] = 1;
		}
	}

	//Making sure ncellz = 1 for a 2D simulation-------------------------------------------------------------------------
	if(AMREX_SPACEDIM==2 && ncellz!=1){
		ncellz=1;
	}

    //Define BoxArray and Geometry-------------------------------------------------------------------------
	BoxArray ba;
	Geometry geom;
		//Lower Boundary of the domain-------------------------------------------------------------------------
		IntVect dom_lo(AMREX_D_DECL(0,0,0));

		//Upper Boundary of the domain-------------------------------------------------------------------------
		IntVect dom_high(AMREX_D_DECL(ncellx-1,ncelly-1,ncellz-1));
		
        //Initialise the domain-------------------------------------------------------------------------
		Box domain(dom_lo,dom_high);
		
		//Real size of the box-------------------------------------------------------------------------
		RealBox real_box({AMREX_D_DECL(Real(0.0),Real(0.0),Real(0.0))},{AMREX_D_DECL(Real(ncellx*dx), Real(ncelly*dy), Real(ncellz*dz))});
		
        //Define the doamin, box, coordinate sytem and boundary condition of the geometry-------------------------------------------------------------------------
	 	geom.define(domain, &real_box, CoordSys::cartesian, is_periodic.data()); 


    //Using pre-defined BoxArray from the checkpoint file if it exists else defining a new BoxArray-------------------------------------------------------------------------
		if(restart == 1 && n>0){
            ba = grids[0];
    	}
    	else{
		    ba.define(domain);
		    ba.maxSize(maxgrid);
		}

	//Printing the box array
	Print()<<ba<<"\n";

	//Number of ghost cells for the Multifab-------------------------------------------------------------------------
	int ghost=1;

	//Mapping processor to the partitions in box array------------------------------------------------------------------------- 
	 DistributionMapping dm(ba, ParallelDescriptor::NProcs());

	//Declaring Multifabs-------------------------------------------------------------------------
	 MultiFab phi_old(ba, dm, nump, ghost);	
	 MultiFab phi_new(ba, dm, nump, ghost);
	 MultiFab mu_old(ba, dm, numcom-1, ghost);
	 MultiFab mu_new(ba, dm, numcom-1, ghost);
	 MultiFab comp_old(ba, dm, numcom-1, ghost);
	 MultiFab comp_new(ba, dm, numcom-1, ghost);
	 MultiFab term1(ba, dm, nump, ghost);
	 MultiFab term2(ba, dm, nump, ghost);
	 MultiFab term3(ba, dm, nump, ghost);
	 MultiFab psi(ba, dm, nump, ghost);
	 MultiFab lambad(ba, dm, nump, ghost);
	 //MultiFab print(ba,dm,nump+2*(numcom-1),0);
	MultiFab print(ba,dm,nump+numcom-1,0);

	//Initialsing MultiFabs-------------------------------------------------------------------------
	 phi_old.setVal(0.0);
	 phi_new.setVal(0.0);
	 mu_old.setVal(0.0);
	 mu_new.setVal(0.0);
	 comp_old.setVal(0.0);
	 comp_new.setVal(0.0);
	 term1.setVal(0.0);
	 term2.setVal(1.0);
	 term3.setVal(0.0);
	 print.setVal(0.0);


    //Define a vector to store boundary condition-------------------------------------------------------------------------
	Vector<BCRec> bc_phi(phi_old.nComp());
	Vector<BCRec> bc_mu(mu_old.nComp());
	//Vector<BCRec> bc_comp(comp_old.nComp());
	
    //Implementing boundary condtion as read from the infile-------------------------------------------------------------------------
	bound_cond(phi_old, bc_phi,bc_hi_phi,bc_lo_phi);
	bound_cond(mu_old, bc_mu,bc_hi_mu,bc_lo_mu);
	//bound_cond(comp_old, bc_comp,bc_hi_comp,bc_lo_comp);

	//Initialise the function pointers and fill phi-------------------------------------------------------------------------
	init_functions(phi_new);
	
	//Calling the appropriate functions-------------------------------------------------------------------------
	function_A();					//Function A points to "function_F_04_function_A" which reads the data from csv file and updates the value of A
	function_B();					//Function B points to "function_F_04_function_B" and calculates the value of B
	function_C();					//Function C points to "function_F_04_function_C" and calculates the value of C
	Mu(mu_new);						//Function Mu initialises the chemical potential
	dc_dmu();						//Calculates dc/dmu
	//Init_comp(phi_new,comp_new);	//Initialise the composition
	Calculate_Tau(phi_new);			//Calculate Tau
	
    //Read the parameters from Checkpoint file if restarting a simulation----------------------------------------------------------------------------------
	if(restart == 1 && restart_chkfile != ""){
		Readchkfile(phi_new, mu_new);
	} 
	
    //Copy phi_new and mu_new to print so that the solution can be exported----------------------------------------------------------------------------------
	MultiFab::Copy(print, phi_new, 0, 0, nump, 0);
	MultiFab::Copy(print, mu_new, 0, nump, numcom-1, 0);
	//MultiFab::Copy(print, comp_new, 0, nump+numcom-1, numcom-1, 0);
	//Print()<<"1\n";
    //This vector is used to name the components of print multifab---------------------------------------------------------------------------------- 
	for(int m=0; m<numcom-1; m++){
		phase.push_back("mu_"+comp[m]);
	}

	// for(int m=0; m<numcom-1; m++){
	// 	phase.push_back("comp_"+comp[m]);
	// }


    //Plot the initial file---------------------------------------------------------------------------------- 
	if(trackprog>0 && restart_chkfile == "")
	 {
	 	const std::string& pltfile  = amrex::Concatenate("plt",0,1);
	 	WriteSingleLevelPlotfile(pltfile, print, phase,geom,timee,0);
		//Writechkfile(phi_new, mu_new,n,chk_file);
	 }
	
    //Iterating loop to calculate the order parameters----------------------------------------------------------------------------------  
	for(n=stepnum+1; n<=nsteps; ++n)
	 {  
        //Copy phi_new to phi_old and use phi_old for further iterations----------------------------------------------------------------------------------   
	 	MultiFab::Copy(phi_old, phi_new, 0,0,nump,0);

        //Copy mu_new to mu_old and use mu_old for further iterations----------------------------------------------------------------------------------
	 	MultiFab::Copy(mu_old, mu_new, 0,0,numcom-1,0);

		//Copy mu_new to comp_old and use comp_old for further iterations----------------------------------------------------------------------------------
		//MultiFab::Copy(comp_old, comp_new, 0,0,numcom-1,0);

        //Filling the boundary cells for phi_old and mu_old----------------------------------------------------------------------------------
	 	FillDomainBoundary(phi_old, geom, bc_phi);
		FillDomainBoundary(mu_old, geom, bc_mu);
		//FillDomainBoundary(comp_old, geom, bc_comp);

		//Print()<<"2\n";
        //Advance the solution to calculate field variables for next time step----------------------------------------------------------------------------------
	 	advance(phi_old, phi_new, mu_old, mu_new, comp_old, comp_new, term1, term2, term3, psi,lambad,geom);
	 	
        //Update the time----------------------------------------------------------------------------------
	 	timee=timee+dt;

        //Print the step number---------------------------------------------------------------------------------- 
		if(n%trackprog==0){
		 	amrex::Print()<<"Advanced step"<<n<<"\n";
			Print()<<"\n";
			Print()<<"anisotropy_max = "<<term1.max(0,0,0)*eps<<"\n";
			Print()<<"anisotropy_min = "<<term1.min(0,0,0)*eps<<"\n";
			Print()<<"doublewell_max = "<<term2.max(0,0,0)/eps<<"\n";
			Print()<<"doublewell_min = "<<term2.min(0,0,0)/eps<<"\n";
			Print()<<"GrandPotFunc_max = "<<term3.max(0,0,0)/Vm<<"\n";
			Print()<<"GrandPotFunc_min = "<<term3.min(0,0,0)/Vm<<"\n";
			Print()<<"mu_max = "<<mu_old.max(0,0,0)<<"\n";
			Print()<<"mu_min = "<<mu_old.min(0,0,0)<<"\n";
			Print()<<"------------------------------------------------\n";
		}
		
        //Copy the updated data to print----------------------------------------------------------------------------------
		MultiFab::Copy(print, phi_new, 0, 0, nump, 0);
		MultiFab::Copy(print, mu_new, 0, nump, numcom-1, 0);
		//MultiFab::Copy(print, comp_new, 0, nump+numcom-1, numcom-1, 0);
    
        //Plot the updated data----------------------------------------------------------------------------------
	 	if(trackprog>0 && n%trackprog==0)
	 	{
	 		const std::string& pltfile = amrex::Concatenate("plt",n,1);
	 		WriteSingleLevelPlotfile( pltfile, print, phase,geom,timee,n);
	 	}

        //Save the checkpoint file----------------------------------------------------------------------------------
		if(n%savet == 0){
            Writechkfile(phi_new, mu_new, n, chk_file);
        }
	 	
	 }

}
