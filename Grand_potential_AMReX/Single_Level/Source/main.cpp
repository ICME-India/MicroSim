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
#include "Stress_solver.h"
#include "Merge_init.h"


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
	Vector<int> bc_lo_comp(AMREX_SPACEDIM,0);			//Boundary condition for composition

	Vector<int> bc_lo_disp(AMREX_SPACEDIM,0);
	Vector<int> bc_lo_str(AMREX_SPACEDIM,0);

	//bc_hi is the boundary condition at the upper edge(eg. X+, Y+, Z+)-------------------------------------------------------------------------
    Vector<int> bc_hi_phi(AMREX_SPACEDIM,0);			//Boundary condition for phi
	Vector<int> bc_hi_mu(AMREX_SPACEDIM,0);				//Boundary condition for chemical potential
	Vector<int> bc_hi_comp(AMREX_SPACEDIM,0);			//Boundary condition for composition

	Vector<int> bc_hi_disp(AMREX_SPACEDIM,0);
	Vector<int> bc_hi_str(AMREX_SPACEDIM,0);

	//Reading the boundary condition from the infile into bc_lo and bc_hi for phi, mu and comp-------------------------------------------------------------------------

	for(int w=0; w<bound.size(); w++){
		if(bound[w][0]=="phi"){
			bc_hi_phi[X] = stod(bound[w][1]);
			bc_hi_phi[Y] = stod(bound[w][3]);
			bc_lo_phi[X] = stod(bound[w][2]);
			bc_lo_phi[Y] = stod(bound[w][4]);

			#if(AMREX_SPACEDIM>2)
			bc_hi_phi[Z] = stod(bound[w][5]);
			bc_lo_phi[Z] = stod(bound[w][6]);
			#endif

			if (bc_hi_phi[X] == 3){
				PERIODIC_BC=1;
			}
		}

		else if(bound[w][0]=="mu"){
			bc_hi_mu[X] = stod(bound[w][1]);
			bc_hi_mu[Y] = stod(bound[w][3]);
			bc_lo_mu[X] = stod(bound[w][2]);
			bc_lo_mu[Y] = stod(bound[w][4]);

			#if(AMREX_SPACEDIM>2)
			bc_hi_mu[Z] = stod(bound[w][5]);
			bc_lo_mu[Z] = stod(bound[w][6]);
			#endif
		}

		else if(bound[w][0]=="c"){
			bc_hi_comp[X] = stod(bound[w][1]);
			bc_hi_comp[Y] = stod(bound[w][3]);
			bc_lo_comp[X] = stod(bound[w][2]);
			bc_lo_comp[Y] = stod(bound[w][4]);

			#if(AMREX_SPACEDIM>2)
			bc_hi_comp[Z] = stod(bound[w][5]);
			bc_lo_comp[Z] = stod(bound[w][6]);
			#endif
		}

		if(ELASTICITY){							//ONLY SUPPORTS PERIODIC BOUNDARY CONDITION
			bc_hi_disp[X] = 3;
			bc_hi_disp[Y] = 3;
			bc_lo_disp[X] = 3;
			bc_lo_disp[Y] = 3;
			#if(AMREX_SPACEDIM>2)
			bc_hi_disp[Z] = 3;
			bc_lo_disp[Z] = 3;
			#endif

			bc_hi_str[X] = 3;
			bc_hi_str[Y] = 3;
			bc_lo_str[X] = 3;
			bc_lo_str[Y] = 3;
			#if(AMREX_SPACEDIM>2)
			bc_hi_str[Z] = 3;
			bc_lo_str[Z] = 3;
			#endif
		}
	}
	
	//Parameter for periodic boundary condition-------------------------------------------------------------------------
	Vector<int> is_periodic(AMREX_SPACEDIM,0);
	for (int idim=0; idim < AMREX_SPACEDIM; ++idim)
	{
		if (PERIODIC_BC)
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
	 MultiFab term4(ba, dm, nump, ghost);
	 MultiFab psi(ba, dm, nump, ghost);
	 MultiFab lambad(ba, dm, nump, ghost);
	 MultiFab disp_X(ba,dm,3,ghost);
	 MultiFab disp_Y(ba,dm,3,ghost);
	 MultiFab disp_Z(ba,dm,3,ghost);
	 MultiFab strain_X(ba,dm,6,ghost);
	 MultiFab strain_Y(ba,dm,6,ghost);
	 MultiFab strain_Z(ba,dm,6,ghost);
	 MultiFab merge_phi(ba,dm,1,1);
	 //MultiFab print(ba,dm,nump+2*(numcom-1),0);
	MultiFab print(ba,dm,nump+2*numcom+2,0);

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
	 term4.setVal(0.0);
	 disp_X.setVal(0.0);
	 disp_Y.setVal(0.0);
	 disp_Z.setVal(0.0);
	 strain_X.setVal(0.0);
	 strain_Y.setVal(0.0);
	 strain_Z.setVal(0.0);
	 print.setVal(0.0);
	 merge_phi.setVal(0.0);


    //Define a vector to store boundary condition-------------------------------------------------------------------------
	Vector<BCRec> bc_phi(phi_old.nComp());
	Vector<BCRec> bc_mu(mu_old.nComp());
	Vector<BCRec> bc_comp(comp_old.nComp());
	Vector<BCRec> bc_disp(disp_X.nComp());
	Vector<BCRec> bc_str(strain_X.nComp());
	
    //Implementing boundary condtion as read from the infile-------------------------------------------------------------------------
	bound_cond(phi_old, bc_phi,bc_hi_phi,bc_lo_phi);
	bound_cond(mu_old, bc_mu,bc_hi_mu,bc_lo_mu);
	bound_cond(comp_old, bc_comp,bc_hi_comp,bc_lo_comp);
	

	if(ELASTICITY){
	bound_cond(disp_X, bc_disp,bc_hi_disp,bc_lo_disp);
	//bound_cond(disp_Y, bc_disp,bc_hi_disp,bc_lo_disp);

	bound_cond(strain_X, bc_str,bc_hi_str,bc_lo_str);
	//bound_cond(strain_Y, bc_str,bc_hi_str,bc_lo_str);

	#if(AMREX_SPACEDIM>2)
		//bound_cond(disp_Z, bc_disp,bc_hi_disp,bc_lo_disp);
		//bound_cond(strain_Z, bc_str,bc_hi_str,bc_lo_str);
	#endif
	}

	//Initialise the function pointers and fill phi-------------------------------------------------------------------------
	init_functions(phi_new);
	
	//Calling the appropriate functions-------------------------------------------------------------------------
	function_A();					//Function A points to "function_F_04_function_A" which reads the data from csv file and updates the value of A
	function_B();					//Function B points to "function_F_04_function_B" and calculates the value of B
	function_C();					//Function C points to "function_F_04_function_C" and calculates the value of C
	Init_mu(mu_new,phi_new);						//Function Mu initialises the chemical potential
	dc_dmu();						//Calculates dc/dmu
	Init_comp(phi_new,comp_new);	//Initialise the composition
	Calculate_Tau(phi_new);			//Calculate Tau
	
    //Read the parameters from Checkpoint file if restarting a simulation----------------------------------------------------------------------------------
	if(restart == 1 && restart_chkfile != ""){
		Readchkfile(phi_new, mu_new, comp_new);
	} 
	
	if(nump>2){
		merge_initial(merge_phi,phi_new,nump);
	}

    //Copy phi_new and mu_new to print so that the solution can be exported----------------------------------------------------------------------------------
	MultiFab::Copy(print, phi_new, 0, 0, nump, 0);
	MultiFab::Copy(print, mu_new, 0, nump, numcom-1, 0);
	MultiFab::Copy(print, comp_new, 0, nump+numcom-1, numcom-1, 0);
	
	if(ELASTICITY){
	 MultiFab::Copy(print,disp_X,2,nump+2*numcom-2,1,0);
	 MultiFab::Copy(print,disp_Y,2,nump+2*numcom-1,1,0);
	#if(AMREX_SPACEDIM>2)
		MultiFab::Copy(print,disp_Z,2,nump+2*numcom,1,0);
	#endif
	}
	if(ELASTICITY==1 && nump>2){
		MultiFab::Copy(print, merge_phi, 0, nump+2*numcom, 1, 0);
	}
	if(ELASTICITY==0 && nump>2){
		MultiFab::Copy(print, merge_phi, 0, nump+2*numcom-2, 1, 0);
	}

	//Print()<<"1\n";
    //This vector is used to name the components of print multifab---------------------------------------------------------------------------------- 
	for(int m=0; m<numcom-1; m++){
		phase.push_back("mu_"+comp[m]);
	}
	
	for(int m=0; m<numcom-1; m++){
		phase.push_back("comp_"+comp[m]);
	}

	if(ELASTICITY){
		phase.push_back("Ux");
		phase.push_back("Uy");
	#if(AMREX_SPACEDIM>2)
		phase.push_back("Uz");
	#endif
	}

	if(nump>2){
		phase.push_back("Merged");
	}


	stiffness_n = Vector<Vector<Real>>(nump,Vector<Real>(voigiso[0].size(),0.0));

	if(ELASTICITY){
		for(int a=0; a<nump; a++){
			stiffness_n[a][0] = a;
			stiffness_n[a][1] = voigiso[a][1]/voigiso[nump-1][3];
			stiffness_n[a][2] = voigiso[a][2]/voigiso[nump-1][3];
			stiffness_n[a][3] = voigiso[a][3]/voigiso[nump-1][3]; 
		}
	}


    //Plot the initial file---------------------------------------------------------------------------------- 
	if(savet>0 && restart_chkfile == "")
	 {
	 	const std::string& pltfile  = amrex::Concatenate("plt",0,1);
	 	WriteSingleLevelPlotfile(pltfile, print, phase,geom,timee,0);

	 }
	
    //Iterating loop to calculate the order parameters----------------------------------------------------------------------------------  
	for(n=stepnum+1; n<=nsteps; ++n)
	 {  
        //Copy phi_new to phi_old and use phi_old for further iterations----------------------------------------------------------------------------------   
	 	MultiFab::Copy(phi_old, phi_new, 0,0,nump,0);

        //Copy mu_new to mu_old and use mu_old for further iterations----------------------------------------------------------------------------------
	 	MultiFab::Copy(mu_old, mu_new, 0,0,numcom-1,0);

		//Copy mu_new to comp_old and use comp_old for further iterations----------------------------------------------------------------------------------
		MultiFab::Copy(comp_old, comp_new, 0,0,numcom-1,0);

        //Filling the boundary cells for phi_old and mu_old----------------------------------------------------------------------------------
	 	FillDomainBoundary(phi_old, geom, bc_phi);
		FillDomainBoundary(mu_old, geom, bc_mu);
		FillDomainBoundary(comp_old, geom, bc_comp);
		
		if(ELASTICITY == 1){
			FillDomainBoundary(disp_X, geom, bc_disp);
			FillDomainBoundary(disp_Y, geom, bc_disp);
		#if(AMREX_SPACEDIM>2)
			FillDomainBoundary(disp_Z, geom, bc_disp);
		#endif
		}


        //Advance the solution to calculate field variables for next time step----------------------------------------------------------------------------------
	 	advance(phi_old, phi_new, mu_old, mu_new, comp_old, comp_new, term1, term2, term3, term4, disp_X, disp_Y , disp_Z, psi, lambad, geom);
	 	

		if(ELASTICITY==1){
		for(int u=0; u<MAX_ITERATIONS; u++){

		FillDomainBoundary(strain_X, geom, bc_str);
		FillDomainBoundary(strain_Y, geom, bc_str);
		FillDomainBoundary(disp_X, geom, bc_disp);
		FillDomainBoundary(disp_Y, geom, bc_disp);
		#if(AMREX_SPACEDIM>2)
			FillDomainBoundary(strain_Z, geom, bc_disp);
			FillDomainBoundary(disp_Z, geom, bc_disp);
		#endif
		
		Iterative_stress_solver(phi_old,disp_X,disp_Y,disp_Z,strain_X, strain_Y,strain_Z,geom);
		
		}
		}

        //Update the time----------------------------------------------------------------------------------
	 	timee=timee+dt;

        //Print the step number---------------------------------------------------------------------------------- 
		if(n%trackprog==0){

		 	amrex::Print()<<"Current step: "<<n<<" , Current time: "<<timee<<"\n";
			Print()<<"\n";
			
			for(int w=0; w<nump; w++){
				Print()<<"Phi "<<phase[w]<<": MAX = "<<phi_new.max(w,0,0)<<" , MIN = "<<phi_new.min(w,0,0)<<"\n";
			}
			for(int w=0; w<numcom-1; w++){
				Print()<<"Mu "<<phase[nump+w]<<": MAX = "<<mu_new.max(w,0,0)<<" , MIN = "<<mu_new.min(w,0,0)<<"\n";
			}
			for(int w=0; w<numcom-1; w++){
				Print()<<"Comp "<<phase[nump+numcom-1+w]<<": MAX = "<<comp_new.max(w,0,0)<<" , MIN = "<<comp_new.min(w,0,0)<<"\n";
			}
			
			Print()<<"------------------------------------------------\n";

		}
    
        //Plot the updated data----------------------------------------------------------------------------------
	 	if(savet>0 && n%savet==0)
	 	{	
			merge_phi.setVal(0.0);
			merge_initial(merge_phi,phi_new,nump);
			//Copy the updated data to print----------------------------------------------------------------------------------
			MultiFab::Copy(print, phi_new, 0, 0, nump, 0);
			MultiFab::Copy(print, mu_new, 0, nump, numcom-1, 0);
			MultiFab::Copy(print, comp_new, 0, nump+numcom-1, numcom-1, 0);
			
			if(ELASTICITY){
				MultiFab::Copy(print,disp_X,2,nump+2*numcom-2,1,0);
				MultiFab::Copy(print,disp_Y,2,nump+2*numcom-1,1,0);
				#if(AMREX_SPACEDIM>2)
					MultiFab::Copy(print,disp_Z,2,nump+2*numcom,1,0);
				#endif
			}

			if(ELASTICITY==1 && nump>2){
				MultiFab::Copy(print, merge_phi, 0, nump+2*numcom, 1, 0);
			}
			if(ELASTICITY==0 && nump>2){
				MultiFab::Copy(print, merge_phi, 0, nump+2*numcom-2, 1, 0);
			}

	 		const std::string& pltfile = amrex::Concatenate("plt",n,1);
	 		WriteSingleLevelPlotfile( pltfile, print, phase,geom,timee,n);
	 	}

        //Save the checkpoint file----------------------------------------------------------------------------------
		//if(n%(50*savet) == 0){
            //Writechkfile(phi_new, mu_new, comp_new, n, chk_file);
        //}
	 	
	 }

}
