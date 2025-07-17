#include <iostream>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

#include<AmrCoreGP.H>

using namespace std;
using namespace amrex; 

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
	auto strt_time = ParallelDescriptor::second();
    AmrCoreGP amr_core_gp;
    amr_core_gp.InitData();    
    amr_core_gp.Evolve();
	auto stop_time = ParallelDescriptor::second()-strt_time;
	const int IOProc = ParallelDescriptor::IOProcessorNumber();
	ParallelDescriptor::ReduceRealMax(stop_time, IOProc);
	amrex::Print()<<"Run time = "<<stop_time<<"\n";
    amrex::Finalize();
}
