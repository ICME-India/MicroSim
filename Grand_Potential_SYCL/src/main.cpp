#include <iostream>
#include <mpi.h>
#include <sycl/sycl.hpp>
#include <hdf5.h>

#include "utilities/Logger.hpp"
#include "comm/device_info.hpp"
#include "utilities/IO.hpp"
#include "utilities/microsim.hpp"
#include "comm/MPI_Info.hpp"
#include "physics/filling.hpp"
#include "utilities/helper_functions.hpp"
#include "physics/PF_model.hpp"
#include "utilities/runtime_pointers.hpp"


using std::atoi, std::cout;

int main(int argc, char*argv[]){

    if(argc < 5){
        std::cerr<<"Number of arguments provided is less than the required !>\n";
        return 1;
    }

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(rank == 0) mkdir("DATA", 0777);

    // Define the logging level 
    LogLevel opt_level = LogLevel::INFO;
    bool cmd_log = false;
    IO::Logger log("output.log", opt_level, cmd_log);

    std::string filename = argv[3];

    // Backend Device selection - Define at compile time
    Device_Info back_end(log);
    back_end.LogInfo();

    // Assigning queue to the corresponding device
    sycl::queue dev_queue(back_end.selected_device);

    microsim::Info *global_mem = new microsim::Info(dev_queue);
    IO::Read_input *input = new IO::Read_input(argv[1], argv[2], global_mem, dev_queue, filename, log);

    // Reading the infile and filling file
    input->Read_infile();
    input->Read_filling_file();


    mpi::mpi_comm *comm = new mpi::mpi_comm(argc, argv, global_mem, dev_queue, log);

    global_mem->Allocate_field();

    microsim::filling *Init = new microsim::filling(global_mem, dev_queue, log);

    std::cout<<"Starting domain filling ... \n";

    Init->fill_domain();

    input->print_infile();

    // Initiate the seed


    //input->print_filling_file();


    
    {

        microsim::Solver solver(global_mem, comm, dev_queue);

        int Dimension = global_mem->domainInfo->DIMENSION;

        microsim::RuntimePointers ptrs(Dimension);

        std::cout<<"reading runtime pointers done\n";

        for(int t=0; t<=global_mem->parameters->ntime_steps; ++t){
            //std::cout<<"Time step: " << t << "\n";

            if(t%global_mem->parameters->saveT == 0){
                (input->*ptrs.write_hdf5)(t, rank, comm->comm_cart);
                (input->*ptrs.write_vtk)(t, rank, comm->comm_cart);
            }

            //std::cout<<"Evolving system ... \n";
            (solver.*ptrs.evolve_system)();

        }
    }

    std::cout<<"Simulation completed successfully !\n";
    
    MPI_Finalize();
    return 0;
}