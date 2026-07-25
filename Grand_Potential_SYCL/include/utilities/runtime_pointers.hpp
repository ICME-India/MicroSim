#pragma once
#include <mpi.h>
#include "utilities/IO.hpp"
#include "physics/PF_model.hpp"

namespace microsim {
    using write_func_ptr = void (IO::Read_input::*)(int, int, MPI_Comm);
    using evolve_func_ptr = void (microsim::Solver::*)();

    struct RuntimePointers {
        write_func_ptr write_hdf5;
        write_func_ptr write_vtk;
        evolve_func_ptr evolve_system;

        RuntimePointers(int Dim) {
             if (Dim == 2) {
                write_hdf5 = &IO::Read_input::write_hdf5<2>;
                write_vtk = &IO::Read_input::write_vtk<2>;
                evolve_system = &microsim::Solver::evolve_system<2>;
            } else if (Dim == 3) {
                write_hdf5 = &IO::Read_input::write_hdf5<3>;
                write_vtk = &IO::Read_input::write_vtk<3>;
                evolve_system = &microsim::Solver::evolve_system<3>;
            } else {
                write_hdf5 = nullptr;
                write_vtk = nullptr;
                evolve_system = nullptr;
            }
        }
    };
}
