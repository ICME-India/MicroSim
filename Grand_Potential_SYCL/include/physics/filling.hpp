#pragma once

#include <mpi.h>
#include <map>

#include "utilities/microsim.hpp"
#include "utilities/Input_parameters.hpp"
#include "utilities/Logger.hpp"

namespace microsim{
    class filling{
        private:
            microsim::Info *global_fd;
            IO::Logger &logger;
            sycl::queue dev_q;

            int nPhase;
            int nComp;
            int Dimension;

            template<int Dim>
            void Initialize_zero();

            template<int Dim>
            void fill_cube(int phase, int count);

            template<int Dim>
            void fill_sphere(int phase, int count);

            template<int Dim>
            void fill_composition();

            template<int Dim>
            void fill_temperature();

        public :
            filling(microsim::Info *global_param, sycl::queue &q, IO::Logger &logg);

            void fill_domain();

            ~filling();
    };
}