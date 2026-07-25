#pragma once

#include <mpi.h>
#include <sycl/sycl.hpp>

#include "utilities/Input_parameters.hpp"
#include "utilities/Logger.hpp"
#include "utilities/microsim.hpp"

namespace mpi{

    class mpi_comm{
        private:
            int argc;
            char **argv;

            microsim::Info *global_param;
            sycl::queue &dev_q;
            IO::Logger &logger;

            int rank;
            int size;
            int shm_rank;
            int shm_size;

            MPI_Comm comm_world;
            MPI_Comm comm_shm;

            MPI_Win win_send;

            double *my_send_buf_base = nullptr;
            std::vector<double*> neighbor_ptr_cache;

            // Helper functions
            bool is_local(int t);

            double* get_my_send_ptr(size_t offset);

            void copy_from_neighbor(int t, size_t off, double* d, size_t c);

            void node_barrier();

            template<int Dim>
            void apply_physical_bcs(double *data, int offset_fields, int num_fields);
        
        public:
            // constructor
            mpi_comm(int num_arg, char** args, microsim ::Info *g_parameters, sycl::queue &q, IO::Logger &log);

            MPI_Comm comm_cart;

            template<int Dim>
            void exchange_boundary(std::vector<FieldRange> ranges);

            // distructor
            ~mpi_comm();
    };

}