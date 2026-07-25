#pragma once

#include <cmath>
#include "utilities/Input_parameters.hpp"
#include "utilities/microsim.hpp"
#include "comm/MPI_Info.hpp"


namespace microsim{
    class Solver{
        private:
            microsim::Info *global_fd;

            mpi::mpi_comm *comm;

            int nPhase;
            int nComp;
            int current_step = 0;

            // ================================================================
            // Multi-pass intermediate buffers
            // ================================================================
            double *term1_buffer = nullptr;      // Gradient energy term
            double *term2_buffer = nullptr;      // Double-well potential term
            double *term3_buffer = nullptr;      // Thermodynamic driving force
            double *divphi_buffer = nullptr;     // Divergence of phi (Laplacian)
            double *phase_comp_buffer = nullptr; // Pre-computed phase compositions

            template <int Dim>
            void solveCsClEq();

            template <int Dim>
            void solve_flux();

            template <int Dim>
            void solve_Phasefield();

            template <int Dim>
            void solve_Catr();

            template <int Dim>
            void solve_diffusion();

            template <int Dim>
            void add_noise(int time_step);

            template <int Dim>
            void solve_temperature();

            sycl::queue &q;

        public:
            Solver(microsim::Info *global_fd_, mpi::mpi_comm *comm_, sycl::queue &q_);

            template <int Dim>
            void evolve_system();


            ~Solver();
    };
}
