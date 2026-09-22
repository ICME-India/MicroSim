#include "MPIContext.hpp"
#include "Grid.hpp"
#include "FreeEnergy.hpp"
#include "Mobility.hpp"
#include "Solver.hpp"
#include "IO.hpp"
#include "Config.hpp"

#include <iostream>
#include <iomanip>
#include <memory>
#include <chrono>
#include <sys/stat.h>

void print_banner(int rank) {
    if (rank != 0) return;
    std::cout << "======================================================================\n";
    std::cout << "  Generalized Multi-Component Cahn-Hilliard (MCCH) Solver with MPI   \n";
    std::cout << "  Slab Domain Decomposition (1D Pencil Parallelization)               \n";
    std::cout << "======================================================================\n";
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    print_banner(rank);

    // Parse config
    mcch::SimulationConfig config;
    std::string config_path = "";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-c" || arg == "--config") {
            if (i + 1 < argc) {
                config_path = argv[++i];
            }
        } else if (arg == "--gpu" || arg == "-g") {
            config.device = "gpu";
        } else if (arg == "--device") {
            if (i + 1 < argc) {
                config.device = argv[++i];
            }
        } else if (arg == "-h" || arg == "--help") {
            if (rank == 0) {
                std::cout << "Usage: mpirun -np <ranks> " << argv[0] << " [-c config.json] [--gpu] [--device cpu|gpu]\n";
            }
            MPI_Finalize();
            return 0;
        } else if (arg[0] != '-') {
            config_path = arg;
        }
    }

    if (!config_path.empty()) {
        if (rank == 0) {
            std::cout << "Loading configuration from: " << config_path << "\n";
        }
        try {
            config.load_from_json(config_path);
        } catch (const std::exception& e) {
            if (rank == 0) {
                std::cerr << "Error loading configuration: " << e.what() << "\n";
            }
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    } else {
        if (rank == 0) {
            std::cout << "No config specified. Using default 3-component spinodal benchmark.\n";
        }
        config.finalize_defaults();
    }

    // Set up MPI Context and Slab Partitioning
    bool is_semi_implicit = config.is_semi_implicit();
    mcch::MPIContext mpi(MPI_COMM_WORLD, config.nx, config.ny, config.nz, config.get_bc_x(), 1, is_semi_implicit);

    // Set up Grid
    mcch::Grid grid(config.dim, config.nx, config.ny, config.nz,
                    config.dx, config.dy, config.dz, mpi,
                    config.get_bc_y(), config.get_bc_z(), 1);

    if (rank == 0) {
        std::cout << "\n[MPI Slab Decomposition]\n";
        std::cout << "  Total MPI Ranks: " << size << "\n";
        std::cout << "  Global Dimensions: " << grid.nx << " x " << grid.ny;
        if (grid.dim == 3) std::cout << " x " << grid.nz;
        std::cout << " (dim = " << grid.dim << "D)\n";
        std::cout << "  Boundary Conditions: X=" << config.bc_x 
                  << ", Y=" << config.bc_y;
        if (grid.dim == 3) std::cout << ", Z=" << config.bc_z;
        std::cout << "\n";
    }

    // Print slab assignment per rank
    for (int r = 0; r < size; ++r) {
        if (rank == r) {
            std::cout << "  Rank " << std::setw(2) << rank 
                      << " owns slab x in [" << std::setw(4) << mpi.local_x_start 
                      << " .. " << std::setw(4) << mpi.local_x_end - 1 
                      << "] (local_nx = " << mpi.local_nx 
                      << ", left=" << mpi.left_rank << ", right=" << mpi.right_rank << ")\n";
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Set up Thermodynamics (Free Energy)
    std::shared_ptr<mcch::FreeEnergy> free_energy;
    if (config.free_energy_type == "regular_solution") {
        free_energy = std::make_shared<mcch::RegularSolution>(config.num_components, config.RT, config.omega);
    } else if (config.free_energy_type == "binary_double_well") {
        free_energy = std::make_shared<mcch::BinaryDoubleWell>(1.0);
    } else if (config.free_energy_type == "Redlich_Kister" || config.free_energy_type == "redlich_kister") {
       free_energy  = std::make_shared<mcch::Redlich_Kister>(config.num_components, config.T);
    } else { // polynomial_multiwell
        free_energy = std::make_shared<mcch::PolynomialMultiWell>(config.num_components, config.W, config.A);
    }

    // Set up Mobility
    std::shared_ptr<mcch::MobilityModel> mobility;
    if (config.mobility_type == "degenerate") {
        mobility = std::make_shared<mcch::DegenerateMobility>(config.num_components, config.mobility_val);
    } else {
        if (config.num_components > 2 || !config.mobility_matrix.empty()) {
            mobility = std::make_shared<mcch::ConstantMobility>(config.num_components, config.mobility_matrix);
        } else {
            mobility = std::make_shared<mcch::ConstantMobility>(config.num_components, config.mobility_val);
        }
    }

    // Determine execution device (CPU or GPU)
    mcch::DeviceType dev_type = config.get_device_type();
    if (dev_type == mcch::DeviceType::GPU) {
#ifdef MCCH_ENABLE_CUDA
        int device_count = 0;
        cudaError_t err = cudaGetDeviceCount(&device_count);
        if (err != cudaSuccess || device_count == 0) {
            if (rank == 0) {
                std::cerr << "Warning: GPU requested, but no CUDA-capable GPU found. Falling back to CPU.\n";
            }
            dev_type = mcch::DeviceType::CPU;
            config.device = "cpu";
        } else {
            int gpu_id = rank % device_count;
            cudaSetDevice(gpu_id);
            cudaDeviceProp prop;
            cudaGetDeviceProperties(&prop, gpu_id);
            if (rank == 0) {
                std::cout << "\n[GPU Acceleration Enabled]\n";
                std::cout << "  CUDA Devices Available: " << device_count << "\n";
                std::cout << "  Rank " << rank << " mapped to GPU " << gpu_id << ": " << prop.name
                          << " (Compute " << prop.major << "." << prop.minor 
                          << ", " << (prop.totalGlobalMem / (1024 * 1024)) << " MB VRAM)\n";
            }
        }
#else
        if (rank == 0) {
            std::cerr << "Warning: GPU execution requested, but MCCH solver was compiled without CUDA. Falling back to CPU.\n";
        }
        dev_type = mcch::DeviceType::CPU;
        config.device = "cpu";
#endif
    }

    // Create Solver
    mcch::Solver solver(mpi, grid, config.num_components, free_energy, mobility, config.kappa);
    if (dev_type == mcch::DeviceType::GPU) {
        solver.set_device(mcch::DeviceType::GPU);
    }

    // Initial Condition
    if (config.initial_condition == "droplet") {
        solver.initialize_droplet(config.droplet_comp, 
                                 0.5 * grid.nx * grid.dx, 
                                 0.5 * grid.ny * grid.dy, 
                                 (grid.dim == 3 ? 0.5 * grid.nz * grid.dz : 0.0),
                                 config.droplet_radius, 
                                 config.droplet_diffuse_width);
    } else {
        solver.initialize_random(config.c_mean, config.noise_amp, config.seed);
    }

    // Time step determination
    mcch::TimeIntegratorType integrator_type = config.get_integrator_type();
    double dt = config.dt;
    if (dt <= 0.0) {
        dt = solver.estimate_stable_dt(integrator_type);
        if (rank == 0) {
            std::cout << "\n[Time Step] Auto-estimated stable dt = " << dt << "\n";
        }
    }

    // Ensure output directory exists
    if (rank == 0) {
        mkdir(config.output_dir.c_str(), 0777);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    std::string csv_file = config.output_dir + "/history.csv";
    mcch::IO::init_history_csv(csv_file, config.num_components, rank);

    // Initial diagnostics and VTK
    mcch::Diagnostics initial_diag = solver.compute_diagnostics(0.0);
    mcch::IO::append_history_csv(csv_file, initial_diag, rank);
    mcch::IO::write_vtk_gathered(config.output_dir + "/solution_00000.vtk", solver);

    if (rank == 0) {
        std::cout << "\n[Simulation Setup]\n";
        std::cout << "  Components: " << config.num_components << "\n";
        std::cout << "  Free Energy: " << config.free_energy_type << "\n";
        if (config.mobility_type == "degenerate") {
            std::cout << "  Mobility: " << config.mobility_type << " (M0 = " << config.mobility_val << ")\n";
        } else if (config.num_components > 2 || !config.mobility_matrix.empty()) {
            std::cout << "  Mobility: constant matrix (" << config.num_components << "x" << config.num_components << ")\n";
        } else {
            std::cout << "  Mobility: " << config.mobility_type << " (M0 = " << config.mobility_val << ")\n";
        }
        std::cout << "  Integrator: " << config.integrator << ", dt = " << dt 
                  << ", total steps = " << config.total_steps << "\n";
        std::cout << "  Execution Device: " << (dev_type == mcch::DeviceType::GPU ? "GPU (CUDA)" : "CPU (Host)") << "\n";
        if (config.is_semi_implicit() && config.stabilization > 0.0) {
            std::cout << "  Stabilization parameter S: " << config.stabilization << "\n";
        }
        std::cout << "  Output dir: " << config.output_dir << "\n";
        std::cout << "  Initial Energy F0 = " << std::setprecision(6) << initial_diag.total_free_energy << "\n\n";

        std::cout << "------------------------------------------------------------------------------------------------------\n";
        std::cout << std::setw(8) << "Step" 
                  << std::setw(12) << "Time" 
                  << std::setw(16) << "Total Energy" 
                  << std::setw(14) << "dF/dt" 
                  << std::setw(14) << "Max |sum(c)-1|" 
                  << std::setw(12) << "c_min" 
                  << std::setw(12) << "c_max" 
                  << std::setw(14) << "MLUPS\n";
        std::cout << "------------------------------------------------------------------------------------------------------\n";
    }

    auto start_time = std::chrono::high_resolution_clock::now();
    auto loop_start = start_time;

    // Main Time Stepping Loop
    for (int step = 1; step <= config.total_steps; ++step) {
        solver.step(dt, integrator_type, config.stabilization);

        if (step % config.diag_interval == 0 || step == config.total_steps) {
            auto now = std::chrono::high_resolution_clock::now();
            double elapsed_sec = std::chrono::duration<double>(now - loop_start).count();
            loop_start = now;

            double steps_in_interval = config.diag_interval;
            double total_updates = steps_in_interval * grid.global_total_size();
            double mlups = (total_updates / elapsed_sec) / 1.0e6;

            mcch::Diagnostics diag = solver.compute_diagnostics(dt * steps_in_interval);
            mcch::IO::append_history_csv(csv_file, diag, rank);

            if (rank == 0) {
                std::cout << std::setw(8) << step 
                          << std::setw(12) << std::fixed << std::setprecision(4) << solver.current_time 
                          << std::setw(16) << std::scientific << std::setprecision(5) << diag.total_free_energy 
                          << std::setw(14) << diag.dF_dt 
                          << std::setw(14) << diag.max_unity_deviation 
                          << std::setw(12) << std::fixed << std::setprecision(4) << diag.min_composition 
                          << std::setw(12) << diag.max_composition 
                          << std::setw(14) << std::fixed << std::setprecision(2) << mlups << "\n";
            }
        }

        if (step % config.output_interval == 0 || step == config.total_steps) {
            std::stringstream ss;
            ss << config.output_dir << "/solution_" << std::setfill('0') << std::setw(5) << step << ".vtk";
            mcch::IO::write_vtk_gathered(ss.str(), solver);
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    double total_wall_time = std::chrono::duration<double>(end_time - start_time).count();
    double total_mlups = (double)config.total_steps * grid.global_total_size() / total_wall_time / 1.0e6;

    mcch::Diagnostics final_diag = solver.compute_diagnostics(dt);

    if (rank == 0) {
        std::cout << "------------------------------------------------------------------------------------------------------\n";
        std::cout << "\n[Simulation Summary]\n";
        std::cout << "  Completed " << config.total_steps << " steps in " 
                  << std::fixed << std::setprecision(2) << total_wall_time << " seconds.\n";
        std::cout << "  Performance: " << std::setprecision(2) << total_mlups << " MLUPS (Mega Lattice Updates / sec)\n";
        std::cout << "  Initial Free Energy: " << initial_diag.total_free_energy << "\n";
        std::cout << "  Final Free Energy:   " << final_diag.total_free_energy << " (decrease: " 
                  << (initial_diag.total_free_energy - final_diag.total_free_energy) << ")\n";
        std::cout << "  Conservation Check:\n";
        for (int m = 0; m < config.num_components; ++m) {
            double initial_mass = initial_diag.average_composition[m];
            double final_mass = final_diag.average_composition[m];
            double drift = std::abs(final_mass - initial_mass);
            std::cout << "    Component " << m << ": initial=" << std::setprecision(6) << initial_mass 
                      << ", final=" << final_mass << ", drift=" << std::scientific << drift << "\n";
        }
        std::cout << "  Partition of unity error: " << std::scientific << final_diag.max_unity_deviation << "\n";
        std::cout << "  Outputs saved to: " << config.output_dir << "/\n";
        std::cout << "======================================================================\n";
    }

    solver.fftw_solver.reset();
    MPI_Finalize();
    return 0;
}
