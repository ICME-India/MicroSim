#pragma once

#include "MPIContext.hpp"
#include "Solver.hpp"
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace mcch {

struct SimulationConfig {
    int dim = 2;
    int nx = 128;
    int ny = 128;
    int nz = 1;
    double dx = 1.0;
    double dy = 1.0;
    double dz = 1.0;

    std::string bc_x = "periodic";
    std::string bc_y = "periodic";
    std::string bc_z = "periodic";

    int num_components = 3;
    std::string free_energy_type = "polynomial_multiwell"; // "regular_solution", "polynomial_multiwell", "binary_double_well"
    double RT = 1.0;
    double T = 1.0;
    std::vector<std::vector<double>> omega; // interaction matrix
    std::vector<std::vector<double>> W;     // multi-well barriers
    std::vector<double> A;                 // self-well barriers

    std::vector<std::vector<double>> kappa; // gradient energy coefficients
    std::string mobility_type = "constant"; // "constant", "degenerate"
    double mobility_val = 1.0;
    std::vector<std::vector<double>> mobility_matrix; // constant mobility matrix for multicomponent systems

    std::string integrator = "rk2"; // "euler", "rk2", "rk4", "semi_implicit", "semi_implicit_bdf2"
    std::string device = "cpu";     // "cpu", "gpu", "cuda"
    double dt = -1.0; // <= 0 means auto-estimate
    double stabilization = 0.0; // Eyre-type convex-concave stabilization parameter S >= 0
    int total_steps = 1000;
    int output_interval = 100;
    int diag_interval = 10;

    std::string initial_condition = "random"; // "random", "droplet"
    std::vector<double> c_mean;
    double noise_amp = 0.02;
    unsigned int seed = 42;

    // Droplet parameters
    int droplet_comp = 0;
    double droplet_radius = 20.0;
    double droplet_diffuse_width = 2.0;

    std::string output_dir = "results";

    // Parse JSON file
    void load_from_json(const std::string& filepath);
    
    // Set sensible defaults if fields are unpopulated
    void finalize_defaults();

    BoundaryType get_bc_x() const {
        return (bc_x == "neumann" || bc_x == "no_flux") ? BoundaryType::NEUMANN : BoundaryType::PERIODIC;
    }
    BoundaryType get_bc_y() const {
        return (bc_y == "neumann" || bc_y == "no_flux") ? BoundaryType::NEUMANN : BoundaryType::PERIODIC;
    }
    BoundaryType get_bc_z() const {
        return (bc_z == "neumann" || bc_z == "no_flux") ? BoundaryType::NEUMANN : BoundaryType::PERIODIC;
    }
    TimeIntegratorType get_integrator_type() const {
        if (integrator == "euler") return TimeIntegratorType::EULER;
        if (integrator == "rk4") return TimeIntegratorType::RK4;
        if (integrator == "semi_implicit" || integrator == "semi_implicit_euler" || integrator == "imex") {
            return TimeIntegratorType::SEMI_IMPLICIT;
        }
        if (integrator == "semi_implicit_bdf2" || integrator == "sbdf2") {
            return TimeIntegratorType::SEMI_IMPLICIT_BDF2;
        }
        return TimeIntegratorType::RK2;
    }
    bool is_semi_implicit() const {
        TimeIntegratorType t = get_integrator_type();
        return (t == TimeIntegratorType::SEMI_IMPLICIT || t == TimeIntegratorType::SEMI_IMPLICIT_BDF2);
    }
    DeviceType get_device_type() const {
        return (device == "gpu" || device == "cuda") ? DeviceType::GPU : DeviceType::CPU;
    }
};

} // namespace mcch
