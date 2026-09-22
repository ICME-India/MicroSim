#include "MPIContext.hpp"
#include "Grid.hpp"
#include "FreeEnergy.hpp"
#include "Mobility.hpp"
#include "Solver.hpp"

#include <iostream>
#include <cassert>
#include <cmath>
#include <stdexcept>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // ========================================================
    // Part 1: Unit tests for Redlich_Kister interface & math
    // ========================================================
    if (rank == 0) {
        std::cout << "[Test 1] Instantiating Redlich_Kister...\n";
    }

    // Verify invalid component count throws exception
    bool caught_invalid_comp = false;
    try {
        mcch::Redlich_Kister rk_bad(3, 800.0);
    } catch (const std::invalid_argument& e) {
        caught_invalid_comp = true;
    }
    assert(caught_invalid_comp);
    (void)caught_invalid_comp;

    // Create valid 2-component Redlich-Kister model at T = 800 K
    mcch::Redlich_Kister fe_rk(2, 800.0);
    assert(fe_rk.num_components() == 2);
    assert(std::abs(fe_rk.temperature() - 800.0) < 1e-9);

    fe_rk.set_temperature(750.0);
    assert(std::abs(fe_rk.temperature() - 750.0) < 1e-9);
    fe_rk.set_temperature(800.0);

    // Verify density and derivative evaluations
    if (rank == 0) {
        std::cout << "[Test 2] Validating thermodynamic derivatives against finite difference...\n";
    }

    double h = 1e-6;
    for (double c0 = 0.1; c0 <= 0.9; c0 += 0.1) {
        double c1 = 1.0 - c0;
        std::vector<double> c = {c0, c1};
        std::vector<double> c_plus = {c0 + h, c1 - h};
        std::vector<double> c_minus = {c0 - h, c1 + h};

        double G = fe_rk.density(c);
        double G_plus = fe_rk.density(c_plus);
        double G_minus = fe_rk.density(c_minus);
        assert(!std::isnan(G));
        (void)G;

        std::vector<double> df_dc;
        fe_rk.chemical_derivatives(c, df_dc);
        assert(df_dc.size() == 2);
        assert(!std::isnan(df_dc[0]) && !std::isnan(df_dc[1]));

        double num_deriv = (G_plus - G_minus) / (2.0 * h);
        double ana_deriv = df_dc[0] - df_dc[1];
        double rel_err = std::abs((num_deriv - ana_deriv) / ana_deriv);
        assert(rel_err < 1e-5);
        (void)rel_err;
    }

    // ========================================================
    // Part 2: Multicomponent Cahn-Hilliard Solver integration
    // ========================================================
    if (rank == 0) {
        std::cout << "[Test 3] Running Cahn-Hilliard solver with Redlich_Kister...\n";
    }

    int nx = 32, ny = 32, nz = 1;
    mcch::MPIContext ctx(MPI_COMM_WORLD, nx, ny, nz, mcch::BoundaryType::PERIODIC, 1);
    mcch::Grid grid(2, nx, ny, nz, 1.0, 1.0, 1.0, ctx, mcch::BoundaryType::PERIODIC, mcch::BoundaryType::PERIODIC, 1);

    int n_comp = 2;
    auto fe_ptr = std::make_shared<mcch::Redlich_Kister>(n_comp, 800.0);
    // Use physical mobility scaling
    double mob_val = 1e-4;
    auto mob_ptr = std::make_shared<mcch::ConstantMobility>(n_comp, mob_val);

    std::vector<std::vector<double>> kappa = {
        {10.0, 0.0},
        {0.0, 10.0}
    };

    mcch::Solver solver(ctx, grid, n_comp, fe_ptr, mob_ptr, kappa);
    solver.initialize_random({0.5, 0.5}, 0.02, 42);

    mcch::Diagnostics d_initial = solver.compute_diagnostics(0.0);
    double initial_energy = d_initial.total_free_energy;
    double current_energy = initial_energy;

    double dt = 0.001;
    int energy_increases = 0;

    for (int step = 1; step <= 50; ++step) {
        solver.step(dt, mcch::TimeIntegratorType::RK2);
        double new_energy = solver.compute_total_free_energy();
        if (new_energy > current_energy + 1e-8) {
            energy_increases++;
        }
        current_energy = new_energy;
    }

    mcch::Diagnostics d_final = solver.compute_diagnostics(dt * 50);

    // Verify mass conservation
    for (int m = 0; m < n_comp; ++m) {
        double drift = std::abs(d_final.average_composition[m] - d_initial.average_composition[m]);
        if (rank == 0) {
            std::cout << "Component " << m << " mass drift: " << drift << "\n";
        }
        assert(drift < 1e-12);
    }

    // Verify partition of unity
    assert(d_final.max_unity_deviation < 1e-12);

    // Verify thermodynamic consistency (free energy dissipation)
    if (rank == 0) {
        std::cout << "Initial Free Energy: " << initial_energy << "\n";
        std::cout << "Final Free Energy:   " << current_energy << "\n";
        std::cout << "Energy Reduction:    " << (initial_energy - current_energy) << "\n";
        std::cout << "Energy increases:    " << energy_increases << "\n";
    }
    assert(energy_increases == 0);
    assert(current_energy <= initial_energy);

    if (rank == 0) {
        std::cout << "[PASS] Redlich_Kister thermodynamic function verified successfully!\n";
    }

    MPI_Finalize();
    return 0;
}
