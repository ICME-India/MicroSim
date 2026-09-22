#include "MPIContext.hpp"
#include "Grid.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int nx = 64, ny = 32, nz = 1;

    // Test 1: Partitioning completeness
    {
        mcch::MPIContext ctx(MPI_COMM_WORLD, nx, ny, nz, mcch::BoundaryType::PERIODIC, 1);
        int total_local = 0;
        MPI_Reduce(&ctx.local_nx, &total_local, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            assert(total_local == nx);
            std::cout << "[PASS] Slab partition coverage: sum(local_nx) == " << total_local << "\n";
        }
    }

    // Test 2: Ghost exchange accuracy (Periodic)
    {
        mcch::MPIContext ctx(MPI_COMM_WORLD, nx, ny, nz, mcch::BoundaryType::PERIODIC, 1);
        mcch::Grid grid(2, nx, ny, nz, 1.0, 1.0, 1.0, ctx, mcch::BoundaryType::PERIODIC, mcch::BoundaryType::PERIODIC, 1);

        std::vector<double> field(grid.local_total_size(), -999.0);

        // Fill interior with unique global indices: 1000 * i_glob + j
        for (int i_loc = 0; i_loc < ctx.local_nx; ++i_loc) {
            int i_glob = ctx.local_x_start + i_loc;
            for (int j = 0; j < ny; ++j) {
                size_t idx = grid.idx(i_loc, j, 0);
                field[idx] = 1000.0 * i_glob + j;
            }
        }

        ctx.exchange_ghosts(field.data());

        // Check left ghost (i_alloc = 0): should match (ctx.local_x_start - 1 + nx) % nx
        int expected_left_glob = (ctx.local_x_start - 1 + nx) % nx;
        for (int j = 0; j < ny; ++j) {
            size_t left_ghost_idx = grid.raw_idx(0, j, 0);
            double expected = 1000.0 * expected_left_glob + j;
            if (std::abs(field[left_ghost_idx] - expected) > 1e-9) {
                std::cerr << "Rank " << rank << " Left ghost mismatch at j=" << j 
                          << ": got " << field[left_ghost_idx] << " expected " << expected << "\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        // Check right ghost (i_alloc = ctx.local_nx + 1): should match (ctx.local_x_end) % nx
        int expected_right_glob = (ctx.local_x_end) % nx;
        for (int j = 0; j < ny; ++j) {
            size_t right_ghost_idx = grid.raw_idx(ctx.local_nx + 1, j, 0);
            double expected = 1000.0 * expected_right_glob + j;
            if (std::abs(field[right_ghost_idx] - expected) > 1e-9) {
                std::cerr << "Rank " << rank << " Right ghost mismatch at j=" << j 
                          << ": got " << field[right_ghost_idx] << " expected " << expected << "\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        if (rank == 0) {
            std::cout << "[PASS] Periodic ghost cell exchange across " << size << " ranks.\n";
        }
    }

    // Test 3: Ghost exchange accuracy (Neumann / No-Flux)
    {
        mcch::MPIContext ctx(MPI_COMM_WORLD, nx, ny, nz, mcch::BoundaryType::NEUMANN, 1);
        mcch::Grid grid(2, nx, ny, nz, 1.0, 1.0, 1.0, ctx, mcch::BoundaryType::NEUMANN, mcch::BoundaryType::NEUMANN, 1);

        std::vector<double> field(grid.local_total_size(), -999.0);
        for (int i_loc = 0; i_loc < ctx.local_nx; ++i_loc) {
            int i_glob = ctx.local_x_start + i_loc;
            for (int j = 0; j < ny; ++j) {
                size_t idx = grid.idx(i_loc, j, 0);
                field[idx] = 1000.0 * i_glob + j;
            }
        }

        ctx.exchange_ghosts(field.data());

        // For rank 0, left ghost mirrors interior i_alloc = 1 (i_loc = 0)
        if (rank == 0) {
            for (int j = 0; j < ny; ++j) {
                size_t ghost_idx = grid.raw_idx(0, j, 0);
                size_t interior_idx = grid.raw_idx(1, j, 0);
                (void)ghost_idx;
                (void)interior_idx;
                assert(field[ghost_idx] == field[interior_idx]);
            }
        }

        // For rank size-1, right ghost mirrors interior i_alloc = local_nx
        if (rank == size - 1) {
            for (int j = 0; j < ny; ++j) {
                size_t ghost_idx = grid.raw_idx(ctx.local_nx + 1, j, 0);
                size_t interior_idx = grid.raw_idx(ctx.local_nx, j, 0);
                (void)ghost_idx;
                (void)interior_idx;
                assert(field[ghost_idx] == field[interior_idx]);
            }
        }

        if (rank == 0) {
            std::cout << "[PASS] Neumann boundary ghost cell exchange across " << size << " ranks.\n";
        }
    }

    MPI_Finalize();
    return 0;
}
