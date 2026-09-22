#pragma once

#include "MPIContext.hpp"
#include <cstddef>
#include <vector>
#include <cmath>

namespace mcch {

class Grid {
public:
    int dim; // 2 or 3
    int nx, ny, nz; // Global dimensions
    double dx, dy, dz; // Grid spacings
    int ghost; // Ghost layer depth (default 1)

    // Local slab decomposition info (along x)
    int local_nx;
    int local_x_start;
    int local_x_end;
    int alloc_nx; // local_nx + 2 * ghost

    // Boundary conditions
    BoundaryType bc_x;
    BoundaryType bc_y;
    BoundaryType bc_z;

    Grid(int dimension, int nx_glob, int ny_glob, int nz_glob,
         double spacing_x, double spacing_y, double spacing_z,
         const MPIContext& mpi_ctx,
         BoundaryType boundary_y = BoundaryType::PERIODIC,
         BoundaryType boundary_z = BoundaryType::PERIODIC,
         int ghost_depth = 1)
        : dim(dimension), nx(nx_glob), ny(ny_glob), nz(nz_glob),
          dx(spacing_x), dy(spacing_y), dz(spacing_z), ghost(ghost_depth),
          bc_x(mpi_ctx.bc_x), bc_y(boundary_y), bc_z(boundary_z) {

        local_nx = mpi_ctx.local_nx;
        local_x_start = mpi_ctx.local_x_start;
        local_x_end = mpi_ctx.local_x_end;
        alloc_nx = local_nx + 2 * ghost;

        if (dim == 2) {
            nz = 1;
            dz = 1.0;
        }
    }

    // Strides & sizes
    inline size_t slice_size() const {
        return (size_t)ny * nz;
    }

    inline size_t local_interior_size() const {
        return (size_t)local_nx * ny * nz;
    }

    inline size_t local_total_size() const {
        return (size_t)alloc_nx * ny * nz;
    }

    inline size_t global_total_size() const {
        return (size_t)nx * ny * nz;
    }

    // Index calculation with ghost offset:
    // i_loc is in [0, local_nx - 1], mapping to local allocated array index i_alloc = i_loc + ghost
    inline size_t idx(int i_loc, int j, int k) const {
        return ((size_t)(i_loc + ghost) * ny + j) * nz + k;
    }

    // Direct access to allocated array index (including ghost cells)
    // i_alloc is in [0, alloc_nx - 1]
    inline size_t raw_idx(int i_alloc, int j, int k) const {
        return ((size_t)i_alloc * ny + j) * nz + k;
    }

    // 2D indexing convenience (k = 0)
    inline size_t idx2d(int i_loc, int j) const {
        return (size_t)(i_loc + ghost) * ny + j;
    }

    inline size_t raw_idx2d(int i_alloc, int j) const {
        return (size_t)i_alloc * ny + j;
    }

    // Physical coordinates
    inline double coord_x(int i_loc) const {
        return (local_x_start + i_loc) * dx;
    }

    inline double coord_y(int j) const {
        return j * dy;
    }

    inline double coord_z(int k) const {
        return k * dz;
    }

    // Total domain volume / area
    inline double volume() const {
        if (dim == 2) {
            return (nx * dx) * (ny * dy);
        } else {
            return (nx * dx) * (ny * dy) * (nz * dz);
        }
    }

    inline double cell_volume() const {
        if (dim == 2) {
            return dx * dy;
        } else {
            return dx * dy * dz;
        }
    }

    // Neighbor index helper along y with boundary conditions
    inline void get_y_neighbors(int j, int& j_m, int& j_p) const {
        if (bc_y == BoundaryType::PERIODIC) {
            j_m = (j == 0) ? ny - 1 : j - 1;
            j_p = (j == ny - 1) ? 0 : j + 1;
        } else { // NEUMANN / Reflection
            j_m = (j == 0) ? 1 : j - 1;
            j_p = (j == ny - 1) ? ny - 2 : j + 1;
        }
    }

    // Neighbor index helper along z with boundary conditions
    inline void get_z_neighbors(int k, int& k_m, int& k_p) const {
        if (dim == 2) {
            k_m = 0;
            k_p = 0;
            return;
        }
        if (bc_z == BoundaryType::PERIODIC) {
            k_m = (k == 0) ? nz - 1 : k - 1;
            k_p = (k == nz - 1) ? 0 : k + 1;
        } else { // NEUMANN / Reflection
            k_m = (k == 0) ? 1 : k - 1;
            k_p = (k == nz - 1) ? nz - 2 : k + 1;
        }
    }
};

} // namespace mcch
