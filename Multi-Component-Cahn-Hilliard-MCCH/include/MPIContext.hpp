#pragma once

#include <mpi.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>
#include <algorithm>

namespace mcch {

enum class BoundaryType {
    PERIODIC,
    NEUMANN // No-flux
};

class MPIContext {
public:
    MPI_Comm comm;
    int rank;
    int size;

    // Global and local domain decomposition (slab along x-axis)
    int nx_global;
    int ny_global;
    int nz_global;

    int local_nx;
    int local_x_start;
    int local_x_end; // exclusive: local_x_start + local_nx

    std::vector<int> all_local_nx;
    std::vector<int> all_local_x_start;

    int left_rank;
    int right_rank;

    BoundaryType bc_x;
    int ghost; // depth of ghost layers (usually 1)

    MPIContext(MPI_Comm communicator, int nx, int ny, int nz, 
               BoundaryType bc_x_type = BoundaryType::PERIODIC, int ghost_depth = 1,
               bool use_fftw_decomp = false)
        : comm(communicator), nx_global(nx), ny_global(ny), nz_global(nz),
          bc_x(bc_x_type), ghost(ghost_depth) {
        
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &size);

        if (size > nx_global) {
            if (rank == 0) {
                std::cerr << "Error: Number of MPI ranks (" << size 
                          << ") cannot exceed global nx (" << nx_global << ") for slab decomposition." << std::endl;
            }
            MPI_Abort(comm, 1);
        }

        if (use_fftw_decomp) {
            fftw_mpi_init();
            ptrdiff_t loc_nx = 0, loc_start = 0;
            if (nz_global > 1) {
                fftw_mpi_local_size_3d(nx_global, ny_global, nz_global, comm, &loc_nx, &loc_start);
            } else {
                fftw_mpi_local_size_2d(nx_global, ny_global, comm, &loc_nx, &loc_start);
            }
            local_nx = (int)loc_nx;
            local_x_start = (int)loc_start;
            local_x_end = local_x_start + local_nx;
        } else {
            // Compute slab distribution along x
            int base_nx = nx_global / size;
            int remainder = nx_global % size;

            local_nx = base_nx + (rank < remainder ? 1 : 0);
            local_x_start = rank * base_nx + std::min(rank, remainder);
            local_x_end = local_x_start + local_nx;
        }

        all_local_nx.resize(size);
        all_local_x_start.resize(size);
        MPI_Allgather(&local_nx, 1, MPI_INT, all_local_nx.data(), 1, MPI_INT, comm);
        MPI_Allgather(&local_x_start, 1, MPI_INT, all_local_x_start.data(), 1, MPI_INT, comm);

        // Determine neighbor ranks along x
        if (bc_x == BoundaryType::PERIODIC) {
            left_rank = (rank - 1 + size) % size;
            right_rank = (rank + 1) % size;
        } else { // NEUMANN / NO-FLUX
            left_rank = (rank == 0) ? MPI_PROC_NULL : rank - 1;
            right_rank = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;
        }
    }

    // Exchange ghost cells for a single 3D/2D scalar field
    // Local array layout: (local_nx + 2*ghost) * ny * nz
    // Interior x indices: [ghost, ghost + local_nx - 1]
    // Left ghost: [0, ghost - 1], Right ghost: [ghost + local_nx, ghost + local_nx + ghost - 1]
    void exchange_ghosts(double* field) const {
        size_t slice_size = (size_t)ny_global * nz_global;
        size_t ghost_data_size = (size_t)ghost * slice_size;

        if (size == 1) {
            // Single rank case
            if (bc_x == BoundaryType::PERIODIC) {
                // Periodic: copy interior right to ghost left, and interior left to ghost right
                std::copy(field + (size_t)ghost * slice_size,
                          field + (size_t)(ghost + ghost) * slice_size,
                          field + (size_t)(ghost + local_nx) * slice_size);
                std::copy(field + (size_t)local_nx * slice_size,
                          field + (size_t)(local_nx + ghost) * slice_size,
                          field);
            } else {
                // Neumann: ghost mirrors adjacent interior slice
                // Left boundary: ghost[0..ghost-1] = interior[ghost..2*ghost-1] (reversed or direct)
                for (int g = 0; g < ghost; ++g) {
                    std::copy(field + (size_t)(ghost + g) * slice_size,
                              field + (size_t)(ghost + g + 1) * slice_size,
                              field + (size_t)(ghost - 1 - g) * slice_size);
                    std::copy(field + (size_t)(ghost + local_nx - 1 - g) * slice_size,
                              field + (size_t)(ghost + local_nx - g) * slice_size,
                              field + (size_t)(ghost + local_nx + g) * slice_size);
                }
            }
            return;
        }

        // Multi-rank MPI exchange
        MPI_Request requests[4];
        int num_requests = 0;

        // Pointers
        double* send_left_buf = field + (size_t)ghost * slice_size;
        double* recv_left_buf = field;
        double* send_right_buf = field + (size_t)local_nx * slice_size;
        double* recv_right_buf = field + (size_t)(ghost + local_nx) * slice_size;

        // Post receives first
        if (left_rank != MPI_PROC_NULL) {
            MPI_Irecv(recv_left_buf, (int)ghost_data_size, MPI_DOUBLE, left_rank, 100, comm, &requests[num_requests++]);
        }
        if (right_rank != MPI_PROC_NULL) {
            MPI_Irecv(recv_right_buf, (int)ghost_data_size, MPI_DOUBLE, right_rank, 200, comm, &requests[num_requests++]);
        }

        // Post sends
        if (right_rank != MPI_PROC_NULL) {
            MPI_Isend(send_right_buf, (int)ghost_data_size, MPI_DOUBLE, right_rank, 100, comm, &requests[num_requests++]);
        }
        if (left_rank != MPI_PROC_NULL) {
            MPI_Isend(send_left_buf, (int)ghost_data_size, MPI_DOUBLE, left_rank, 200, comm, &requests[num_requests++]);
        }

        if (num_requests > 0) {
            MPI_Waitall(num_requests, requests, MPI_STATUSES_IGNORE);
        }

        // Physical boundary condition handling for non-periodic boundaries at domain ends
        if (bc_x == BoundaryType::NEUMANN) {
            if (rank == 0) {
                // Mirror left boundary
                for (int g = 0; g < ghost; ++g) {
                    std::copy(field + (size_t)(ghost + g) * slice_size,
                              field + (size_t)(ghost + g + 1) * slice_size,
                              field + (size_t)(ghost - 1 - g) * slice_size);
                }
            }
            if (rank == size - 1) {
                // Mirror right boundary
                for (int g = 0; g < ghost; ++g) {
                    std::copy(field + (size_t)(ghost + local_nx - 1 - g) * slice_size,
                              field + (size_t)(ghost + local_nx - g) * slice_size,
                              field + (size_t)(ghost + local_nx + g) * slice_size);
                }
            }
        }
    }

    // Exchange ghosts for multiple fields simultaneously
    void exchange_ghosts_multiple(std::vector<double*>& fields) const {
        for (double* f : fields) {
            exchange_ghosts(f);
        }
    }

    // Global reductions
    double allreduce_sum(double val) const {
        double result = 0.0;
        MPI_Allreduce(&val, &result, 1, MPI_DOUBLE, MPI_SUM, comm);
        return result;
    }

    double allreduce_max(double val) const {
        double result = 0.0;
        MPI_Allreduce(&val, &result, 1, MPI_DOUBLE, MPI_MAX, comm);
        return result;
    }

    double allreduce_min(double val) const {
        double result = 0.0;
        MPI_Allreduce(&val, &result, 1, MPI_DOUBLE, MPI_MIN, comm);
        return result;
    }

    // Gather local interior field from all ranks to root rank
    void gather_field(const double* local_field, double* global_field, int root = 0) const {
        size_t slice_size = (size_t)ny_global * nz_global;
        const double* send_buf = local_field + (size_t)ghost * slice_size;
        int send_count = local_nx * (int)slice_size;

        std::vector<int> recv_counts(size);
        std::vector<int> displs(size);

        for (int r = 0; r < size; ++r) {
            recv_counts[r] = all_local_nx[r] * (int)slice_size;
            displs[r] = all_local_x_start[r] * (int)slice_size;
        }

        MPI_Gatherv(send_buf, send_count, MPI_DOUBLE,
                    global_field, recv_counts.data(), displs.data(), MPI_DOUBLE,
                    root, comm);
    }
};

} // namespace mcch
