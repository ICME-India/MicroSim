#pragma once

#include "Solver.hpp"
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <sys/stat.h>

namespace mcch {

class IO {
public:
    // Write legacy VTK Structured Points file (gathered on rank 0)
    static void write_vtk_gathered(const std::string& filename, const Solver& solver) {
        const_cast<Solver&>(solver).sync_from_device();
        std::vector<double> global_buffer;
        if (solver.mpi.rank == 0) {
            global_buffer.resize(solver.grid.global_total_size());
        }

        // Rank 0 opens file and writes header
        std::ofstream vtk_file;
        if (solver.mpi.rank == 0) {
            vtk_file.open(filename);
            if (!vtk_file.is_open()) {
                std::cerr << "Error: Could not open " << filename << " for writing VTK." << std::endl;
                return;
            }

            vtk_file << "# vtk DataFile Version 3.0\n";
            vtk_file << "MCCH Simulation Step " << solver.current_step << " Time " << solver.current_time << "\n";
            vtk_file << "ASCII\n";
            vtk_file << "DATASET STRUCTURED_POINTS\n";
            vtk_file << "DIMENSIONS " << solver.grid.nx << " " << solver.grid.ny << " " << solver.grid.nz << "\n";
            vtk_file << "ORIGIN 0.0 0.0 0.0\n";
            vtk_file << "SPACING " << solver.grid.dx << " " << solver.grid.dy << " " << solver.grid.dz << "\n";
            vtk_file << "POINT_DATA " << solver.grid.global_total_size() << "\n";
        }

        // Write each component field
        for (int m = 0; m < solver.n_comp; ++m) {
            solver.mpi.gather_field(solver.c[m].data(), global_buffer.data(), 0);

            if (solver.mpi.rank == 0) {
                vtk_file << "SCALARS c" << m << " double 1\n";
                vtk_file << "LOOKUP_TABLE default\n";

                // VTK point ordering is x fastest or z fastest?
                // In standard VTK STRUCTURED_POINTS:
                // point index = x + y * nx + z * nx * ny
                // But our C++ array index is: ((i_glob * ny) + j) * nz + k
                // Let's write in standard VTK order (iterate z, then y, then x):
                for (int k = 0; k < solver.grid.nz; ++k) {
                    for (int j = 0; j < solver.grid.ny; ++j) {
                        for (int i = 0; i < solver.grid.nx; ++i) {
                            size_t cpp_idx = ((size_t)i * solver.grid.ny + j) * solver.grid.nz + k;
                            vtk_file << std::setprecision(6) << global_buffer[cpp_idx] << " ";
                        }
                        vtk_file << "\n";
                    }
                }
            }
        }

        // Also write chemical potentials mu_m
        for (int m = 0; m < solver.n_comp; ++m) {
            solver.mpi.gather_field(solver.mu[m].data(), global_buffer.data(), 0);

            if (solver.mpi.rank == 0) {
                vtk_file << "SCALARS mu" << m << " double 1\n";
                vtk_file << "LOOKUP_TABLE default\n";
                for (int k = 0; k < solver.grid.nz; ++k) {
                    for (int j = 0; j < solver.grid.ny; ++j) {
                        for (int i = 0; i < solver.grid.nx; ++i) {
                            size_t cpp_idx = ((size_t)i * solver.grid.ny + j) * solver.grid.nz + k;
                            vtk_file << std::setprecision(6) << global_buffer[cpp_idx] << " ";
                        }
                        vtk_file << "\n";
                    }
                }
            }
        }

        if (solver.mpi.rank == 0) {
            vtk_file.close();
        }
    }

    // Initialize CSV history file with header
    static void init_history_csv(const std::string& filename, int n_comp, int rank) {
        if (rank == 0) {
            std::ofstream f(filename);
            if (f.is_open()) {
                f << "step,time,total_free_energy,dF_dt,min_c,max_c,max_unity_dev";
                for (int m = 0; m < n_comp; ++m) {
                    f << ",avg_c" << m;
                }
                f << "\n";
                f.close();
            }
        }
    }

    // Append diagnostics to CSV history file
    static void append_history_csv(const std::string& filename, const Diagnostics& diag, int rank) {
        if (rank == 0) {
            std::ofstream f(filename, std::ios::app);
            if (f.is_open()) {
                f << diag.step << ","
                  << std::scientific << std::setprecision(8) << diag.time << ","
                  << diag.total_free_energy << ","
                  << diag.dF_dt << ","
                  << diag.min_composition << ","
                  << diag.max_composition << ","
                  << diag.max_unity_deviation;
                for (size_t m = 0; m < diag.average_composition.size(); ++m) {
                    f << "," << diag.average_composition[m];
                }
                f << "\n";
                f.close();
            }
        }
    }

    // Binary checkpoint output per rank
    static void write_checkpoint(const std::string& base_path, int step, const Solver& solver) {
        std::stringstream ss;
        ss << base_path << "_step" << step << "_rank" << solver.mpi.rank << ".bin";
        std::string filename = ss.str();

        std::ofstream f(filename, std::ios::binary);
        if (!f.is_open()) return;

        // Write header
        int n_c = solver.n_comp;
        int loc_nx = solver.grid.local_nx;
        int ny = solver.grid.ny;
        int nz = solver.grid.nz;
        double t = solver.current_time;

        f.write(reinterpret_cast<const char*>(&step), sizeof(int));
        f.write(reinterpret_cast<const char*>(&t), sizeof(double));
        f.write(reinterpret_cast<const char*>(&n_c), sizeof(int));
        f.write(reinterpret_cast<const char*>(&loc_nx), sizeof(int));
        f.write(reinterpret_cast<const char*>(&ny), sizeof(int));
        f.write(reinterpret_cast<const char*>(&nz), sizeof(int));

        // Write interior slab data
        for (int m = 0; m < n_c; ++m) {
            for (int i_loc = 0; i_loc < loc_nx; ++i_loc) {
                size_t src_offset = solver.grid.idx(i_loc, 0, 0);
                size_t count = solver.grid.slice_size();
                f.write(reinterpret_cast<const char*>(&solver.c[m][src_offset]), count * sizeof(double));
            }
        }
        f.close();
    }
};

} // namespace mcch
