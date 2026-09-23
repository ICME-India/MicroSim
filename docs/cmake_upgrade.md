# CMake Build System Upgrade

This document describes the migration of the **MicroSim** build system from a
collection of hand-written `Makefile`s to a unified [CMake](https://cmake.org/)
build system.  The original `Makefile`s are retained alongside the new
`CMakeLists.txt` files to ease the transition.

---

## Table of Contents

1. [Motivation](#motivation)
2. [Prerequisites](#prerequisites)
3. [Repository Layout](#repository-layout)
4. [Quick Start](#quick-start)
5. [Building Individual Solvers](#building-individual-solvers)
6. [Configuration Options](#configuration-options)
7. [Solver-Specific Notes](#solver-specific-notes)
   - [Grand_Potential_Serial](#grand_potential_serial)
   - [Grand_Potential_MPI](#grand_potential_mpi)
   - [KKS_CuFFT](#kks_cufft)
   - [KKS_FD_CUDA_MPI](#kks_fd_cuda_mpi)
   - [KKS_OpenCl](#kks_opencl)
   - [Cahn_Hilliard_FFT_2D](#cahn_hilliard_fft_2d)
8. [Out-of-Scope Modules](#out-of-scope-modules)
9. [Troubleshooting](#troubleshooting)

---

## Motivation

The previous build system consisted of independent, hand-crafted `Makefile`s
inside each solver subdirectory.  This caused several pain points:

| Problem | CMake Solution |
|---|---|
| Hard-coded library paths (`/share/apps/gsl/lib`) | `find_package()` locates libraries portably |
| No out-of-source build support | CMake separates source and build trees by default |
| Compiler flags duplicated across solvers | Shared CMake logic and variables |
| No way to build all solvers in one step | Top-level `CMakeLists.txt` with `add_subdirectory()` |
| No install support | `install()` targets for all executables |

---

## Prerequisites

| Dependency | Minimum version | Required by |
|---|---|---|
| CMake | 3.18 | all |
| C compiler (GCC, Clang, NVHPC) | C11 | all |
| C++ compiler | C++17 | KKS_FD_CUDA_MPI |
| CUDA Toolkit | 11.0 | KKS_CuFFT, KKS_FD_CUDA_MPI |
| MPI (OpenMPI / MPICH) | any | Grand_Potential_MPI, KKS_FD_CUDA_MPI, KKS_OpenCl |
| HDF5 (parallel build) | 1.10+ | Grand_Potential_MPI, KKS_FD_CUDA_MPI, KKS_OpenCl |
| GSL | 2.0+ | Grand_Potential_MPI, KKS_FD_CUDA_MPI, KKS_OpenCl |
| OpenCL | 1.2+ | KKS_OpenCl |
| FFTW3 | 3.x | Cahn_Hilliard_FFT_2D |

---

## Repository Layout

```
MicroSim/
├── CMakeLists.txt                  ← NEW top-level CMake entry point
├── cmake/
│   └── FindFFTW3.cmake             ← NEW custom FFTW3 find module
├── docs/
│   └── cmake_upgrade.md            ← this document
│
├── Grand_Potential_Serial/
│   ├── CMakeLists.txt              ← NEW
│   └── Makefile                    ← original (retained)
│
├── Grand_Potential_MPI/
│   ├── CMakeLists.txt              ← NEW
│   └── Makefile                    ← original (retained)
│
├── KKS_CuFFT/
│   ├── CMakeLists.txt              ← NEW
│   └── Makefile                    ← original (retained)
│
├── KKS_FD_CUDA_MPI/
│   ├── CMakeLists.txt              ← NEW
│   └── Makefile                    ← original (retained)
│
├── KKS_OpenCl/
│   ├── CMakeLists.txt              ← NEW
│   └── Makefile                    ← original (retained)
│
├── Cahn_Hilliard_FFT_2D/
│   ├── CMakeLists.txt              ← NEW
│   └── Makefile                    ← original (retained)
│
├── Grand_potential_AMReX/          ← uses AMReX GNUmake (unchanged)
└── Grand_potential_OpenFOAM/       ← uses OpenFOAM wmake (unchanged)
```

---

## Quick Start

### Build all solvers

```bash
# 1. Clone the repository (or enter the existing clone)
cd MicroSim

# 2. Create an out-of-source build directory
cmake -S . -B build

# 3. Compile (replace 8 with the number of CPU cores available)
cmake --build build -j8

# 4. (Optional) Install executables to /usr/local/bin
cmake --install build
```

To change the install prefix:

```bash
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/.local
cmake --build build -j8
cmake --install build
```

---

## Building Individual Solvers

You can configure CMake to build only the solvers you need:

```bash
# Build only the serial Grand-Potential solver
cmake -S . -B build \
      -DBUILD_GRAND_POTENTIAL_SERIAL=ON  \
      -DBUILD_GRAND_POTENTIAL_MPI=OFF    \
      -DBUILD_KKS_CUFFT=OFF             \
      -DBUILD_KKS_FD_CUDA_MPI=OFF       \
      -DBUILD_KKS_OPENCL=OFF            \
      -DBUILD_CAHN_HILLIARD_FFT_2D=OFF
cmake --build build -j8
```

Alternatively, change into a specific subdirectory and invoke CMake
independently:

```bash
cd Grand_Potential_Serial
cmake -S . -B build
cmake --build build -j8
```

---

## Configuration Options

| CMake option | Default | Description |
|---|---|---|
| `BUILD_GRAND_POTENTIAL_SERIAL` | `ON` | Build the serial Grand-Potential solver |
| `BUILD_GRAND_POTENTIAL_MPI` | `ON` | Build the MPI Grand-Potential solver |
| `BUILD_KKS_CUFFT` | `ON` | Build the KKS cuFFT solver |
| `BUILD_KKS_FD_CUDA_MPI` | `ON` | Build the KKS FD CUDA MPI solver |
| `BUILD_KKS_OPENCL` | `ON` | Build the KKS OpenCL solver |
| `BUILD_CAHN_HILLIARD_FFT_2D` | `ON` | Build the Cahn-Hilliard FFT 2D solver |
| `CUDA_ARCH` | `70` | CUDA GPU architecture (e.g. `70`, `80`, `90`) |
| `ENABLE_HDF5` *(KKS_FD_CUDA_MPI)* | `ON` | Enable HDF5 output |
| `ENABLE_CUFFTMP` *(KKS_FD_CUDA_MPI)* | `OFF` | Enable multi-GPU cuFFTMp support |
| `CUFFTMP_INCLUDE_DIR` | *(unset)* | Path to cuFFTMp headers (needed when `ENABLE_CUFFTMP=ON`) |
| `CUFFTMP_LIB_DIR` | *(unset)* | Path to cuFFTMp library (needed when `ENABLE_CUFFTMP=ON`) |
| `FFTW3_ROOT` | *(unset)* | Hint path for FFTW3 installation prefix |

Options are passed to CMake using the `-D` flag, for example:

```bash
cmake -S . -B build -DCUDA_ARCH=80 -DBUILD_KKS_FD_CUDA_MPI=ON -DENABLE_CUFFTMP=OFF
```

---

## Solver-Specific Notes

### Grand_Potential_Serial

- **Languages**: C
- **Dependencies**: libm (standard library)
- **Executables produced**: `microsim_gp`
- The solver is a single-translation-unit program; all `functions/` and
  `solverloop/` files are header-only includes.

### Grand_Potential_MPI

- **Languages**: C
- **Dependencies**: MPI, HDF5 (parallel), GSL
- **Executables produced**: `microsim_gp`, `write_xdmf` (if source present)
- On systems where HDF5 was built with MPI support (the `h5pcc` wrapper), CMake
  will find the library automatically via `find_package(HDF5 REQUIRED)`.
- If your GSL is installed in a non-standard location set `GSL_ROOT_DIR`:
  ```bash
  cmake -S . -B build -DGSL_ROOT_DIR=/path/to/gsl
  ```

### KKS_CuFFT

- **Languages**: CUDA
- **Dependencies**: CUDA Toolkit (cuFFT)
- **Executables produced**: `microsim_kks_cufft`
- A `DATA/` directory is created in the build folder at configure time to match
  the solver's expected output path.
- Set the target GPU architecture with `-DCUDA_ARCH=<sm>` (default: `70`).

### KKS_FD_CUDA_MPI

- **Languages**: C, C++, CUDA
- **Dependencies**: MPI, HDF5, GSL, CUDA Toolkit (cuFFT or cuFFTMp)
- **Executables produced**: `microsim_kks_fd_cuda_mpi`, `write_xdmf_kks`,
  `reconstruct_kks`
- CUDA **separable compilation** (`-dc` in the original Makefile) is handled
  automatically by setting `CUDA_SEPARABLE_COMPILATION ON`.
- To enable multi-GPU cuFFTMp:
  ```bash
  cmake -S . -B build \
        -DENABLE_CUFFTMP=ON \
        -DCUFFTMP_INCLUDE_DIR=/opt/nvidia/hpc_sdk/.../cufftmp \
        -DCUFFTMP_LIB_DIR=/opt/nvidia/hpc_sdk/.../math_libs/lib64
  ```

### KKS_OpenCl

- **Languages**: C
- **Dependencies**: MPI, OpenCL, HDF5, GSL
- **Executables produced**: `microsim_kks_opencl`, `reconstruct_opencl`
- OpenCL is found via the standard CMake `FindOpenCL` module.  If CMake cannot
  locate it automatically, set `OpenCL_ROOT`:
  ```bash
  cmake -S . -B build -DOpenCL_ROOT=/usr/local/cuda
  ```

### Cahn_Hilliard_FFT_2D

- **Languages**: C
- **Dependencies**: FFTW3
- **Executables produced**: `microsim_ch_fft`
- FFTW3 is found via the bundled `cmake/FindFFTW3.cmake` module.  If FFTW3 is
  in a non-standard location set `FFTW3_ROOT`:
  ```bash
  cmake -S . -B build -DFFTW3_ROOT=/path/to/fftw3
  ```

---

## Out-of-Scope Modules

The following modules use build systems tightly coupled to their respective
frameworks and are **not** migrated to CMake:

| Module | Build system | Reason |
|---|---|---|
| `Grand_potential_AMReX/` | AMReX GNUmake | AMReX provides its own Make infrastructure; a CMake port requires [AMReX cmake support](https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#cmake) and is a larger effort |
| `Grand_potential_OpenFOAM/` | OpenFOAM `wmake` | OpenFOAM solvers must be built with `wmake`; CMake is not officially supported |
| `Bridgman/` | OpenFOAM `wmake` | Same reason as above |

---

## Troubleshooting

### CMake cannot find MPI

```
CMake Error: Could not find MPI
```

Ensure that the MPI compilers are on `PATH` (`mpicc`, `mpicxx`) or set
`MPI_ROOT`:

```bash
cmake -S . -B build -DMPI_ROOT=/path/to/mpi
```

### CMake cannot find HDF5

On systems that provide a parallel HDF5 wrapper (e.g., `h5pcc`), tell CMake to
prefer the parallel variant:

```bash
HDF5_PREFER_PARALLEL=ON cmake -S . -B build
```

Or set the hint manually:

```bash
cmake -S . -B build -DHDF5_ROOT=/path/to/hdf5
```

### CMake cannot find CUDA

Ensure `nvcc` is on `PATH` or set `CUDAToolkit_ROOT`:

```bash
cmake -S . -B build -DCUDAToolkit_ROOT=/usr/local/cuda
```

### CUDA architecture mismatch

If you see PTX/SASS compatibility warnings or errors at runtime, update
`CUDA_ARCH` to match your GPU:

```bash
cmake -S . -B build -DCUDA_ARCH=80   # A100, RTX 30xx
cmake -S . -B build -DCUDA_ARCH=90   # H100
```

### CMake cannot find FFTW3

```
CMake Error: Could not find FFTW3
```

Set `FFTW3_ROOT` to the FFTW3 installation prefix:

```bash
cmake -S . -B build -DFFTW3_ROOT=/usr/local
```

Or load the appropriate module on HPC systems:

```bash
module load fftw3
cmake -S . -B build
```
