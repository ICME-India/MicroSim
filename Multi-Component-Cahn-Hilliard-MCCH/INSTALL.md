# Generalized Installation & Portability Guide

This document provides complete instructions for building, installing, and deploying the **Generalized Multi-Component Cahn-Hilliard (MCCH)** solver on **any machine**—including Linux workstations, macOS laptops (Intel & Apple Silicon), Conda environments, HPC supercomputers, and containers.

---

## 1. Supported Platforms & Environments

| Category | Supported Platforms |
|---|---|
| **Operating Systems** | Linux (Ubuntu/Debian, Fedora, RHEL/CentOS/Rocky/Alma, Arch, openSUSE, Alpine), macOS (11+ Big Sur through Sequoia) |
| **CPU Architectures** | `x86_64` (AMD/Intel), `aarch64` / `arm64` (Apple Silicon M1/M2/M3/M4, AWS Graviton, Nvidia Jetson, ARM HPC) |
| **C++ Compilers** | GCC 9+, Clang 10+, AppleClang, Intel oneAPI (`icpx`/`icx`), NVHPC |
| **MPI Implementations**| OpenMPI 3+, MPICH 3+, Intel MPI, MS-MPI, MVAPICH |
| **Acceleration** | NVIDIA CUDA Toolkit 11.0+ (optional; automatic fallback to CPU-only on systems without GPU) |
| **FFT Library** | FFTW 3.3+ with MPI support (`libfftw3` and `libfftw3_mpi`) |
| **Python** | Python 3.8+ (`numpy`, `scipy`, `matplotlib`, `pandas`, `flask`, `pillow`, `mpi4py`) |

---

## 2. Quick Start: Automated Universal Installer

The repository includes a self-configuring installer script ([`install.sh`](install.sh)) that auto-detects your operating system, CPU architecture, available compilers, package managers, and hardware accelerators:

```bash
# Standard user installation (installs to ~/.local or active Conda env)
./install.sh

# Run verification test suite during installation
./install.sh --test
```

### Common Installation Presets:

```bash
# 1. Non-root user on shared machine or cluster (using Conda):
./install.sh --conda

# 2. Automated system install with dependency resolution (requires sudo):
./install.sh --system-deps --test

# 3. Isolated local Python virtual environment:
./install.sh --python-venv --test

# 4. Systems with PEP 668 externally managed environment (force host install):
./install.sh --break-system-packages

# 5. High-performance host-tuned build (-march=native / -mcpu=native):
./install.sh --native --test

# 6. Custom installation prefix:
./install.sh --prefix=/opt/mcch

# 7. Build strictly CPU-only (skip CUDA even if nvcc is present):
./install.sh --no-cuda
```

---

## 3. Installation via Conda / Mamba (Recommended for Non-Root Users & macOS)

For shared systems, university clusters, or macOS where you don't have administrative (`sudo`) privileges, **Conda-Forge** packages the entire toolchain (compilers, MPI, FFTW with MPI, Python, and libraries) in an isolated, user-space environment:

```bash
# 1. Create the environment from environment.yml
conda env create -f environment.yml

# 2. Activate the environment
conda activate mcch-env

# 3. Build and install directly into the Conda environment
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX
cmake --build build -j$(nproc)
ctest --test-dir build --output-on-failure
cmake --install build
```

---

## 4. Native Linux Package Managers

If you prefer installing dependencies via your native Linux distribution package manager:

### Ubuntu / Debian / Linux Mint / Pop!_OS
```bash
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    cmake \
    pkg-config \
    libopenmpi-dev \
    openmpi-bin \
    libfftw3-dev \
    libfftw3-mpi-dev \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv

# Build and install:
./install.sh --test
```

### Fedora / RHEL / CentOS / Rocky Linux / AlmaLinux
```bash
sudo dnf install -y \
    gcc-c++ \
    make \
    cmake \
    pkgconf-pkg-config \
    openmpi-devel \
    fftw-devel \
    python3-devel \
    python3-pip

# Load OpenMPI module (standard on RHEL/Fedora):
module load mpi/openmpi-x86_64

# Build and install:
./install.sh --test
```

### Arch Linux / Manjaro
```bash
sudo pacman -Sy --noconfirm \
    base-devel \
    cmake \
    pkg-config \
    openmpi \
    fftw \
    python \
    python-pip

# Build and install:
./install.sh --test
```

### openSUSE / SLES
```bash
sudo zypper install -y \
    gcc-c++ \
    make \
    cmake \
    pkg-config \
    openmpi-devel \
    fftw3-devel \
    python3-devel \
    python3-pip

./install.sh --test
```

---

## 5. macOS Installation (Apple Silicon M1/M2/M3/M4 and Intel)

On macOS, install the dependencies using **Homebrew**:

```bash
# 1. Install prerequisites via Homebrew
brew install cmake open-mpi fftw python pkg-config

# 2. Run the automated installer
./install.sh --test
```

> [!NOTE]
> On Apple Silicon (`arm64`), Homebrew packages are located in `/opt/homebrew`. The CMake build system automatically detects `/opt/homebrew/include` and `/opt/homebrew/lib`.

---

## 6. Supercomputers & HPC Clusters (Environment Modules, Spack, SLURM)

HPC supercomputers typically manage software through **Lmod** or **Environment Modules**:

### Using Environment Modules (e.g. OLCF, NERSC, ALCF, Universities)
```bash
# Load cluster modules
module purge
module load gcc/11.2.0    # or intel-oneapi / clang
module load openmpi/4.1.2 # or cray-mpich / intel-mpi
module load fftw3/3.3.10  # or fftw
module load cmake/3.24.0
module load cuda/12.2     # optional for GPU nodes

# Configure with custom prefix (e.g., in your scratch/project directory):
./install.sh --prefix=$HOME/mcch_install --test
```

### Using Spack
```bash
spack env create mcch
spack env activate mcch
spack add cmake gcc openmpi fftw+mpi
spack install
./install.sh --prefix=$SPACK_ENV/view --test
```

### Sample SLURM Batch Job Script (`job.slurm`)
```bash
#!/bin/bash
#SBATCH --job-name=mcch_spinodal
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --time=01:00:00
#SBATCH --partition=compute

# Run 8-rank MPI slab solver
srun mcch_solver -c examples/ternary_spinodal.json
```

---

## 7. Containerized Deployment (Docker & Apptainer / Singularity)

Containers provide 100% reproducible execution across heterogeneous operating systems without installing host dependencies:

### Docker (Workstations & Cloud)
```bash
# Build Docker image (CPU)
docker build -t mcch-solver:latest .

# Build Docker image with CUDA GPU acceleration:
docker build --build-arg BASE_IMAGE=nvidia/cuda:12.2.2-devel-ubuntu22.04 -t mcch-solver:gpu .

# Run the Phase-Field Studio Web Dashboard:
docker run -p 5000:5000 --rm -it mcch-solver:latest

# Run a CLI simulation with MPI:
docker run --rm -it mcch-solver:latest mpirun -np 4 mcch_solver -c examples/ternary_spinodal.json
```

### Apptainer / Singularity (HPC Supercomputers)
HPC facilities prohibit Docker due to root daemon security restrictions. Use **Apptainer**:
```bash
# Build Apptainer image
apptainer build mcch.sif Apptainer.def

# Run across cluster nodes with SLURM and GPU support:
srun apptainer exec --nv mcch.sif mcch_solver -c examples/gpu_ternary_spinodal.json
```

---

## 8. CMake Build Configuration Reference

For fine-grained control, invoke CMake directly:

```bash
mkdir -p build && cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$HOME/.local \
    -DENABLE_CUDA=AUTO \
    -DENABLE_NATIVE_TUNING=OFF \
    -DBUILD_TESTING=ON

make -j$(nproc)
ctest --output-on-failure
make install
```

### Configuration Options:

| Option | Values | Default | Description |
|---|---|---|---|
| `CMAKE_INSTALL_PREFIX` | Path | `~/.local` | Target directory for binaries (`bin/`), libraries (`lib/`), headers (`include/mcch/`), and assets |
| `CMAKE_BUILD_TYPE` | `Release` / `Debug` / `RelWithDebInfo` | `Release` | Optimization and debugging level |
| `ENABLE_CUDA` | `ON` / `OFF` / `AUTO` | `AUTO` | Enables NVIDIA CUDA kernels if `nvcc` is detected; falls back to CPU cleanly if absent |
| `REQUIRE_CUDA` | `ON` / `OFF` | `OFF` | If `ON`, fails configuration if CUDA compiler is not present |
| `ENABLE_NATIVE_TUNING` | `ON` / `OFF` | `OFF` | Enables host-specific CPU microarchitecture tuning (`-march=native` / `-mcpu=native`). Set to `OFF` for generic portable binaries |
| `BUILD_TESTING` | `ON` / `OFF` | `ON` | Compiles the 12-test MPI/thermodynamics verification suite |
| `FFTW3_ROOT` | Path | Auto | Explicit path to custom FFTW3 installation root directory |

---

## 9. Verifying the Installation

### 1. Verify Command-Line Solver
```bash
# Ensure mcch_solver is in your PATH
which mcch_solver

# Display help and supported parameters
mcch_solver --help

# Run a 4-rank MPI simulation
mpirun -np 4 mcch_solver -c examples/ternary_spinodal.json
```

### 2. Verify Interactive Web App Studio
```bash
python3 run_app.py --port 5000
```
Navigate to `http://localhost:5000` in your web browser.

### 3. Run Verification Tests
```bash
cd build
ctest --output-on-failure
```

---

## 10. Troubleshooting & FAQ

### Q1: `FFTW3 or FFTW3_MPI library not found`
- **Debian / Ubuntu**: Install both packages:
  ```bash
  sudo apt-get install libfftw3-dev libfftw3-mpi-dev
  ```
- **Fedora / RHEL**:
  ```bash
  sudo dnf install fftw-devel
  ```
- **macOS**:
  ```bash
  brew install fftw
  ```
- **Non-root**: Specify the FFTW path manually:
  ```bash
  cmake -B build -DFFTW3_ROOT=/path/to/fftw ...
  ```
  or use Conda: `conda install -c conda-forge fftw`.

### Q2: OpenMPI reports `There are not enough slots available in the system`
- When running 4 MPI ranks on machines with fewer physical CPU cores (e.g. 1-2 core virtual machines or CI runners), pass `--oversubscribe`:
  ```bash
  mpirun -np 4 --oversubscribe mcch_solver -c examples/ternary_spinodal.json
  ```
  The CMake test suite automatically includes this flag whenever supported.

### Q3: `mcch_solver: command not found` after installation
- Ensure your installation `bin` directory is in your `PATH`:
  ```bash
  export PATH="$HOME/.local/bin:$PATH"
  ```
  Add this line to your `~/.bashrc` or `~/.zshrc`.

### Q4: `CMake Error: The source ... does not match the source ... used to generate cache`
- **Cause**: `CMakeCache.txt` stores absolute filesystem paths to the source tree, build directory, compilers, and dependencies. If the repository directory is moved, copied, extracted elsewhere, or compiled on a different host, CMake halts on path mismatches.
- **Solution**:
  - Run `./install.sh`: The installer automatically detects path mismatches or missing compilers in stale caches and purges the build directory automatically.
  - Or manually clean the build directory before re-configuring:
    ```bash
    ./install.sh --distclean    # Or: rm -rf build
    ./install.sh                # Fresh build
    ```

### Q5: `error: externally-managed-environment` during Python dependency installation
- **Cause**: Modern Linux distributions (Ubuntu 23.04+, Ubuntu 24.04 LTS, Debian 12 Bookworm, Fedora 38+, Arch Linux) and Homebrew Python 3.12+ enforce **PEP 668**. This marks system Python as "externally managed" so that `pip` cannot accidentally alter system packages managed by `apt`/`dnf`/`pacman`.
- **Solutions**:
  1. **Automatic Virtual Environment** (Recommended):
     The `./install.sh` script automatically detects PEP 668 and creates a local virtual environment in `.venv`. You can also trigger it explicitly:
     ```bash
     ./install.sh --python-venv
     ```
     To use the Phase-Field Studio web app or Python CLI tools later:
     ```bash
     source .venv/bin/activate
     python3 run_app.py --port 5000
     ```
  2. **Conda / Mamba** (Recommended for non-root / HPC users):
     ```bash
     ./install.sh --conda
     ```
  3. **Allow System/User Package Installation**:
     If you prefer installing Python packages into `~/.local` without a virtualenv:
     ```bash
     ./install.sh --break-system-packages
     # Or manual pip:
     pip install --user --break-system-packages -e .
     ```
  4. **Manual Virtual Environment**:
     ```bash
     python3 -m venv .venv
     source .venv/bin/activate
     pip install -e .
     ```
