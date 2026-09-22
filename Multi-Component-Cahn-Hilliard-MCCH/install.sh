#!/usr/bin/env bash
# =============================================================================
# Generalized Universal Installer for MCCH Solver & Phase-Field Studio
# Portable across Linux (Ubuntu, Debian, Fedora, RHEL, Arch, openSUSE), macOS,
# Conda environments, and HPC Supercomputers.
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Color output helpers
# -----------------------------------------------------------------------------
if [[ -t 1 ]]; then
    BOLD="$(printf '\033[1m')"
    RED="$(printf '\033[31m')"
    GREEN="$(printf '\033[32m')"
    YELLOW="$(printf '\033[33m')"
    BLUE="$(printf '\033[34m')"
    CYAN="$(printf '\033[36m')"
    RESET="$(printf '\033[0m')"
else
    BOLD=""
    RED=""
    GREEN=""
    YELLOW=""
    BLUE=""
    CYAN=""
    RESET=""
fi

log_info()    { echo -e "${BLUE}${BOLD}[INFO]${RESET} $*"; }
log_success() { echo -e "${GREEN}${BOLD}[SUCCESS]${RESET} $*"; }
log_warn()    { echo -e "${YELLOW}${BOLD}[WARNING]${RESET} $*"; }
log_error()   { echo -e "${RED}${BOLD}[ERROR]${RESET} $*" >&2; }
log_step()    { echo -e "\n${CYAN}${BOLD}==>${RESET} ${BOLD}$*${RESET}"; }

# -----------------------------------------------------------------------------
# Default Options
# -----------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSPACE_ROOT="${SCRIPT_DIR}"
BUILD_DIR="${WORKSPACE_ROOT}/build"
BUILD_TYPE="Release"
INSTALL_PREFIX=""
ENABLE_CUDA="AUTO"
ENABLE_NATIVE="OFF"
INSTALL_SYSTEM_DEPS=false
USE_CONDA=false
USE_VENV=false
BREAK_SYSTEM_PACKAGES=false
RUN_TESTS=false
CLEAN_BUILD=false
DISTCLEAN_ONLY=false
PARALLEL_JOBS=""

# Check if the host Python environment is externally managed (PEP 668)
is_pep668_managed() {
    # If already inside an active virtualenv or conda env, it is not externally managed
    if [[ -n "${VIRTUAL_ENV:-}" || -n "${CONDA_PREFIX:-}" ]]; then
        return 1
    fi
    if command -v python3 &>/dev/null; then
        if python3 -c "
import sys, sysconfig, os
if hasattr(sys, 'real_prefix') or (hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix):
    sys.exit(1)
for p in [sysconfig.get_path('stdlib'), sysconfig.get_path('platstdlib')]:
    if p and os.path.isfile(os.path.join(p, 'EXTERNALLY-MANAGED')):
        sys.exit(0)
sys.exit(1)
" 2>/dev/null; then
            return 0
        fi
    fi
    return 1
}

# Detect default parallel jobs
if command -v nproc &>/dev/null; then
    PARALLEL_JOBS="$(nproc)"
elif command -v sysctl &>/dev/null; then
    PARALLEL_JOBS="$(sysctl -n hw.ncpu 2>/dev/null || echo 4)"
else
    PARALLEL_JOBS=4
fi

# -----------------------------------------------------------------------------
# Usage & Help
# -----------------------------------------------------------------------------
usage() {
    cat << EOF
${BOLD}Usage:${RESET} ./install.sh [OPTIONS]

${BOLD}Options:${RESET}
  --prefix=<DIR>          Installation prefix directory.
                          Defaults:
                            • In active Conda environment: \$CONDA_PREFIX
                            • As root: /usr/local
                            • As standard user: \$HOME/.local
  --system-deps           Automatically install missing dependencies using the
                          system package manager (apt, dnf, pacman, brew).
                          Requires sudo access on Linux.
  --conda                 Create or use an isolated Conda environment ('mcch-env')
                          using conda-forge (does NOT require sudo or root!).
  --python-venv           Create and use a Python virtual environment in .venv.
  --break-system-packages Allow pip to install to system/user Python even if PEP 668
                          (externally-managed-environment) is active.
  --cuda                  Force enable CUDA GPU acceleration (fails if nvcc missing).
  --no-cuda               Disable CUDA acceleration (build CPU-only solver).
  --native                Enable CPU native architecture optimization (-march=native).
  --test                  Run the verification test suite (ctest) after building.
  --clean                 Remove build directory before compiling.
  --distclean             Purge all build artifacts, CMake cache, and temporary files
                          to restore clean repository portability, then exit.
  -j, --jobs <N>          Number of parallel compile jobs (default: ${PARALLEL_JOBS}).
  -h, --help              Show this help message and exit.

${BOLD}Examples:${RESET}
  # Quick installation for current user:
  ./install.sh

  # Clean non-portable build caches and reset workspace:
  ./install.sh --distclean

  # Non-root installation using Conda:
  ./install.sh --conda

  # Automated system install with dependency installation & testing:
  ./install.sh --system-deps --test

  # High-performance tuned build for current machine:
  ./install.sh --native --test
EOF
    exit 0
}

# -----------------------------------------------------------------------------
# Parse CLI arguments
# -----------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --prefix=*)
            INSTALL_PREFIX="${1#*=}"
            shift
            ;;
        --prefix)
            INSTALL_PREFIX="$2"
            shift 2
            ;;
        --system-deps)
            INSTALL_SYSTEM_DEPS=true
            shift
            ;;
        --conda)
            USE_CONDA=true
            shift
            ;;
        --python-venv)
            USE_VENV=true
            shift
            ;;
        --break-system-packages)
            BREAK_SYSTEM_PACKAGES=true
            shift
            ;;
        --cuda)
            ENABLE_CUDA="ON"
            shift
            ;;
        --no-cuda)
            ENABLE_CUDA="OFF"
            shift
            ;;
        --native)
            ENABLE_NATIVE="ON"
            shift
            ;;
        --test)
            RUN_TESTS=true
            shift
            ;;
        --clean)
            CLEAN_BUILD=true
            shift
            ;;
        --distclean)
            DISTCLEAN_ONLY=true
            shift
            ;;
        -j|--jobs)
            PARALLEL_JOBS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            log_error "Unknown option: $1"
            usage
            ;;
    esac
done

# -----------------------------------------------------------------------------
# Standalone Distclean Action
# -----------------------------------------------------------------------------
if [[ "${DISTCLEAN_ONLY}" == true ]]; then
    log_step "Purging Build Artifacts & CMake Cache for Maximum Portability"
    if [[ -d "${BUILD_DIR}" ]]; then
        log_info "Removing build directory: ${BUILD_DIR}"
        rm -rf "${BUILD_DIR}"
    fi
    log_info "Cleaning Python and packaging artifacts..."
    rm -rf "${WORKSPACE_ROOT}"/*.egg-info "${WORKSPACE_ROOT}/.venv"
    find "${WORKSPACE_ROOT}" -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
    find "${WORKSPACE_ROOT}" -type f -name "*.pyc" -delete 2>/dev/null || true
    log_success "Solver repository is now completely clean and portable!"
    exit 0
fi

# -----------------------------------------------------------------------------
# OS and Architecture Detection
# -----------------------------------------------------------------------------
log_step "Detecting System Environment"

OS_NAME="$(uname -s)"
ARCH_NAME="$(uname -m)"
DISTRO=""

if [[ "${OS_NAME}" == "Linux" ]]; then
    if [[ -f /etc/os-release ]]; then
        # shellcheck disable=SC1091
        source /etc/os-release
        DISTRO="${ID:-linux}"
    elif [[ -f /etc/redhat-release ]]; then
        DISTRO="rhel"
    elif [[ -f /etc/debian_version ]]; then
        DISTRO="debian"
    else
        DISTRO="generic-linux"
    fi
elif [[ "${OS_NAME}" == "Darwin" ]]; then
    DISTRO="macos"
fi

log_info "Operating System: ${OS_NAME} (${DISTRO})"
log_info "Architecture:     ${ARCH_NAME}"

# Determine default prefix if not supplied
if [[ -z "${INSTALL_PREFIX}" ]]; then
    if [[ -n "${CONDA_PREFIX:-}" ]]; then
        INSTALL_PREFIX="${CONDA_PREFIX}"
        log_info "Detected active Conda environment: ${INSTALL_PREFIX}"
    elif [[ $EUID -eq 0 ]]; then
        INSTALL_PREFIX="/usr/local"
    else
        INSTALL_PREFIX="${HOME}/.local"
    fi
fi
log_info "Install Prefix:   ${INSTALL_PREFIX}"

# -----------------------------------------------------------------------------
# Conda Environment Setup (if requested)
# -----------------------------------------------------------------------------
if [[ "${USE_CONDA}" == true ]]; then
    log_step "Setting up Conda Environment ('mcch-env')"
    CONDA_CMD=""
    if command -v mamba &>/dev/null; then
        CONDA_CMD="mamba"
    elif command -v conda &>/dev/null; then
        CONDA_CMD="conda"
    elif command -v micromamba &>/dev/null; then
        CONDA_CMD="micromamba"
    else
        log_error "Conda/Mamba executable not found in PATH."
        log_error "Please install Miniforge or Miniconda first, or run without --conda."
        exit 1
    fi

    log_info "Using package manager: ${CONDA_CMD}"
    if ${CONDA_CMD} env list | grep -q "^mcch-env "; then
        log_info "Updating existing conda environment 'mcch-env'..."
        ${CONDA_CMD} env update -n mcch-env -f "${WORKSPACE_ROOT}/environment.yml" --prune
    else
        log_info "Creating new conda environment 'mcch-env'..."
        ${CONDA_CMD} env create -f "${WORKSPACE_ROOT}/environment.yml"
    fi

    log_warn "To activate this environment in your shell, run:"
    log_warn "  conda activate mcch-env"
    
    # Try activating for the current script
    # shellcheck disable=SC1091
    CONDA_BASE="$(${CONDA_CMD} info --base 2>/dev/null || echo "")"
    if [[ -n "${CONDA_BASE}" && -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
        source "${CONDA_BASE}/etc/profile.d/conda.sh"
        conda activate mcch-env
        INSTALL_PREFIX="${CONDA_PREFIX}"
        log_info "Active Conda prefix set to: ${INSTALL_PREFIX}"
    fi
fi

# -----------------------------------------------------------------------------
# Dependency Detection & Automated System Package Installation
# -----------------------------------------------------------------------------
log_step "Checking System Dependencies"

MISSING_DEPS=()

# Check C++ compiler
CXX_FOUND=false
for cxx in "${CXX:-}" g++ clang++ icpx c++; do
    if [[ -n "${cxx}" ]] && command -v "${cxx}" &>/dev/null; then
        CXX_FOUND=true
        CXX_BIN="$(command -v "${cxx}")"
        break
    fi
done
if [[ "${CXX_FOUND}" == false ]]; then
    MISSING_DEPS+=("c++-compiler")
else
    log_info "C++ Compiler:     ${CXX_BIN}"
fi

# Check CMake
if command -v cmake &>/dev/null; then
    CMAKE_VERSION="$(cmake --version | head -n1 | awk '{print $3}')"
    log_info "CMake:            $(command -v cmake) (v${CMAKE_VERSION})"
else
    MISSING_DEPS+=("cmake")
fi

# Check MPI
if command -v mpicxx &>/dev/null || command -v mpic++ &>/dev/null; then
    MPI_BIN="$(command -v mpicxx 2>/dev/null || command -v mpic++)"
    log_info "MPI C++ Wrapper:  ${MPI_BIN}"
else
    MISSING_DEPS+=("mpi")
fi

# Check FFTW3
FFTW3_DETECTED=false
if pkg-config --exists fftw3 2>/dev/null || [[ -f /usr/include/fftw3.h ]] || [[ -f "${INSTALL_PREFIX}/include/fftw3.h" ]] || [[ -f /opt/homebrew/include/fftw3.h ]]; then
    FFTW3_DETECTED=true
    log_info "FFTW3:            Detected in system/prefix include path"
else
    # Check common library locations
    for p in /usr/lib /usr/lib64 /usr/lib/*-linux-gnu /usr/local/lib /opt/homebrew/lib; do
        if compgen -G "${p}/libfftw3*.so*" >/dev/null 2>&1 || compgen -G "${p}/libfftw3*.dylib" >/dev/null 2>&1; then
            FFTW3_DETECTED=true
            log_info "FFTW3:            Detected in ${p}"
            break
        fi
    done
fi
if [[ "${FFTW3_DETECTED}" == false ]]; then
    MISSING_DEPS+=("fftw3")
fi

# Check Python3 & Pip
if command -v python3 &>/dev/null; then
    PY_VERSION="$(python3 --version 2>&1 | awk '{print $2}')"
    log_info "Python:           $(command -v python3) (v${PY_VERSION})"
else
    MISSING_DEPS+=("python3")
fi

# Check CUDA
if command -v nvcc &>/dev/null; then
    CUDA_VERSION="$(nvcc --version | grep "release" | awk '{print $5}' | tr -d ',')"
    log_info "CUDA Compiler:    $(command -v nvcc) (v${CUDA_VERSION})"
    HAS_CUDA=true
else
    log_info "CUDA Compiler:    Not found in PATH"
    HAS_CUDA=false
fi

# Handle Missing Dependencies
if [[ ${#MISSING_DEPS[@]} -gt 0 ]]; then
    log_warn "Missing required dependencies: ${MISSING_DEPS[*]}"

    if [[ "${INSTALL_SYSTEM_DEPS}" == true ]]; then
        log_info "Attempting automated system package installation..."
        case "${DISTRO}" in
            ubuntu|debian|linuxmint|pop)
                sudo apt-get update -y
                sudo apt-get install -y build-essential cmake libopenmpi-dev \
                    libfftw3-dev libfftw3-mpi-dev python3-dev python3-pip python3-venv pkg-config
                ;;
            fedora|rhel|centos|rocky|almalinux)
                sudo dnf install -y gcc-c++ make cmake openmpi-devel \
                    fftw-devel python3-devel python3-pip pkgconf-pkg-config
                ;;
            arch|manjaro)
                sudo pacman -Sy --noconfirm base-devel cmake openmpi fftw python python-pip
                ;;
            opensuse*|sles)
                sudo zypper install -y gcc-c++ make cmake openmpi-devel \
                    fftw3-devel python3-devel python3-pip pkg-config
                ;;
            macos)
                if command -v brew &>/dev/null; then
                    brew install cmake open-mpi fftw python pkg-config
                else
                    log_error "Homebrew not found. Please install Homebrew or dependencies manually."
                    exit 1
                fi
                ;;
            *)
                log_error "Automatic package installation not supported for '${DISTRO}'."
                log_error "Please install the missing tools manually, or use './install.sh --conda'."
                exit 1
                ;;
        esac
    else
        echo ""
        log_warn "To install dependencies automatically with your package manager, run:"
        log_warn "  ./install.sh --system-deps"
        echo ""
        log_info "Alternatively, to install EVERYTHING without root/sudo access using Conda:"
        log_info "  ./install.sh --conda"
        echo ""
        log_error "Cannot proceed without required dependencies. Aborting."
        exit 1
    fi
fi

# -----------------------------------------------------------------------------
# Python Environment Setup
# -----------------------------------------------------------------------------
AUTO_VENV_ACTIVATED=false

if [[ "${USE_VENV}" == true ]] || ( [[ "${USE_CONDA}" == false ]] && is_pep668_managed && [[ "${BREAK_SYSTEM_PACKAGES}" == false ]] ); then
    log_step "Configuring Python Virtual Environment (.venv)"
    if [[ -d "${WORKSPACE_ROOT}/.venv" && -f "${WORKSPACE_ROOT}/.venv/bin/activate" ]]; then
        # shellcheck disable=SC1091
        source "${WORKSPACE_ROOT}/.venv/bin/activate"
        log_info "Reusing existing virtual environment: ${VIRTUAL_ENV}"
        AUTO_VENV_ACTIVATED=true
    elif command -v python3 &>/dev/null && python3 -m venv "${WORKSPACE_ROOT}/.venv" 2>/dev/null; then
        # shellcheck disable=SC1091
        source "${WORKSPACE_ROOT}/.venv/bin/activate"
        log_info "Created and activated isolated virtual environment: ${VIRTUAL_ENV}"
        AUTO_VENV_ACTIVATED=true
    else
        if [[ "${USE_VENV}" == true ]]; then
            log_error "Failed to create virtual environment (.venv)."
            log_error "On Debian/Ubuntu, please install python3-venv: 'sudo apt install python3-venv'"
            exit 1
        else
            log_warn "Detected externally managed environment (PEP 668), but python3 -m venv failed."
            log_info "Falling back to system environment with --break-system-packages flag."
            BREAK_SYSTEM_PACKAGES=true
        fi
    fi
fi

# -----------------------------------------------------------------------------
# CMake Build and Install
# -----------------------------------------------------------------------------
log_step "Configuring and Building C++ MCCH Solver"

if [[ "${CLEAN_BUILD}" == true && -d "${BUILD_DIR}" ]]; then
    log_info "Cleaning previous build directory: ${BUILD_DIR}"
    rm -rf "${BUILD_DIR}"
fi

# Detect non-portable or stale CMakeCache.txt (e.g. copied from another machine or path)
if [[ -f "${BUILD_DIR}/CMakeCache.txt" ]]; then
    CACHED_SOURCE_DIR="$(grep "^CMAKE_HOME_DIRECTORY:INTERNAL=" "${BUILD_DIR}/CMakeCache.txt" 2>/dev/null | cut -d= -f2- || true)"
    CACHED_BUILD_DIR="$(grep "^CMAKE_CACHEFILE_DIR:INTERNAL=" "${BUILD_DIR}/CMakeCache.txt" 2>/dev/null | cut -d= -f2- || true)"
    CACHED_CXX="$(grep "^CMAKE_CXX_COMPILER:FILEPATH=" "${BUILD_DIR}/CMakeCache.txt" 2>/dev/null | cut -d= -f2- || true)"

    STALE_REASON=""
    if [[ -n "${CACHED_SOURCE_DIR}" && "${CACHED_SOURCE_DIR}" != "${WORKSPACE_ROOT}" ]]; then
        STALE_REASON="source directory moved from '${CACHED_SOURCE_DIR}' to '${WORKSPACE_ROOT}'"
    elif [[ -n "${CACHED_BUILD_DIR}" && "${CACHED_BUILD_DIR}" != "${BUILD_DIR}" ]]; then
        STALE_REASON="build directory moved from '${CACHED_BUILD_DIR}' to '${BUILD_DIR}'"
    elif [[ -n "${CACHED_CXX}" && ! -x "${CACHED_CXX}" ]]; then
        STALE_REASON="cached C++ compiler ('${CACHED_CXX}') is missing or not executable on this machine"
    fi

    if [[ -n "${STALE_REASON}" ]]; then
        log_warn "Non-portable or stale CMakeCache.txt detected: ${STALE_REASON}."
        log_info "Automatically resetting build directory '${BUILD_DIR}' for clean portable build..."
        rm -rf "${BUILD_DIR}"
    fi
fi

mkdir -p "${BUILD_DIR}"

CMAKE_ARGS=(
    "-DCMAKE_BUILD_TYPE=${BUILD_TYPE}"
    "-DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}"
    "-DENABLE_NATIVE_TUNING=${ENABLE_NATIVE}"
)

if [[ "${ENABLE_CUDA}" == "ON" ]]; then
    CMAKE_ARGS+=("-DENABLE_CUDA=ON" "-DREQUIRE_CUDA=ON")
elif [[ "${ENABLE_CUDA}" == "OFF" ]]; then
    CMAKE_ARGS+=("-DENABLE_CUDA=OFF")
else
    CMAKE_ARGS+=("-DENABLE_CUDA=ON")
fi

log_info "Running CMake configuration:"
log_info "  cmake -B \"${BUILD_DIR}\" ${CMAKE_ARGS[*]}"

cmake -B "${BUILD_DIR}" "${CMAKE_ARGS[@]}" "${WORKSPACE_ROOT}"

log_step "Compiling Solver (-j${PARALLEL_JOBS})"
cmake --build "${BUILD_DIR}" -j"${PARALLEL_JOBS}"

# -----------------------------------------------------------------------------
# Testing
# -----------------------------------------------------------------------------
if [[ "${RUN_TESTS}" == true ]]; then
    log_step "Running CTest Verification Suite"
    (
        cd "${BUILD_DIR}"
        ctest --output-on-failure
    )
    log_success "All test suites passed successfully!"
fi

# -----------------------------------------------------------------------------
# Target Installation
# -----------------------------------------------------------------------------
log_step "Installing to: ${INSTALL_PREFIX}"

mkdir -p "${INSTALL_PREFIX}"

if [[ -w "${INSTALL_PREFIX}" ]]; then
    cmake --install "${BUILD_DIR}"
else
    log_info "Elevated permissions required for ${INSTALL_PREFIX}. Running with sudo..."
    sudo cmake --install "${BUILD_DIR}"
fi

# -----------------------------------------------------------------------------
# Install Python Package & CLI tools
# -----------------------------------------------------------------------------
log_step "Installing Python Dependencies & CLI Tools"

if command -v pip3 &>/dev/null || command -v pip &>/dev/null; then
    PIP_CMD="$(command -v pip3 2>/dev/null || command -v pip)"
    log_info "Using pip: ${PIP_CMD}"

    PIP_FLAGS=("--no-warn-script-location")

    # If externally managed or --break-system-packages requested, add flag if pip supports it
    if [[ "${BREAK_SYSTEM_PACKAGES}" == true ]] || is_pep668_managed; then
        if ${PIP_CMD} install --help 2>&1 | grep -q -- "--break-system-packages"; then
            PIP_FLAGS+=("--break-system-packages")
            log_info "Applied --break-system-packages to satisfy PEP 668"
        fi
        # If installing as non-root outside a venv/conda, default to --user
        if [[ -z "${VIRTUAL_ENV:-}" && -z "${CONDA_PREFIX:-}" && "${EUID}" -ne 0 ]]; then
            PIP_FLAGS+=("--user")
        fi
    fi

    log_info "Installing MCCH Python components..."
    PY_INSTALL_SUCCESS=false
    if ${PIP_CMD} install "${PIP_FLAGS[@]}" -e "${WORKSPACE_ROOT}"; then
        PY_INSTALL_SUCCESS=true
        log_success "MCCH Python package and CLI tools installed in editable mode."
    else
        log_warn "Editable install failed, falling back to requirements.txt..."
        if ${PIP_CMD} install "${PIP_FLAGS[@]}" -r "${WORKSPACE_ROOT}/requirements.txt"; then
            PY_INSTALL_SUCCESS=true
            log_success "MCCH Python dependencies installed from requirements.txt."
        fi
    fi

    if [[ "${PY_INSTALL_SUCCESS}" == false ]]; then
        log_warn "Python package installation could not be completed automatically."
        if is_pep668_managed; then
            log_warn "Detected externally managed Python environment (PEP 668)."
            echo -e "  To resolve and install Python packages, choose one of:"
            echo -e "    1. Use a virtualenv:  ${BOLD}./install.sh --python-venv${RESET}"
            echo -e "    2. Use Conda:         ${BOLD}./install.sh --conda${RESET}"
            echo -e "    3. Force host pip:    ${BOLD}./install.sh --break-system-packages${RESET}"
        else
            echo -e "  Install dependencies manually with: ${BOLD}${PIP_CMD} install -r requirements.txt${RESET}"
        fi
    fi
else
    log_warn "pip not found. Please install Python packages manually: pip install -r requirements.txt"
fi

# -----------------------------------------------------------------------------
# Post-Install Verification & Summary
# -----------------------------------------------------------------------------
SOLVER_BIN="${INSTALL_PREFIX}/bin/mcch_solver"
if [[ ! -f "${SOLVER_BIN}" ]]; then
    SOLVER_BIN="${BUILD_DIR}/mcch_solver"
fi

echo ""
echo -e "${GREEN}${BOLD}======================================================================${RESET}"
echo -e "${GREEN}${BOLD}  MCCH Solver Installation Completed Successfully!                     ${RESET}"
echo -e "${GREEN}${BOLD}======================================================================${RESET}"
echo -e "  • Solver Binary:       ${BOLD}${SOLVER_BIN}${RESET}"
echo -e "  • Install Prefix:      ${INSTALL_PREFIX}"
echo -e "  • Architecture:        ${ARCH_NAME} (${DISTRO})"
echo -e "  • CUDA Acceleration:   $([[ "${HAS_CUDA}" == true && "${ENABLE_CUDA}" != "OFF" ]] && echo "Enabled" || echo "Disabled (CPU)")"
echo -e "  • Native Tuning:       ${ENABLE_NATIVE}"
if [[ -n "${VIRTUAL_ENV:-}" ]]; then
    echo -e "  • Python Virtualenv:   ${BOLD}${VIRTUAL_ENV}${RESET}"
elif [[ -n "${CONDA_PREFIX:-}" ]]; then
    echo -e "  • Python Conda Prefix: ${BOLD}${CONDA_PREFIX}${RESET}"
fi

# Check PATH
if [[ ":$PATH:" != *":${INSTALL_PREFIX}/bin:"* ]]; then
    echo ""
    echo -e "${YELLOW}${BOLD}[NOTICE]${RESET} '${INSTALL_PREFIX}/bin' is not in your current PATH."
    echo -e "Add it to your shell configuration (e.g. ~/.bashrc or ~/.zshrc):"
    echo -e "  ${BOLD}export PATH=\"${INSTALL_PREFIX}/bin:\$PATH\"${RESET}"
fi

echo ""
echo -e "${BOLD}Next Steps:${RESET}"
echo -e "  1. Run a sample MPI simulation:"
echo -e "     ${BOLD}mpirun -np 4 ${SOLVER_BIN} -c examples/ternary_spinodal.json${RESET}"
echo ""
echo -e "  2. Launch the interactive Web App Studio:"
if [[ -n "${VIRTUAL_ENV:-}" ]]; then
    echo -e "     ${BOLD}source ${VIRTUAL_ENV}/bin/activate${RESET}"
    echo -e "     ${BOLD}python3 run_app.py --port 5000${RESET}"
else
    echo -e "     ${BOLD}python3 run_app.py --port 5000${RESET}"
fi
echo -e "     Then navigate to ${CYAN}http://localhost:5000${RESET} in your browser."
echo ""
echo -e "  3. Re-run tests at any time:"
echo -e "     ${BOLD}cd ${BUILD_DIR} && ctest --output-on-failure${RESET}"
echo -e "${GREEN}${BOLD}======================================================================${RESET}\n"
