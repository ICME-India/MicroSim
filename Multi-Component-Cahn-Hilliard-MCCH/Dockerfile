# =============================================================================
# Multi-Stage Dockerfile for MCCH Solver & Phase-Field Studio
# Supports CPU-only (default) or CUDA GPU builds (via --build-arg BASE_IMAGE=...)
# =============================================================================

ARG BASE_IMAGE=ubuntu:22.04
FROM ${BASE_IMAGE} AS builder

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC
ENV PIP_BREAK_SYSTEM_PACKAGES=1

# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    pkg-config \
    libopenmpi-dev \
    openmpi-bin \
    libfftw3-dev \
    libfftw3-mpi-dev \
    python3 \
    python3-pip \
    python3-dev \
    python3-venv \
    git \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . /src

# Build and install solver to /opt/mcch
RUN ./install.sh --prefix=/opt/mcch --no-cuda --clean

# -----------------------------------------------------------------------------
# Runtime Image
# -----------------------------------------------------------------------------
FROM ${BASE_IMAGE} AS runtime

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC
ENV PIP_BREAK_SYSTEM_PACKAGES=1

RUN apt-get update && apt-get install -y --no-install-recommends \
    openmpi-bin \
    libopenmpi3 \
    libfftw3-double3 \
    libfftw3-mpi3 \
    python3 \
    python3-pip \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Copy installed solver, libraries, headers and assets from builder
COPY --from=builder /opt/mcch /opt/mcch
COPY --from=builder /src /workspace

WORKDIR /workspace

# Install Python requirements
RUN pip3 install --no-cache-dir -e /workspace

# Create non-root user
RUN useradd -m -u 1000 -s /bin/bash mcchuser && \
    chown -R mcchuser:mcchuser /workspace /opt/mcch

USER mcchuser

ENV PATH="/opt/mcch/bin:${PATH}"
ENV LD_LIBRARY_PATH="/opt/mcch/lib:${LD_LIBRARY_PATH:-}"
ENV PYTHONPATH="/workspace:${PYTHONPATH:-}"

# Expose Web App Studio Port
EXPOSE 5000

# Default entrypoint: Phase-Field Studio Web App
ENTRYPOINT ["python3", "run_app.py", "--host", "0.0.0.0", "--port", "5000"]
