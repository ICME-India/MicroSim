#!/usr/bin/env python3
from setuptools import setup, find_packages

setup(
    name="mcch-solver",
    version="1.0.0",
    description="Generalized Multi-Component Cahn-Hilliard Phase-Field Solver & Studio",
    author="MCCH Team",
    packages=find_packages(include=["webapp*", "python*"]),
    py_modules=["run_app"],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0",
        "matplotlib>=3.3.0",
        "pandas>=1.2.0",
        "pillow>=8.0.0",
        "flask>=2.0.0",
    ],
    extras_require={
        "mpi": ["mpi4py>=3.0.0"],
        "gpu": ["cupy>=9.0.0"],
    },
    entry_points={
        "console_scripts": [
            "mcch-studio=run_app:main",
            "mcch-plot=python.plot_diagnostics:main",
            "mcch-visualize=python.visualize_vtk:main",
            "mcch-solver-py=python.mcch_solver:main",
        ],
    },
)
