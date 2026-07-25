# Readme
### portable exascale model for Phase Field simulation based on Grand Potential formulation 
---

This is phase field simulation model capable exascale computation

## Features
- Based Grand Potential Formation extended for rapid solidifcation
- Using SYCL - a vendor indepented parallel paradime


## Libraries
The current solver is extended with following libraries

- **[Intel OneAPI HPC toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html)** - uses OneAPI icpx and mpiicpx compilers
- **[Codeplay OneAPI](https://developer.codeplay.com/products/oneapi/nvidia/2023.0.0/guides/get-started-guide-nvidia)** - to use nvidia GPUs in the backent
- **[OneMATH](https://github.com/uxlfoundation/oneMath)** - for BLAS and LAPACK calculation in a device indepent routine
- **[GSL scientific library C++](https://www.gnu.org/software/gsl/)** - for matrix inversion and spline interpolation for hessian and composition information
- **[HDF5 - High-performance data manegment and storage suite](https://www.hdfgroup.org/solutions/hdf5/)** - multiple paralle processes share and operate on the save HDF5 file concurrently


## Usefull URLS
- **[Device Initated communication using RMA - CUDA Backent](https://www.intel.com/content/www/us/en/docs/mpi-library/developer-guide-linux/2021-16/device-initiated-communications.html)** - for device initiated communication using cuda backent