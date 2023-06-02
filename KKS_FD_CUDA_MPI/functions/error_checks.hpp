#ifndef ERROR_CHECKS_HPP_
#define ERROR_CHECKS_HPP_

#ifndef ENABLE_CUFFTMP
#define ENABLE_CUFFTMP 0
#endif

#if ENABLE_CUFFTMP == 1

#include <cufftMp.h>
#include <cstdio>
//#include "box_iterator.hpp"

#define CUDA_CHECK(ans) { gpu_checkAssert((ans), __FILE__, __LINE__); }
void gpu_checkAssert(cudaError_t code, const char *file, int line, bool abort=true);

#define CUFFT_CHECK(ans) { cufft_check((ans), __FILE__, __LINE__); }
void cufft_check(int code, const char *file, int line, bool abort=true);

// template<typename T>
// double compute_error(const T& ref, const T& test, cufftBox3d box);

int assess_error(double error);

#endif //ENABLE_CUFFTMP
#endif //ifndef ERROR_CHECKS_HPP_
