#include "error_checks.hpp"

#if ENABLE_CUFFTMP == 1

void gpu_checkAssert(cudaError_t code, const char *file, int line, bool abort)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr,"CUDA_CHECK: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

void cufft_check(int code, const char *file, int line, bool abort)
{
    if (code != CUFFT_SUCCESS)
    {
        fprintf(stderr,"CUFFT_CHECK: %d %s %d\n", code, file, line);
        if (abort) exit(code);
    }
}

// template<typename T>
// double compute_error(const T& ref, const T& test, cufftBox3d box) {
//     auto[begin, end] = BoxIterators(box, ref.data());
//     double max_diff = 0;
//     double max_norm = 0;
//     for(auto it = begin; it != end; it++) {
//         max_norm = std::max<double>(max_norm, std::abs(ref[it.i()]));
//         max_diff = std::max<double>(max_diff, std::abs(ref[it.i()] - test[it.i()]));
//     }
//
//     return max_diff / max_norm;
// }

int assess_error(double error) {
    if(error > 1e-6) {
        printf("FAILED with error %e\n", error);
        return 1;
    } else {
        printf("PASSED with error %e\n", error);
        return 0;
    }
}

#endif
