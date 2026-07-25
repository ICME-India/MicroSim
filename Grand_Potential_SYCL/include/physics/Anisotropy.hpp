#pragma once
#include <sycl/sycl.hpp>

struct AnisoVars{
    float ac;
    sycl::double3 dadn;
};


inline AnisoVars cubic_anisotropy(sycl::double3 ni, double delta){
    AnisoVars res;

    float nx2 = ni.x() * ni.x();
    float ny2 = ni.y() * ni.y();
    float nz2 = ni.z() * ni.z();

    float nx4 = nx2 * nx2;
    float ny4 = ny2 * ny2;
    float nz4 = nz2 * nz2;

    // Calculate scalat ac
    res.ac = 1.0f - (float)delta * (3.0f - 4.0f * (nx4 + ny4 + nz4));

    // calculate vector v
    res.dadn.x() = 16.0f * (float)delta * ni.x() * nx2;
    res.dadn.y() = 16.0f * (float)delta * ni.y() * ny2;
    res.dadn.z() = 16.0f * (float)delta * ni.z() * nz2;

    return res;
}