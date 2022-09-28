#ifndef CUDA_ATOMICS_HPP
#define CUDA_ATOMICS_HPP

#include "cuda_keywords.hpp"

namespace gmx_gpu_fmm{

DEVICE __forceinline__
void
__atomicAdd(float* address, float val)
{
    atomicAdd(address,val);
}

DEVICE __forceinline__
void
__atomicAdd(double* address, double val)
{

    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val +__longlong_as_double(assumed)));
    } while (assumed != old);

    //atomicAdd(address,val);
}

}//namespace end

#endif // CUDA_ATOMICS_HPP
