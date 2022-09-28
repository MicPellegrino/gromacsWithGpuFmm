#ifndef CUDA_DIPOLE_COMPENSATION_HPP
#define CUDA_DIPOLE_COMPENSATION_HPP

#include "cuda_P2L.hpp"
#include "cuda_P2M.hpp"

namespace gmx_gpu_fmm{

template <typename Real4, typename CoeffMatrix>
__global__
void
__P2L_dipole_corr(Real4* fake_particles, CoeffMatrix** mu, size_t p, size_t fake_parts_size)
{
    const size_t i = threadIdx.x;

    if(i < fake_parts_size)
    {
        __P2L_cubereduce(fake_particles[i], *mu[0], p);
    }
}

template <typename Real4, typename CoeffMatrix>
__global__
void
__P2M_dipole_corr(Real4* q0abc, CoeffMatrix** omega, size_t p)
{
    const size_t i = threadIdx.x;

    if(i < 4)
    {
        __P2M_cubreduce(q0abc[i].q, q0abc[i], *omega[0], p);
    }
}

template <typename Real4, typename CoeffMatrix>
__global__
void
__P2M_P2L_dipole_corr(Real4* q0abc, CoeffMatrix** omega, Real4* fake_particles, CoeffMatrix** mu, size_t p, size_t fake_parts_size)
{
    const size_t i = threadIdx.x;

    if(i < fake_parts_size)
    {
        __P2L_cubereduce(fake_particles[i], *mu[0], p);
    }

    if(i < 4)
    {
        __P2M_cubreduce(q0abc[i].q, q0abc[i], *omega[0], p);
    }
}

template <typename Real4>
__global__
void
__P2M_P2L_dipole_corr_flush(Real4* q0abc, Real4* fake_particles, size_t fake_parts_size)
{
    const size_t i = threadIdx.x;

    if(i < fake_parts_size)
    {
        fake_particles[i] = Real4(0.,0.,0.,0.);
    }

    if(i < 4)
    {
        q0abc[i] = Real4(0.,0.,0.,0.);
    }
}

}//namespace end

#endif // CUDA_DIPOLE_COMPENSATION_HPP
