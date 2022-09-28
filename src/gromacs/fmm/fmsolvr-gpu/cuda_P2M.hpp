#ifndef _BK_FMM_cuda_P2M_hpp
#define _BK_FMM_cuda_P2M_hpp

#include "cuda_keywords.hpp"
#include "cub_h.h"

namespace gmx_gpu_fmm{

namespace p2m{

template <typename Real>
DEVICE
inline Real __reciprocal(Real a)
{
    return (Real) 1. / a;
}

}

//reference cuda solution
template <typename Real3, typename CoefficientMatrix>
__device__
void __P2M(
        const typename Real3::value_type & q_,
        const Real3 & xyz,
        CoefficientMatrix & omega,
        const size_t p)
{
    typedef typename Real3::value_type real_type;
    typedef typename CoefficientMatrix::value_type complex_type;
    //typedef double real_type;
    //typedef fmm::complex<real_type> complex_type;

    typedef typename Vec_traits<real_type>::scalar_type scalar_type;
    Vec_traits<real_type> VT;

    const real_type x = xyz.x;
    const real_type y = xyz.y;
    const real_type z = xyz.z;
    const real_type q = q_;

    real_type dist_squared = x * x + y * y + z * z;
    real_type twice_z = z + z;

    complex_type complex_x_y(x, -y);
    complex_type complex_x_y_m(VT.same(1));     // iterative computation of ((x - iy)^m)

    // omega_0_0  (for the first iteration)
    real_type omega_mplus1_mplus1 = q;
    real_type e_m = z + twice_z;                // iterative computation of ((2*m + 3)*z)

    // omega_0_0 upto omega_p_p-1
    for (size_t m = 0; m < p; ++m) {
        // omega_m_m  (from previous iteration)
        real_type omega_m_m = omega_mplus1_mplus1;
        omega(m, m) %= complex_x_y_m * omega_m_m;

        // omega_m+1_m+1  (for the next iteration)
        omega_mplus1_mplus1 = VT.same(p2m::__reciprocal(scalar_type(2 * (m + 1)))) * omega_m_m;

        // omega_m+1_m
        real_type omega_mplus1_m = z * omega_m_m;
        omega(m + 1, m) %= complex_x_y_m * omega_mplus1_m;

        // omega_m+2_m upto omega_p_m
        real_type omega_lminus2_m = omega_m_m;
        real_type omega_lminus1_m = omega_mplus1_m;
        real_type f_l = e_m;                    // iterative computation of ((2*l - 1)*z)
        for (size_t l = m + 2; l <= p; ++l) {
            // omega_l_m
            real_type omega_l_m = VT.same(p2m::__reciprocal(scalar_type(l * l - m * m))) *
                (f_l * omega_lminus1_m - dist_squared * omega_lminus2_m);
            omega(l, m) %= complex_x_y_m * omega_l_m;

            // (for the next iteration)
            omega_lminus2_m = omega_lminus1_m;
            omega_lminus1_m = omega_l_m;
            f_l += twice_z;
        }

        // (for the next iteration)
        complex_x_y_m *= complex_x_y;
        e_m += twice_z;
    }

    // omega_p_p
    omega(p, p) %= complex_x_y_m * omega_mplus1_mplus1;
}

//uses warp reduction
template <typename Real3, typename CoefficientMatrix>
__device__
void __P2M_cub_block_reduce(
        const typename Real3::value_type & q_,
        const Real3 & xyz,
        CoefficientMatrix & omega,
        const size_t p)
{
    typedef typename Real3::value_type real_type;
    typedef typename CoefficientMatrix::value_type complex_type;

    //typedef double real_type;
    //typedef fmm::complex<real_type> complex_type;
    typedef typename Vec_traits<real_type>::scalar_type scalar_type;
    Vec_traits<real_type> VT;

    typedef cub::BlockReduce<complex_type, 224> BlockReduce;
    __shared__ typename BlockReduce::TempStorage temp_storage;

    const real_type x = xyz.x;
    const real_type y = xyz.y;
    const real_type z = xyz.z;
    const real_type q = q_;

    real_type dist_squared = x * x + y * y + z * z;
    real_type twice_z = z + z;

    complex_type complex_x_y(x, -y);
    complex_type complex_x_y_m(VT.same(1));     // iterative computation of ((x - iy)^m)

    // omega_0_0  (for the first iteration)
    real_type omega_mplus1_mplus1 = q;
    real_type e_m = z + twice_z;
    // iterative computation of ((2*m + 3)*z)

    //size_t lane_id = threadIdx.x%32;

    // omega_0_0 upto omega_p_p-1
    for (size_t m = 0; m < p; ++m)
    {
        // omega_m_m  (from previous iteration)
        real_type omega_m_m = omega_mplus1_mplus1;

        complex_type out_omega_m_m = complex_x_y_m * omega_m_m;
        complex_type sum_out_omega_m_m = BlockReduce(temp_storage).Sum(out_omega_m_m);
        if(threadIdx.x == 0)
            omega(m, m) %=sum_out_omega_m_m;

        // omega_m+1_m+1  (for the next iteration)
        omega_mplus1_mplus1 = VT.same(p2m::__reciprocal(scalar_type(2 * (m + 1)))) * omega_m_m;

        // omega_m+1_m
        real_type omega_mplus1_m = z * omega_m_m;
        complex_type out_omega_mplus1_m = complex_x_y_m * omega_mplus1_m;
        complex_type sum_out_omega_mplus1_m = BlockReduce(temp_storage).Sum(out_omega_mplus1_m);
        if(threadIdx.x == 0)
            omega(m + 1, m) %= sum_out_omega_mplus1_m;

        // omega_m+2_m upto omega_p_m
        real_type omega_lminus2_m = omega_m_m;
        real_type omega_lminus1_m = omega_mplus1_m;
        real_type f_l = e_m;                    // iterative computation of ((2*l - 1)*z)
        for (size_t l = m + 2; l <= p; ++l) {
            // omega_l_m
            real_type omega_l_m = VT.same(p2m::__reciprocal(scalar_type(l * l - m * m))) *
                (f_l * omega_lminus1_m - dist_squared * omega_lminus2_m);

            complex_type out_omega_l_m = complex_x_y_m * omega_l_m;
            complex_type sum_out_omega_l_m = BlockReduce(temp_storage).Sum(out_omega_l_m);
            if(threadIdx.x == 0)
                omega(l, m) %= sum_out_omega_l_m;

            // (for the next iteration)
            omega_lminus2_m = omega_lminus1_m;
            omega_lminus1_m = omega_l_m;
            f_l += twice_z;
        }

        // (for the next iteration)
        complex_x_y_m *= complex_x_y;
        e_m += twice_z;
    }

    // omega_p_p
    complex_type out_omega_p_p = complex_x_y_m * omega_mplus1_mplus1;
    complex_type sum_out_omega_p_p = BlockReduce(temp_storage).Sum(out_omega_p_p);
    if(threadIdx.x == 0)
        omega(p, p) %=sum_out_omega_p_p;
}


//uses warp reduction
template <typename Real3, typename CoefficientMatrix>
__device__
void __P2M_cubreduce(
        const typename Real3::value_type & q_,
        const Real3 & xyz,
        CoefficientMatrix & omega,
        const size_t p)
{
    typedef typename Real3::value_type real_type;
    typedef typename CoefficientMatrix::value_type complex_type;

    //typedef double real_type;
    //typedef fmm::complex<real_type> complex_type;
    typedef typename Vec_traits<real_type>::scalar_type scalar_type;
    Vec_traits<real_type> VT;

    typedef cub::WarpReduce<complex_type> WarpReducer;
    __shared__ typename WarpReducer::TempStorage temp[128];

    const real_type x = xyz.x;
    const real_type y = xyz.y;
    const real_type z = xyz.z;
    const real_type q = q_;

    real_type dist_squared = x * x + y * y + z * z;
    real_type twice_z = z + z;

    complex_type complex_x_y(x, -y);
    complex_type complex_x_y_m(VT.same(1));     // iterative computation of ((x - iy)^m)

    // omega_0_0  (for the first iteration)
    real_type omega_mplus1_mplus1 = q;
    real_type e_m = z + twice_z;
    // iterative computation of ((2*m + 3)*z)

    size_t lane_id = threadIdx.x%32;

    // omega_0_0 upto omega_p_p-1
    for (size_t m = 0; m < p; ++m)
    {
        // omega_m_m  (from previous iteration)
        real_type omega_m_m = omega_mplus1_mplus1;

        complex_type out_omega_m_m = complex_x_y_m * omega_m_m;
        complex_type sum_out_omega_m_m = WarpReducer(temp[lane_id]).Sum(out_omega_m_m);
        if(lane_id == 0)
            omega(m, m) %=sum_out_omega_m_m;

        // omega_m+1_m+1  (for the next iteration)
        omega_mplus1_mplus1 = VT.same(p2m::__reciprocal(scalar_type(2 * (m + 1)))) * omega_m_m;

        // omega_m+1_m
        real_type omega_mplus1_m = z * omega_m_m;
        complex_type out_omega_mplus1_m = complex_x_y_m * omega_mplus1_m;
        complex_type sum_out_omega_mplus1_m = WarpReducer(temp[lane_id]).Sum(out_omega_mplus1_m);
        if(lane_id == 0)
            omega(m + 1, m) %= sum_out_omega_mplus1_m;

        // omega_m+2_m upto omega_p_m
        real_type omega_lminus2_m = omega_m_m;
        real_type omega_lminus1_m = omega_mplus1_m;
        real_type f_l = e_m;                    // iterative computation of ((2*l - 1)*z)
        for (size_t l = m + 2; l <= p; ++l) {
            // omega_l_m
            real_type omega_l_m = VT.same(p2m::__reciprocal(scalar_type(l * l - m * m))) *
                (f_l * omega_lminus1_m - dist_squared * omega_lminus2_m);

            complex_type out_omega_l_m = complex_x_y_m * omega_l_m;
            complex_type sum_out_omega_l_m = WarpReducer(temp[lane_id]).Sum(out_omega_l_m);
            if(lane_id == 0)
                omega(l, m) %= sum_out_omega_l_m;

            // (for the next iteration)
            omega_lminus2_m = omega_lminus1_m;
            omega_lminus1_m = omega_l_m;
            f_l += twice_z;
        }

        // (for the next iteration)
        complex_x_y_m *= complex_x_y;
        e_m += twice_z;
    }

    // omega_p_p
    complex_type out_omega_p_p = complex_x_y_m * omega_mplus1_mplus1;
    complex_type sum_out_omega_p_p = WarpReducer(temp[lane_id]).Sum(out_omega_p_p);
    if(lane_id == 0)
        omega(p, p) %=sum_out_omega_p_p;
}

//reference cuda solution
template <typename Real3, typename CoefficientMatrix>
__device__
void __P2M_reduce(
        const typename Real3::value_type & q_,
        const Real3 & xyz,
        CoefficientMatrix & omega,
        const size_t p)
{
    typedef typename Real3::value_type real_type;
    typedef typename CoefficientMatrix::value_type complex_type;

    //typedef double real_type;
    //typedef fmm::complex<real_type> complex_type;
    typedef typename Vec_traits<real_type>::scalar_type scalar_type;
    Vec_traits<real_type> VT;

    const real_type x = xyz.x;
    const real_type y = xyz.y;
    const real_type z = xyz.z;
    const real_type q = q_;

    real_type dist_squared = x * x + y * y + z * z;
    real_type twice_z = z + z;

    complex_type complex_x_y(x, -y);
    complex_type complex_x_y_m(VT.same(1));     // iterative computation of ((x - iy)^m)

    // omega_0_0  (for the first iteration)
    real_type omega_mplus1_mplus1 = q;
    real_type e_m = z + twice_z;                // iterative computation of ((2*m + 3)*z)

    // omega_0_0 upto omega_p_p-1
    for (size_t m = 0; m < p; ++m) {
        // omega_m_m  (from previous iteration)
        real_type omega_m_m = omega_mplus1_mplus1;
        omega(m, m) %= complex_x_y_m * omega_m_m;

        // omega_m+1_m+1  (for the next iteration)
        omega_mplus1_mplus1 = VT.same(p2m::__reciprocal(scalar_type(2 * (m + 1)))) * omega_m_m;

        // omega_m+1_m
        real_type omega_mplus1_m = z * omega_m_m;
        omega(m + 1, m) %= complex_x_y_m * omega_mplus1_m;

        // omega_m+2_m upto omega_p_m
        real_type omega_lminus2_m = omega_m_m;
        real_type omega_lminus1_m = omega_mplus1_m;
        real_type f_l = e_m;                    // iterative computation of ((2*l - 1)*z)
        for (size_t l = m + 2; l <= p; ++l) {
            // omega_l_m
            real_type omega_l_m = VT.same(p2m::__reciprocal(scalar_type(l * l - m * m))) *
                (f_l * omega_lminus1_m - dist_squared * omega_lminus2_m);
            omega(l, m) %= complex_x_y_m * omega_l_m;

            // (for the next iteration)
            omega_lminus2_m = omega_lminus1_m;
            omega_lminus1_m = omega_l_m;
            f_l += twice_z;
        }

        // (for the next iteration)
        complex_x_y_m *= complex_x_y;
        e_m += twice_z;
    }

    // omega_p_p
    omega(p, p) %= complex_x_y_m * omega_mplus1_mplus1;
}

template <typename CoefficientMatrix, typename Real3, typename Real4>
__global__
void
__P2M_tree(size_t p,
           size_t offset,
           size_t *box_particle_offset,
           Real3 *expansion_point,
           Real4 *ordered_particles,
           size_t *block_map,
           size_t *offset_map,
           CoefficientMatrix ** omega
           )
{
    //maps particles blocks to cuda threadblocks
    const size_t id = block_map[blockIdx.x*blockDim.y + threadIdx.y];
    //threadblocks contain some index overhead as they have constant size
    //particles are enumerated coninously so the thread index (particle index) has to be corrected by offset computed on host
    const size_t i = (blockIdx.x*blockDim.y+threadIdx.y)*blockDim.x + threadIdx.x - offset_map[blockIdx.x*blockDim.y+threadIdx.y];

    typedef typename Real3::value_type Real;

    Real3 xyz = Real3(0.,0.,0.);
    Real q = 0.;
    size_t idx = offset + id;
    if(i < box_particle_offset[id+1]) 
    {

        Real3 expn_point = expansion_point[idx];
        Real4 particle = ordered_particles[i];
        xyz = Real3(particle.x,particle.y,particle.z) - expn_point;
        q = particle.q;

        __P2M_cubreduce(q, make_xyzq<Real>(xyz, q), *omega[idx], p);
        //__P2M_reduce(q, make_xyzq<Real>(xyz, q), *omega[idx], p);
        // far field Field
    }
}

}//namespace end

#endif
