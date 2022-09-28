#ifndef _BK_cuda_FE
#define _BK_cuda_FE

#include "cuda_lib.hpp"
#include <vector>
#include "local2local_derivative.hpp"
#include "cub_h.h"

namespace gmx_gpu_fmm{

namespace fe{

DEVICE
__forceinline__
float __rsqrt(float x)
{
    return rsqrtf(x);
}

DEVICE
__forceinline__
double __rsqrt(double x)
{
    return 1.0/sqrt(x);
}

template<typename Real, typename Real3>
DEVICE
__forceinline__
Real __squaredlength(Real3 x)
{
    return x.x * x.x  +  x.y * x.y  +  x.z * x.z;
}

template <typename Real>
DEVICE
__forceinline__
Real __rcplength(const XYZ<Real>& a)
{
    return __rsqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

template <typename Real>
DEVICE
__forceinline__
Real __reciprocal(Real a)
{
    return (Real) 1./ a ;
}

template <typename Real, typename Real3>
CUDA
__forceinline__
void
one_coulomb(Real3& x_i, Real3& x_j, Real& q_j, Real3 &efield, Real &potential)
{
    //coulomb
    Real3 diff = x_i - x_j;                 // (x_i - x_j)
    Real rlen = __rcplength(diff);            // 1 / |x_i - x_j|
    Real q_j_rlen = q_j * rlen;             // q_j / |x_i - x_j|
    Real rlen2 = rlen * rlen;               // 1 / |x_i - x_j|^2
    Real q_j_rlen3 = q_j_rlen * rlen2;      // q_j / |x_i - x_j|^3
    efield += diff * q_j_rlen3;             // (x_i - x_j) * q_j / |x_i - x_j|^3
    potential += q_j_rlen;
}

template <typename Real, typename Real3>
CUDA
__forceinline__
void
one_coulomb(Real3& diff, Real& q_j, Real3 &efield, Real &potential)
{
    //coulomb
    Real rlen = __rcplength(diff);            // 1 / |x_i - x_j|
    Real q_j_rlen = q_j * rlen;             // q_j / |x_i - x_j|
    Real rlen2 = rlen * rlen;               // 1 / |x_i - x_j|^2
    Real q_j_rlen3 = q_j_rlen * rlen2;      // q_j / |x_i - x_j|^3
    efield += diff * q_j_rlen3;             // (x_i - x_j) * q_j / |x_i - x_j|^3
    potential += q_j_rlen;
}

}

template <typename CoefficientMatrix,typename Real, typename Real3, typename Real4>
__global__
void
__MU_derivative_ref(size_t p,
     size_t offset,
     size_t num_boxes_lowest,
     CoefficientMatrix ** mu,
     CoefficientMatrix ** dmux,
     CoefficientMatrix ** dmuy,
     CoefficientMatrix ** dmuz
     )
{
    typedef typename CoefficientMatrix::complex_type complex_type;
    typedef typename complex_type::value_type scalar_type;
    size_t id = blockIdx.x * blockDim.x + threadIdx.x;
    size_t idx = offset + id;


    if(id < num_boxes_lowest)
    {
        CoefficientMatrix & mu_(*mu[idx]);
        DxL2L_reference(mu_, *dmux[id], p);
        DyL2L_reference(mu_, *dmuy[id], p);
        DzL2L_reference(mu_, *dmuz[id], p);
    }
}

template <typename CoefficientMatrix,typename Real, typename Real3, typename Real4>
__global__
void
__MU_derivative(size_t p,
     size_t offset,
     size_t num_boxes_lowest,
     CoefficientMatrix ** mu,
     CoefficientMatrix ** dmux,
     CoefficientMatrix ** dmuy,
     CoefficientMatrix ** dmuz
     )
{
    typedef typename CoefficientMatrix::complex_type complex_type;
    typedef typename complex_type::value_type scalar_type;
    size_t id = blockIdx.x;
    size_t idx = offset + id;
    CoefficientMatrix & mu_(*mu[idx]);

    if(id < num_boxes_lowest)
    {
        ssize_t l = threadIdx.x;
        ssize_t m = threadIdx.y;

        if(m < l+1 && l < p)
        {
            complex_type tmp1 = mu_.get_vectorized(l + 1, m + 1);
            complex_type tmp2 = mu_.get_vectorized(l + 1, m - 1);
            (*dmux[id])(l, m) = (tmp1 - tmp2) * scalar_type(-0.5);
            complex_type tmp  = (tmp1 + tmp2) * scalar_type(0.5);
            (*dmuy[id])(l, m) = complex_type(-tmp.imag(), tmp.real());
            (*dmuz[id])(l, m) = -mu_.get_vectorized(l + 1, m);
        }

        //(*dmux[id])(p, m) = complex_type(0.);
        //(*dmuy[id])(p, m) = complex_type(0.);
        //(*dmuz[id])(p, m) = complex_type(0.);
    }
}

template <typename CoefficientMatrix,typename Real, typename Real3, typename Real4>
__global__
void
__MU_derivative2(size_t p,
     size_t offset,
     size_t num_boxes_lowest,
     CoefficientMatrix ** mu,
     CoefficientMatrix ** dmux,
     CoefficientMatrix ** dmuy,
     CoefficientMatrix ** dmuz
     )
{
    typedef typename CoefficientMatrix::complex_type complex_type;
    typedef typename complex_type::value_type scalar_type;
    size_t id = blockIdx.x;
    size_t idx = offset + id;
    CoefficientMatrix & mu_(*mu[idx]);

    if(id < num_boxes_lowest)
    {
        ssize_t l = blockIdx.y;
        ssize_t m = blockIdx.z;

        if(m < l+1 && l < p)
        {
            complex_type tmp1 = mu_.get_vectorized(l + 1, m + 1);
            complex_type tmp2 = mu_.get_vectorized(l + 1, m - 1);
            (*dmux[id])(l, m) = (tmp1 - tmp2) * scalar_type(-0.5);
            complex_type tmp  = (tmp1 + tmp2) * scalar_type(0.5);
            (*dmuy[id])(l, m) = complex_type(-tmp.imag(), tmp.real());
            (*dmuz[id])(l, m) = -mu_.get_vectorized(l + 1, m);
        }

        //(*dmux[id])(p, m) = complex_type(0.);
        //(*dmuy[id])(p, m) = complex_type(0.);
        //(*dmuz[id])(p, m) = complex_type(0.);
    }
}

template <typename outputadapter_type, typename CoefficientMatrix,typename VCoefficientMatrix,typename Real, typename Real3, typename Real4>
__global__
void
__FE(size_t p,
     size_t offset,
     outputadapter_type *result,
     size_t *box_particle_offset,
     CoefficientMatrix ** mu,
     Real3 *expansion_point,
     CoefficientMatrix ** dmux,
     CoefficientMatrix ** dmuy,
     CoefficientMatrix ** dmuz,
     Real4 *ordered_particles,
     size_t *block_map,
     size_t *offset_map,
     size_t stream_offset
     )
{
    /*
    //maps particles blocks to cuda treadblocks
    const size_t id = block_map[blockIdx.x+stream_offset];
    //threadblocks contain some index overhead as they have constant size
    //particles are enumerated coninuously so the thread index (particle index) has to be corrected by offset computed on host
    const size_t i = (blockIdx.x+stream_offset) * blockDim.x + threadIdx.x - offset_map[blockIdx.x+stream_offset];
    */

    //maps particles blocks to cuda treadblocks
    const size_t id = block_map[blockIdx.x*blockDim.y + threadIdx.y + stream_offset];
    //threadblocks contain some index overhead as they have constant size
    //particles are enumerated coninuously so the thread index (particle index) has to be corrected by offset computed on host
    const size_t i = (blockIdx.x*blockDim.y+threadIdx.y + stream_offset)*blockDim.x + threadIdx.x - offset_map[blockIdx.x*blockDim.y+threadIdx.y + stream_offset];

    if(i < box_particle_offset[id+1])
    {
        Real q = -1.;  // Field computation only
        typedef typename CoefficientMatrix::complex_type complex_type;
        {
            //shared memory?
            size_t idx = offset + id;
            CoefficientMatrix & mu_(*mu[idx]);  // needed for each thread
            Real3 expansion_point_mu_ = expansion_point[idx]; // compute locally

            CoefficientMatrix *dmux_id = *(dmux+id);
            CoefficientMatrix *dmuy_id = *(dmuy+id);
            CoefficientMatrix *dmuz_id = *(dmuz+id);

            {
                //vectorizing does not help a lot
                Real3 xyz = Real3(ordered_particles[i].x, ordered_particles[i].y, ordered_particles[i].z) - expansion_point_mu_;
                //Real3 xyz = ordered_particles[i] - expansion_point_mu_;

                // far field Field
                // FIXME: optimize multiplication: (re, 0) = (re, im) * (re, im)
                complex_type Phi_i(0.);
                complex_type Fx(0.);
                complex_type Fy(0.);
                complex_type Fz(0.);
                {
                    //size_t index = i;
                    typedef typename Real3::value_type real_type;
                    //typedef typename VCoefficientMatrix::value_type complex_type;
                    typedef typename Vec_traits<real_type>::scalar_type scalar_type;
                    Vec_traits<real_type> VT;

                    const real_type x = xyz.x;
                    const real_type y = xyz.y;
                    const real_type z = xyz.z;
                    //const real_type q = q_;

                    real_type dist_squared = x * x + y * y + z * z;
                    real_type twice_z = z + z;
                    complex_type complex_x_y(x, -y);
                    complex_type complex_x_y_m(VT.same(1));     // iterative computation of ((x - iy)^m)
                    complex_type omega_lm(0.0);
                    complex_type omega_l_minus_m(0.0);

                    // omega_0_0  (for the first iteration)
                    real_type omega_mplus1_mplus1 = q;
                    real_type e_m = z + twice_z;                // iterative computation of ((2*m + 3)*z)

                    //if m==0
                    real_type scale(0.);
                    complex_type mu_m(0.0);

                    complex_type dmux_(0.0);
                    complex_type dmuy_(0.0);
                    complex_type dmuz_(0.0);

                    // omega_0_0 upto omega_p_p-1
                    for (size_t m = 0; m < p; ++m)
                    {
                        // omega_m_m  (from previous iteration)
                        real_type omega_m_m = omega_mplus1_mplus1;
                        
                        omega_lm = complex_x_y_m * omega_m_m;
                        omega_l_minus_m = cuda_toggle_sign_if_odd(m, conj(omega_lm)) * scale;
                        mu_m = mu_.get_vectorized(m, m);

                        Phi_i -= mu_m * omega_lm + cuda_toggle_sign_if_odd(m, conj(mu_m)) * omega_l_minus_m;

                        dmux_ = dmux_id->get_vectorized(m, m);
                        dmuy_ = dmuy_id->get_vectorized(m, m);
                        dmuz_ = dmuz_id->get_vectorized(m, m);

                        Fx += dmux_ * omega_lm + cuda_toggle_sign_if_odd(m, conj(dmux_)) * omega_l_minus_m;
                        Fy += dmuy_ * omega_lm + cuda_toggle_sign_if_odd(m, conj(dmuy_)) * omega_l_minus_m;
                        Fz += dmuz_ * omega_lm + cuda_toggle_sign_if_odd(m, conj(dmuz_)) * omega_l_minus_m;
                        
                        // omega_m+1_m+1  (for the next iteration)
                        omega_mplus1_mplus1 = VT.same(fe::__reciprocal(scalar_type(2 * (m + 1)))) * omega_m_m;

                        // omega_m+1_m
                        real_type omega_mplus1_m = z * omega_m_m;
                        
                        omega_lm = complex_x_y_m * omega_mplus1_m;
                        omega_l_minus_m = cuda_toggle_sign_if_odd(m, conj(omega_lm)) * scale;

                        mu_m = mu_.get_vectorized(m+1, m);

                        dmux_ = dmux_id->get_vectorized(m+1, m);
                        dmuy_ = dmuy_id->get_vectorized(m+1, m);
                        dmuz_ = dmuz_id->get_vectorized(m+1, m);
                        
                        Phi_i -=mu_m * omega_lm + cuda_toggle_sign_if_odd(m, conj(mu_m)) * omega_l_minus_m;
                        Fx += dmux_ * omega_lm + cuda_toggle_sign_if_odd(m, conj(dmux_)) * omega_l_minus_m;;
                        Fy += dmuy_ * omega_lm + cuda_toggle_sign_if_odd(m, conj(dmuy_)) * omega_l_minus_m;
                        Fz += dmuz_ * omega_lm + cuda_toggle_sign_if_odd(m, conj(dmuz_)) * omega_l_minus_m;

                        // omega_m+2_m upto omega_p_m
                        real_type omega_lminus2_m = omega_m_m;
                        real_type omega_lminus1_m = omega_mplus1_m;
                        real_type f_l = e_m;                    // iterative computation of ((2*l - 1)*z)
                        for (size_t l = m + 2; l <= p; ++l)
                        {
                            // omega_l_m
                            real_type omega_l_m = VT.same(fe::__reciprocal(scalar_type(l * l - m * m))) *
                                    (f_l * omega_lminus1_m - dist_squared * omega_lminus2_m);
                             
                            omega_lm = complex_x_y_m * omega_l_m;
                        
                            omega_l_minus_m = cuda_toggle_sign_if_odd(m, conj(omega_lm)) * scale;

                            mu_m = mu_.get_vectorized(l, m);

                            dmux_ = dmux_id->get_vectorized(l, m);
                            dmuy_ = dmuy_id->get_vectorized(l, m);
                            dmuz_ = dmuz_id->get_vectorized(l, m);

                            Phi_i -=mu_m * omega_lm + cuda_toggle_sign_if_odd(m, conj(mu_m)) * omega_l_minus_m;
                            Fx += dmux_ * omega_lm + cuda_toggle_sign_if_odd(m, conj(dmux_)) * omega_l_minus_m;;
                            Fy += dmuy_ * omega_lm + cuda_toggle_sign_if_odd(m, conj(dmuy_)) * omega_l_minus_m;
                            Fz += dmuz_ * omega_lm + cuda_toggle_sign_if_odd(m, conj(dmuz_)) * omega_l_minus_m;

                            // (for the next iteration)
                            omega_lminus2_m = omega_lminus1_m;
                            omega_lminus1_m = omega_l_m;
                            f_l += twice_z;
                        }

                        // (for the next iteration)
                        complex_x_y_m *= complex_x_y;
                        e_m += twice_z;
                        scale = 1.0;
                    }                     
                    omega_lm = complex_x_y_m * omega_mplus1_mplus1;
                    omega_l_minus_m = cuda_toggle_sign_if_odd(p, conj(omega_lm));

                    mu_m = mu_.get_vectorized(p, p);

                    dmux_ = dmux_id->get_vectorized(p, p);
                    dmuy_ = dmuy_id->get_vectorized(p, p);
                    dmuz_ = dmuz_id->get_vectorized(p, p);

                    Phi_i -=mu_m * omega_lm + cuda_toggle_sign_if_odd(p, conj(mu_m)) * omega_l_minus_m;
                    Fx += dmux_ * omega_lm + cuda_toggle_sign_if_odd(p, conj(dmux_)) * omega_l_minus_m;;
                    Fy += dmuy_ * omega_lm + cuda_toggle_sign_if_odd(p, conj(dmuy_)) * omega_l_minus_m;
                    Fz += dmuz_ * omega_lm + cuda_toggle_sign_if_odd(p, conj(dmuz_)) * omega_l_minus_m;
                }

                //atomic not needed but faster in this case because vectorized
                //no blocking, each thread is writing on his own memory address
                q = ordered_particles[i].q;
                result->atomic_reduce_pf(i, Phi_i.real(), -Real3(Fx.real(), Fy.real(), Fz.real())*q);
            }
        }
    }
}

__global__
void
map_fmm_ids_to_original(size_t* orig_ids, size_t* fmm_ids, size_t n)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= n)
        return;

    fmm_ids[orig_ids[i]] = i;
}

template <typename Real4, typename outputadapter_type, typename Real33>
__global__
void
compute_exclusionsv3(Real4* unordered_particles, outputadapter_type* result, size_t* fmm_ids, int2* exclusion_pairs, Real33 abc, Real33 halfabc, size_t n)
{
    typedef typename outputadapter_type::Real Real;
    typedef typename outputadapter_type::Real3 Real3;
    size_t ti = blockIdx.x * blockDim.x + threadIdx.x;
    if(ti >= n)
        return;

    Real3 efield(0.,0.,0.);
    Real potential = 0;

    size_t i = exclusion_pairs[ti].x;
    size_t j = exclusion_pairs[ti].y;

    REAL4 tmp  = *reinterpret_cast<REAL4*>(&unordered_particles[i]);
    Real3 x_i  = Real3(tmp.x, tmp.y, tmp.z);
    Real q_i   = tmp.w;

    tmp = *reinterpret_cast<REAL4*>(&unordered_particles[j]);
    Real3 diff  = x_i - Real3(tmp.x, tmp.y, tmp.z);
    Real q_j   = tmp.w;

    if (diff.x > halfabc.a.x)
    {
        diff.x -= abc.a.x;
    }
    if(diff.x <= -halfabc.a.x)
    {
        diff.x += abc.a.x;
    }
    if (diff.y > halfabc.b.y)
    {
        diff.y -= abc.b.y;
    }
    if(diff.y <= -halfabc.b.y)
    {
        diff.y += abc.b.y;
    }
    if (diff.z > halfabc.c.z)
    {
        diff.z -= abc.c.z;
    }
    if(diff.z <= -halfabc.c.z)
    {
        diff.z += abc.c.z;
    }
    fe::one_coulomb(diff, tmp.w, efield, potential);

    result->atomic_reduce_pf(fmm_ids[i], -potential, -efield*q_i);
    result->atomic_reduce_pf(fmm_ids[j], -potential,  efield*q_i);
}

template <typename Real4, typename outputadapter_type, typename Real33>
__global__
void
compute_exclusionsv2(Real4* unordered_particles, outputadapter_type* result, size_t* fmm_ids, int** excl, int* excl_sizes, Real33 abc, Real33 halfabc, size_t n)
{
    typedef typename outputadapter_type::Real Real;
    typedef typename outputadapter_type::Real3 Real3;

    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= n)
        return;

    Real3 efield(0.,0.,0.);
    Real potential = 0;

    REAL4 tmp  = *reinterpret_cast<REAL4*>(&unordered_particles[i]);
    Real3 x_i  = Real3(tmp.x, tmp.y, tmp.z);
    Real q_i   = tmp.w;

    int exclsizes = excl_sizes[i];
    int* exclptr  = excl[i];

    for(size_t excl_i = 0; excl_i < exclsizes; ++excl_i)
    {
        int j = exclptr[excl_i];

        tmp = *reinterpret_cast<REAL4*>(&unordered_particles[j]);

        Real3 diff  = x_i - Real3(tmp.x, tmp.y, tmp.z);

        if (diff.x > halfabc.a.x)
        {
            diff.x -= abc.a.x;
        }
        if(diff.x <= -halfabc.a.x)
        {
            diff.x += abc.a.x;
        }
        if (diff.y > halfabc.b.y)
        {
            diff.y -= abc.b.y;
        }
        if(diff.y <= -halfabc.b.y)
        {
            diff.y += abc.b.y;
        }
        if (diff.z > halfabc.c.z)
        {
            diff.z -= abc.c.z;
        }
        if(diff.z <= -halfabc.c.z)
        {
            diff.z += abc.c.z;
        }
        fe::one_coulomb(diff, tmp.w, efield, potential);
    }
    if (exclsizes > 0)
        result->atomic_reduce_pf(fmm_ids[i], -potential, -efield*q_i);
}

template <typename Real4, typename outputadapter_type, typename Real33>
__global__
void
compute_exclusions(Real4* unordered_particles, outputadapter_type* result, size_t* orig_ids, int** excl, int* excl_sizes, Real33 abc, size_t n)
{

    typedef typename outputadapter_type::Real Real;
    typedef typename outputadapter_type::Real3 Real3;
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= n)
        return;
    extern __shared__ Real3 shifts[];

    if(threadIdx.x == 0)
    {
        int index = 0;
        for (int a = -1; a <= +1; ++a)
        {
            for (int b = -1; b <= +1; ++b)
            {
                for (int c = -1; c <= +1; ++c)
                {
                    shifts[index++] = abc * Real3(a, b, c);
                }
            }
        }
    }
    __syncthreads();

    Real3 efield(0.,0.,0.);
    Real potential = 0;

    int orig_i = orig_ids[i];
    REAL4 tmp = *reinterpret_cast<REAL4*>(&unordered_particles[orig_i]);
    Real3 x_i = Real3(tmp.x, tmp.y, tmp.z);
    Real q_i = tmp.w;

    int exclsizes = excl_sizes[orig_i];
    int* exclptr = excl[orig_i];

    for(size_t j = 0; j < exclsizes; ++j)
    {
        int orig_j = exclptr[j];

        if(orig_i != orig_j)
        {
            tmp = *reinterpret_cast<REAL4*>(&unordered_particles[orig_j]);

            int offsetindex = 13;

            Real3 diff  = Real3(tmp.x, tmp.y, tmp.z) - x_i;

            Real  dist2 = fe::__squaredlength<Real,Real3>(diff);
            Real  d2;
            if(dist2 > shifts[26].x)
            {
                for (int s = 0; s < 27; ++s)
                {
                    d2 = fe::__squaredlength<Real, Real3>(diff - shifts[s]);
                    if (d2 < dist2)
                    {
                        dist2  = d2;
                        offsetindex = s;
                    }
                }
            }
            diff += x_i - shifts[offsetindex];
#if 0
            Real  dd    = fe::__squaredlength<Real, Real3>(diff);
            if(dist2 != dd)
            {
                printf("orig distance %e, corrected distance %e, index pair(%lu, %lu)\n",sqrt(dd)/abc.a.x, sqrt(dist2)/abc.a.x, orig_i, orig_j);
            }
#endif

            //std::cout<<orig_i<<"---"<<unordered_particles[orig_i]<<" - "<<orig_j<<"---"<<unordered_particles[orig_j]<<x_jj<<dist2<<std::endl;
            fe::one_coulomb(x_i, diff, tmp.w, efield, potential);
        }
    }
    if (exclsizes > 1)
        result->atomic_reduce_pf(i, -potential, -efield*q_i);
}

template <typename Real, int BLOCK_THREADS, int ITEMS_PER_THREAD,  cub::BlockReduceAlgorithm ALGORITHM>
__global__
void __energy_kernel(Real* potential, Real *Ec, size_t n)
{
    // Specialize BlockReduce type for our thread block
    typedef cub::BlockReduce<Real, BLOCK_THREADS, ALGORITHM> BlockReduceT;

    // Shared memory
    __shared__ typename BlockReduceT::TempStorage temp_storage;

    // Per-thread tile data
    Real data[ITEMS_PER_THREAD];
    for(int i = 0; i < ITEMS_PER_THREAD; ++i)
    {
        int index = blockIdx.x*BLOCK_THREADS*ITEMS_PER_THREAD + threadIdx.x*ITEMS_PER_THREAD + i;
        if(index >= n)
            data[i] = 0.0;
        else
            data[i] = potential[index];
    }
    __syncthreads();
    //cub::LoadDirectStriped<BLOCK_THREADS>(threadIdx.x, &potential[blockIdx.x*BLOCK_THREADS*ITEMS_PER_THREAD], data);
    Real aggregate = BlockReduceT(temp_storage).Sum(data);
    if (threadIdx.x == 0)
    {
        Ec[blockIdx.x] = aggregate;
    }
}

template <typename Real4, typename Real>
__global__
void __charge_potential(Real* potential, Real4* particles, size_t n)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= n)
        return;
    potential[i] = potential[i] * particles[i].q;
}

template <typename Real>
__global__
void __reduce_last_chunk(Real* potential, Real *Ec, int targetindex, size_t begin, size_t end)
{
    Real E = 0.;
    for (int i = begin; i < end; ++i)
    {
        E += potential[i];
    }
    Ec[targetindex] = E;
}

template <typename Real3, typename Real>
__global__
void
__ordered_to_unordered_forces(Real3* forces_input, Real3* forces_output, size_t* orig_ids, Real boxeps, size_t n)
{

    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= n)
        return;

    int orig_i    = orig_ids[i];
    REAL3 tmp     = *reinterpret_cast<REAL3*>(&forces_input[i]);
    Real3 force_i = Real3(tmp.x, tmp.y, tmp.z) * boxeps;
    forces_output[orig_i] = force_i;
}

template <typename Real>
__global__
void
__ordered_to_unordered_potential(Real* input, Real* output, size_t* orig_ids, Real boxeps, size_t n)
{

    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= n)
        return;

    int orig_i     = orig_ids[i];
    output[orig_i] = input[i] * boxeps;
}

template <typename Real3>
__global__ void
__set_fmm_to_gmx_forces(Real3* fmm_forces, float3* gmx_forces, int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= n)
        return;

    //printf("gpu gmx %d %e %e %e\n",i, gmx_forces[i].x, gmx_forces[i].y, gmx_forces[i].z);

    gmx_forces[i].x += fmm_forces[i].x;
    gmx_forces[i].y += fmm_forces[i].y;
    gmx_forces[i].z += fmm_forces[i].z;

    //printf("gpu fmm %d %e %e %e\n",i, fmm_forces[i].x, fmm_forces[i].y, fmm_forces[i].z);
}

template <typename outputadapter_type, typename Real, typename Real4>
__global__
void
__energy_dump_kernel(outputadapter_type *result_ptr, Real4 *ordered_particles, size_t n, Real scale)
{
    Real Ec = 0.;
    for (size_t i = 0; i < n; ++i)
    {
        Real q = ordered_particles[i].q;
        printf("%g   %.20e\n",q, result_ptr->vPhi_ptr[i]);
        Ec += q * result_ptr->vPhi_ptr[i];
    }
    printf("Energy        %.20e\n",Ec*0.5*scale);
}

template <typename outputadapter_type, typename Real, typename Real3>
__global__
void
__force_dump_kernel(outputadapter_type *result_ptr, size_t n, size_t type)
{
    Real3 Norm(0.0,0.0,0.0);

    if(type == 0)
    {
        for (size_t i = 0; i < n; ++i)
        {
            Real3 tmp = result_ptr->vF_ptr[i];
            printf("Forcex         %.20e\n", tmp.x);
            printf("Forcey         %.20e\n", tmp.y);
            printf("Forcez         %.20e\n", tmp.z);
            tmp *= tmp;
            Norm += tmp;
        }
    }
    else
    {
        for (size_t i = 0; i < n; ++i)
        {
            Real3 tmp = result_ptr->vF_ptr[i];
            tmp *= tmp;
            Norm += tmp;
        }
    }

    Norm /=n;

    if(type == 1)
        printf("Forcex         %.20e\n", Norm.xyz_sqrt().x);
    if(type == 2)
        printf("Forcey         %.20e\n", Norm.xyz_sqrt().y);
    if(type == 3)
        printf("Forcez         %.20e\n", Norm.xyz_sqrt().z);
}

}//namespace end
#endif
