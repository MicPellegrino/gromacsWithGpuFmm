#ifndef _BK_FMM_cuda_M2M_hpp
#define _BK_FMM_cuda_M2M_hpp

#include "cub_h.h"
#include "cuda_keywords.hpp"

namespace gmx_gpu_fmm{

namespace m2m{

__device__
size_t __make_boxid(size_t x, size_t y, size_t z, unsigned depth)
{
    return (z << (depth * 2)) + (y << depth) + x;
}

}

template <typename CoefficientMatrixSoA, typename Real3>
__global__
void __get_dipole(CoefficientMatrixSoA* omegaSoA, Real3* dipole_device)
{

    dipole_device->x =  omegaSoA->get_lin_SoA_ptr(threadIdx.x + 2,0)->real() * 2.0;
    dipole_device->y = -omegaSoA->get_lin_SoA_ptr(threadIdx.x + 2,0)->imag() * 2.0;
    dipole_device->z =  omegaSoA->get_lin_SoA_ptr(threadIdx.x + 1,0)->real();
}


template <typename CoefficientMatrix>
__device__
void m2m_kmin_kmax(
        const CoefficientMatrix& omega_in,
        const M2M_Operator<CoefficientMatrix>  &A,
        CoefficientMatrix  &omega_out,
        const int p, int l, int m, int j, int k_min, int k_max)
{
   

    typedef typename CoefficientMatrix::complex_type complex_type;
    complex_type omega_l_m(0.);
    for (int k = k_min; k <= k_max; ++k)
    {
        omega_l_m += A.get(l - j, m - k) *  omega_in.get(j, k);
    }
    omega_out(l, m) %=omega_l_m;
}


template <typename CoefficientMatrix>
__global__ void __M2Mkernel_sequential_p_reduction(size_t depth_offset_p,
                             size_t depth_offset_c,
                             size_t d_p, size_t d_c,
                             size_t d_delta,
                             CoefficientMatrix ** omega,
                             M2M_Operator<CoefficientMatrix> *** A,
                             size_t p,
                             size_t dim_p)
{
    int x =  blockIdx.x * blockDim.x + threadIdx.x;
    int y =  blockIdx.y * blockDim.y + threadIdx.y;
    int z =  blockIdx.z * blockDim.z + threadIdx.z;

    int z_range = gridDim.z*blockDim.z;
    int y_range = gridDim.y*blockDim.y;

    int q = x*y_range*z_range + y*z_range + z;

    int p_out = p+1;
    int p_out2 = p_out*p_out;
    int p_out3 = p_out*p_out2;

    int i  = (q/(dim_p*dim_p*8*p_out3))%dim_p;
    int j  = (q/(dim_p*8*p_out3))%dim_p;
    int k  = (q/(8*p_out3))%dim_p;
    int ii = (q/(4*p_out3))%2;
    int jj = (q/(2*p_out3))%2;
    int kk = (q/p_out3)%2;
    int l  = (q/p_out2)%p_out;
    int m  = (q/p_out)%p_out;
    int j_max = q%p_out;

    if(m>l)
        return;
    if(j_max>l)
        return;

    const int idx_p = depth_offset_p + m2m::__make_boxid(i,j,k,d_p);

    const int idx_c = depth_offset_c + m2m::__make_boxid(2*i+ii,2 *j+jj,2*k+kk, d_c);
       
    const int idx_op = m2m::__make_boxid(ii,jj,kk,d_delta);

    const int k_min = static_cast <int>(max(static_cast<long long>(-j_max), static_cast<long long>(m - (l - j_max))));
    const int k_max = static_cast <int>(min(static_cast<long long>(+j_max), static_cast<long long>(m + (l - j_max))));

    m2m_kmin_kmax(*omega[idx_c], *A[d_p][idx_op], *omega[idx_p], (int)p, l,m,j_max, k_min,k_max);
}

template <typename Box>
__global__ void
__activate_boxes_above(Box *box, const size_t depth_offset)
{
    size_t box_id = depth_offset + blockDim.x*blockIdx.x + threadIdx.x;

    if(box[box_id].active == 1)
    {
        size_t target_id = box[box_id].a_target_ids[0];
        box[target_id].active = 1;
    }
}

template <typename CoefficientMatrix,typename CoefficientMatrixSoA, typename Box, typename Real, typename Real3, typename complex_type, typename inttype, int depthtype>
__global__ void
__M2M_one(Box *box,
         CoefficientMatrixSoA *omegaSoA,
         complex_type*** targets_SoA,
         size_t child_box_id,
         const size_t num_boxes_above,
         const size_t num_boxes_tree,
         const size_t depth_offset,
         const size_t p,
         const size_t p1,
         const size_t p1xx2,
         const size_t pxp1_2,
         const size_t p1xp2_2)
{

    typedef typename Box::M2M_Operator Operator;
    extern __shared__ complex_type shared_op[];
    size_t l = blockIdx.x;
    size_t m = blockIdx.y;
    if(m>l)
        return;

    size_t box_id, tx;

    if (depthtype == 0)
    {
        box_id = blockIdx.z * num_boxes_above + depth_offset;
        tx = threadIdx.x;
    }
    else
    {
        box_id = child_box_id * num_boxes_above + depth_offset;
        tx = blockIdx.z * blockDim.x + threadIdx.x;
    }

    Operator **A_ptr = box[box_id].a_operators;
    Operator *A2S = A_ptr[0];

    if(threadIdx.x < l+1)
    {
        ssize_t j = threadIdx.x;//%(l+1);
        ssize_t k_min;
        ssize_t k_max;
        size_t index = 0;
        for(ssize_t j_ = 0; j_ < j; ++j_)
        {
            k_min = static_cast <ssize_t>(max(static_cast<long long>(-j_), static_cast<long long>(m - (l - j_))));
            k_max = static_cast <ssize_t>(min(static_cast<long long>(+j_), static_cast<long long>(m + (l - j_))));
            index += k_max - k_min + 1;
        }
        k_min = static_cast <ssize_t>(max(static_cast<long long>(-j), static_cast<long long>(m - (l - j))));
        k_max = static_cast <ssize_t>(min(static_cast<long long>(+j), static_cast<long long>(m + (l - j))));

        for (ssize_t k = k_min; k <= k_max; ++k)
        {
            shared_op[index++] = A2S->get_vectorized(l-j, m-k);
        }
    }
    __syncthreads();

    if(tx < num_boxes_above)
    {
        complex_type omega_l_m(0.0);
        complex_type om;
        complex_type op;
        size_t index = 0;
        for (ssize_t j = 0; j <= l; ++j)
        {
            const ssize_t k_min = static_cast <ssize_t>(max(static_cast<long long>(-j), static_cast<long long>(m - (l - j))));
            const ssize_t k_max = static_cast <ssize_t>(min(static_cast<long long>(+j), static_cast<long long>(m + (l - j))));

            for (ssize_t k = k_min; k <= k_max; ++k)
            {
                om = omegaSoA->getSoA_full(j, k, box_id + tx);
                //om = (*box[box_id + threadIdx.x].omega).get(j,k);
                //printf("%e %e\n",(om-om1).real(), (om-om1).imag());
                //op = A2S->get_vectorized(l-j, m-k);
                op = shared_op[index++];
                omega_l_m += om*op;
            }
        }
#ifdef M2M_SOA_OPTIMIZATION
        const size_t lm_index = (l * (l + 1)) / 2 + m;
        size_t target_box_index = box_id - 1 + tx;
        complex_type** omega_ptr = targets_SoA[lm_index];
        omega_ptr[target_box_index]->atomic_adder(omega_l_m);
#else
        size_t target_box_id = box[box_id + tx].a_target_ids[0];
        *(omegaSoA->get_SoA_ptr(l, m, target_box_id)) %= omega_l_m;
#endif
    }
}

}//namespace end

#endif
