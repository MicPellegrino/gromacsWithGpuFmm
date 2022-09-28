#ifndef _BK_FMM_cuda_L2L_hpp
#define _BK_FMM_cuda_L2L_hpp

#include "cub_h.h"

namespace gmx_gpu_fmm{

namespace l2l{

/*!
 * \brief __make_boxid
 * \param x
 * \param y
 * \param z
 * \param depth
 * \return
 */
__device__
size_t __make_boxid(size_t x, size_t y, size_t z, unsigned depth)
{
    return (z << (depth * 2)) + (y << depth) + x;
}

}

template <typename CoefficientMatrix>
__device__
void l2l_kmin_kmax_(
        const CoefficientMatrix  &mu_in,
        const L2L_Operator<CoefficientMatrix>  &C,
        CoefficientMatrix  &mu_out,
        const int p, int l, int m, int j, int k_min, int k_max)
{
   

    typedef typename CoefficientMatrix::complex_type complex_type;
    complex_type mu_l_m(0.);
    for (int k = k_min; k <= k_max; ++k)
    {
        mu_l_m += C.get(j-l, k-m) * mu_in.get(j, k);
    }
    mu_out(l, m) %=mu_l_m;
}


//last for loop not parallel
template <typename CoefficientMatrix>
__global__ void __L2Lkernel_sequential_p_reduction(size_t depth_offset_p,
                             size_t depth_offset_c,
                             size_t d_p, size_t d_c,
                             size_t d_delta,
                             CoefficientMatrix ** mu,
                             L2L_Operator<CoefficientMatrix> *** C,
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
    if(j_max>p || j_max< l)
        return;

    const int idx_p = depth_offset_p + l2l::__make_boxid(i,j,k,d_p);

    const int idx_c = depth_offset_c + l2l::__make_boxid(2*i+ii,2 *j+jj,2*k+kk, d_c);
       
    const int idx_op = l2l::__make_boxid(ii,jj,kk,d_delta);

    const int k_min = m - (j_max-l);
    const int k_max = m + (j_max-l);
    l2l_kmin_kmax_(*mu[idx_p], *C[d_p][idx_op], *mu[idx_c], (int)p, l,m,j_max, k_min,k_max);
}


template <typename CoefficientMatrix,typename CoefficientMatrixSoA, typename Box, typename Real, typename Real3, typename complex_type, typename inttype, int depthtype>
__global__ void
__L2L_one(Box *box,
         CoefficientMatrixSoA *muSoA,
         complex_type*** targets_SoA,
         size_t op_id,
         const size_t num_boxes_above,
         const size_t num_boxes_tree,
         const size_t depth_offset,
         const size_t p,
         const size_t p1,
         const size_t p1xx2,
         const size_t pxp1_2,
         const size_t p1xp2_2)
{

    typedef typename Box::L2L_Operator Operator;
    extern __shared__ complex_type shared_op[];
    ssize_t l = blockIdx.x;
    ssize_t m =  blockIdx.y;
    if(m>l)
        return;

    size_t box_id = depth_offset;
    size_t tx;

    Operator **C_ptr = box[box_id].c_operators;
    Operator *C2S;
    size_t target_box_index;

    if (depthtype == 0)
    {
        C2S = C_ptr[blockIdx.z];
        tx = threadIdx.x;
        target_box_index = box_id + threadIdx.x + blockIdx.z * (num_boxes_tree-1);
    }
    else
    {
        C2S = C_ptr[op_id];
        tx = blockIdx.z * blockDim.x + threadIdx.x;
        target_box_index = box_id + tx + op_id * (num_boxes_tree-1);
    }

    if(threadIdx.x >= l && threadIdx.x < p1)
    {
        ssize_t j = threadIdx.x;
        size_t index = 0;
        ssize_t k_min;
        ssize_t k_max;

        for(ssize_t j_ = l; j_ < j; ++j_)
        {
            k_min = m - (j_ - l);
            k_max = m + (j_ - l);
            index += k_max - k_min + 1;
        }
        k_min = m - (j - l);
        k_max = m + (j - l);
        for (ssize_t k = k_min; k <= k_max; ++k)
        {
            shared_op[index++] = C2S->get_vectorized(j-l, k-m);
        }
    }
    __syncthreads();

    if(tx < num_boxes_above)
    {
        complex_type mu_l_m(0.0);
        complex_type om;
        complex_type op;
        size_t index = 0;
        for (ssize_t j = l; j <= p; ++j)
        {
            const ssize_t k_min = m - (j - l);
            const ssize_t k_max = m + (j - l);

            for (ssize_t k = k_min; k <= k_max; ++k)
            {
                om = muSoA->getSoA_full(j, k, box_id + tx);
                //om = (*box[box_id + threadIdx.x].mu).get(j,k);
                //printf("%e %e\n",(om-om1).real(), (om-om1).imag());
                //op = C2S->get_vectorized(j-l, k-m);
                op = shared_op[index++];
                mu_l_m += om*op;
            }
        }
#ifdef L2L_SOA_OPTIMIZATION
        const size_t lm_index = (l * (l + 1)) / 2 + m;
        complex_type** mu_ptr = targets_SoA[lm_index];
        mu_ptr[target_box_index]->atomic_adder(mu_l_m);
#else
        if (depthtype == 0)
        {
            target_box_index = box[box_id + tx].c_target_ids[blockIdx.z];
            *(muSoA->get_SoA_ptr(l, m, target_box_index)) %= mu_l_m;
        }
        else
        {
            target_box_index = box[box_id + tx].c_target_ids[op_id];
            *(muSoA->get_SoA_ptr(l, m, target_box_index)) %= mu_l_m;
        }
#endif
    }
}

}//namespace end
#endif
