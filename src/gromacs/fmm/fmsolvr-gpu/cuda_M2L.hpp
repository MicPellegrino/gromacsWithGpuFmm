#ifndef _BK_FMM_cuda_M2L_hpp
#define _BK_FMM_cuda_M2L_hpp

#include "cuda_keywords.hpp"

namespace gmx_gpu_fmm{

namespace m2l{

__device__
size_t __make_boxid(size_t x, size_t y, size_t z, unsigned depth)
{
    return (z << (depth * 2)) + (y << depth) + x;
}

}

#ifndef GMX_FMM_DOUBLE
DEVICE int2 index_lm_map[11000];
#else
DEVICE int2 index_lm_map[1000];
#endif


__global__ void
__prepare_lm_map(size_t p)
{
    int index = 0;
    for (size_t l = 0; l <= p; ++l)
    {
        for (size_t m = 0; m <= l; ++m)
        {
            index_lm_map[index++] = make_int2(l,m);
        }
    }
    for (ssize_t l = 1; l <= p; ++l)
    {
        for (ssize_t m = 1; m <= l; ++m)
        {
            index_lm_map[index++] = make_int2(l,-m);
        }
    }
}

template <typename Box, typename CoefficientMatrix>
CUDA
__forceinline__
CoefficientMatrix** get_target_ptr_ptr(Box *box, size_t boxid)
{
    return (box+boxid)->b_targets;
}

template <typename Box, typename CoefficientMatrix>
CUDA
__forceinline__
CoefficientMatrix** get_operator_ptr_ptr(Box *box, size_t boxid)
{
    return (box+boxid)->b_operators;
}

template <typename CoefficientMatrix, typename VCoefficientMatrix, typename Real, typename Real3,typename complex_type>
__global__ void __M2Lkernel_dynamic_child(CoefficientMatrix **  omega,
                                          CoefficientMatrix ** mu,
                                          M2L_Operators<VCoefficientMatrix> **  Bv,
                                          size_t omega_id,
                                          size_t mu_id,
                                          size_t op_id,
                                          size_t p,
                                          size_t p_out,
                                          size_t p_out2,
                                          const bool open_boundary_conditions,
                                          const size_t dim,
                                          const size_t dim2,
                                          const size_t i_,
                                          const size_t j_,
                                          const size_t k_,
                                          const ssize_t ii,
                                          const ssize_t jj,
                                          const ssize_t kk,
                                          const size_t mu_size,
                                          const size_t op_size,
                                          const size_t d,
                                          size_t op_vec_id
                                          )

{
    const ssize_t x = blockIdx.x;
    const ssize_t y = blockIdx.y;
    const ssize_t z = blockIdx.z;

    const ssize_t ii_x = ii+x-1;
    const ssize_t jj_y = jj+y-1;
    const ssize_t kk_z = kk+z-1;
    if (open_boundary_conditions && !(i_ < ii_x || i_ >= ii_x+3 || j_ < jj_y || j_ >= jj_y+3 || k_ < kk_z || k_ >= kk_z+3))
        return;
    
    const ssize_t l = threadIdx.x;
    const ssize_t m = threadIdx.y;
    const ssize_t tx = l*p_out + m;
    //

    size_t mu_id_l;
    mu_id_l = mu_id + x + y*dim + z*dim2;

    size_t op_vec_id_l = op_vec_id + z*4 + y*2 + x;
    extern __shared__ complex_type cache[];
    complex_type* op_jkml = cache;
    complex_type* omega_jk =(&cache[op_size]);
    
    if(tx<p_out2)
    {
        if(l>0 && m>0)
        {
            if(m <= l)
            {

                //self
                complex_type opp(Bv[op_vec_id_l]->get_vectorized(l,m, op_id));
                ssize_t l2 = l*l;
                op_jkml[l2+l+m] = opp;
                complex_type omm( omega[omega_id]->get_vectorized(l,m));
                omega_jk[l2+l+m] = omm;

                ssize_t lm = -m;

                op_jkml[l2+l+lm] = cuda_toggle_sign_if_odd(m, conj(opp));
                omega_jk[l2+l+lm] = cuda_toggle_sign_if_odd(m, conj(omm));

                //bottom right
                ssize_t lp = l+p;
                ssize_t rbl2 = lp*lp;
                complex_type opp_br(Bv[op_vec_id_l]->get_vectorized(lp,m, op_id));
                op_jkml[rbl2+lp+m] = opp_br;

                //bottom left
                op_jkml[rbl2+lp+lm] = cuda_toggle_sign_if_odd(m, conj(opp_br));

                //corner bottom right
                ssize_t rbm = m+p;
                complex_type opp_rbm(Bv[op_vec_id_l]->get_vectorized(lp,rbm, op_id));
                op_jkml[rbl2+lp+rbm] = opp_rbm;

                //corner bottom left
                ssize_t lbm = -rbm;
                 op_jkml[rbl2+lp+lbm] = cuda_toggle_sign_if_odd(m, conj(opp_rbm));;
            }
            if(m > l)
            {

                //bottom right
                ssize_t lp = l+p;
                ssize_t rbl2 = lp*lp;
                complex_type opp_br(Bv[op_vec_id_l]->get_vectorized(lp,m, op_id));
                op_jkml[rbl2+lp+m] = opp_br;

                //bottom left
                ssize_t lm = -m;
                op_jkml[rbl2+lp+lm] = cuda_toggle_sign_if_odd(m, conj(opp_br));
            }
        }
        if(m==0)
        {
            {

                //self
                ssize_t l2 = l*l;
                op_jkml[l2+l] = Bv[op_vec_id_l]->get_vectorized(l,0, op_id);
                omega_jk[l2+l] = omega[omega_id]->get_vectorized(l,0);
                //bottom
                ssize_t lp = l+p;
                ssize_t rbl2 = lp*lp;
                op_jkml[rbl2+lp+0] = Bv[op_vec_id_l]->get_vectorized(lp,0, op_id);

            }
        }
        
    } 
__syncthreads();

#if 1
    ssize_t l_,t1,ll,mm,f,f1;
    ssize_t k_start,k_end;
    t1 = tx-l;
    l_ = t1/p_out;
    f = 1-(l-l_);
    f1 = 1-f;
    ll = f1*l_ + f*m;
    mm = f1*m + f*l_;

    complex_type mu_l_m(0.);

    for (ssize_t j=0; j<=p; j++)
    {
        k_start = f*j;
        k_end = f*j + j + f;

        int j_ = static_cast<int>(j);
        ssize_t jj2 = static_cast<ssize_t>(j_*j_);

        ssize_t jl = j+ll;
        ssize_t jl2 = jl*jl;
        ssize_t lk= 0;

        complex_type mu_l_m_j(0.);
        ssize_t jl2_jl_m = jl2+jl+mm;
        ssize_t jj2_j = jj2+j;


        for (ssize_t k_ = k_start; k_ < k_end; k_++)
        {
            lk = k_-j;
            complex_type op(op_jkml[jl2_jl_m+lk]);
            complex_type om(omega_jk[jj2_j+lk]);
            mu_l_m_j += op*om;
        }
        mu_l_m +=cuda_toggle_sign_if_odd(j, mu_l_m_j);
    }
    (*mu[mu_id_l])(ll, mm) %= mu_l_m;

#endif
}



template <typename CoefficientMatrix, typename VCoefficientMatrix, typename Real, typename Real3>
__global__ void __M2Lkernel_dynamic_parent_open(const size_t depth_offset,
                                    const size_t d,
                                    CoefficientMatrix ** omega,
                                    CoefficientMatrix ** mu,
                                    M2L_Operators<VCoefficientMatrix> ** Bv,
                                    const size_t p,
                                    const size_t p_out,
                                    const size_t p_out2,
                                    const size_t dim,
                                    const size_t dim2,
                                    const bool open_boundary_conditions,
                                    const ssize_t ws,
                                    const ssize_t ws_dim,
                                    const ssize_t ws_dim2,
                                    ssize_t io_beg,
                                    ssize_t jo_beg,
                                    ssize_t ko_beg,
                                    const size_t ii,
                                    const size_t jj,
                                    const size_t kk,
                                    const size_t iip,
                                    const size_t jjp,
                                    const size_t kkp,
                                    const size_t omega_id,
                                    size_t mu_size,
                                    size_t op_size,
                                    size_t mu_and_op
                                    )
{

    const ssize_t x = threadIdx.x*2;
    const ssize_t y = threadIdx.y*2;
    const ssize_t z = threadIdx.z*2;

    const ssize_t i  = io_beg + x;
    if (i < 0 || i >= dim)
    {
        return;
    }
    const ssize_t j  = jo_beg + y;
    if (j < 0 || j >= dim)
    {
        return;
    }
    const ssize_t k  = ko_beg + z;
    if (k < 0 || k >= dim)
    {
        return;
    }
    const size_t ip = i/2;
    const size_t jp = j/2;
    const size_t kp = k/2;
    if(ip != (iip) || jp != (jjp) || kp != (kkp))
    {
        typedef typename CoefficientMatrix::complex_type complex_type;
        size_t mu_id = depth_offset + m2l::__make_boxid(i, j, k, d);
        size_t op_vec_id = (d-1)*8;

        size_t iii = (ii-i) + static_cast<size_t>(__sad((int)(ip-ws)*2, (int)i,0));
        size_t jjj = (jj-j) + static_cast<size_t>(__sad((int)(jp-ws)*2, (int)j,0));
        size_t kkk = (kk-k) + static_cast<size_t>(__sad((int)(kp-ws)*2, (int)k,0));
        size_t op_id = iii*ws_dim2+jjj*ws_dim+kkk;

        dim3 b(p_out+1,p_out,1);
        dim3 g(2,2,2);
        
        __M2Lkernel_dynamic_child<CoefficientMatrix,VCoefficientMatrix,Real,Real3,complex_type><<<g,b,sizeof(complex_type)*(mu_and_op)>>>
        (omega, mu, Bv, omega_id, mu_id, op_id, p,p_out,p_out2,open_boundary_conditions,dim,dim2,ii,jj,kk,i,j,k,mu_size,op_size,d,op_vec_id);
    }
}


template <typename CoefficientMatrix, typename VCoefficientMatrix, typename Real, typename Real3>
__global__ void __M2Lkernel_dynamic_parent_periodic(const size_t depth_offset,
                                    const size_t d,
                                    CoefficientMatrix ** omega,
                                    CoefficientMatrix ** mu,
                                    M2L_Operators<VCoefficientMatrix> ** Bv,
                                    const size_t p,
                                    const size_t p_out,
                                    const size_t p_out2,
                                    const size_t dim,
                                    const size_t dim2,
                                    const bool open_boundary_conditions,
                                    const ssize_t ws,
                                    const size_t ws_dim,
                                    const size_t ws_dim2,
                                    ssize_t io_beg,
                                    ssize_t jo_beg,
                                    ssize_t ko_beg,
                                    const size_t ii,
                                    const size_t jj,
                                    const size_t kk,
                                    const size_t iip,
                                    const size_t jjp,
                                    const size_t kkp,
                                    const size_t omega_id,
                                    size_t mu_size,
                                    size_t op_size,
                                    size_t mu_and_op
                                    )
{

    const ssize_t x = threadIdx.x*2;
    const ssize_t y = threadIdx.y*2;
    const ssize_t z = threadIdx.z*2;

    const ssize_t i  = io_beg + x;
    const ssize_t j  = jo_beg + y;
    const ssize_t k  = ko_beg + z;

    ssize_t iff = i;
    while (iff < 0)
    {
        iff += dim;
    }
    while (iff >= dim)
    {
        iff -= dim;
    }
    ssize_t jf = j;
    while (jf < 0)
    {
        jf += dim;
    }
    while(jf >= dim)
    {
        jf -= dim;
    }
    ssize_t kf = k;
    while(kf < 0)
    {
        kf += dim;
    }while(kf >= dim)
    {
        kf -= dim;
    }

    const size_t ip = i/2;
    const size_t jp = j/2;
    const size_t kp = k/2;
    if(ip != (iip) || jp != (jjp) || kp != (kkp))
    {
        typedef typename CoefficientMatrix::complex_type complex_type;
        size_t mu_id = depth_offset + m2l::__make_boxid(iff, jf, kf, d);
        size_t op_vec_id = (d-1)*8;
        size_t iii = (ii-i) + static_cast<size_t>(__sad((int)(ip-ws)*2, (int)i,0));
        size_t jjj = (jj-j) + static_cast<size_t>(__sad((int)(jp-ws)*2, (int)j,0));
        size_t kkk = (kk-k) + static_cast<size_t>(__sad((int)(kp-ws)*2, (int)k,0));
        size_t op_id = iii*ws_dim2+jjj*ws_dim+kkk;

        dim3 b(p_out+1,p_out,1);
        dim3 g(2,2,2);

        __M2Lkernel_dynamic_child<CoefficientMatrix,VCoefficientMatrix,Real,Real3,complex_type><<<g,b,sizeof(complex_type)*(mu_and_op)>>>
        (omega, mu, Bv, omega_id, mu_id, op_id, p,p_out,p_out2,open_boundary_conditions,dim,dim2,ii,jj,kk,i,j,k,mu_size,op_size,d,op_vec_id);
    }
}

template <typename CoefficientMatrix,typename Operator, typename Box, typename Real, typename Real3, typename complex_type>
__global__ void __M2L_one_p2_no_shared(Box* box,
                             CoefficientMatrix **omega,
                             const size_t depth_offset,
                             const size_t num_boxes_tree,
                             const size_t p,
                             const size_t p1,
                             const size_t p1xx2,
                             const size_t op_p1xx2
                             )

{

    size_t first_box = threadIdx.y * gridDim.x + depth_offset;
    size_t op_id = blockIdx.y;
    size_t box_id = first_box + blockIdx.x;

    if(box[box_id].active == 0)
        return;

    size_t target_id = box[box_id].b_target_ids[op_id];

    if(box[target_id].active == 0)
      return;

    Operator* B =  box[box_id].b_operators[op_id];
    CoefficientMatrix* Omega = box[box_id].omega;

    const ssize_t l = threadIdx.x;
    const ssize_t m = blockIdx.z;

    if(m>l)
        return;

    complex_type mu_l_m(0.);

    for (ssize_t j=0; j<=p; j++)
    {
        complex_type mu_l_m_j(0.);
        for (ssize_t k = -j; k <= j; k++)
        {
            mu_l_m_j += Omega->get_vectorized(j,k) * B->get_vectorized(l+j,m+k);
        }
        mu_l_m +=cuda_toggle_sign_if_odd(j, mu_l_m_j);
    }

    CoefficientMatrix* mu =  box[box_id].b_targets[op_id];
    (*mu)(l, m) %= mu_l_m;
}

template <typename CoefficientMatrix,typename Operator, typename Box, typename Real, typename Real3, typename complex_type>
__global__ void __M2L_one_p2(Box* box,
                             CoefficientMatrix **omega,
                             const size_t depth_offset,
                             const size_t num_boxes_tree,
                             const size_t p,
                             const size_t p1,
                             const size_t p1xx2,
                             const size_t op_p1xx2
                             )

{
    extern __shared__ complex_type cache[];
    complex_type* op_jkml = cache;
    complex_type* omega_jk =(&cache[op_p1xx2]);

    size_t first_box = blockIdx.z * gridDim.x + depth_offset;
    size_t op_id     = blockIdx.y;
    size_t box_id    = first_box + blockIdx.x;

    Operator* B2S =  box[box_id].b_operators[op_id];
    CoefficientMatrix* Omega = box[box_id].omega;

    const ssize_t l = threadIdx.x;
    const ssize_t m = threadIdx.y;
    const ssize_t tx = l*p1 + m;

    if(tx<p1xx2)
    {
        if(l>0 && m>0)
        {
            if(m <= l)
            {
                //self
                complex_type opp(B2S->get_vectorized(l,m));
                ssize_t l2 = l*l;
                op_jkml[l2+l+m] = opp;
                complex_type omm(Omega->get(l,m));
                omega_jk[l2+l+m] = omm;

                ssize_t lm = -m;
                op_jkml[l2+l+lm] = cuda_toggle_sign_if_odd(m, conj(opp));
                omega_jk[l2+l+lm] = cuda_toggle_sign_if_odd(m, conj(omm));

                //bottom right
                ssize_t lp = l+p;
                ssize_t rbl2 = lp*lp;
                complex_type opp_br(B2S->get_vectorized(lp,m));
                op_jkml[rbl2+lp+m] = opp_br;

                //bottom left
                op_jkml[rbl2+lp+lm] = cuda_toggle_sign_if_odd(m, conj(opp_br));

                //corner bottom right
                ssize_t rbm = m+p;
                complex_type opp_rbm(B2S->get_vectorized(lp,rbm));
                op_jkml[rbl2+lp+rbm] = opp_rbm;

                //corner bottom left
                ssize_t lbm = -rbm;
                op_jkml[rbl2+lp+lbm] = cuda_toggle_sign_if_odd(m, conj(opp_rbm));
            }
            if(m > l)
            {
                //bottom right
                ssize_t lp = l+p;
                ssize_t rbl2 = lp*lp;
                complex_type opp_br(B2S->get_vectorized(lp,m));
                op_jkml[rbl2+lp+m] = opp_br;

                //bottom left
                ssize_t lm = -m;
                op_jkml[rbl2+lp+lm] = cuda_toggle_sign_if_odd(m, conj(opp_br));
            }
        }
        if(m==0)
        {
            //self
            ssize_t l2 = l*l;
            op_jkml[l2+l] = B2S->get_vectorized(l,0);
            omega_jk[l2+l] = Omega->get(l,0);
            //bottom
            ssize_t lp = l+p;
            ssize_t rbl2 = lp*lp;
            op_jkml[rbl2+lp+0] = B2S->get_vectorized(lp,0);
        }
    }
    __syncthreads();

    ssize_t l_,t1,ll,mm,f,f1;
    ssize_t k_start,k_end;
    t1 = tx-l;
    l_ = t1/p1;
    f = 1-(l-l_);
    f1 = 1-f;
    ll = f1*l_ + f*m;
    mm = f1*m + f*l_;

    complex_type mu_l_m(0.);

    for (ssize_t j=0; j<=p; j++)
    {
        k_start = f*j;
        k_end = f*j + j + f;

        int j_ = static_cast<int>(j);
        ssize_t jj2 = static_cast<ssize_t>(j_*j_);

        ssize_t jl = j+ll;
        ssize_t jl2 = jl*jl;
        ssize_t lk= 0;

        complex_type mu_l_m_j(0.);
        ssize_t jl2_jl_m = jl2+jl+mm;
        ssize_t jj2_j = jj2+j;

        for (ssize_t k_ = k_start; k_ < k_end; k_++)
        {
            lk = k_-j;
            complex_type op(op_jkml[jl2_jl_m+lk]);
            complex_type om(omega_jk[jj2_j+lk]);
            mu_l_m_j += op*om;
        }
        mu_l_m +=cuda_toggle_sign_if_odd(j, mu_l_m_j);
    }

    CoefficientMatrix* mu =  box[box_id].b_targets[op_id];
    (*mu)(ll, mm) %= mu_l_m;
}

template <typename CoefficientMatrix,typename CoefficientMatrixSoA, typename Box, typename Real, typename Real3, typename complex_type, bool SPARSE>
__global__ void __M2L_one_p2_SoA(Box *box,
                             CoefficientMatrixSoA *omegaSoA,
                             CoefficientMatrixSoA *muSoA,
                             complex_type*** targets_SoA,
                             const size_t depth_offset,
                             const size_t num_boxes_tree,
                             const size_t p,
                             const size_t p1,
                             const size_t p1xx2,
                             const size_t op_p1xx2
                             )

{
#if 1
    extern __shared__ complex_type cache[];
    complex_type* op_jkml = cache;
    complex_type* omega_jk =(&cache[op_p1xx2]);

    size_t first_box = blockIdx.z * gridDim.x + depth_offset;
    size_t op_id = blockIdx.y;
    size_t box_id = first_box + blockIdx.x;

    if(SPARSE)
    {
        if(box[box_id].active == 0)
            return;
        size_t target_id = box[box_id].b_target_ids[op_id];

        if(box[target_id].active == 0)
            return;
    }

    CoefficientMatrix* B2S =  box[box_id].b_operators[op_id];

    const ssize_t l = threadIdx.y;
    const ssize_t m = threadIdx.x;
    const ssize_t tx = l*p1 + m;

    if(tx<p1xx2)
    {
        if(l>0 && m>0)
        {
            if(m <= l)
            {
                //self
                complex_type opp(B2S->get_vectorized(l,m));
                ssize_t l2 = l*l;
                op_jkml[l2+l+m] = opp;
                complex_type omm(omegaSoA->getSoA(l,m,box_id));
                omega_jk[l2+l+m] = omm;

                ssize_t lm = -m;
                op_jkml[l2+l+lm] = cuda_toggle_sign_if_odd(m, conj(opp));
                omega_jk[l2+l+lm] = cuda_toggle_sign_if_odd(m, conj(omm));

                //bottom right
                ssize_t lp = l+p;
                ssize_t rbl2 = lp*lp;
                complex_type opp_br(B2S->get_vectorized(lp,m));
                op_jkml[rbl2+lp+m] = opp_br;

                //bottom left
                op_jkml[rbl2+lp+lm] = cuda_toggle_sign_if_odd(m, conj(opp_br));

                //corner bottom right
                ssize_t rbm = m+p;
                complex_type opp_rbm(B2S->get_vectorized(lp,rbm));
                op_jkml[rbl2+lp+rbm] = opp_rbm;

                //corner bottom left
                ssize_t lbm = -rbm;
                op_jkml[rbl2+lp+lbm] = cuda_toggle_sign_if_odd(m, conj(opp_rbm));
            }
            if(m > l)
            {
                //bottom right
                ssize_t lp = l+p;
                ssize_t rbl2 = lp*lp;
                complex_type opp_br(B2S->get_vectorized(lp,m));
                op_jkml[rbl2+lp+m] = opp_br;

                //bottom left
                ssize_t lm = -m;
                op_jkml[rbl2+lp+lm] = cuda_toggle_sign_if_odd(m, conj(opp_br));
            }
        }
        if(m==0)
        {
            //self
            ssize_t l2 = l*l;
            op_jkml[l2+l] = B2S->get_vectorized(l,0);
            omega_jk[l2+l] = omegaSoA->getSoA(l,0,box_id);
            //bottom
            ssize_t lp = l+p;
            ssize_t rbl2 = lp*lp;
            op_jkml[rbl2+lp+0] = B2S->get_vectorized(lp,0);
        }
    }
    __syncthreads();

    ssize_t l_,t1,ll,mm,f,f1;
    ssize_t k_start,k_end;
    t1 = tx-l;
    l_ = t1/p1;
    f = 1-(l-l_);
    f1 = 1-f;
    ll = f1*l_ + f*m;
    mm = f1*m + f*l_;

    complex_type mu_l_m(0.);

    for (ssize_t j=0; j<=p; j++)
    {
        k_start = f*j;
        k_end = f*j + j + f;

        int j_ = static_cast<int>(j);
        ssize_t jj2 = static_cast<ssize_t>(j_*j_);

        ssize_t jl = j+ll;
        ssize_t jl2 = jl*jl;
        ssize_t lk= 0;

        complex_type mu_l_m_j(0.);
        ssize_t jl2_jl_m = jl2+jl+mm;
        ssize_t jj2_j = jj2+j;

        for (ssize_t k_ = k_start; k_ < k_end; k_++)
        {
            lk = k_-j;
            complex_type op(op_jkml[jl2_jl_m+lk]);
            complex_type om(omega_jk[jj2_j+lk]);
            mu_l_m_j += op*om;
        }
        mu_l_m +=cuda_toggle_sign_if_odd(j, mu_l_m_j);
    }

#ifdef M2L_SOA_OPTIMIZATION
    const size_t lm_index = (ll * (ll + 1)) / 2 + mm;
    size_t target_box_index = box_id - 1 + op_id * (num_boxes_tree-1);
    complex_type** mu_ptr = targets_SoA[lm_index];
    mu_ptr[target_box_index]->atomic_adder(mu_l_m);
#else
    size_t target_box_index = box[box_id].b_target_ids[op_id];
    *(muSoA->get_SoA_ptr(ll, mm, target_box_index)) %= mu_l_m;
#endif

#else

    extern __shared__ Real real_cache[];
    Real* op_jkml_real  = real_cache;
    Real* op_jkml_imag  = &real_cache[op_p1xx2];
    Real* omega_jk_real = (&real_cache[op_p1xx2*2]);
    Real* omega_jk_imag = (&real_cache[op_p1xx2*2 + p1xx2]);

    size_t first_box = blockIdx.z * gridDim.x + depth_offset;
    size_t op_id = blockIdx.y;
    size_t box_id = first_box + blockIdx.x;

    if(box[box_id].active == 0)
        return;

    size_t target_id = box[box_id].b_target_ids[op_id];

    if(box[target_id].active == 0)
        return;

    CoefficientMatrix* B =  box[box_id].b_operators[op_id];

    const ssize_t l = threadIdx.y;
    const ssize_t m = threadIdx.x;
    const ssize_t tx = l*p1 + m;

    int ll, mm, l1, t1, f, f1;

    t1 = tx - l;
    l1 = t1/p1;
    f = 1 -(l - l1);
    f1 = 1 - f;
    ll = f1 * l1 + f * m;
    mm = f1 * m  + f * l1;

    if(f1)
    {
        //self
        int l2 = ll*ll;
        int index = l2+ll+mm;
        complex_type omm( omegaSoA->getSoA(ll,mm,box_id) );
        omega_jk_real[index] = omm.real();
        omega_jk_imag[index] = omm.imag();

        complex_type opp(B->get_vectorized(ll,mm));
        op_jkml_real[index] = opp.real();
        op_jkml_imag[index] = opp.imag();

        //left
        index = l2+ll-mm;
        omm = cuda_toggle_sign_if_odd((size_t)mm, conj(omm));
        omega_jk_real[index] = omm.real();
        omega_jk_imag[index] = omm.imag();

        opp = cuda_toggle_sign_if_odd((size_t)mm, conj(opp));
        op_jkml_real[index] = opp.real();
        op_jkml_imag[index] = opp.imag();

        if(mm==0)
        {
            //self
            index = l2+ll;
            complex_type l0(B->get_vectorized(ll,0));
            op_jkml_real[index] = l0.real();
            op_jkml_imag[index] = l0.imag();
            complex_type m0 = omegaSoA->getSoA(ll,0,box_id);
            omega_jk_real[index] = m0.real();
            omega_jk_imag[index] = m0.imag();
        }
    }
    else
    {
        //bottom right
        int lp = ll+p;
        int rbl2 = lp*lp;
        int index = rbl2+lp+mm;
        complex_type opp_br(B->get_vectorized(lp,mm));
        op_jkml_real[index] = opp_br.real();
        op_jkml_imag[index] = opp_br.imag();

        //bottom left
        index = rbl2+lp-mm;
        opp_br = cuda_toggle_sign_if_odd((size_t)mm, conj(opp_br));
        op_jkml_real[index] = opp_br.real();
        op_jkml_imag[index] = opp_br.imag();

        //corner bottom right
        index = rbl2+lp+mm+p;
        complex_type opp_rbm(B->get_vectorized(lp,mm+p));
        op_jkml_real[index] = opp_rbm.real();
        op_jkml_imag[index] = opp_rbm.imag();

        //corner bottom left (never used)
        //index = rbl2+lp-mm-p;
        //opp_rbm = cuda_toggle_sign_if_odd(mm, fmm::conj(opp_rbm));
        //op_jkml_real[index] = opp_rbm.real();
        //op_jkml_imag[index] = opp_rbm.imag();

        if(mm==0)
        {
            //bottom
            index = rbl2+lp;
            complex_type l0 = B->get_vectorized(lp,0);
            op_jkml_real[index] = l0.real();
            op_jkml_imag[index] = l0.imag();
        }
    }
    if(tx < p1xx2 && m > l)
    {
        //bottom right
        int lp = l+p;
        int rbl2 = lp*lp;
        int index = rbl2+lp+m;
        complex_type opp_br(B->get_vectorized(lp,m));
        op_jkml_real[index] = opp_br.real();
        op_jkml_imag[index] = opp_br.imag();

        //bottom left
        index = rbl2+lp-m;
        opp_br = cuda_toggle_sign_if_odd((size_t)m, conj(opp_br));
        op_jkml_real[index] = opp_br.real();
        op_jkml_imag[index] = opp_br.imag();
    }
__syncthreads();


#if 1

    int k_start, k_end;
    complex_type mu_l_m(0.);

    for (int j = 0; j <= p; ++j )
    {
        k_start = f * j;
        k_end = k_start + j + f;

        if(k_start >= k_end)
            continue;

        int jl = j + ll;
        int lk;

        complex_type mu_l_m_j(0.);
        int jl2_jl_m = jl * jl + jl + mm;
        int jj2_j = j * j + j;
        complex_type op,om;
        for (int k = k_start; k < k_end; ++k)
        {
            lk = k - j;
            op.real(op_jkml_real[jl2_jl_m + lk]);
            op.imag(op_jkml_imag[jl2_jl_m + lk]);
            om.real(omega_jk_real[jj2_j + lk]);
            om.imag(omega_jk_imag[jj2_j + lk]);
            mu_l_m_j += op*om;
        }

        mu_l_m += cuda_toggle_sign_if_odd((size_t)j, mu_l_m_j);
    }

    //const size_t lm_index = (ll * (ll + 1)) / 2 + mm;
    //size_t target_box_index = box_id - 1 + op_id * (num_boxes_tree-1);
    //complex_type** mu_ptr = targets_SoA[lm_index];
    //mu_ptr[target_box_index]->atomic_adder(mu_l_m);

    CoefficientMatrix* mu =  box[box_id].b_targets[op_id];
    (*mu)(ll, mm) %= mu_l_m;

#endif

#endif
}

template <typename Box>
__global__ void __set_M2L_sparse_list(Box *box, int2* interactions_list, size_t depth_offset, int* sparse_list_size)
{
    size_t first_box = blockIdx.y * gridDim.x + depth_offset;
    size_t op_id = threadIdx.x;
    size_t box_id = first_box + blockIdx.x;
    size_t target_id = box[box_id].b_target_ids[op_id];

    if(box[box_id].active == 1 && box[target_id].active == 1)
        interactions_list[atomicAdd(&sparse_list_size[0],1)] = make_int2(box_id,op_id);
}

template <typename CoefficientMatrix,typename CoefficientMatrixSoA, typename Box, typename Real, typename Real3, typename complex_type>
__global__ void __M2L_one_p2_SoA_interactions_list(Box *box,
                                                   CoefficientMatrixSoA *omegaSoA,
                                                   complex_type*** targets_SoA,
                                                   const size_t depth_offset,
                                                   const size_t num_boxes_tree,
                                                   const size_t p,
                                                   const size_t p1,
                                                   const size_t p1xx2,
                                                   const size_t op_p1xx2,
                                                   int2* box_interactions_list
                                                   )

{

    extern __shared__ Real real_cache[];
    Real* op_jkml_real  = real_cache;
    Real* op_jkml_imag  = &real_cache[op_p1xx2];
    Real* omega_jk_real = (&real_cache[op_p1xx2*2]);
    Real* omega_jk_imag = (&real_cache[op_p1xx2*2 + p1xx2]);

    size_t box_id        = box_interactions_list[blockIdx.x].x;
    size_t op_id         = box_interactions_list[blockIdx.x].y;

    CoefficientMatrix* B =  box[box_id].b_operators[op_id];

    const ssize_t l = threadIdx.y;
    const ssize_t m = threadIdx.x;
    const ssize_t tx = l*p1 + m;

    int ll, mm, l1, t1, f, f1;

    t1 = tx - l;
    l1 = t1/p1;
    f = 1 -(l - l1);
    f1 = 1 - f;
    ll = f1 * l1 + f * m;
    mm = f1 * m  + f * l1;

    if(f1)
    {
        //self
        int l2 = ll*ll;
        int index = l2+ll+mm;
        complex_type omm( omegaSoA->getSoA(ll,mm,box_id) );
        omega_jk_real[index] = omm.real();
        omega_jk_imag[index] = omm.imag();

        complex_type opp(B->get_vectorized(ll,mm));
        op_jkml_real[index] = opp.real();
        op_jkml_imag[index] = opp.imag();

        //left
        index = l2+ll-mm;
        omm = cuda_toggle_sign_if_odd((size_t)mm, conj(omm));
        omega_jk_real[index] = omm.real();
        omega_jk_imag[index] = omm.imag();

        opp = cuda_toggle_sign_if_odd((size_t)mm, conj(opp));
        op_jkml_real[index] = opp.real();
        op_jkml_imag[index] = opp.imag();

        if(mm==0)
        {
            //self
            index = l2+ll;
            complex_type l0(B->get_vectorized(ll,0));
            op_jkml_real[index] = l0.real();
            op_jkml_imag[index] = l0.imag();
            complex_type m0 = omegaSoA->getSoA(ll,0,box_id);
            omega_jk_real[index] = m0.real();
            omega_jk_imag[index] = m0.imag();
        }
    }
    else
    {
        //bottom right
        int lp = ll+p;
        int rbl2 = lp*lp;
        int index = rbl2+lp+mm;
        complex_type opp_br(B->get_vectorized(lp,mm));
        op_jkml_real[index] = opp_br.real();
        op_jkml_imag[index] = opp_br.imag();

        //bottom left
        index = rbl2+lp-mm;
        opp_br = cuda_toggle_sign_if_odd((size_t)mm, conj(opp_br));
        op_jkml_real[index] = opp_br.real();
        op_jkml_imag[index] = opp_br.imag();

        //corner bottom right
        index = rbl2+lp+mm+p;
        complex_type opp_rbm(B->get_vectorized(lp,mm+p));
        op_jkml_real[index] = opp_rbm.real();
        op_jkml_imag[index] = opp_rbm.imag();

        //corner bottom left (never used)
        //index = rbl2+lp-mm-p;
        //opp_rbm = cuda_toggle_sign_if_odd(mm, fmm::conj(opp_rbm));
        //op_jkml_real[index] = opp_rbm.real();
        //op_jkml_imag[index] = opp_rbm.imag();

        if(mm==0)
        {
            //bottom
            index = rbl2+lp;
            complex_type l0 = B->get_vectorized(lp,0);
            op_jkml_real[index] = l0.real();
            op_jkml_imag[index] = l0.imag();
        }
    }
    if(tx < p1xx2 && m > l)
    {
        //bottom right
        int lp = l+p;
        int rbl2 = lp*lp;
        int index = rbl2+lp+m;
        complex_type opp_br(B->get_vectorized(lp,m));
        op_jkml_real[index] = opp_br.real();
        op_jkml_imag[index] = opp_br.imag();

        //bottom left
        index = rbl2+lp-m;
        opp_br = cuda_toggle_sign_if_odd((size_t)m, conj(opp_br));
        op_jkml_real[index] = opp_br.real();
        op_jkml_imag[index] = opp_br.imag();
    }
    __syncthreads();

    int k_start, k_end;
    complex_type mu_l_m(0.);

    for (int j = 0; j <= p; ++j )
    {
        k_start = f * j;
        k_end = k_start + j + f;

        if(k_start >= k_end)
            continue;

        int jl = j + ll;
        int lk;

        complex_type mu_l_m_j(0.);
        int jl2_jl_m = jl * jl + jl + mm;
        int jj2_j = j * j + j;
        complex_type op,om;
        for (int k = k_start; k < k_end; ++k)
        {
            lk = k - j;
            op.real(op_jkml_real[jl2_jl_m + lk]);
            op.imag(op_jkml_imag[jl2_jl_m + lk]);
            om.real(omega_jk_real[jj2_j + lk]);
            om.imag(omega_jk_imag[jj2_j + lk]);
            mu_l_m_j += op*om;
        }

        mu_l_m += cuda_toggle_sign_if_odd((size_t)j, mu_l_m_j);
    }
    const size_t lm_index = (ll * (ll + 1)) / 2 + mm;
    size_t target_box_index = box_id - 1 + op_id * (num_boxes_tree-1);
    complex_type** mu_ptr = targets_SoA[lm_index];
    mu_ptr[target_box_index]->atomic_adder(mu_l_m);
}

CUDA
__forceinline__
unsigned int bitset_offset(unsigned int pos)
{
    return ( pos + (( pos >> 31) & (1<<5 + ~0))) >> 5;
}

CUDA
__forceinline__
unsigned int get_sign(unsigned int pos, unsigned int &bits)
{
    unsigned int sign = bits<<(pos&31);
    return (sign & 0x80000000);

    //unsigned int sign =  bits << (pos);
    //unsigned int sign =  bits & pos;
    //return (0x80000000);
}

CUDA
__forceinline__
uint2 get_signs(unsigned int pos, unsigned int &bits)
{
    unsigned int sign1 = bits<<(pos&31);
    return make_uint2(sign1 & 0x80000000, (sign1<<1) & 0x80000000);
}

template <typename T>
__device__
__forceinline__
T cuda_toggle_sign_if_one(unsigned int b, T &x)
{
#ifndef GMX_FMM_DOUBLE
    return __int_as_float( __float_as_int(x)^b);
    //return 1.0;
#else
    return __hiloint2double(__double2hiint(x)^b, __double2loint(x));

#endif
}

template <typename T>
__device__
__forceinline__
T cuda_toggle_sign_and_subtract_if_one(unsigned int &a, unsigned int &b, T &x, T &y)
{
#ifndef GMX_FMM_DOUBLE
    return __int_as_float( __float_as_int(x)^a) - __int_as_float( __float_as_int(y)^b);
#else
    return __hiloint2double(__double2hiint(x)^a, __double2loint(x)) - __hiloint2double(__double2hiint(y)^b, __double2loint(y));

#endif
}

template <typename T>
__device__
__forceinline__
T cuda_toggle_sign_and_add_if_one(unsigned int &a, unsigned int &b, T &x, T &y)
{
#ifndef GMX_FMM_DOUBLE
    return __int_as_float( __float_as_int(x)^a) + __int_as_float( __float_as_int(y)^b);
#else
    return __hiloint2double(__double2hiint(x)^a, __double2loint(x)) + __hiloint2double(__double2hiint(y)^b, __double2loint(y));

#endif
}

CUDA
__forceinline__ size_t _boxes_on_depth(size_t d)
{
    // 8^d
    return size_t(1) << (3 * d);
}
CUDA
__forceinline__ size_t _boxes_above_depth(size_t d)
{
    // sum(0, d-1, 8^i)
    return ((size_t(1) << (3 * d)) - 1 ) / 7;
}

template <typename CoefficientMatrix,typename CoefficientMatrixSoA, typename Box, typename Real, typename Real3, typename complex_type, size_t offset, size_t type, typename inttype>
__global__ void
__M2L_one_shared_v8_all(Box *box,
                         CoefficientMatrixSoA *omegaSoA,
                         CoefficientMatrixSoA *muSoA,
                         complex_type*** targets_SoA,
                         size_t child_box_id,
                         const size_t num_boxes_above,
                         const size_t num_boxes_tree,
                         const size_t depth_offset,
                         const size_t p,
                         const size_t p1,
                         const size_t p1xx2,
                         const size_t pxp1_2,
                         const size_t p1xp2_2,
                         const unsigned int num_of_bitset_elems)
{

    typedef typename Box::M2L_Operator Operator;
    extern __shared__ complex_type s[];
    size_t l = (size_t)index_lm_map[blockIdx.x].x;
    size_t m = (size_t)index_lm_map[blockIdx.x].y;
    size_t op_id = offset + blockIdx.y*type;

    complex_type* right = s;
    complex_type* left  = &right[p1xp2_2];

    size_t box_id = child_box_id * num_boxes_above + depth_offset;

    size_t thx = blockIdx.z * blockDim.x + threadIdx.x;

    Operator **B_ptr = get_operator_ptr_ptr<Box,Operator>(box, box_id);

    Operator *B2S = B_ptr[op_id];

    Operator *B2S1;
    Operator *B2S2;
    Operator *B2S3;
    Operator *B2S4;
    Operator *B2S5;
    Operator *B2S6;
    Operator *B2S7;

    if(type > 1)
    {
        B2S1 = B_ptr[op_id+1];
    }

    if(type > 2)
    {
        B2S2 = B_ptr[op_id+2];
        B2S3 = B_ptr[op_id+3];
    }

    if(type > 4)
    {
        B2S4 = B_ptr[op_id+4];
        B2S5 = B_ptr[op_id+5];
        B2S6 = B_ptr[op_id+6];
        B2S7 = B_ptr[op_id+7];
    }

    unsigned int* bits1;
    unsigned int* bits2;
    unsigned int* bits3;
    unsigned int* bits4;
    unsigned int* bits5;
    unsigned int* bits6;
    unsigned int* bits7;

    if(type > 1)
    {
        bits1 = reinterpret_cast<unsigned int*>(&left[pxp1_2]);
    }

    if(type > 2)
    {
        bits2 = &bits1[num_of_bitset_elems];
        bits3 = &bits2[num_of_bitset_elems];
    }

    if(type > 4)
    {
        bits4 = &bits3[num_of_bitset_elems];
        bits5 = &bits4[num_of_bitset_elems];
        bits6 = &bits5[num_of_bitset_elems];
        bits7 = &bits6[num_of_bitset_elems];
    }

    unsigned int bits;
    const size_t lm_index = (l * (l + 1)) / 2 + m;

    if(threadIdx.x < p1xx2)
    {
        if(threadIdx.x < num_of_bitset_elems)
        {
            if(type > 1)
            {
                bits = B2S->bitset[lm_index].bits[threadIdx.x];
                bits1[threadIdx.x] = B2S1->bitset[lm_index].bits[threadIdx.x]^bits;
            }
            if(type > 2)
            {
                bits2[threadIdx.x] = B2S2->bitset[lm_index].bits[threadIdx.x]^bits;
                bits3[threadIdx.x] = B2S3->bitset[lm_index].bits[threadIdx.x]^bits;
            }
            if(type > 4)
            {
                bits4[threadIdx.x] = B2S4->bitset[lm_index].bits[threadIdx.x]^bits;
                bits5[threadIdx.x] = B2S5->bitset[lm_index].bits[threadIdx.x]^bits;
                bits6[threadIdx.x] = B2S6->bitset[lm_index].bits[threadIdx.x]^bits;
                bits7[threadIdx.x] = B2S7->bitset[lm_index].bits[threadIdx.x]^bits;
            }
        }

        //writes into both right and left shared array
        right[threadIdx.x] = B2S->get_vectorized(l+index_lm_map[threadIdx.x].x, m+index_lm_map[threadIdx.x].y);
    }
    __syncthreads();

    if(thx < num_boxes_above)
    {
        complex_type mu_l_m0(0.);
        complex_type mu_l_m1(0.);
        complex_type mu_l_m2(0.);
        complex_type mu_l_m3(0.);
        complex_type mu_l_m4(0.);
        complex_type mu_l_m5(0.);
        complex_type mu_l_m6(0.);
        complex_type mu_l_m7(0.);

        size_t index_r = 0;
        size_t index_l = 0;
        complex_type om;
        complex_type op_r, op_l, val;
        Real a,b,c,d;
        Real ac,bd,ad,bc;

        unsigned int bit_index_r;
        unsigned int bit_index_l;

        //unsigned int a_sgn_change;
        //unsigned int b_sgn_change;
        uint2 sgn_change;
        unsigned int bitsetoffset;

        complex_type mu_l_m_j0;
        complex_type mu_l_m_j1;
        complex_type mu_l_m_j2;
        complex_type mu_l_m_j3;
        complex_type mu_l_m_j4;
        complex_type mu_l_m_j5;
        complex_type mu_l_m_j6;
        complex_type mu_l_m_j7;

        for (size_t j = 0; j <= p; ++j)
        {
            mu_l_m_j0 = 0.0;

            bit_index_r = j*(j+1)*2;
            bit_index_l = bit_index_r;

            om = omegaSoA->get_lin_SoA(index_r, box_id + thx);
            op_r = right[index_r];

            a = op_r.real();
            b = op_r.imag();
            c = om.real();
            d = om.imag();

            ac = a*c;
            //bd = b*d;
            //ad = a*d;
            bc = b*c;

            val.real(ac);
            val.imag(bc);
            mu_l_m_j0 += val;

            if(type > 1)
            {
                mu_l_m_j1 = 0.0;

                bitsetoffset = bitset_offset(bit_index_r);
                bits = bits1[bitsetoffset];
                sgn_change = get_signs(bit_index_r,bits);
                //sgn_change.x = get_sign(bit_index_r,bits);
                //sgn_change.y = get_sign(bit_index_r+1,bits);

                //printf("%u %u\n",a_sgn_change,b_sgn_change );

                val.real(cuda_toggle_sign_if_one(sgn_change.x,ac));// - cuda_toggle_sign_if_one(b_sgn_change,bd));
                val.imag(/*cuda_toggle_sign_if_one(a_sgn_change,ad) + */cuda_toggle_sign_if_one(sgn_change.y,bc));
                mu_l_m_j1 += val;
            }

            if(type > 2)
            {
                mu_l_m_j2 = 0.0;
                mu_l_m_j3 = 0.0;

                bits = bits2[bitsetoffset];
                sgn_change = get_signs(bit_index_r,bits);
                //sgn_change.x = get_sign(bit_index_r,bits);
                //sgn_change.y = get_sign(bit_index_r+1,bits);

                val.real(cuda_toggle_sign_if_one(sgn_change.x,ac));// - cuda_toggle_sign_if_one(b_sgn_change,bd));
                val.imag(/*cuda_toggle_sign_if_one(a_sgn_change,ad) + */cuda_toggle_sign_if_one(sgn_change.y,bc));
                mu_l_m_j2 += val;

                bits = bits3[bitsetoffset];
                sgn_change = get_signs(bit_index_r,bits);
                //sgn_change.x = get_sign(bit_index_r,bits);
                //sgn_change.y = get_sign(bit_index_r+1,bits);

                val.real(cuda_toggle_sign_if_one(sgn_change.x,ac));// - cuda_toggle_sign_if_one(b_sgn_change,bd));
                val.imag(/*cuda_toggle_sign_if_one(a_sgn_change,ad) + */cuda_toggle_sign_if_one(sgn_change.y,bc));
                mu_l_m_j3 += val;
            }

            if(type > 4)
            {
                mu_l_m_j4 = 0.0;
                mu_l_m_j5 = 0.0;
                mu_l_m_j6 = 0.0;
                mu_l_m_j7 = 0.0;

                bits = bits4[bitsetoffset];
                sgn_change = get_signs(bit_index_r,bits);
                //sgn_change.x = get_sign(bit_index_r,bits);
                //sgn_change.y = get_sign(bit_index_r+1,bits);

                val.real(cuda_toggle_sign_if_one(sgn_change.x,ac));// - cuda_toggle_sign_if_one(b_sgn_change,bd));
                val.imag(/*cuda_toggle_sign_if_one(a_sgn_change,ad) + */cuda_toggle_sign_if_one(sgn_change.y,bc));
                mu_l_m_j4 += val;

                bits = bits5[bitsetoffset];
                sgn_change = get_signs(bit_index_r,bits);
                //sgn_change.x = get_sign(bit_index_r,bits);
                //sgn_change.y = get_sign(bit_index_r+1,bits);

                val.real(cuda_toggle_sign_if_one(sgn_change.x,ac));// - cuda_toggle_sign_if_one(b_sgn_change,bd));
                val.imag(/*cuda_toggle_sign_if_one(a_sgn_change,ad) + */cuda_toggle_sign_if_one(sgn_change.y,bc));
                mu_l_m_j5 += val;

                bits = bits6[bitsetoffset];
                sgn_change = get_signs(bit_index_r,bits);
                //sgn_change.x = get_sign(bit_index_r,bits);
                //sgn_change.y = get_sign(bit_index_r+1,bits);

                val.real(cuda_toggle_sign_if_one(sgn_change.x,ac));// - cuda_toggle_sign_if_one(b_sgn_change,bd));
                val.imag(/*cuda_toggle_sign_if_one(a_sgn_change,ad) + */cuda_toggle_sign_if_one(sgn_change.y,bc));
                mu_l_m_j6 += val;

                bits = bits7[bitsetoffset];
                sgn_change = get_signs(bit_index_r,bits);
                //sgn_change.x = get_sign(bit_index_r,bits);
                //sgn_change.y = get_sign(bit_index_r+1,bits);

                val.real(cuda_toggle_sign_if_one(sgn_change.x,ac));// - cuda_toggle_sign_if_one(b_sgn_change,bd));
                val.imag(/*cuda_toggle_sign_if_one(a_sgn_change,ad) + */cuda_toggle_sign_if_one(sgn_change.y,bc));
                mu_l_m_j7 += val;
            }

            ++index_r;

            for (size_t k = 1; k <= j; ++k)
            {
                bit_index_r+=2;
                bit_index_l-=2;

                om  = omegaSoA->get_lin_SoA(index_r, box_id + thx);

                op_r = right[index_r];
                op_l = left[index_l];

                a = op_r.real();
                b = op_r.imag();
                c = om.real();
                d = om.imag();

                ac = a*c;
                bd = b*d;
                ad = a*d;
                bc = b*c;

                complex_type val(ac-bd,ad+bc);

                mu_l_m_j0 += val;

                if(type > 1)
                {
                    bitsetoffset = bitset_offset(bit_index_r);
                    bits = bits1[bitsetoffset];
                    sgn_change = get_signs(bit_index_r,bits);
                    //sgn_change.x = get_sign(bit_index_r,bits);
                    //sgn_change.y = get_sign(bit_index_r+1,bits);

                    val.real(cuda_toggle_sign_and_subtract_if_one(sgn_change.x, sgn_change.y, ac, bd));
                    //val.real(cuda_toggle_sign_if_one(a_sgn_change,ac) - cuda_toggle_sign_if_one(b_sgn_change,bd));
                    val.imag(cuda_toggle_sign_and_add_if_one(sgn_change.x, sgn_change.y, ad, bc));
                    //val.imag(cuda_toggle_sign_if_one(a_sgn_change,ad) + cuda_toggle_sign_if_one(b_sgn_change,bc));
                    mu_l_m_j1 += val;
                }

                if(type > 2)
                {
                    bits = bits2[bitsetoffset];
                    sgn_change = get_signs(bit_index_r,bits);
                    //sgn_change.x = get_sign(bit_index_r,bits);
                    //sgn_change.y = get_sign(bit_index_r+1,bits);

                    val.real(cuda_toggle_sign_and_subtract_if_one(sgn_change.x, sgn_change.y, ac, bd));
                    //val.real(cuda_toggle_sign_if_one(a_sgn_change,ac) - cuda_toggle_sign_if_one(b_sgn_change,bd));
                    val.imag(cuda_toggle_sign_and_add_if_one(sgn_change.x, sgn_change.y, ad, bc));
                    //val.imag(cuda_toggle_sign_if_one(a_sgn_change,ad) + cuda_toggle_sign_if_one(b_sgn_change,bc));
                    mu_l_m_j2 += val;

                    bits = bits3[bitsetoffset];
                    sgn_change = get_signs(bit_index_r,bits);
                    //sgn_change.x = get_sign(bit_index_r,bits);
                    //sgn_change.y = get_sign(bit_index_r+1,bits);

                    val.real(cuda_toggle_sign_and_subtract_if_one(sgn_change.x, sgn_change.y, ac, bd));
                    //val.real(cuda_toggle_sign_if_one(a_sgn_change,ac) - cuda_toggle_sign_if_one(b_sgn_change,bd));
                    val.imag(cuda_toggle_sign_and_add_if_one(sgn_change.x, sgn_change.y, ad, bc));
                    //val.imag(cuda_toggle_sign_if_one(a_sgn_change,ad) + cuda_toggle_sign_if_one(b_sgn_change,bc));
                    mu_l_m_j3 += val;
                }

                if(type > 4)
                {
                    bits = bits4[bitsetoffset];
                    sgn_change = get_signs(bit_index_r,bits);
                    //sgn_change.x = get_sign(bit_index_r,bits);
                    //sgn_change.y = get_sign(bit_index_r+1,bits);

                    val.real(cuda_toggle_sign_and_subtract_if_one(sgn_change.x, sgn_change.y, ac, bd));
                    //val.real(cuda_toggle_sign_if_one(a_sgn_change,ac) - cuda_toggle_sign_if_one(b_sgn_change,bd));
                    val.imag(cuda_toggle_sign_and_add_if_one(sgn_change.x, sgn_change.y, ad, bc));
                    //val.imag(cuda_toggle_sign_if_one(a_sgn_change,ad) + cuda_toggle_sign_if_one(b_sgn_change,bc));
                    mu_l_m_j4 += val;

                    bits = bits5[bitsetoffset];
                    sgn_change = get_signs(bit_index_r,bits);
                    //sgn_change.x = get_sign(bit_index_r,bits);
                    //sgn_change.y = get_sign(bit_index_r+1,bits);

                    val.real(cuda_toggle_sign_and_subtract_if_one(sgn_change.x, sgn_change.y, ac, bd));
                    //val.real(cuda_toggle_sign_if_one(a_sgn_change,ac) - cuda_toggle_sign_if_one(b_sgn_change,bd));
                    val.imag(cuda_toggle_sign_and_add_if_one(sgn_change.x, sgn_change.y, ad, bc));
                    //val.imag(cuda_toggle_sign_if_one(a_sgn_change,ad) + cuda_toggle_sign_if_one(b_sgn_change,bc));
                    mu_l_m_j5 += val;

                    bits = bits6[bitsetoffset];
                    sgn_change = get_signs(bit_index_r,bits);
                    //sgn_change.x = get_sign(bit_index_r,bits);
                    //sgn_change.y = get_sign(bit_index_r+1,bits);

                    val.real(cuda_toggle_sign_and_subtract_if_one(sgn_change.x, sgn_change.y, ac, bd));
                    //val.real(cuda_toggle_sign_if_one(a_sgn_change,ac) - cuda_toggle_sign_if_one(b_sgn_change,bd));
                    val.imag(cuda_toggle_sign_and_add_if_one(sgn_change.x, sgn_change.y, ad, bc));
                    //val.imag(cuda_toggle_sign_if_one(a_sgn_change,ad) + cuda_toggle_sign_if_one(b_sgn_change,bc));
                    mu_l_m_j6 += val;

                    bits = bits7[bitsetoffset];
                    sgn_change = get_signs(bit_index_r,bits);
                    //sgn_change.x = get_sign(bit_index_r,bits);
                    //sgn_change.y = get_sign(bit_index_r+1,bits);

                    val.real(cuda_toggle_sign_and_subtract_if_one(sgn_change.x, sgn_change.y, ac, bd));
                    //val.real(cuda_toggle_sign_if_one(a_sgn_change,ac) - cuda_toggle_sign_if_one(b_sgn_change,bd));
                    val.imag(cuda_toggle_sign_and_add_if_one(sgn_change.x, sgn_change.y, ad, bc));
                    //val.imag(cuda_toggle_sign_if_one(a_sgn_change,ad) + cuda_toggle_sign_if_one(b_sgn_change,bc));
                    mu_l_m_j7 += val;
                }

                om = cuda_toggle_sign_if_odd(k, conj(om));
                a = op_l.real();
                b = op_l.imag();
                c = om.real();
                d = om.imag();

                ac = a*c;
                bd = b*d;
                ad = a*d;
                bc = b*c;

                val.real(ac - bd);
                val.imag(ad + bc);
                mu_l_m_j0 += val;

                if(type > 1)
                {
                    bitsetoffset = bitset_offset(bit_index_l);
                    bits = bits1[bitsetoffset];
                    sgn_change = get_signs(bit_index_l,bits);
                    //sgn_change.x = get_sign(bit_index_l,bits);
                    //sgn_change.y = get_sign(bit_index_l+1,bits);

                    val.real(cuda_toggle_sign_and_subtract_if_one(sgn_change.x, sgn_change.y, ac, bd));
                    //val.real(cuda_toggle_sign_if_one(a_sgn_change,ac) - cuda_toggle_sign_if_one(b_sgn_change,bd));
                    val.imag(cuda_toggle_sign_and_add_if_one(sgn_change.x, sgn_change.y, ad, bc));
                    //val.imag(cuda_toggle_sign_if_one(a_sgn_change,ad) + cuda_toggle_sign_if_one(b_sgn_change,bc));
                    mu_l_m_j1 += val;
                }

                if(type > 2)
                {
                    bits = bits2[bitsetoffset];
                    sgn_change = get_signs(bit_index_l,bits);
                    //sgn_change.x = get_sign(bit_index_l,bits);
                    //sgn_change.y = get_sign(bit_index_l+1,bits);

                    val.real(cuda_toggle_sign_and_subtract_if_one(sgn_change.x, sgn_change.y, ac, bd));
                    //val.real(cuda_toggle_sign_if_one(a_sgn_change,ac) - cuda_toggle_sign_if_one(b_sgn_change,bd));
                    val.imag(cuda_toggle_sign_and_add_if_one(sgn_change.x, sgn_change.y, ad, bc));
                    //val.imag(cuda_toggle_sign_if_one(a_sgn_change,ad) + cuda_toggle_sign_if_one(b_sgn_change,bc));
                    mu_l_m_j2 += val;

                    bits = bits3[bitsetoffset];
                    sgn_change = get_signs(bit_index_l,bits);
                    //sgn_change.x = get_sign(bit_index_l,bits);
                    //sgn_change.y = get_sign(bit_index_l+1,bits);

                    val.real(cuda_toggle_sign_and_subtract_if_one(sgn_change.x, sgn_change.y, ac, bd));
                    //val.real(cuda_toggle_sign_if_one(a_sgn_change,ac) - cuda_toggle_sign_if_one(b_sgn_change,bd));
                    val.imag(cuda_toggle_sign_and_add_if_one(sgn_change.x, sgn_change.y, ad, bc));
                    //val.imag(cuda_toggle_sign_if_one(a_sgn_change,ad) + cuda_toggle_sign_if_one(b_sgn_change,bc));
                    mu_l_m_j3 += val;
                }

                if(type > 4)
                {
                    bits = bits4[bitsetoffset];
                    sgn_change = get_signs(bit_index_l,bits);
                    //sgn_change.x = get_sign(bit_index_l,bits);
                    //sgn_change.y = get_sign(bit_index_l+1,bits);

                    val.real(cuda_toggle_sign_and_subtract_if_one(sgn_change.x, sgn_change.y, ac, bd));
                    //val.real(cuda_toggle_sign_if_one(a_sgn_change,ac) - cuda_toggle_sign_if_one(b_sgn_change,bd));
                    val.imag(cuda_toggle_sign_and_add_if_one(sgn_change.x, sgn_change.y, ad, bc));
                    //val.imag(cuda_toggle_sign_if_one(a_sgn_change,ad) + cuda_toggle_sign_if_one(b_sgn_change,bc));
                    mu_l_m_j4 += val;

                    bits = bits5[bitsetoffset];
                    sgn_change = get_signs(bit_index_l,bits);
                    //sgn_change.x = get_sign(bit_index_l,bits);
                    //sgn_change.y = get_sign(bit_index_l+1,bits);

                    val.real(cuda_toggle_sign_and_subtract_if_one(sgn_change.x, sgn_change.y, ac, bd));
                    //val.real(cuda_toggle_sign_if_one(a_sgn_change,ac) - cuda_toggle_sign_if_one(b_sgn_change,bd));
                    val.imag(cuda_toggle_sign_and_add_if_one(sgn_change.x, sgn_change.y, ad, bc));
                    //val.imag(cuda_toggle_sign_if_one(a_sgn_change,ad) + cuda_toggle_sign_if_one(b_sgn_change,bc));
                    mu_l_m_j5 += val;

                    bits = bits6[bitsetoffset];
                    sgn_change = get_signs(bit_index_l,bits);
                    //sgn_change.x = get_sign(bit_index_l,bits);
                    //sgn_change.y = get_sign(bit_index_l+1,bits);

                    val.real(cuda_toggle_sign_and_subtract_if_one(sgn_change.x, sgn_change.y, ac, bd));
                    //val.real(cuda_toggle_sign_if_one(a_sgn_change,ac) - cuda_toggle_sign_if_one(b_sgn_change,bd));
                    val.imag(cuda_toggle_sign_and_add_if_one(sgn_change.x, sgn_change.y, ad, bc));
                    //val.imag(cuda_toggle_sign_if_one(a_sgn_change,ad) + cuda_toggle_sign_if_one(b_sgn_change,bc));
                    mu_l_m_j6 += val;

                    bits = bits7[bitsetoffset];
                    sgn_change = get_signs(bit_index_l,bits);
                    //sgn_change.x = get_sign(bit_index_l,bits);
                    //sgn_change.y = get_sign(bit_index_l+1,bits);

                    val.real(cuda_toggle_sign_and_subtract_if_one(sgn_change.x, sgn_change.y, ac, bd));
                    //val.real(cuda_toggle_sign_if_one(a_sgn_change,ac) - cuda_toggle_sign_if_one(b_sgn_change,bd));
                    val.imag(cuda_toggle_sign_and_add_if_one(sgn_change.x, sgn_change.y, ad, bc));
                    //val.imag(cuda_toggle_sign_if_one(a_sgn_change,ad) + cuda_toggle_sign_if_one(b_sgn_change,bc));
                    mu_l_m_j7 += val;
                }
                ++index_l;
                ++index_r;
            }

            inttype mask = bitmask(j);
            mu_l_m0 += cuda_toggle_sign_mask(mask, mu_l_m_j0);
            if(type > 1)
            {
                mu_l_m1 += cuda_toggle_sign_mask(mask, mu_l_m_j1);
            }
            if(type > 2)
            {
                mu_l_m2 += cuda_toggle_sign_mask(mask, mu_l_m_j2);
                mu_l_m3 += cuda_toggle_sign_mask(mask, mu_l_m_j3);
            }
            if(type > 4)
            {
                mu_l_m4 += cuda_toggle_sign_mask(mask, mu_l_m_j4);
                mu_l_m5 += cuda_toggle_sign_mask(mask, mu_l_m_j5);
                mu_l_m6 += cuda_toggle_sign_mask(mask, mu_l_m_j6);
                mu_l_m7 += cuda_toggle_sign_mask(mask, mu_l_m_j7);
            }
        }

#ifdef M2L_SOA_OPTIMIZATION
        size_t target_box_index = box_id - 1 + thx + op_id * (num_boxes_tree-1);
        complex_type** mu_ptr = targets_SoA[lm_index];
        mu_ptr[target_box_index]->atomic_adder(mu_l_m0);

        if(type > 1)
        {
            mu_ptr[target_box_index + (num_boxes_tree-1)]->atomic_adder(mu_l_m1);
        }
        if(type > 2)
        {
            mu_ptr[target_box_index + 2*(num_boxes_tree-1)]->atomic_adder(mu_l_m2);
            mu_ptr[target_box_index + 3*(num_boxes_tree-1)]->atomic_adder(mu_l_m3);
        }
        if(type > 4)
        {
            mu_ptr[target_box_index + 4*(num_boxes_tree-1)]->atomic_adder(mu_l_m4);
            mu_ptr[target_box_index + 5*(num_boxes_tree-1)]->atomic_adder(mu_l_m5);
            mu_ptr[target_box_index + 6*(num_boxes_tree-1)]->atomic_adder(mu_l_m6);
            mu_ptr[target_box_index + 7*(num_boxes_tree-1)]->atomic_adder(mu_l_m7);
        }
#else
        size_t target_box_index =  box[box_id+thx].b_target_ids[op_id];
        *(muSoA->get_SoA_ptr(l, m, target_box_index)) %= mu_l_m0;
        if(type > 1)
        {
            target_box_index =  box[box_id+thx].b_target_ids[op_id+1];
            *(muSoA->get_SoA_ptr(l, m, target_box_index)) %= mu_l_m1;
        }
        if(type > 2)
        {
            target_box_index =  box[box_id+thx].b_target_ids[op_id+2];
            *(muSoA->get_SoA_ptr(l, m, target_box_index)) %= mu_l_m2;
            target_box_index =  box[box_id+thx].b_target_ids[op_id+3];
            *(muSoA->get_SoA_ptr(l, m, target_box_index)) %= mu_l_m3;
        }
        if(type > 4)
        {
            target_box_index =  box[box_id+thx].b_target_ids[op_id+4];
            *(muSoA->get_SoA_ptr(l, m, target_box_index)) %= mu_l_m4;
            target_box_index =  box[box_id+thx].b_target_ids[op_id+5];
            *(muSoA->get_SoA_ptr(l, m, target_box_index)) %= mu_l_m5;
            target_box_index =  box[box_id+thx].b_target_ids[op_id+6];
            *(muSoA->get_SoA_ptr(l, m, target_box_index)) %= mu_l_m6;
            target_box_index =  box[box_id+thx].b_target_ids[op_id+7];
            *(muSoA->get_SoA_ptr(l, m, target_box_index)) %= mu_l_m7;
        }
#endif
    }
}

}//namespace end

#endif
