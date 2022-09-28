#ifndef _BK_FMM_cuda_P2P_hpp
#define _BK_FMM_cuda_P2P_hpp

#include "xyz.hpp"
#include "floatdouble234.hpp"
#include "cuda_lib.hpp"

#define FULL_MASK 0xffffffff

namespace gmx_gpu_fmm{

namespace p2p{

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

template <typename Real, typename Real3>
DEVICE
__forceinline__
Real __rcplength(const Real3& a)
{
    return __rsqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

template <typename Real, typename Real3>
DEVICE __forceinline__
void
one_coulomb(Real3& x_i, Real3& x_j, Real& q_j, Real3 &efield, Real &potential)
{
    //coulomb
    Real3 diff = x_i - x_j;
    Real q_j_rlen = __rcplength<Real, Real3>(diff);
    diff  *= q_j_rlen * q_j_rlen;
    q_j_rlen *= q_j;

    potential += q_j_rlen;
    efield += diff * q_j_rlen;
}

template <typename Real, typename Real3>
CUDA __forceinline__
void
one_coulomb_ref(Real3& x_i, Real3& x_j, Real& q_j, Real3 &efield, Real &potential)
{
    //coulomb
    Real3 diff = x_i - x_j;                          // (x_i - x_j)
    Real rlen = __rcplength<Real, Real3>(diff);            // 1 / |x_i - x_j|
    Real q_j_rlen = q_j * rlen;                      // q_j / |x_i - x_j|
    Real rlen2 = rlen * rlen;                        // 1 / |x_i - x_j|^2
    Real q_j_rlen3 = q_j_rlen * rlen2;               // q_j / |x_i - x_j|^3
    efield += diff * q_j_rlen3;                      // (x_i - x_j) * q_j / |x_i - x_j|^3
    potential += q_j_rlen;
}

template <typename Real, typename Real3>
DEVICE __forceinline__
void
one_coulomb_forceonly(Real3& x_i, Real3& x_j, Real& q_j, Real3 &efield)
{
    Real3 diff = x_i - x_j;
    Real rlen = __rcplength<Real, Real3>(diff);
    efield   += diff * rlen * rlen * rlen * q_j;
}

}

template <typename Box, typename Real, typename Real3, typename Real4, typename outputadapter_type>
__global__
void
__P2P_self(Box* box, Real4* ordered_particles, outputadapter_type* result, size_t global_offset, size_t* box_particle_offset)
{
    size_t target_box_id = global_offset + blockIdx.y;
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;

    size_t global_p_offset_id = box[target_box_id].particle_offset_ids[0];
    size_t pi_start = box_particle_offset[global_p_offset_id];
    size_t pi_end = box_particle_offset[global_p_offset_id + 1];

    size_t pi = i + pi_start;

    if(pi >= pi_end)
        return;

    Real3 x_i = Real3(ordered_particles[pi].x, ordered_particles[pi].y, ordered_particles[pi].z );
    Real3 x_j;
    Real q_j;

    Real3 efield(0.,0.,0.);
    Real potential = 0.0;

    size_t pj_end = min((int)pi_end, (int)pi);
    for(size_t pj = pi_start; pj < pj_end ; ++pj)
    {
        q_j = ordered_particles[pj].q;
        x_j = Real3(ordered_particles[pj].x, ordered_particles[pj].y, ordered_particles[pj].z);
        p2p::one_coulomb(x_i, x_j, q_j, efield, potential);
    }

    size_t pj_start = max((int)pi_start, (int)pi+1);
    for(size_t pj = pj_start; pj < pi_end; ++pj)
    {

        q_j = ordered_particles[pj].q;
        x_j = Real3(ordered_particles[pj].x, ordered_particles[pj].y, ordered_particles[pj].z);
        p2p::one_coulomb(x_i, x_j, q_j, efield, potential);
    }

    result->atomic_reduce_pf(pi, potential, efield);
}

template <typename Box, typename Real, typename Real3, typename Real4, typename outputadapter_type>
__global__
void
__P2P_self_v2(Box* box, Real4* ordered_particles, outputadapter_type* result, size_t global_offset, size_t* box_particle_offset)
{
    size_t target_box_id = global_offset + blockIdx.y;
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;

    size_t global_p_offset_id = box[target_box_id].particle_offset_ids[0];
    size_t pi_start = box_particle_offset[global_p_offset_id];
    size_t pi_end = box_particle_offset[global_p_offset_id + 1];

    size_t pi = i + pi_start;


    if(pi >= pi_end)
        return;

    //printf("i = %lu, pi =%lu, pi_start = %lu, pi_end = %lu\n",i,pi,pi_start, pi_end);

    REAL4 tmp = *reinterpret_cast<REAL4*>(&ordered_particles[pi]);
    REAL3 x_i = make_real3(tmp.x,tmp.y,tmp.z);
    Real q_i = tmp.w;
    REAL3 x_j;
    REAL3 efield = make_real3(0.,0.,0.);
    Real potential = 0.0;

    size_t pj_end = min((int)pi_end, (int)pi);
    for(size_t pj = pi_start; pj < pj_end ; ++pj)
    {
        tmp = *reinterpret_cast<REAL4*>(&ordered_particles[pj]);
        x_j = make_real3(tmp.x,tmp.y,tmp.z);
        p2p::one_coulomb_ref(x_i, x_j, tmp.w, efield, potential);
    }

    size_t pj_start = max((int)pi_start, (int)pi+1);
    for(size_t pj = pj_start; pj < pi_end; ++pj)
    {
        tmp = *reinterpret_cast<REAL4*>(&ordered_particles[pj]);
        x_j = make_real3(tmp.x,tmp.y,tmp.z);
        p2p::one_coulomb_ref(x_i, x_j, tmp.w, efield, potential);
    }



    result->atomic_reduce_pf(pi, potential, efield*q_i);
}

template <typename Box, typename Real, typename Real3, typename Real4, typename outputadapter_type>
__global__
void
__P2P_self_v3_to_xi(Box* box, Real4* ordered_particles, outputadapter_type* result, size_t global_offset, size_t* box_particle_offset)
{
    size_t target_box_id = global_offset + blockIdx.y;
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;

    size_t global_p_offset_id = box[target_box_id].particle_offset_ids[0];
    size_t pi_start = box_particle_offset[global_p_offset_id];
    size_t pi_end = box_particle_offset[global_p_offset_id + 1];

    size_t pi = i + pi_start;

    if(pi >= pi_end)
        return;

    REAL4* ordered_particles_ptr = reinterpret_cast<REAL4*>(&ordered_particles[pi]);
    REAL4 tmp = (*ordered_particles_ptr);
    REAL3 x_i = make_real3(tmp.x,tmp.y,tmp.z);
    REAL3 x_j;
    Real q_j;
    Real q_i = tmp.w;

    REAL3 efield = make_real3(0.,0.,0.);
    Real potential = 0.0;

    size_t pj_end = min((int)pi_end, (int)pi);
    for(size_t pj = pi_start; pj < pj_end ; ++pj)
    {
        ordered_particles_ptr = reinterpret_cast<REAL4*>(&ordered_particles[pj]);
        tmp = (*ordered_particles_ptr);
        q_j = tmp.w;
        x_j = make_real3(tmp.x,tmp.y,tmp.z);
        p2p::one_coulomb(x_i, x_j, q_j, efield, potential);
    }

    result->atomic_reduce_pf(pi, potential, efield*q_i);
}

template <typename Box, typename Real, typename Real3, typename Real4, typename outputadapter_type>
__global__
void
__P2P_self_v3_from_xi(Box* box, Real4* ordered_particles, outputadapter_type* result, size_t global_offset, size_t* box_particle_offset)
{
    size_t target_box_id = global_offset + blockIdx.y;
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;

    size_t global_p_offset_id = box[target_box_id].particle_offset_ids[0];
    size_t pi_start = box_particle_offset[global_p_offset_id];
    size_t pi_end = box_particle_offset[global_p_offset_id + 1];

    size_t pi = i + pi_start;

    if(pi >= pi_end)
        return;

    REAL4* ordered_particles_ptr = reinterpret_cast<REAL4*>(&ordered_particles[pi]);
    REAL4 tmp = (*ordered_particles_ptr);
    REAL3 x_i = make_real3(tmp.x,tmp.y,tmp.z);
    REAL3 x_j;
    Real q_j;

    REAL3 efield = make_real3(0.,0.,0.);
    Real potential = 0.0;

    size_t pj_start = max((int)pi_start, (int)pi+1);
    for(size_t pj = pj_start; pj < pi_end; ++pj)
    {

        ordered_particles_ptr = reinterpret_cast<REAL4*>(&ordered_particles[pj]);
        tmp = (*ordered_particles_ptr);
        q_j = tmp.w;
        x_j = make_real3(tmp.x,tmp.y,tmp.z);
        p2p::one_coulomb(x_i, x_j, q_j, efield, potential);
    }

    result->atomic_reduce_pf(pi, potential, efield);
}


template <typename Box, typename Real, typename Real3, typename Real4, typename outputadapter_type>
__global__
void
__P2P_self_v4(Real4* ordered_particles, outputadapter_type* result, size_t* box_particle_offset, size_t id_)
{

    size_t id = blockIdx.y * gridDim.z + blockIdx.z;
    //size_t id = id_;
    size_t pi_start = box_particle_offset[id];
    size_t pi_end =   box_particle_offset[id + 1];

    size_t i = blockIdx.x * blockDim.x + threadIdx.x;

    size_t pi = i + pi_start;

    if(pi >= pi_end)
    {
        return;
    }

    REAL4 tmp = *reinterpret_cast<REAL4*>(&ordered_particles[pi]);
    REAL3 x_i = make_real3(tmp.x,tmp.y,tmp.z);
    REAL3 x_j;
    REAL q_i = tmp.w;

    REAL3 efield = make_real3(0.,0.,0.);
    Real potential = 0.0;

#if 0
    for(size_t pj = pi_start; pj < pi_end ; ++pj)
    {
        if(pi != pj)
        {
            tmp = *reinterpret_cast<REAL4*>(&ordered_particles[pj]);
            x_j = make_real3(tmp.x,tmp.y,tmp.z);
            p2p::one_coulomb(x_i, x_j, tmp.w, efield, potential);
        }
    }
#else

    size_t pj_end = min((int)pi_end, (int)pi);
    for(size_t pj = pi_start; pj < pj_end ; ++pj)
    {
        tmp = *reinterpret_cast<REAL4*>(&ordered_particles[pj]);
        x_j = make_real3(tmp.x,tmp.y,tmp.z);
        p2p::one_coulomb(x_i, x_j, tmp.w, efield, potential);
    }

    size_t pj_start = max((int)pi_start, (int)pi+1);
    for(size_t pj = pj_start; pj < pi_end; ++pj)
    {

        tmp = *reinterpret_cast<REAL4*>(&ordered_particles[pj]);
        x_j = make_real3(tmp.x,tmp.y,tmp.z);
        p2p::one_coulomb(x_i, x_j, tmp.w, efield, potential);
    }
#endif

    result->atomic_reduce_pf(pi, potential, efield*q_i);
}

template <typename Box, typename Real, typename Real3, typename Real4, typename outputadapter_type>
__global__
void
__P2P_self_v4_parent(Real4* ordered_particles, outputadapter_type* result, size_t* box_particle_offset, size_t num_boxes_lowest)
{
    size_t id       = blockIdx.x * blockDim.x + threadIdx.x;
    size_t pi_start = box_particle_offset[id];
    size_t pi_end   = box_particle_offset[id + 1];

    if(id >= num_boxes_lowest)
        return;

    dim3 P2P_block(64,1,1);
    dim3 P2P_grid_self((pi_end - pi_start - 1)/P2P_block.x + 1, 1, 1);
    //waiting for copy particles
    __P2P_self_v4<Box, Real, Real3, Real4, outputadapter_type><<<P2P_grid_self, P2P_block>>>
    (ordered_particles, result, box_particle_offset, id);
}

template <typename Box, typename Real, typename Real3, typename Real4, typename outputadapter_type>
__global__
void
__P2P_self_v4_prep(Real4* ordered_particles, outputadapter_type* result, size_t* box_particle_offset, size_t id_)
{

    size_t id = blockIdx.y;
    size_t pi_start = box_particle_offset[id];
    size_t pi_end =   box_particle_offset[id + 1];

    size_t i = blockIdx.x * blockDim.x + threadIdx.x;

    size_t pi = i + pi_start;

    if(pi >= pi_end)
    {
        return;
    }
}

template <typename Box, typename Real, typename Real3, typename Real4, typename outputadapter_type>
__global__
void
__P2P_self_v4_forceonly(Real4* ordered_particles, outputadapter_type* result, size_t* box_particle_offset)
{

    size_t id = blockIdx.y * gridDim.z + blockIdx.z;
    //size_t id = id_;
    size_t pi_start = box_particle_offset[id];
    size_t pi_end =   box_particle_offset[id + 1];

    size_t i = blockIdx.x * blockDim.x + threadIdx.x;

    size_t pi = i + pi_start;

    if(pi >= pi_end)
    {
        return;
    }

    REAL4 tmp = *reinterpret_cast<REAL4*>(&ordered_particles[pi]);
    REAL3 x_i = make_real3(tmp.x,tmp.y,tmp.z);
    REAL3 x_j;
    REAL q_i = tmp.w;

    REAL3 efield = make_real3(0.,0.,0.);
    //Real potential = 0.0;

#if 0
    for(size_t pj = pi_start; pj < pi_end ; ++pj)
    {
        if(pi != pj)
        {
            tmp = *reinterpret_cast<REAL4*>(&ordered_particles[pj]);
            x_j = make_real3(tmp.x,tmp.y,tmp.z);
            p2p::one_coulomb(x_i, x_j, tmp.w, efield, potential);
        }
    }
#else

    size_t pj_end = min((int)pi_end, (int)pi);
    for(size_t pj = pi_start; pj < pj_end ; ++pj)
    {
        tmp = *reinterpret_cast<REAL4*>(&ordered_particles[pj]);
        x_j = make_real3(tmp.x,tmp.y,tmp.z);
        //p2p::one_coulomb(x_i, x_j, tmp.w, efield, potential);
        p2p::one_coulomb_forceonly(x_i, x_j, tmp.w, efield);
    }

    size_t pj_start = max((int)pi_start, (int)pi+1);
    for(size_t pj = pj_start; pj < pi_end; ++pj)
    {

        tmp = *reinterpret_cast<REAL4*>(&ordered_particles[pj]);
        x_j = make_real3(tmp.x,tmp.y,tmp.z);
        //p2p::one_coulomb(x_i, x_j, tmp.w, efield, potential);
        p2p::one_coulomb_forceonly(x_i, x_j, tmp.w, efield);
    }
#endif

    result->atomic_reduce_f(pi, efield*q_i);
}

template <typename Box, typename Real, typename Real3, typename Real4, typename outputadapter_type, bool open>
__global__
void
__P2P_off(Box* box, Real4* ordered_particles, outputadapter_type* result, size_t global_offset, size_t* box_particle_offset)
{
    typedef REAL3 Real3_;
    size_t source_box_local_index = blockIdx.z + 1;
    size_t target_box_id = global_offset + blockIdx.y;
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;

    if(box[target_box_id].active == 0)
        return;

    size_t global_p_offset_id = box[target_box_id].particle_offset_ids[0];
    size_t pi_start = box_particle_offset[global_p_offset_id];
    size_t pi_end = box_particle_offset[global_p_offset_id + 1];

    size_t pi = i + pi_start;

    if(pi >= pi_end)
        return;

    size_t global_sp_offset_id = box[target_box_id].particle_offset_ids[source_box_local_index];
    size_t pj_start = box_particle_offset[global_sp_offset_id];
    size_t pj_end = box_particle_offset[global_sp_offset_id + 1];

    Real3 periodic_offset = box[target_box_id].particle_periodic_shifts[source_box_local_index];
    Real3_ offset = make_real3(periodic_offset.x, periodic_offset.y, periodic_offset.z);

    if(open)
    {
        if( offset.x != 0.0 || offset.y != 0.0 || offset.z != 0.0)
            return;
    }

    const REAL4* ordered_particles_ptr = reinterpret_cast<const REAL4*>(&ordered_particles[pi]);
    const REAL4 tmp = (*ordered_particles_ptr);

    Real3_ x_i = make_real3(tmp.x,tmp.y,tmp.z) + offset;
    Real3_ x_j;
    Real q_j;
    Real q_i = tmp.w;

    Real3_ efield = make_real3(0.,0.,0.);
    Real potential = 0.0;

    for(size_t pj = pj_start; pj < pj_end; ++pj)
    {
        const REAL4* ordered_particles_ptr = reinterpret_cast<const REAL4*>(&ordered_particles[pj]);
        const REAL4 tmp = (*ordered_particles_ptr);
        q_j = tmp.w;
        x_j = make_real3(tmp.x,tmp.y,tmp.z);
        p2p::one_coulomb(x_i, x_j, q_j, efield, potential);
    }

    result->atomic_reduce_pf(pi, potential, efield*q_i);
}

template <typename Box, typename Real, typename Real3, typename Real4, typename outputadapter_type, bool open>
__global__
void
__P2P_off_streams(Box* box, Real4* ordered_particles, outputadapter_type* result, size_t global_offset, size_t* box_particle_offset, size_t source_box_local_index)
{
    typedef REAL3 Real3_;
    size_t target_box_id = global_offset + blockIdx.y * gridDim.z + blockIdx.z;
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;

    size_t global_p_offset_id = box[target_box_id].particle_offset_ids[0];
    size_t pi_start = box_particle_offset[global_p_offset_id];
    size_t pi_end = box_particle_offset[global_p_offset_id + 1];

    size_t pi = i + pi_start;

    if(pi >= pi_end)
        return;

    size_t global_sp_offset_id = box[target_box_id].particle_offset_ids[source_box_local_index];
    size_t pj_start = box_particle_offset[global_sp_offset_id];
    size_t pj_end = box_particle_offset[global_sp_offset_id + 1];

    Real3 periodic_offset = box[target_box_id].particle_periodic_shifts[source_box_local_index];
    Real3_ offset = make_real3(periodic_offset.x, periodic_offset.y, periodic_offset.z);

    if(open)
    {
        if( offset.x != 0.0 || offset.y != 0.0 || offset.z != 0.0)
            return;
    }

    const REAL4* ordered_particles_ptr = reinterpret_cast<const REAL4*>(&ordered_particles[pi]);
    const REAL4 tmp = (*ordered_particles_ptr);

    Real3_ x_i = make_real3(tmp.x,tmp.y,tmp.z) + offset;
    Real3_ x_j;
    Real q_j;

    Real3_ efield = make_real3(0.,0.,0.);
    Real potential = 0.0;

    for(size_t pj = pj_start; pj < pj_end; ++pj)
    {
        const REAL4* ordered_particles_ptr = reinterpret_cast<const REAL4*>(&ordered_particles[pj]);
        const REAL4 tmp = (*ordered_particles_ptr);
        q_j = tmp.w;
        x_j = make_real3(tmp.x,tmp.y,tmp.z);
        p2p::one_coulomb(x_i, x_j, q_j, efield, potential);
    }

    result->atomic_reduce_pf(pi, potential, efield);
}

template <typename outputadapter_type, typename particle_type, typename Real, typename Real3>
__global__
void
__lambda_P2P_self(outputadapter_type** results,particle_type** targets, particle_type** sources, size_t* targets_range, size_t* sources_range, Real* lambdas)
{
    size_t group = blockIdx.x;
    size_t pi = blockIdx.y * blockDim.y + threadIdx.x;

    if(pi>=targets_range[group])
        return;

    Real3 x_i = Real3(targets[group][pi]);
    Real3 x_j;
    Real q_j;

    Real3 efield(0.,0.,0.);
    Real potential = 0;

    for(size_t pj = 0; pj < sources_range[group]; ++pj)
    {
        if(pi!=pj)
        {
            q_j = sources[group][pj].s;
            x_j = Real3(sources[group][pj]);
            p2p::one_coulomb(x_i, x_j, q_j, efield, potential);
        }
    }
    Real lambda = lambdas[group];
    results[group]->atomic_reduce_pf(pi, potential*lambda, efield*lambda);
}

template <typename outputadapter_type, typename particle_type, typename Real, typename Real3>
__global__
void
__lambda_P2P(outputadapter_type** results,particle_type** targets, particle_type** sources, size_t* targets_range, size_t* sources_range, Real* lambdas, size_t offset)
{
    size_t group = blockIdx.x + offset;
    size_t pi = blockIdx.y * blockDim.y + threadIdx.x;

    if(pi>=targets_range[group])
        return;

    Real3 x_i = Real3(targets[group][pi]);
    Real3 x_j;
    Real q_j;

    Real3 efield(0.,0.,0.);
    Real potential = 0;

    for(size_t pj = 0; pj < sources_range[group]; ++pj)
    {
        q_j = sources[group][pj].s;
        x_j = Real3(sources[group][pj]);
        p2p::one_coulomb(x_i, x_j, q_j, efield, potential);
    }

    Real lambda = lambdas[group];
    results[group]->atomic_reduce_pf(pi, potential*lambda, efield*lambda);
}

static __forceinline__ __device__
void reduce_force_j_warp_shfl(REAL4 f, REAL4* fout, int tidxi, int aidx)
{
    f.x += __shfl_down_sync(FULL_MASK, f.x, 1);
    f.y += __shfl_up_sync  (FULL_MASK, f.y, 1);

    f.z += __shfl_down_sync(FULL_MASK, f.z, 1);
    f.w += __shfl_up_sync  (FULL_MASK, f.w, 1);

    if (tidxi & 1)
    {
        f.x = f.y;
        f.z = f.w;
    }

    f.x += __shfl_down_sync(FULL_MASK, f.x, 2);
    f.z += __shfl_up_sync  (FULL_MASK, f.z, 2);

    if (tidxi & 2)
    {
        f.x = f.z;
    }

    f.x += __shfl_down_sync(FULL_MASK, f.x, 4);

    if (tidxi <  4)
    {
        __atomicAdd((&fout[aidx].x) + tidxi, f.x);
    }
}

static __forceinline__ __device__
void reduce_force_j_warp_shfl(REAL3 f, REAL4* fout, int tidxi, int aidx)
{
    f.x += __shfl_down_sync(FULL_MASK,f.x, 1);
    f.y += __shfl_up_sync  (FULL_MASK,f.y, 1);
    f.z += __shfl_down_sync(FULL_MASK,f.z, 1);

    if (tidxi & 1)
    {
        f.x = f.y;
    }

    f.x += __shfl_down_sync(FULL_MASK,f.x, 2);
    f.z += __shfl_up_sync(FULL_MASK,f.z, 2);

    if (tidxi & 2)
    {
        f.x = f.z;
    }

    f.x += __shfl_down_sync(FULL_MASK,f.x, 4);

    if (tidxi <  3)
    {
        __atomicAdd((&fout[aidx].x) + tidxi, f.x);
    }
}

static __forceinline__ __device__
void reduce_force_i_warp_shfl(REAL4 fin, REAL4 *fout, int tidxj, int aidx)
{
    fin.x += __shfl_down_sync(FULL_MASK, fin.x, 8);
    fin.y += __shfl_up_sync  (FULL_MASK, fin.y, 8);

    fin.z += __shfl_down_sync(FULL_MASK, fin.z, 8);
    fin.w += __shfl_up_sync  (FULL_MASK, fin.w, 8);

    if (tidxj & 1)
    {
        fin.x = fin.y;
        fin.z = fin.w;
    }

    fin.x += __shfl_down_sync(FULL_MASK, fin.x, 16);
    fin.z += __shfl_up_sync(FULL_MASK, fin.z, 16);

    if (tidxj & 2)
    {
        fin.x = fin.z;
    }

    /* Threads 0,1,2 and 4,5,6 increment x,y,z for their warp */
    if ((tidxj & 3) < 4)
    {
        __atomicAdd(&fout[aidx].x + (tidxj & 3), fin.x);
    }
}

static __forceinline__ __device__
void reduce_force_i_warp_shfl(REAL3 fin, REAL4 *fout, int tidxj, int aidx)
{
    fin.x += __shfl_down_sync(FULL_MASK, fin.x, 8);
    fin.y += __shfl_up_sync  (FULL_MASK, fin.y, 8);
    fin.z += __shfl_down_sync(FULL_MASK, fin.z, 8);

    if (tidxj & 1)
    {
        fin.x = fin.y;
    }

    fin.x += __shfl_down_sync(FULL_MASK, fin.x, 16);
    fin.z += __shfl_up_sync(FULL_MASK, fin.z, 16);

    if (tidxj & 2)
    {
        fin.x = fin.z;
    }

    /* Threads 0,1,2 and 4,5,6 increment x,y,z for their warp */
    if ((tidxj & 3) < 3)
    {
        __atomicAdd(&fout[aidx].x + (tidxj & 3), fin.x);
    }
}

template <typename REAL4, typename REAL3, int bdx_z = 100>
__global__
void p2p_half_stencilv2(REAL4** particles, int* p_sizes, int* p_sizes_div8, int2** pairs, REAL4** results, REAL3** periodic_shifts)
{

    size_t blockIdxz, blockIdxy;
    if(bdx_z == 100)
    {
        blockIdxz = blockIdx.z;
        blockIdxy = blockIdx.y;
    }
    else
    {
        blockIdxz = bdx_z;
        blockIdxy = blockIdx.y * gridDim.z + blockIdx.z;
    }

    __shared__ REAL4 xqib[64];

    REAL4 fcj_buf, xqbuf, fci_buf[8];

    REAL3 x_j, diff;

    REAL  q_j, rlen, q_j_rlen, q_i_rlen;

    int target_index = threadIdx.y * 8 + threadIdx.x;
    int ai = blockIdx.x*64 + target_index;

    int2 pair = pairs[blockIdxz][blockIdxy];

    if(ai >= p_sizes[pair.x])
        return;

    REAL4* partsx = particles[pair.x];
    xqib[target_index] = partsx[ai] + make_real4(periodic_shifts[blockIdxz][blockIdxy].x, periodic_shifts[blockIdxz][blockIdxy].y, periodic_shifts[blockIdxz][blockIdxy].z, 0.0);

    size_t y_size = p_sizes_div8[pair.y];

    if(y_size == 0)
    {
        return;
    }

    for (int i = 0; i < 8; i++)
    {
        fci_buf[i] = make_real4(0.0, 0.0, 0.0, 0.0);
    }
    __syncthreads();

    REAL4* partsy = particles[pair.y];

    for(int j4 = 0; j4 < y_size; ++j4)
    {
        int aj = j4 * 8 + threadIdx.y;
        x_j = make_real3(partsy[aj].x, partsy[aj].y, partsy[aj].z);
        q_j = partsy[aj].w;
        fcj_buf = make_real4(0.0, 0.0, 0.0, 0.0);

        for (int i = 0; i < 8; ++i)
        {
            xqbuf = xqib[i * 8 + threadIdx.x];

            diff = make_real3(xqbuf.x, xqbuf.y, xqbuf.z) - x_j;

            rlen = p2p::__rcplength<REAL,REAL3>(diff);
            //if(isnan(rlen) || isinf(rlen))
            //printf("(%d %d %d) --  %d<->%d : %d %d (%e %e %e) - (%e %e %e) = (%e %e %e), qi = %f, qj = %f\n",target_index,threadIdx.x, threadIdx.y, pair.x,pair.y,ai,aj, xqbuf.x,xqbuf.y,xqbuf.z ,x_j.x,x_j.y,x_j.z, diff.x, diff.y, diff.z, xqbuf.w, q_j);

            q_j_rlen =  q_j * rlen;       //phi
            q_i_rlen =  xqbuf.w * rlen;   //phi

            fci_buf[i].w +=  q_j_rlen;
            fcj_buf.w    +=  q_i_rlen;

            rlen *= q_i_rlen * q_j_rlen;
            diff *= rlen;

            fci_buf[i] += diff;
            fcj_buf    -= diff;
        }

#ifndef GMX_FMM_DOUBLE
        reduce_force_j_warp_shfl(fcj_buf, results[pair.y], threadIdx.x, aj);
#else
        __atomicAdd(&(results[pair.y][aj].x), fcj_buf.x);
        __atomicAdd(&(results[pair.y][aj].y), fcj_buf.y);
        __atomicAdd(&(results[pair.y][aj].z), fcj_buf.z);
        __atomicAdd(&(results[pair.y][aj].w), fcj_buf.w);
#endif

    }

    for (int i = 0; i < 8; i++)
    {

        ai = blockIdx.x*64 + i*8 + threadIdx.x;

#ifndef GMX_FMM_DOUBLE
        reduce_force_i_warp_shfl(fci_buf[i], results[pair.x], threadIdx.y, ai);
#else
        __atomicAdd(&(results[pair.x][ai].x), fci_buf[i].x);
        __atomicAdd(&(results[pair.x][ai].y), fci_buf[i].y);
        __atomicAdd(&(results[pair.x][ai].z), fci_buf[i].z);
        __atomicAdd(&(results[pair.x][ai].w), fci_buf[i].w);
#endif
    }
}

template <typename REAL4, typename REAL3, int bdx_z = 100>
__global__
void p2p_half_stencilv2_forcesonly(REAL4** particles, int* p_sizes, int* p_sizes_div8, int2** pairs, REAL4** results, REAL3** periodic_shifts)
{
    size_t blockIdxz, blockIdxy;
    if(bdx_z == 100)
    {
        blockIdxz = blockIdx.z;
        blockIdxy = blockIdx.y;
    }
    else
    {
        blockIdxz = bdx_z;
        blockIdxy = blockIdx.y * gridDim.z + blockIdx.z;
    }

    __shared__ REAL4 xqib[64];

    REAL4 xqbuf;
    REAL3 fcj_buf, fci_buf[8];

    REAL3 x_j, diff;

    REAL  q_j, rlen;

    int target_index = threadIdx.y * 8 + threadIdx.x;
    int ai = blockIdx.x*64 + target_index;

    int2 pair = pairs[blockIdxz][blockIdxy];

    if(ai >= p_sizes[pair.x])
        return;

    REAL4* partsx = particles[pair.x];
    xqib[target_index] = partsx[ai] + make_real4(periodic_shifts[blockIdxz][blockIdxy].x, periodic_shifts[blockIdxz][blockIdxy].y, periodic_shifts[blockIdxz][blockIdxy].z, 0.0);

    size_t y_size = p_sizes_div8[pair.y];

    if(y_size == 0)
    {
        return;
    }

    for (int i = 0; i < 8; i++)
    {
        fci_buf[i] = make_real3(0.0, 0.0, 0.0);
    }
    __syncthreads();

    REAL4* partsy = particles[pair.y];

    for(int j4 = 0; j4 < y_size; ++j4)
    {
        int aj = j4 * 8 + threadIdx.y;
        x_j = make_real3(partsy[aj].x, partsy[aj].y, partsy[aj].z);
        q_j = partsy[aj].w;
        fcj_buf = make_real3(0.0, 0.0, 0.0);

        for (int i = 0; i < 8; ++i)
        {
            xqbuf = xqib[i * 8 + threadIdx.x];

            diff = make_real3(xqbuf.x, xqbuf.y, xqbuf.z) - x_j;

            rlen = p2p::__rcplength<REAL,REAL3>(diff);

            diff *= rlen * rlen * rlen * q_j * xqbuf.w;

            fci_buf[i] += diff;
            fcj_buf    -= diff;
        }

#ifndef GMX_FMM_DOUBLE
        reduce_force_j_warp_shfl(fcj_buf, results[pair.y], threadIdx.x, aj);
#else
        __atomicAdd(&(results[pair.y][aj].x), fcj_buf.x);
        __atomicAdd(&(results[pair.y][aj].y), fcj_buf.y);
        __atomicAdd(&(results[pair.y][aj].z), fcj_buf.z);
#endif

    }

    for (int i = 0; i < 8; i++)
    {

        ai = blockIdx.x*64 + i*8 + threadIdx.x;

#ifndef GMX_FMM_DOUBLE
        reduce_force_i_warp_shfl(fci_buf[i], results[pair.x], threadIdx.y, ai);
#else
        __atomicAdd(&(results[pair.x][ai].x), fci_buf[i].x);
        __atomicAdd(&(results[pair.x][ai].y), fci_buf[i].y);
        __atomicAdd(&(results[pair.x][ai].z), fci_buf[i].z);
#endif
    }
}

//particles, x_pairs, output_gmx_v3, SIZE_PER_BOX
template <typename REAL3, typename REAL4, int Splitsize>
__global__
void p2p_half_stencilv2_streams_splitted(REAL4** particles, size_t* p_sizes,size_t* p_sizes_div8, int2** pairs, REAL4** results, REAL3** periodic_shifts, int pair_index = 0)
{
    __shared__ REAL3 xqib[64];
    __shared__ REAL  xqib_q[64];

    REAL3 fcj_buf, fci_buf[Splitsize];
    REAL  phi_buf_i[Splitsize], phi_buf_j;

    REAL  q_j, q_j_rlen, q_i_rlen;
    REAL3 x_j,diff;

    int target_index = threadIdx.y * 8 + threadIdx.x;
    int ai = blockIdx.z*64 + target_index;
    int i, i_x;
    i_x = Splitsize * blockIdx.x;

    int2 pair = pairs[pair_index][blockIdx.y];
    REAL4* partsx =  particles[pair.x];

    if(ai >= p_sizes[pair.x])
        return;

    REAL4* partsy =  particles[pair.y];
    size_t y_size = p_sizes_div8[pair.y];

    xqib[target_index] = make_real3(partsx[ai].x, partsx[ai].y, partsx[ai].z) + periodic_shifts[pair_index][blockIdx.y];
    xqib_q[target_index] = partsx[ai].w;

    for (i = 0; i < Splitsize; i++)
    {
        fci_buf[i] = make_real3(0.0, 0.0, 0.0);
        phi_buf_i[i] = 0.0;
    }

    __syncthreads();

    for (int j4 = 0; j4 < y_size; ++j4)
    {
        int aj = j4 * 8 + threadIdx.y;
        x_j = make_real3(partsy[aj].x, partsy[aj].y, partsy[aj].z);
        q_j = partsy[aj].w;
        fcj_buf = make_real3(0.0, 0.0, 0.0);
        phi_buf_j = 0.0;

        //#pragma unroll
        for (i = 0; i < Splitsize; ++i)
        {
            diff = xqib[(i+i_x) * 8 + threadIdx.x] - x_j;

            q_j_rlen = p2p::__rcplength<REAL,REAL3>(diff);

            diff *= q_j_rlen * q_j_rlen;

            q_i_rlen = q_j_rlen * xqib_q[(i+i_x) * 8 + threadIdx.x];   //phi at x_j
            q_j_rlen = q_j_rlen * q_j;                                 //phi at x_i

            phi_buf_i[i] +=  q_j_rlen;
            phi_buf_j    +=  q_i_rlen;

            fci_buf[i] += diff * q_j_rlen;
            fcj_buf    -= diff * q_i_rlen;
        }

#ifndef GMX_FMM_DOUBLE
        reduce_force_j_warp_shfl( make_real4(fcj_buf.x, fcj_buf.y, fcj_buf.z, phi_buf_j), results[pair.y], threadIdx.x, aj);
#else
        __atomicAdd(&(results[pair.y][aj].x), fcj_buf.x);
        __atomicAdd(&(results[pair.y][aj].y), fcj_buf.y);
        __atomicAdd(&(results[pair.y][aj].z), fcj_buf.z);
        __atomicAdd(&(results[pair.y][aj].w), phi_buf_j);
#endif

    }

    for (i = 0; i < Splitsize; i++)
    {

        ai = blockIdx.z*64 + (i+i_x)*8 + threadIdx.x;

#ifndef GMX_FMM_DOUBLE
        reduce_force_i_warp_shfl(make_real4(fci_buf[i].x, fci_buf[i].y, fci_buf[i].z, phi_buf_i[i]), results[pair.x], threadIdx.y, ai);
#else
        __atomicAdd(&(results[pair.x][ai].x), fci_buf[i].x);
        __atomicAdd(&(results[pair.x][ai].y), fci_buf[i].y);
        __atomicAdd(&(results[pair.x][ai].z), fci_buf[i].z);
        __atomicAdd(&(results[pair.x][ai].w), phi_buf_i[i]);
#endif
    }
}


}//namespace end

#endif
