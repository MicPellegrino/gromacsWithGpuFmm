#include "fmm.hpp"
#include "cuda_P2P.hpp"
#include "cuda_lib.hpp"
//#include "vec_RIO.h"
#include <omp.h>


#define FARFARAWAY 1E+16

namespace gmx_gpu_fmm{

size_t periodic_remapped_index(ssize_t i, ssize_t j, ssize_t k, ssize_t ii, ssize_t jj, ssize_t kk, ssize_t dim, size_t depth)
{
    ssize_t iif = ii+i;
    if(iif < 0)
    {
        iif+=dim;
    }
    if(iif >= dim)
    {
        iif-=dim;
    }

    ssize_t jjf = jj+j;
    if(jjf < 0)
    {
        jjf+=dim;
    }
    if(jjf >= dim)
    {
        jjf-=dim;
    }

    ssize_t kkf = kk+k;
    if(kkf < 0)
    {
        kkf+=dim;
    }
    if(kkf >= dim)
    {
        kkf-=dim;
    }

    return make_boxid(iif, jjf, kkf, depth);
}

template<typename Real33, typename Real3>
Real3 periodic_shift(Real33 &abc, ssize_t i, ssize_t j, ssize_t k, ssize_t ii, ssize_t jj, ssize_t kk, ssize_t dim)
{
    Real3 shift(0., 0., 0.);
    ssize_t iif = ii+i;
    while (iif < 0) {
        iif += dim;
        shift -= abc.a;
    }
    while (iif >= dim) {
        iif -= dim;
        shift += abc.a;
    }
    ssize_t jjf = jj+j;
    while (jjf < 0) {
        jjf += dim;
        shift -= abc.b;
    }
    while (jjf >= dim) {
        jjf -= dim;
        shift += abc.b;
    }
    ssize_t kkf = kk+k;
    while(kkf < 0) {
        kkf += dim;
        shift -= abc.c;
    }
    while(kkf >= dim) {
        kkf -= dim;
        shift += abc.c;
    }

    return shift;
}

void fmm_algorithm::prepare_p2p_impl(int version){

    p2p_version = version;
if(version == 1)
{
    Device::custom_alloc(p2p_particles_box_pairs, 13*sizeof(int2*));
    Device::custom_alloc(p2p_particles_periodic_shifts, 13*sizeof(REAL3*));

    Host::custom_alloc(p2p_particles_box_pairs_host, 13*sizeof(int2*));
    Host::custom_alloc(p2p_particles_periodic_shifts_host, 13*sizeof(REAL3*));

    for(int i = 0; i < 13; ++i)
    {
        Device::custom_alloc(p2p_particles_box_pairs[i], num_boxes_lowest*sizeof(int2));
        Device::custom_alloc(p2p_particles_periodic_shifts[i], num_boxes_lowest*sizeof(REAL3));

        Host::custom_alloc(p2p_particles_box_pairs_host[i], num_boxes_lowest*sizeof(int2));
        Host::custom_alloc(p2p_particles_periodic_shifts_host[i], num_boxes_lowest*sizeof(REAL3));
    }

    Device::custom_alloc(p2p_particles, num_boxes_lowest*sizeof(REAL4*));
    Device::custom_alloc(p2p_results, num_boxes_lowest*sizeof(REAL4*));

    Host::custom_alloc(p2p_particles_host, num_boxes_lowest*sizeof(Real4*));
    Host::custom_alloc(p2p_results_host, num_boxes_lowest*sizeof(Real4*));

    Device::custom_alloc(p2p_particles_sizes_div8, num_boxes_lowest*sizeof(size_t));
    Device::custom_alloc(p2p_rounded_particles_sizes, num_boxes_lowest*sizeof(size_t));

    ssize_t dim = ssize_t(1) << depth;
    int index = 0;
    //x dir 1
    for (ssize_t k = 0; k < dim ; ++k)
    {
        for (ssize_t j = 0; j < dim ; ++j)
        {
            for (ssize_t i = -1; i < dim-1 ; ++i)
            {
                Real3 shift = periodic_shift<Real33, Real3>(io->abc,i,j,k,0,0,0,dim);

                p2p_particles_periodic_shifts[0][index].x = shift.x;
                p2p_particles_periodic_shifts[0][index].y = shift.y;
                p2p_particles_periodic_shifts[0][index].z = shift.z;

                p2p_particles_box_pairs[0][index].x = periodic_remapped_index(i,j,k,0,0,0,dim,depth);
                p2p_particles_box_pairs[0][index++].y = periodic_remapped_index(i+1,j,k,0,0,0,dim,depth);

            }
        }
    }
    index = 0;
    //y dir 2
    for (ssize_t k = 0; k < dim ; ++k)
    {
        for (ssize_t i = 0; i < dim ; ++i)
        {
            for (ssize_t j = -1; j < dim-1 ; ++j)
            {
                Real3 shift = periodic_shift<Real33, Real3>(io->abc,i,j,k,0,0,0,dim);

                p2p_particles_periodic_shifts[1][index].x = shift.x;
                p2p_particles_periodic_shifts[1][index].y = shift.y;
                p2p_particles_periodic_shifts[1][index].z = shift.z;

                p2p_particles_box_pairs[1][index].x = periodic_remapped_index(i,j,k,0,0,0,dim,depth);
                p2p_particles_box_pairs[1][index++].y = periodic_remapped_index(i,j+1,k,0,0,0,dim,depth);
            }
        }
    }
    index = 0;
    //xy dir 3
    for (ssize_t k = 0; k < dim ; ++k)
    {
        for (ssize_t i = -1; i < dim-1 ; ++i)
        {
            for (ssize_t j = -1; j < dim-1 ; ++j)
            {
                Real3 shift = periodic_shift<Real33, Real3>(io->abc,i,j,k,0,0,0,dim);

                p2p_particles_periodic_shifts[2][index].x = shift.x;
                p2p_particles_periodic_shifts[2][index].y = shift.y;
                p2p_particles_periodic_shifts[2][index].z = shift.z;

                p2p_particles_box_pairs[2][index].x = periodic_remapped_index(i,j,k,0,0,0,dim,depth);
                p2p_particles_box_pairs[2][index++].y = periodic_remapped_index(i+1,j+1,k,0,0,0,dim,depth);
            }
        }
    }
    index = 0;
    //x-y dir 4
    for (ssize_t k = 0; k < dim ; ++k)
    {
        for (ssize_t i = -1; i < dim-1 ; ++i)
        {
            for (ssize_t j = dim; j > 0 ; --j)
            {
                Real3 shift = periodic_shift<Real33, Real3>(io->abc,i,j,k,0,0,0,dim);

                p2p_particles_periodic_shifts[3][index].x = shift.x;
                p2p_particles_periodic_shifts[3][index].y = shift.y;
                p2p_particles_periodic_shifts[3][index].z = shift.z;

                p2p_particles_box_pairs[3][index].x = periodic_remapped_index(i,j,k,0,0,0,dim,depth);
                p2p_particles_box_pairs[3][index++].y = periodic_remapped_index(i+1,j-1,k,0,0,0,dim,depth);
            }
        }
    }
    index = 0;
    //z dir 5
    for (ssize_t i = 0; i < dim ; ++i)
    {
        for (ssize_t j = 0; j < dim ; ++j)
        {
            for (ssize_t k = -1; k < dim-1 ; ++k)
            {
                Real3 shift = periodic_shift<Real33, Real3>(io->abc,i,j,k,0,0,0,dim);

                p2p_particles_periodic_shifts[4][index].x = shift.x;
                p2p_particles_periodic_shifts[4][index].y = shift.y;
                p2p_particles_periodic_shifts[4][index].z = shift.z;

                p2p_particles_box_pairs[4][index].x = periodic_remapped_index(i,j,k,0,0,0,dim,depth);
                p2p_particles_box_pairs[4][index++].y = periodic_remapped_index(i,j,k+1,0,0,0,dim,depth);
            }
        }
    }
    index = 0;
    //xz dir 6
    for (ssize_t j = 0; j < dim ; ++j)
    {
        for (ssize_t i = -1; i < dim-1 ; ++i)
        {
            for (ssize_t k = -1; k < dim-1 ; ++k)
            {
                Real3 shift = periodic_shift<Real33, Real3>(io->abc,i,j,k,0,0,0,dim);

                p2p_particles_periodic_shifts[5][index].x = shift.x;
                p2p_particles_periodic_shifts[5][index].y = shift.y;
                p2p_particles_periodic_shifts[5][index].z = shift.z;

                p2p_particles_box_pairs[5][index].x = periodic_remapped_index(i,j,k,0,0,0,dim,depth);
                p2p_particles_box_pairs[5][index++].y = periodic_remapped_index(i+1,j,k+1,0,0,0,dim,depth);
            }
        }
    }
    index = 0;
    //x-z dir 7
    for (ssize_t j = 0; j < dim ; ++j)
    {
        for (ssize_t i = -1; i < dim-1 ; ++i)
        {
            for (ssize_t k = dim; k > 0; --k)
            {
                Real3 shift = periodic_shift<Real33, Real3>(io->abc,i,j,k,0,0,0,dim);

                p2p_particles_periodic_shifts[6][index].x = shift.x;
                p2p_particles_periodic_shifts[6][index].y = shift.y;
                p2p_particles_periodic_shifts[6][index].z = shift.z;

                p2p_particles_box_pairs[6][index].x = periodic_remapped_index(i,j,k,0,0,0,dim,depth);
                p2p_particles_box_pairs[6][index++].y = periodic_remapped_index(i+1,j,k-1,0,0,0,dim,depth);
            }
        }
    }
    index = 0;
    //yz dir 8
    for (ssize_t j = -1; j < dim-1 ; ++j)
    {
        for (ssize_t i = 0; i < dim ; ++i)
        {
            for (ssize_t k = -1; k < dim-1; ++k)
            {
                Real3 shift = periodic_shift<Real33, Real3>(io->abc,i,j,k,0,0,0,dim);

                p2p_particles_periodic_shifts[7][index].x = shift.x;
                p2p_particles_periodic_shifts[7][index].y = shift.y;
                p2p_particles_periodic_shifts[7][index].z = shift.z;

                p2p_particles_box_pairs[7][index].x = periodic_remapped_index(i,j,k,0,0,0,dim,depth);
                p2p_particles_box_pairs[7][index++].y = periodic_remapped_index(i,j+1,k+1,0,0,0,dim,depth);
            }
        }
    }
    index = 0;
    //y-z dir 9
    for (ssize_t j = -1; j < dim-1 ; ++j)
    {
        for (ssize_t i = 0; i < dim ; ++i)
        {
            for (ssize_t k = dim; k > 0; --k)
            {
                Real3 shift = periodic_shift<Real33, Real3>(io->abc,i,j,k,0,0,0,dim);

                p2p_particles_periodic_shifts[8][index].x = shift.x;
                p2p_particles_periodic_shifts[8][index].y = shift.y;
                p2p_particles_periodic_shifts[8][index].z = shift.z;

                p2p_particles_box_pairs[8][index].x = periodic_remapped_index(i,j,k,0,0,0,dim,depth);
                p2p_particles_box_pairs[8][index++].y = periodic_remapped_index(i,j+1,k-1,0,0,0,dim,depth);
            }
        }
    }
    index = 0;
    //xyz dir 10
    for (ssize_t j = -1; j < dim-1 ; ++j)
    {
        for (ssize_t i = -1; i < dim-1 ; ++i)
        {
            for (ssize_t k = -1; k < dim-1; ++k)
            {
                Real3 shift = periodic_shift<Real33, Real3>(io->abc,i,j,k,0,0,0,dim);

                p2p_particles_periodic_shifts[9][index].x = shift.x;
                p2p_particles_periodic_shifts[9][index].y = shift.y;
                p2p_particles_periodic_shifts[9][index].z = shift.z;

                p2p_particles_box_pairs[9][index].x = periodic_remapped_index(i,j,k,0,0,0,dim,depth);
                p2p_particles_box_pairs[9][index++].y = periodic_remapped_index(i+1,j+1,k+1,0,0,0,dim,depth);
            }
        }
    }
    index = 0;
    //x-yz dir 11
    for (ssize_t j = dim; j > 0; --j)
    {
        for (ssize_t i = -1; i < dim-1 ; ++i)
        {
            for (ssize_t k = -1; k < dim-1; ++k)
            {
                Real3 shift = periodic_shift<Real33, Real3>(io->abc,i,j,k,0,0,0,dim);

                p2p_particles_periodic_shifts[10][index].x = shift.x;
                p2p_particles_periodic_shifts[10][index].y = shift.y;
                p2p_particles_periodic_shifts[10][index].z = shift.z;

                p2p_particles_box_pairs[10][index].x = periodic_remapped_index(i,j,k,0,0,0,dim,depth);
                p2p_particles_box_pairs[10][index++].y = periodic_remapped_index(i+1,j-1,k+1,0,0,0,dim,depth);
            }
        }
    }
    index = 0;
    //xy-z dir 12
    for (ssize_t j = -1; j < dim-1 ; ++j)
    {
        for (ssize_t i = -1; i < dim-1 ; ++i)
        {
            for (ssize_t k = dim; k > 0; --k)
            {
                Real3 shift = periodic_shift<Real33, Real3>(io->abc,i,j,k,0,0,0,dim);

                p2p_particles_periodic_shifts[11][index].x = shift.x;
                p2p_particles_periodic_shifts[11][index].y = shift.y;
                p2p_particles_periodic_shifts[11][index].z = shift.z;

                p2p_particles_box_pairs[11][index].x = periodic_remapped_index(i,j,k,0,0,0,dim,depth);
                p2p_particles_box_pairs[11][index++].y = periodic_remapped_index(i+1,j+1,k-1,0,0,0,dim,depth);
            }
        }
    }
    index = 0;
    //x-y-z dir 13
    for (ssize_t j = dim; j > 0; --j)
    {
        for (ssize_t i = -1; i < dim-1 ; ++i)
        {
            for (ssize_t k = dim; k > 0; --k)
            {
                Real3 shift = periodic_shift<Real33, Real3>(io->abc,i,j,k,0,0,0,dim);

                p2p_particles_periodic_shifts[12][index].x = shift.x;
                p2p_particles_periodic_shifts[12][index].y = shift.y;
                p2p_particles_periodic_shifts[12][index].z = shift.z;

                p2p_particles_box_pairs[12][index].x = periodic_remapped_index(i,j,k,0,0,0,dim,depth);
                p2p_particles_box_pairs[12][index++].y = periodic_remapped_index(i+1,j-1,k-1,0,0,0,dim,depth);
            }
        }
    }

    //printf("alloc\n");
    for(size_t i = 0; i < num_boxes_lowest; ++i)
    {
        size_t num_of_parts = initial_box_particle_mem_sizes[i];

        if(num_of_parts==0)
        {
            ++num_of_parts;
        }

        size_t rounded_num_of_parts =  ((num_of_parts-1)/64  + 1) * 64;

        //printf("box = %d, particles number = %d, rounded = %d\n", i, initial_box_particle_mem_sizes[i], rounded_num_of_parts);

        Device::custom_alloc(p2p_particles[i], rounded_num_of_parts * sizeof(REAL4));
        Device::custom_alloc(p2p_results[i], rounded_num_of_parts * sizeof(REAL4));

        //CUDA_SAFE_CALL(cudaMemPrefetchAsync(p2p_particles[i], rounded_num_of_parts * sizeof(REAL4), 0,  priority_streams[current_priority_stream]));
        //cudaMemPrefetchAsync(p2p_results[i], rounded_num_of_parts * sizeof(REAL4), 0,  priority_streams[current_priority_stream]);

    }

    for(int pair_index = 0; pair_index < 13; ++pair_index)
    {
        for(int box_index = 0; box_index < num_boxes_lowest; ++box_index)
        {
            p2p_particles_periodic_shifts_host[pair_index][box_index] =  p2p_particles_periodic_shifts[pair_index][box_index];
            p2p_particles_box_pairs_host[pair_index][box_index] =  p2p_particles_box_pairs[pair_index][box_index];
        }
    }

    for(int pair_index = 0; pair_index < 13; ++pair_index)
    {

        //cudaMemPrefetchAsync(&p2p_particles_periodic_shifts[pair_index], num_boxes_lowest*sizeof(int2), 0,  streams[current_stream]);
        //cudaMemPrefetchAsync(&p2p_particles_box_pairs[pair_index], num_boxes_lowest*sizeof(Real3), 0,  streams[current_stream]);
    }

    if(0)
    {
        Device::devSync();
        for (size_t boxid = global_offset; boxid < num_boxes_tree;++boxid)
        {
            for (size_t localid = 0; localid < 27;++localid)
            {
                printf("TargetBoxId %lu-localid %lu, SourceBoxId %lu (%e %e %e)\n",boxid, localid, box[boxid].particle_offset_ids[localid], box[boxid].particle_periodic_shifts[localid].x, box[boxid].particle_periodic_shifts[localid].y, box[boxid].particle_periodic_shifts[localid].z);
            }
        }
    }
}//END P2P VERSION SWITCH

ssize_t dim = ssize_t(1) << depth;
size_t local_index;
for (ssize_t i = 0; i < dim; ++i) {
    for (ssize_t j = 0; j < dim; ++j) {
        for (ssize_t k = 0; k < dim; ++k) {
            //map 3D box id to 1D box id
            size_t idt = make_boxid(i, j, k, depth);
            local_index = 0;

            box[box_id_map[global_offset+idt]].set_offsets(local_index++, idt);

            for (ssize_t ii = i - ws; ii < i + ws + 1; ++ii) {
                for (ssize_t jj = j - ws; jj < j + ws + 1; ++jj) {
                    for (ssize_t kk = k - ws; kk < k + ws + 1; ++kk) {

                        if (ii != i || jj != j || kk != k)
                        {
                            Real3 p_shift = periodic_shift<Real33, Real3>(io->abc,ii,jj,kk,0,0,0,dim);
                            size_t ids = periodic_remapped_index(ii,jj,kk,0,0,0,dim,depth);
                            box[box_id_map[global_offset+idt]].set_offsets(local_index++, ids, -p_shift);
                        }
                    }
                }
            }
        }
    }
}
}

void fmm_algorithm::p2p_host_impl(int host_device, int type)
{

    if(host_device > 0)
    {
    }
    printf("the function requires ordered particles on the host!!!");
    cudaEventSynchronize(copy_particles_to_host_event);
    start_dummy_p2p(1,1);
#if 0

    const int NSIMD = SIMD_BYTES / int(sizeof(Real));
    //printf("SIMDBYTES = %d, NSIMD = %d\n",SIMD_BYTES, NSIMD);
    typedef vec<NSIMD,Real> simdvec;

    typedef typename io_type::host_potential_vector_type potential_type;
    typedef typename io_type::host_field_vector_type force_type;

    int offset = (int)sizeof(Real4);

    for(int pair_index = host_device; pair_index < 13; ++pair_index)
    {
        #pragma omp parallel for
        for(int box_index = 0; box_index < num_boxes_lowest; ++box_index)
        {
            int target_box = p2p_particles_box_pairs_host[pair_index][box_index].x;
            int source_box = p2p_particles_box_pairs_host[pair_index][box_index].y;

            Real3 periodic_shift = Real3(p2p_particles_periodic_shifts_host[pair_index][box_index].x, p2p_particles_periodic_shifts_host[pair_index][box_index].y, p2p_particles_periodic_shifts_host[pair_index][box_index].z);

            int first_ptcl_target = box_particle_offset_host[target_box];
            int first_ptcl_source = box_particle_offset_host[source_box];

            int num_of_targets = box_particle_offset_host[target_box + 1]  - first_ptcl_target;
            int num_of_sources = box_particle_offset_host[source_box + 1]  - first_ptcl_source;

            std::pair<int,int> target_source[2];
            target_source[0] = std::pair<int,int>(first_ptcl_target, first_ptcl_source);
            target_source[1] = std::pair<int,int>(first_ptcl_source, first_ptcl_target);

            potential_type& output_potential = io->potential_host_p2p;
            force_type& output_force = io->efield_host_p2p;

            for(int box_to_box = 0; box_to_box < 2; ++box_to_box)
            {
                int i, j;
            #pragma omp parallel for
                for (i=0; i<num_of_targets/NSIMD; ++i)
                {
                    simdvec zero((Real)0);
                    simdvec rlen(zero);
                    simdvec rlen2(zero);
                    simdvec q_j_rlen(zero);
                    simdvec q_j_rlen3(zero);
                    simdvec diffx(zero);
                    simdvec diffy(zero);
                    simdvec diffz(zero);

                    simdvec phi(zero);
                    simdvec Fx(zero);
                    simdvec Fy(zero);
                    simdvec Fz(zero);

                    Real4* target_i = &ordered_particles_host[target_source[box_to_box].first + i * NSIMD];

                    simdvec xi(&target_i->x, offset);
                    simdvec yi(&target_i->y, offset);
                    simdvec zi(&target_i->z, offset);

                    for (j = 0; j < num_of_sources; j++)
                    {
                        Real4 source_j = ordered_particles_host[target_source[box_to_box].second + j];

                        diffx = xi - source_j.x + periodic_shift.x;
                        diffy = yi - source_j.y + periodic_shift.y;
                        diffz = zi - source_j.z + periodic_shift.z;

                        rlen2 = diffx * diffx + diffy * diffy + diffz * diffz;

                        rlen = rsqrt(rlen2);

                        q_j_rlen = rlen * source_j.q;
                        rlen2 = rlen * rlen;
                        q_j_rlen3 = q_j_rlen * rlen2;

                        Fx  += diffx * q_j_rlen3;
                        Fy  += diffy * q_j_rlen3;
                        Fz  += diffz * q_j_rlen3;
                        phi += q_j_rlen;
                    }

                    Real3* out     = &output_force[target_source[box_to_box].first + i * NSIMD];
                    Real*  phi_out = &output_potential[target_source[box_to_box].first + i * NSIMD];
                    for (int k = 0; k < NSIMD; k++)
                    {
                        *(phi_out + k) += phi[k];
                        (out + k) -> x +=  Fx[k];
                        (out + k) -> y +=  Fy[k];
                        (out + k) -> z +=  Fz[k];
                    }
                }
                periodic_shift *= -1.0;
            }

#if 0
            #pragma omp parallel for private(j)

            for (i=0; i<num_of_sources/NSIMD; ++i)
            {
                int offset = (int)sizeof(Real4);
                simdvec zero((Real)0);
                simdvec rlen(zero);
                simdvec rlen2(zero);
                simdvec q_j_rlen(zero);
                simdvec q_j_rlen3(zero);
                simdvec diffx(zero);
                simdvec diffy(zero);
                simdvec diffz(zero);

                simdvec phi(zero);
                simdvec Fx(zero);
                simdvec Fy(zero);
                simdvec Fz(zero);

                Real4* target_i = &ordered_particles_host[first_ptcl_source + i * NSIMD];

                simdvec xi(&target_i->x, offset);
                simdvec yi(&target_i->y, offset);
                simdvec zi(&target_i->z, offset);

                for (j = 0; j < num_of_targets; j++)
                {
                    Real4 source_j = ordered_particles_host[first_ptcl_target  + j];

                    diffx = xi - source_j.x - periodic_shift.x;
                    diffy = yi - source_j.y - periodic_shift.y;
                    diffz = zi - source_j.z - periodic_shift.z;

                    rlen2 = diffx * diffx + diffy * diffy + diffz * diffz;

                    rlen = rsqrt(rlen2);

                    q_j_rlen = rlen * source_j.q;
                    rlen2 = rlen * rlen;
                    q_j_rlen3 = q_j_rlen * rlen2;

                    Fx  += diffx * q_j_rlen3;
                    Fy  += diffy * q_j_rlen3;
                    Fz  += diffz * q_j_rlen3;
                    phi += q_j_rlen;
                }

                Real3* out     = &output_force[first_ptcl_source + i * NSIMD];
                Real*  phi_out = &output_potential[first_ptcl_source + i * NSIMD];
                for (int k = 0; k < NSIMD; k++)
                {
                    *(phi_out + k) += phi[k];
                    (out + k) -> x +=  Fx[k];
                    (out + k) -> y +=  Fy[k];
                    (out + k) -> z +=  Fz[k];
                }
            }
#endif
        }
    }
    start_dummy_p2p(1, 2);
    cudaEventRecord(copy_particles_to_host_event, NULL);

#endif
}

template <typename Real4, typename Real3>
__global__
void __make_p2p_particles(Real4 *ordered_particles, REAL4** p2p_particles, REAL4** p2p_results, size_t* box_particle_offset, int* p2p_rounded_particles_sizes, int* p2p_particles_sizes_div8, Real3* expansion_points, size_t global_id_offset)
{
    size_t box = blockIdx.y * gridDim.z + blockIdx.z;

    size_t num_of_parts = box_particle_offset[box+1] - box_particle_offset[box];

    size_t i = blockIdx.x * blockDim.x + threadIdx.x;

    int rounded_num_of_parts = ((num_of_parts-1)/64 + 1) * 64;

    if( i < num_of_parts)
    {
        Real4 part = ordered_particles[box_particle_offset[box] + i];

        p2p_particles[box][i].x = part.x;
        p2p_particles[box][i].y = part.y;
        p2p_particles[box][i].z = part.z;
        p2p_particles[box][i].w = part.q;
        p2p_results[box][i]     = make_real4(0.,0.,0.,0.);
    }

    else if(i < rounded_num_of_parts)
    {
        //Real3 box_center = expansion_points[box + global_id_offset];
        p2p_particles[box][i] = make_real4(FARFARAWAY, FARFARAWAY, FARFARAWAY, 0.0);
        p2p_results[box][i] = make_real4(0.,0.,0.,0.);
    }
    p2p_particles_sizes_div8[box] = ((num_of_parts-1)/8 + 1);
    p2p_rounded_particles_sizes[box] = rounded_num_of_parts;
}

template <typename Real, typename Real3, typename Real4, typename outputadapter_type, int forceonly = 0>
__global__
void __copy_result(REAL4** p2p_results, size_t* box_particle_offset, outputadapter_type *result)
{
    size_t box = blockIdx.y * gridDim.z + blockIdx.z;

    size_t num_of_parts = box_particle_offset[box+1] - box_particle_offset[box];

    size_t i = blockIdx.x * blockDim.x + threadIdx.x;

    if( i < num_of_parts)
    {
        //vectorize these loads
        if(forceonly == 1)
        {
            Real3 efield = Real3(p2p_results[box][i].x, p2p_results[box][i].y, p2p_results[box][i].z);
            result->atomic_reduce_f(box_particle_offset[box]+i, efield);
        }
        else
        {
            Real potential = p2p_results[box][i].w;
            Real3 efield = Real3(p2p_results[box][i].x, p2p_results[box][i].y, p2p_results[box][i].z);
            result->atomic_reduce_pf(box_particle_offset[box]+i, potential, efield);
        }
    }
}

void fmm_algorithm::p2p_impl_periodic(outputadapter_type *result, bool calc_energies, int host_device, int version){

    size_t num_of_streams = P2P_STREAMS;

    //printf("max_particles_in_box=%d\n",max_particles_in_box);
    dim3 P2P_block(128,1,1);
    int gy = (num_boxes_lowest-1)/16 + 1;
    int gz = gy > 1 ? 16 : num_boxes_lowest;
    dim3 P2P_grid_self((max_particles_in_box-1)/P2P_block.x + 1, gy, gz);
    //waiting for copy particles
    cudaStreamWaitEvent(streams[current_stream], ordered_particles_set_event, 0);

    if(calc_energies)
    {
        __P2P_self_v4<Box, Real, Real3, Real4, outputadapter_type><<<P2P_grid_self, P2P_block, 0, streams[current_stream]>>>
        (&ordered_particles[0], result, &box_particle_offset[0], 0);
    }
    else
    {
        __P2P_self_v4_forceonly<Box, Real, Real3, Real4, outputadapter_type><<<P2P_grid_self, P2P_block, 0, streams[current_stream] >>>
        (&ordered_particles[0], result, &box_particle_offset[0]);
    }

    //__P2P_self_v4_prep<Box, Real, Real3, Real4, outputadapter_type><<<P2P_grid_self, P2P_block, 0, streams[current_stream] >>>
    //(&ordered_particles[0], result, &box_particle_offset[0], 0);

    /*
    dim3 P2P_block(512,1,1);
    dim3 P2P_grid_self((num_boxes_lowest-1)/P2P_block.x + 1, 1, 1);
    __P2P_self_v4_parent<Box, Real, Real3, Real4, outputadapter_type><<<P2P_grid_self, P2P_block, 0, streams[current_stream] >>>
    (&ordered_particles[0], result, &box_particle_offset[0], num_boxes_lowest);
    */

    CUDA_CHECK_ERROR();

    ssize_t dim = ssize_t(1) << depth;
    if(0)
    {
        size_t index = 0;
        for (size_t k = 0; k < dim ; ++k)
        {
            for (size_t j = 0; j < dim ; ++j)
            {
                for (size_t i = 0; i < dim ; ++i)
                {
                    printf("index %d, id %d\n",index++, make_boxid(i, j, k, depth));
                }
            }
        }
    }

    if(0)
    {
        Device::devSync();
        for(size_t i = 0; i < num_boxes_lowest; ++i)
        {
            size_t num_of_parts = box_particle_offset_host[i+1] - box_particle_offset_host[i];

            for(size_t j = 0; j < num_of_parts; ++j)
            {
                Real4 part = ordered_particles_host[box_particle_offset_host[i]+j];
                p2p_particles[i][j].x =  part.x;
                p2p_particles[i][j].y =  part.y;
                p2p_particles[i][j].z =  part.z;
                p2p_particles[i][j].w =  part.q;
                p2p_results[i][j] = make_real4(0.,0.,0.,0.);
                //std::cout<<i<<" - "<<j<<" "<<part<<std::endl;
            }
            int rounded_num_of_parts = ((num_of_parts-1)/64 + 1) * 64;
            for(size_t j = num_of_parts; j < rounded_num_of_parts; ++j)
            {
                p2p_particles[i][j] = make_real4(expansion_points[i+global_offset].x, expansion_points[i+global_offset].y, expansion_points[i+global_offset].z, 0.0 );
                p2p_results[i][j] = make_real4(0.,0.,0.,0.);
                //std::cout<<i<<" - "<<j<<" "<<"("<<p2p_particles[i][j].x<<", "<<p2p_particles[i][j].x<<", "<<p2p_particles[i][j].z<<": "<<p2p_particles[i][j].w<<")"<<std::endl;
            }
            p2p_rounded_particles_sizes[i] = rounded_num_of_parts;
            p2p_particles_sizes_div8[i] = ((num_of_parts-1)/8 + 1);
            printf("box = %d, particles number = %d, rounded = %d  ysize = %d\n", i, num_of_parts, rounded_num_of_parts, p2p_particles_sizes_div8[i]);
        }
    }

    if(0)
    {
        Device::devSync();
        for(size_t b = 0; b < num_boxes_lowest; ++b)
        {
            for(size_t i = box_particle_offset_host[b]; i < box_particle_offset_host[b+1]; ++i)
            {
                Real4 part = ordered_particles_host[i];
                printf("box %d, particle %f %f %f %f\n", b, part.x,part.y,part.z,part.q);
            }
        }
        for(size_t t = 0; t < 13; ++t)
        {
            printf("type nr %d\n",t);
            for(size_t b = 0; b < num_boxes_lowest; ++b)
            {
                int from_box = p2p_particles_box_pairs[t][b].x;
                int to_box = p2p_particles_box_pairs[t][b].y;
                std::cout<<"box index="<<b<<" "<<from_box<<" -> "<<to_box<<std::endl;
                size_t num_of_parts_from = box_particle_offset_host[from_box+1] - box_particle_offset_host[from_box];
                size_t num_of_parts_to = box_particle_offset_host[to_box+1] - box_particle_offset_host[to_box];

                for(size_t i = 0; i < num_of_parts_from; ++i)
                {
                    std::cout<<"("<<p2p_particles[from_box][i].x<<" "<<p2p_particles[from_box][i].y<<" "<<p2p_particles[from_box][i].z<<")";

                    for(size_t j = 0; j < num_of_parts_to; ++j)
                    {
                        std::cout<<"<----------"<<p2p_particles[to_box][j].x<<" "<<p2p_particles[to_box][j].y<<" "<<p2p_particles[to_box][j].z;
                    }
                    std::cout<<std::endl;
                }
            }
        }
    }

    current_stream = (++current_stream)%num_of_streams;

    if(version == 0)
    {
        dim3 P2P_grid_off((max_particles_in_box-1)/P2P_block.x + 1, num_boxes_lowest, 26);

        if(P2P_STREAMS == 1)
        {

            //waiting for copy particles
            cudaStreamWaitEvent(streams[current_stream], ordered_particles_set_event, 0);
            __P2P_off<Box, Real, Real3, Real4, outputadapter_type, false><<<P2P_grid_off, P2P_block, 0, streams[current_stream]>>>
                 (box, &ordered_particles[0], result, global_offset, &box_particle_offset[0]);
            cudaEventRecord(events[current_stream], streams[current_stream]);
        }
        else
        {
            int gy = (num_boxes_lowest-1)/16 + 1;
            int gz = gy > 1 ? 16 : num_boxes_lowest;
            P2P_grid_off.y = gy;
            P2P_grid_off.z = gz;
            for(size_t i = 1; i < 27; ++i)
            {
                //waiting for copy particles
                cudaStreamWaitEvent(streams[current_stream], ordered_particles_set_event, 0);
                __P2P_off_streams<Box, Real, Real3, Real4, outputadapter_type, false><<<P2P_grid_off, P2P_block, 0, streams[current_stream]>>>
                        (box, &ordered_particles[0], result, global_offset, &box_particle_offset[0], i);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;
            }
            for(int i=0; i < num_of_streams; ++i)
            {
                cudaStreamWaitEvent(streams[current_stream], events[i], 0);
            }
            //dummy for syncing
            start_async_dummy_kernel(1,1,streams[current_stream]);
            cudaEventRecord(events[current_stream], streams[current_stream]);
        }
        CUDA_CHECK_ERROR();
    }
    else
    {
        size_t max_rounded_block = ((max_particles_in_box-1)/64 + 1);

        dim3 prepare_block(64);
        int gy = (num_boxes_lowest-1)/16 + 1;
        int gz = gy > 1 ? 16 : num_boxes_lowest;
        dim3 prepare_grid(max_rounded_block, gy, gz);

        int make_particles_event = current_stream;
        //waiting for copy particles

        cudaStreamWaitEvent(streams[current_stream], ordered_particles_set_event, 0);
        __make_p2p_particles<Real4><<<prepare_grid, prepare_block, 0, streams[current_stream]>>>(&ordered_particles[0], p2p_particles, p2p_results, &box_particle_offset[0], p2p_rounded_particles_sizes, p2p_particles_sizes_div8, &expansion_points[0], global_offset);
        cudaEventRecord(events[make_particles_event], streams[make_particles_event]);
        CUDA_CHECK_ERROR();

        if(depth < 3)
        {
            dim3 block(8,8,1);
            dim3 grid(max_rounded_block,num_boxes_lowest,13);

            cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
            if(host_device > 0)
            {
                if(calc_energies)
                {
                    p2p_half_stencilv2<REAL4, REAL3><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                }
                else
                {
                    p2p_half_stencilv2_forcesonly<REAL4, REAL3><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                }
            }
            CUDA_CHECK_ERROR();
            cudaEventRecord(events[current_stream], streams[current_stream]);
        }
        else
        {
            dim3 block(8,8,1);
            int gy = (num_boxes_lowest-1)/16 + 1;
            int gz = gy > 1 ? 16 : num_boxes_lowest;
            dim3 grid(max_rounded_block, gy, gz);

            if(calc_energies)
            {
                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2<REAL4, REAL3, 0><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2<REAL4, REAL3, 1><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2<REAL4, REAL3, 2><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2<REAL4, REAL3, 3><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2<REAL4, REAL3, 4><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2<REAL4, REAL3, 5><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2<REAL4, REAL3, 6><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2<REAL4, REAL3, 7><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2<REAL4, REAL3, 8><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2<REAL4, REAL3, 9><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2<REAL4, REAL3, 10><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2<REAL4, REAL3, 11><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2<REAL4, REAL3, 12><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                //current_stream = (++current_stream)%num_of_streams;

            }
            else
            {
                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2_forcesonly<REAL4, REAL3, 0><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2_forcesonly<REAL4, REAL3, 1><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2_forcesonly<REAL4, REAL3, 2><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2_forcesonly<REAL4, REAL3, 3><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2_forcesonly<REAL4, REAL3, 4><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2_forcesonly<REAL4, REAL3, 5><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2_forcesonly<REAL4, REAL3, 6><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2_forcesonly<REAL4, REAL3, 7><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2_forcesonly<REAL4, REAL3, 8><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2_forcesonly<REAL4, REAL3, 9><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2_forcesonly<REAL4, REAL3, 10><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2_forcesonly<REAL4, REAL3, 11><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                current_stream = (++current_stream)%num_of_streams;

                cudaStreamWaitEvent(streams[current_stream], events[make_particles_event], 0);
                p2p_half_stencilv2_forcesonly<REAL4, REAL3, 12><<<grid,block, 0, streams[current_stream]>>>(p2p_particles, p2p_rounded_particles_sizes, p2p_particles_sizes_div8, p2p_particles_box_pairs, p2p_results, p2p_particles_periodic_shifts);
                cudaEventRecord(events[current_stream], streams[current_stream]);
                //current_stream = (++current_stream)%num_of_streams;
            }

            CUDA_CHECK_ERROR();

            for(int i=0; i < num_of_streams; ++i)
            {
                cudaStreamWaitEvent(streams[current_stream], events[i], 0);
            }
        }
        if(calc_energies)
        {
            __copy_result<Real, Real3, Real4, outputadapter_type,0><<<prepare_grid, prepare_block, 0, streams[current_stream]>>>(p2p_results, &box_particle_offset[0], result);
        }
        else
        {
            __copy_result<Real, Real3, Real4, outputadapter_type,1><<<prepare_grid, prepare_block, 0, streams[current_stream]>>>(p2p_results, &box_particle_offset[0], result);
        }

        cudaEventRecord(events[current_stream], streams[current_stream]);
        CUDA_CHECK_ERROR();
    }

    if(0)
    {
        Device::devSync();
        for(size_t i = 0; i < num_boxes_lowest; ++i)
        {
            size_t num_of_parts = box_particle_offset_host[i+1] - box_particle_offset_host[i];
            for(size_t j = 0; j < num_of_parts; ++j)
            {
                Real potential = p2p_results[i][j].w;
                Real3 efield = Real3(p2p_results[i][j].x, p2p_results[i][j].y, p2p_results[i][j].z);
                result->not_atomic_reduce_pf(box_particle_offset[i]+j, potential, efield);

            }
        }
    }
}

void fmm_algorithm::p2p_impl_open(outputadapter_type *result){

    const size_t num_of_streams = 1;//P2P_STREAMS;
    dim3 P2P_block(256,1,1);
    dim3 P2P_grid_self((max_particles_in_box-1)/P2P_block.x + 1,num_boxes_lowest, 1);
    dim3 P2P_grid_off((max_particles_in_box-1)/P2P_block.x + 1,num_boxes_lowest, 26);
    //waiting for copy particles
    cudaStreamWaitEvent(streams[current_stream], ordered_particles_set_event, 0);
    __P2P_self_v2<Box, Real, Real3, Real4, outputadapter_type><<<P2P_grid_self, P2P_block, 0, streams[current_stream]>>>
            (box, &ordered_particles[0], result, global_offset, &box_particle_offset[0]);
    cudaEventRecord(events[current_stream], streams[current_stream]);
    CUDA_CHECK_ERROR();

    current_stream = (++current_stream)%num_of_streams;

    if(depth > 0)
    {
        //waiting for copy particles
        cudaStreamWaitEvent(streams[current_stream], ordered_particles_set_event, 0);
        __P2P_off<Box, Real, Real3, Real4, outputadapter_type, true><<<P2P_grid_off, P2P_block, 0, streams[current_stream]>>>
                (box, &ordered_particles[0], result, global_offset, &box_particle_offset[0]);
        cudaEventRecord(events[current_stream], streams[current_stream]);
        CUDA_CHECK_ERROR();
    }


    for(int i = 0; i < num_of_streams; ++i)
    {
        cudaStreamWaitEvent(streams[current_stream], events[i], 0);
    }
    //dummy for syncing
    start_async_dummy_kernel(1,1,streams[current_stream]);
    cudaEventRecord(events[current_stream], streams[current_stream]);
    CUDA_CHECK_ERROR();
}

#if 0
__global__
void
print_lj(REAL2* lj_comb, int n)
{
    for(int i = 0; i < n; ++i)
    {
        printf("%d %e %e\n", i, lj_comb[i].x, lj_comb[i].y);
    }
}

void fmm_algorithm::test_lj()
{
    print_lj<<<1,1>>>(lj_comb, io->n);
}

#endif


}//namespace end


