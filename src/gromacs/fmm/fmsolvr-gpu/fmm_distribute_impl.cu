#include "fmm.hpp"
#include "cuda_DISTRIBUTE.hpp"
#include "cuda_lib.hpp"
#include "cuda_PREPARE_DATA.hpp"

namespace gmx_gpu_fmm{

void fmm_algorithm::gmx_copy_particles_buffer_impl(REAL* gmx_particles, cudaEvent_t* gmx_h2d_ready_event)
{
    size_t num_particles = io->n;
    dim3 block(512,1,1);
    dim3 grid((num_particles - 1)/block.x + 1,1,1);
    cudaStreamWaitEvent(priority_streams[current_priority_stream], *gmx_h2d_ready_event, 0);
    __gmx2fmm_particle_copy<REAL3, Real4, Real><<<grid, block, 0, priority_streams[current_priority_stream]>>>(gmx_particles, &io->unordered_particles[0], io->box_scale, io->n);
}

// its not needed every step, only if particles change positions to different box after integration step
void fmm_algorithm::distribute_particles_impl(bool gmx_does_neighbor_search){

    CUDA_CHECK_ERROR();
    __flush_result<outputadapter_type, Real3><<<(io->n-1)/512+1, 512, 0, priority_streams[current_priority_stream]>>>(io->result_ptr, io->n);
    CUDA_CHECK_ERROR();

    __zero_box_index<<<(num_boxes_tree-1)/256 + 1, 256, 0, priority_streams[current_priority_stream]>>>(box, &particles_per_box_int[0], &box_particle_offset_int[0], num_boxes_tree, num_boxes_lowest);
    cudaEventRecord(priority_events[current_priority_stream], priority_streams[current_priority_stream]);
    CUDA_CHECK_ERROR();

    size_t zero_box_index_event = current_priority_stream;

    //needed for p2m and forces grid dim sizes
    n_blocks = 0;
    bs = 32;
    size_t average_bs = io->n/num_boxes_lowest;

    while(bs < average_bs)
    {
        bs +=32;
    }

    if(bs>128)
        bs = 128;

    Real33 convert_to_abc = change_of_basis_from_standard_basis(io->abc);
    Real scale = 1 << depth;

    dim3 block(512,1,1);
    size_t num_of_streams = gmx_does_neighbor_search ? STREAMS : 1;
    if(io->n < 512)
        num_of_streams = 1;
    size_t particles_tile = io->n/num_of_streams;
    size_t last_tile = io->n%particles_tile;
    size_t particles_offset;

    dim3 grid((particles_tile-1)/block.x+1,1,1);

    if(0)
    {
        cudaDeviceSynchronize();
        for(size_t i = 0; i < io->unordered_particles.size(); ++i)
        {
            std::cout<<i<<"---"<<io->unordered_particles[i]<<std::endl;
        }
    }
    //CUDA_CHECK_ERROR();
    //CUDA_SAFE_CALL(cudaSetDevice(0));
    //int dst = -10;
    //size_t map_size = box_id_map.size();
    //CUDA_SAFE_CALL(cudaMemPrefetchAsync(&box_id_map[0], map_size * sizeof(size_t),dst));

    for(size_t i = 0; i < num_of_streams; ++i)
    {
        current_priority_stream = (++current_priority_stream)%num_of_streams;

        cudaStreamWaitEvent(priority_streams[current_priority_stream], priority_events[zero_box_index_event], 0);
        particles_offset = i*particles_tile;
        if(gmx_does_neighbor_search)
            cudaMemcpyAsync(&io->unordered_particles[particles_offset], &io->unordered_particles_host[particles_offset], particles_tile*sizeof(Real4), cudaMemcpyHostToDevice, priority_streams[current_priority_stream]);
        __distribute_particles<<<grid,block,0,priority_streams[current_priority_stream]>>>
        (io->abc, convert_to_abc, io->reference_corner, &io->unordered_particles[0], &particles_per_box_int[0], depth, global_offset, &box_id_map[0], box, scale, particles_tile, particles_offset);
        cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
        CUDA_CHECK_ERROR();
    }
    if(last_tile != 0)
    {
        current_priority_stream = (++current_priority_stream)%num_of_streams;
        cudaStreamWaitEvent(priority_streams[current_priority_stream], priority_events[zero_box_index_event], 0);
        particles_offset = num_of_streams*particles_tile;
        particles_tile = last_tile;
        if(gmx_does_neighbor_search)
            cudaMemcpyAsync(&io->unordered_particles[particles_offset], &io->unordered_particles_host[particles_offset], particles_tile*sizeof(Real4), cudaMemcpyHostToDevice, priority_streams[current_priority_stream]);
        __distribute_particles<<<grid,block,0,priority_streams[current_priority_stream]>>>
        (io->abc, convert_to_abc, io->reference_corner, &io->unordered_particles[0], &particles_per_box_int[0], depth, global_offset, &box_id_map[0], box, scale, particles_tile, particles_offset);
        cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
        CUDA_CHECK_ERROR();
    }

    dim3 OFFSET_block(256,1,1);
    int gy = (num_boxes_lowest-1)/16 + 1;
    int gz = gy > 1 ? 16 : num_boxes_lowest;
    dim3 OFFSET_grid(num_boxes_lowest/OFFSET_block.x+1,gy,gz);

    //wait till all particles are distributed
    for(size_t i = 0; i < num_of_streams; ++i)
        cudaStreamWaitEvent(priority_streams[current_priority_stream], priority_events[i], 0);

    __compute_offset<<<OFFSET_grid, OFFSET_block, 0, priority_streams[current_priority_stream]>>>(&particles_per_box_int[0], &box_particle_offset_int[0], num_boxes_lowest);
    CUDA_CHECK_ERROR();

    OFFSET_grid.y = 1;
    OFFSET_grid.z = 1;
    __recast_offset<<<OFFSET_grid, OFFSET_block, 0, priority_streams[current_priority_stream]>>>(&box_particle_offset[0], &box_particle_offset_int[0], num_boxes_lowest);
    cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
    CUDA_CHECK_ERROR();

    //only for checking not needed in production
    //__flush_particles<Real4><<<(io->n-1)/512+1, 512, 0, priority_streams[current_priority_stream]>>>(&ordered_particles[0], &orig_ids[0], io->n);
    //CUDA_CHECK_ERROR();

    if(0)
    {
        cudaDeviceSynchronize();
        for(size_t i = 0; i< num_boxes_tree;++i)
        {
            printf ("boxid %d mapped %d\n ", i, box_id_map[i]);
        }
    }

    if(0)
    {
        cudaDeviceSynchronize();
        size_t dim = size_t(1) << depth;
        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j < dim; j++)
            {
                for (size_t k = 0; k < dim; k++)
                {
                    size_t box_id = global_offset + make_boxid(i, j, k, depth);
                    printf ("boxid %lu num of particles %d ", box_id, box[ box_id_map[box_id] ].ptcl_index,particles_per_box_int[box_id - global_offset]);
                    printf (" active = %d\n", box[ box_id_map[box_id] ].active);
                    //for (int particle_index = 0; particle_index < box[ box_id_map[box_id] ].ptcl_index; particle_index++)
                    {
                        //printf ("boxid %lu particle %d ", box_id, particle_index);
                        //std::cout<<"original index "<<box[box_id_map[box_id]].orig_ptcl_ids[particle_index]<<std::endl;
                    }
                }
            }
        }
    }

}

}//namespace end

