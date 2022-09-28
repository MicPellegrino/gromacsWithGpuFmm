#include "fmm.hpp"
#include "cuda_FE.hpp"

namespace gmx_gpu_fmm{

void fmm_algorithm::forces_impl(outputadapter_type *result, Real one4pieps0){

    size_t sync_mu_derivative = current_priority_stream;
    size_t num_of_streams = 1;
    if(io->n < 512)
        num_of_streams = 1;



#ifndef GMX_FMM_DOUBLE
    dim3 MU_d_block(p1,p1,1);
    dim3 MU_d_grid(num_boxes_lowest,1,1);
    __MU_derivative<CoeffMatrix,Real, Real3, Real4><<<MU_d_grid, MU_d_block, 0, priority_streams[current_priority_stream]>>>(p,global_offset,num_boxes_lowest,mu,dmux,dmuy,dmuz);
#else
    if(p<21)
    {
        dim3 MU_d_block(p1,p1,1);
        dim3 MU_d_grid(num_boxes_lowest,1,1);
        __MU_derivative<CoeffMatrix,Real, Real3, Real4><<<MU_d_grid, MU_d_block, 0, priority_streams[current_priority_stream]>>>(p,global_offset,num_boxes_lowest,mu,dmux,dmuy,dmuz);
    }
    else
    {
        dim3 MU_d_block(64,1,1);
        dim3 MU_d_grid((num_boxes_lowest - 1)/MU_d_block.x + 1,1,1);
        __MU_derivative_ref<CoeffMatrix,Real, Real3, Real4><<<MU_d_grid, MU_d_block, 0, priority_streams[current_priority_stream]>>>(p,global_offset,num_boxes_lowest,mu,dmux,dmuy,dmuz);
    }

#endif
    cudaEventRecord(priority_events[sync_mu_derivative],priority_streams[sync_mu_derivative]);

    dim3 EXCL_block(512,1,1);
    dim3 EXCL_grid((io->n - 1)/EXCL_block.x + 1,1,1);

    //dim3 EXCL_blockv3(512,1,1);
    //dim3 EXCL_gridv3((exclusion_pairs_size - 1)/EXCL_block.x + 1,1,1);

    current_priority_stream = (++current_priority_stream)%num_of_streams;
    size_t sync_exclusions = current_priority_stream;

    map_fmm_ids_to_original<<<EXCL_grid, EXCL_block, 0, priority_streams[current_priority_stream]>>>(&orig_ids[0], &fmm_ids[0], io->n);

    EXCL_grid.x = (io->excl_n - 1)/EXCL_block.x + 1;
    //compute_exclusionsv3<<<EXCL_gridv3, EXCL_blockv3, 0, priority_streams[current_priority_stream]>>>(&io->unordered_particles[0], io->result_ptr, &fmm_ids[0], exclusion_pairs, io->abc, io->half_abc, exclusion_pairs_size);
    compute_exclusionsv2<<<EXCL_grid, EXCL_block, 0, priority_streams[current_priority_stream]>>>(&io->unordered_particles[0], io->result_ptr, &fmm_ids[0], exclusions, exclusions_sizes, io->abc, io->half_abc, io->excl_n);
    //compute_exclusions<<<EXCL_grid, EXCL_block, 27 * sizeof(Real3), priority_streams[current_priority_stream]>>>(&io->unordered_particles[0], io->result_ptr, &orig_ids[0], exclusions, exclusions_sizes, io->abc, io->n);
    cudaEventRecord(priority_events[sync_exclusions],priority_streams[sync_exclusions]);

    size_t block_tile = n_blocks/num_of_streams;
    size_t last_block = n_blocks%num_of_streams;

    //printf("n_blocks = %d, block_tile = %d, bs = %d ",n_blocks, block_tile, bs);
    //printf("last_block = %d\n ", last_block);
    size_t split_blocks = 1;
    dim3 FE_block(bs,split_blocks,1);
    dim3 FE_grid(block_tile/split_blocks,1,1);

    cudaEventSynchronize(priority_events[sync_exclusions]);
    cudaEventSynchronize(priority_events[sync_mu_derivative]);
    //wait for p2p
    cudaEventSynchronize(events[current_stream]);

    size_t last_stream = current_priority_stream;
    for(size_t i = 0; i < num_of_streams; ++i)
    {
        size_t stream_offset = i*block_tile;
        size_t start_block_id = block_map_host[stream_offset];
        size_t end_block_id = block_map_host[block_tile + stream_offset];
        if(i == num_of_streams - 1)
        {
            end_block_id = block_map_host[n_blocks-1] + 1;
            FE_grid.x += last_block;
        }

        //why does it need to be synchronized with its self?????????
        //these are completely independent force computations on different chunks of the forcearray
        cudaStreamWaitEvent(priority_streams[current_priority_stream], priority_events[last_stream], 0);
        __FE<outputadapter_type,CoeffMatrix,CoeffMatrixMemory,Real, Real3, Real4><<<FE_grid,FE_block,0, priority_streams[current_priority_stream]>>>
        (p, global_offset, result, &box_particle_offset[0],mu, &expansion_points[0], dmux,dmuy,dmuz, &ordered_particles[0],&block_map[0],&offset_map[0],stream_offset);
        cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);

        //size_t particle_range = box_particle_offset_host[end_block_id] - box_particle_offset_host[start_block_id];
        //size_t start_particle_id = box_particle_offset_host[start_block_id];
        //printf("start_block = (%d)->%d, end_block = (%d)%d\n", stream_offset, start_block_id, block_tile + stream_offset, end_block_id);
        //CUDA_CHECK_ERROR();
        //cudaMemcpyAsync(&io->potential_host[start_particle_id], &io->potential[start_particle_id], particle_range*sizeof(Real), cudaMemcpyDeviceToHost, priority_streams[current_priority_stream]);
        //CUDA_CHECK_ERROR();
        //cudaMemcpyAsync(&io->efield_host[start_particle_id], &io->efield[start_particle_id], particle_range*sizeof(Real3), cudaMemcpyDeviceToHost, priority_streams[current_priority_stream]);
        //CUDA_CHECK_ERROR();

        last_stream = current_priority_stream;
        current_priority_stream = (++current_priority_stream)%num_of_streams;
    }
    size_t num_particles          = io->n;
    Real box_scale                = io->box_scale;
    Real force_box_scaling_factor = box_scale * box_scale;
    Real box_eps                  = force_box_scaling_factor * one4pieps0;

    dim3 block(512,1,1);
    dim3 grid((num_particles - 1)/block.x + 1,1,1);

    __ordered_to_unordered_forces<Real3, Real><<<grid, block, 0, priority_streams[current_priority_stream]>>>(&io->efield[0], &io->forces_orig_order[0], &orig_ids[0], box_eps, num_particles);
    cudaEventRecord(priority_events[current_priority_stream], priority_streams[current_priority_stream]);
}

void fmm_algorithm::gather_results()
{
    cudaEventSynchronize(copy_particles_to_host_event);
    for(size_t i = 0; i < io->n; ++i)
    {
        io->efield_host[i]    += io->efield_host_p2p[i];
        io->potential_host[i] += io->potential_host_p2p[i];

        io->efield_host_p2p[i] = Real3(0.,0.,0.);
        io->potential_host_p2p[i] = 0.0;
    }
}

void fmm_algorithm::get_forces_in_orig_order_impl(float3* gmx_forces, cudaStream_t *gmx_sync_stream)
{
    dim3 block(512,1,1);
    dim3 grid((io->n - 1)/block.x + 1,1,1);

    cudaStreamWaitEvent(*gmx_sync_stream, priority_events[current_priority_stream], 0);
    __set_fmm_to_gmx_forces<Real3><<<grid, block, 0, *gmx_sync_stream>>>(&io->forces_orig_order[0], gmx_forces, io->n);
}

fmm_algorithm::Real fmm_algorithm::energy_impl(){

    Real Ec_host = 0;
    int block = 128;
    int grid = (io->n - 1)/(block * 32) + 1 ;

    int block_c = 512;
    int grid_c = (io->n - 1)/block_c + 1;

    cudaMemsetAsync(Ec, 0, grid*sizeof(Real),  priority_streams[current_priority_stream]);
    CUDA_CHECK_ERROR();
    __charge_potential<<<grid_c, block_c, 0, priority_streams[current_priority_stream]>>>(&io->potential[0], &ordered_particles[0], io->n);
    CUDA_CHECK_ERROR();
    __energy_kernel<Real,128,32,cub::BLOCK_REDUCE_WARP_REDUCTIONS><<<grid,block,0,priority_streams[current_priority_stream]>>>(&io->potential[0], Ec, io->n);
    CUDA_CHECK_ERROR();
    cudaMemcpyAsync(&io->potential_host[0], Ec, grid*sizeof(Real), cudaMemcpyDeviceToHost, priority_streams[current_priority_stream]);
    cudaStreamSynchronize(priority_streams[current_priority_stream]);
    for(int i = 0; i < grid; ++i)
    {
        Ec_host += io->potential_host[i];
    }

    //printf("ECHOST        %.20e\n",  io->potential_host[0]);
    //printf("ECHOST        %.20e\n", Ec_host*0.5*io->box_scale);
    //printf("ECHOST        %.20e\n", io->box_scale);
    //printf("ECHOST        %.20e\n", Ec_host);


    return Ec_host;
}

void fmm_algorithm::energy_dump_impl(outputadapter_type *result){

    __energy_dump_kernel<outputadapter_type, Real, Real4><<<1,1>>>(result, &ordered_particles[0], io->n, io->box_scale);
}

void fmm_algorithm::force_dump_impl(outputadapter_type *result, int type){

    __force_dump_kernel<outputadapter_type, Real, Real3><<<1,1>>>(result, io->n, type);
}

}//namespace end
