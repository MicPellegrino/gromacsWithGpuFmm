#include "fmm.hpp"
#include "cuda_PREPARE_DATA.hpp"
#include "cuda_lib.hpp"
namespace gmx_gpu_fmm{

// its not needed every step, only if particles change positions to different box after integration step
void fmm_algorithm::prepare_data_impl(){

    cudaMemcpyAsync(&box_particle_offset_host[0], &box_particle_offset[0], (num_boxes_lowest+1)*sizeof(size_t), cudaMemcpyDeviceToHost, priority_streams[current_priority_stream] );
    cudaStreamSynchronize(priority_streams[current_priority_stream]);

    if(0)
    {
        cudaDeviceSynchronize();

        for (ssize_t index = 0; index < (ssize_t)box_particle_offset_host.size(); index++)
        {
            printf ("index %ld offset %lu int_offset %d\n", index, box_particle_offset_host[index], box_particle_offset_int[index]);
        }
    }

    n_blocks = 0;
    size_t min_particle = 0;
    size_t max_particle = 0;
    size_t particle_range = 0;

    size_t block_offset = 0;
    size_t interval = 0;
    size_t index = 0;
    for (size_t id = 0; id <num_boxes_lowest; ++id)
    {
        min_particle = max_particle;
        max_particle = box_particle_offset_host[id+1];
        particle_range = max_particle - min_particle;
        interval = 0;
        max_particles_in_box = std::max(particle_range, max_particles_in_box);
        //printf("box = %d, particles = %d\n", id, particle_range);
        for (size_t i = 0; i < particle_range; i+= bs)
        {
            //cumulative remaining offset within each threadblock
            //printf("index = %d, block offset = %d, id = %d\n",index, block_offset,id);
            offset_map_host[index] = block_offset;
            //which box-id belongs to this threadblock
            block_map_host[index] = id;
            index++;
            n_blocks++;
            interval += bs;

        }
        if(interval > 0)
            block_offset += interval - particle_range;

    }

    cudaMemcpyAsync(&offset_map[0], &offset_map_host[0], offset_map_host.size()*sizeof(size_t), cudaMemcpyHostToDevice, priority_streams[current_priority_stream]);
    cudaMemcpyAsync(&block_map[0], &block_map_host[0], block_map_host.size()*sizeof(size_t), cudaMemcpyHostToDevice, priority_streams[current_priority_stream]);

    if(0)
    {
        printf("n_blocks %lu\n",n_blocks);
        cudaDeviceSynchronize();
        for(size_t id=0; id < offset_map_host.size(); id++)
        {
            printf ("id %lu offset %lu block %lu \n",id, offset_map[id],block_map[id]);
        }
    }

    size_t split_blocks = 1;
    dim3 COPY_block(bs,split_blocks,1);
    dim3 COPY_grid(n_blocks/split_blocks,1,1);
    //printf("bs %d n_blocks %d\n",bs,n_blocks);

    if(0)
    {
        size_t id = 0;
        //for(size_t id = 0; id < box_particle_offset.size();id++)
        {
            size_t counter = 0.0;
            size_t idx = global_offset + id;
            for (size_t p = 0; p <box[box_id_map[idx]].ptcl_index; ++p)
            {
                printf("idx %d orig i %d*\n",idx,box[box_id_map[idx]].orig_ptcl_ids[p]);
                counter++;
            }
            printf("num of particles %d\n",counter);
        }
    }

    //cudaMemPrefetchAsync(&box_particle_offset[0], box_particle_offset.size() * sizeof(size_t), 0,  priority_streams[current_priority_stream]);
    //cudaMemPrefetchAsync(&offset_map[0], offset_map_host.size() * sizeof(size_t), 0,  priority_streams[current_priority_stream]);
    //cudaMemPrefetchAsync(&block_map[0], block_map_host.size() * sizeof(size_t), 0,  priority_streams[current_priority_stream]);

    __COPY_parts<<<COPY_grid, COPY_block, 0, priority_streams[current_priority_stream]>>>
    (global_offset, &box_particle_offset[0], &ordered_particles[0], &io->unordered_particles[0], &orig_ids[0], box, &box_id_map[0], &block_map[0], &offset_map[0]);
    //cudaEventRecord(priority_events[current_priority_stream], priority_streams[current_priority_stream]);
    cudaEventRecord(ordered_particles_set_event, priority_streams[current_priority_stream]);

    typedef typename CoeffMatrix::value_type complex_type;
    complex_type init_val(0.0,0.0);
    dim3 SoA_block(512,1,1);
    dim3 SoA_grid((num_boxes_tree-1)/SoA_block.x+1,p1xp2_2,1);
    __initSoA<<<SoA_grid, SoA_block, 0, priority_streams[current_priority_stream]>>>(muSoA, num_boxes_tree, init_val);

    if(0)
    {
        //check if we need it at all, id map for original particles will do the same thing
        cudaMemcpyAsync(&ordered_particles_host[0], &ordered_particles[0], io->n*sizeof(Real4), cudaMemcpyDeviceToHost, priority_streams[current_priority_stream]);
        cudaEventRecord(copy_particles_to_host_event, priority_streams[current_priority_stream]);
        cudaDeviceSynchronize();
        Real netcharge = 0.0;
        for(size_t i = 0; i < io->n; ++i)
        {
            //std::cout<<i<<"+++"<<io->unordered_particles_host[orig_ids[i]]<<std::endl;
            //std::cout<<i<<"+++"<<ordered_particles_host[i]<<std::endl;
            Real3 a = Real3(io->unordered_particles_host[orig_ids[i]]);
            Real3 b = Real3(ordered_particles_host[i]);
            std::cout<<i<<"+++"<<a/io->box_scale<<std::endl;
            //std::cout<<i<<"+++"<<b/io->box_scale<<std::endl;
            Real qa = io->unordered_particles_host[orig_ids[i]].q;
            Real qb = ordered_particles_host[i].q;
            a-=b;
            qa-=qb;
            //std::cout<<i<<"+++"<<orig_ids[i]<<"++++++++++++++++++++++++++++++++"<<a<<"--------"<<qa<<std::endl;
            netcharge +=io->unordered_particles_host[orig_ids[i]].q;
        }
        std::cout<<"NETCHARGE "<<netcharge<<std::endl;
    }

    if(0)
    {
        cudaDeviceSynchronize();
        for(size_t of = 0; of < box_particle_offset_host.size()-1; ++of)
        {
            for(size_t i = box_particle_offset_host[of]; i < box_particle_offset_host[of+1]; ++i)
            {
                std::cout<<of<<"--"<<i<<"--"<<orig_ids[i]<<"+++"<<io->unordered_particles_host[orig_ids[i]]<<std::endl;
                std::cout<<of<<"--"<<i<<"--"<<orig_ids[i]<<"+++"<<ordered_particles_host[i]<<std::endl;
                Real3 a = Real3(io->unordered_particles_host[orig_ids[i]]);
                Real3 b = Real3(ordered_particles_host[i]);
                a-=b;
                std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<a<<std::endl;

            }
        }
    }

    if(0)
    {
        cudaDeviceSynchronize();
        std::vector<size_t>perm_ids(io->n);
        for(size_t i = 0; i < io->n; ++i)
        {
            perm_ids[orig_ids[i]] = i;
        }
        for(size_t i = 0; i <  io->n; ++i)
        {

            Real3 a = Real3(io->unordered_particles_host[i]);
            Real3 b = Real3(ordered_particles_host[perm_ids[i]]);
            a-=b;
            if(a.x != 0.0 || a.y != 0.0 || a.z != 0.0)
            {
                std::cout<<i<<"+++"<<io->unordered_particles_host[i]<<"---"<<perm_ids[i]<<std::endl;
                std::cout<<i<<"+++"<<ordered_particles_host[perm_ids[i]]<<std::endl;
                std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<a<<std::endl;
            }
        }
    }
}

}//namespace end

