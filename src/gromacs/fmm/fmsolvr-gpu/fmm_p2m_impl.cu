#include "fmm.hpp"
#include "cuda_P2M.hpp"
#include "cuda_lib.hpp"
#include "cuda_PREPARE_DATA.hpp"

namespace gmx_gpu_fmm{


void fmm_algorithm::p2m_impl(){

#ifndef GMX_FMM_DOUBLE
    __flush_multipoles<CoeffMatrix><<<num_boxes_tree, p1xp2_2,0,priority_streams[current_priority_stream]>>>(omega, mu);
#else
    if(p<21)
    {
        __flush_multipoles<CoeffMatrix><<<num_boxes_tree, p1xp2_2,0,priority_streams[current_priority_stream]>>>(omega, mu);
    }
    else
    {
        int gy = (num_boxes_tree-1)/16 + 1;
        int gz = gy > 1 ? 16 : num_boxes_tree;
        dim3 flush_block(128,1,1);
        dim3 flush_grid((p1xp2_2 - 1)/flush_block.x + 1 ,gy, gz);
        __flush_multipoles2<CoeffMatrix><<<flush_grid, flush_block,0,priority_streams[current_priority_stream]>>>(omega, mu, p1xp2_2, num_boxes_tree);
    }

#endif

    //should be dynamic if too small block size
    size_t split_blocks = 1;
    dim3 P2M_block(bs,split_blocks,1);
    dim3 P2M_grid(n_blocks/split_blocks,1,1);

    __P2M_tree<CoeffMatrix, Real3, Real4><<<P2M_grid,P2M_block, 0, priority_streams[current_priority_stream]>>>
    (p, global_offset, &box_particle_offset[0], &expansion_points[0], &ordered_particles[0], &block_map[0], &offset_map[0], omega);

    __AoS2SoA_omega__(box, omegaSoA, num_boxes_tree, p1xp2_2, priority_streams[current_priority_stream]);

    if(0)
    {
        Device::devSync();
        for (size_t boxid = 0; boxid < num_boxes_tree;++boxid)
        {
            printf("active = %d\n", box[box_id_map[boxid]].active);
            printf("%d -- expansion point (%f %f %f)\n",boxid,  expansion_points[boxid].x, expansion_points[boxid].y, expansion_points[boxid].z);
            dump(*omega[boxid],p);
        }
    }
}

}//namespace end

