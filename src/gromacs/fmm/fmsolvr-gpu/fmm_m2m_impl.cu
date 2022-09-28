#include "fmm.hpp"
#include "cuda_lib.hpp"
#include "cuda_M2M.hpp"

namespace gmx_gpu_fmm{

template <typename CoefficientMatrix,typename CoefficientMatrixSoA, typename Box, typename complex_type>
__global__ void
__prepare_a_SoA(Box *box,CoefficientMatrixSoA *omegas, complex_type*** targets_SoA, const size_t num_boxes_tree, size_t p1xp2_2,size_t p, size_t num_of_operators)
{

    size_t op_id = blockIdx.x * blockDim.x + threadIdx.x;
    if(op_id >= num_of_operators)
        return;

    size_t target_box_id = blockIdx.y;
    size_t source_box_id = blockIdx.z;

    if(source_box_id == 0)
        return;

    CoefficientMatrix *ptr = box[source_box_id].a_targets[op_id];

    if( ptr != box[target_box_id].omega )
    {
        return;
    }
    size_t soa_index = op_id*(num_boxes_tree-1) + source_box_id - 1;
    for(size_t index = 0; index < p1xp2_2; ++index)
        targets_SoA[index][soa_index] = omegas->get_lin_SoA_ptr(index,target_box_id);

}

void fmm_algorithm::set_m2m_operator(){


    for (size_t d_c = depth; d_c > 0; --d_c)
    {
        size_t d_p = d_c - 1;
        size_t depth_offset_p = boxes_above_depth(d_p);
        size_t depth_offset_c = boxes_above_depth(d_c);

        size_t d_delta = d_c - d_p;

        //size_t num_operators = boxes_on_depth(d_delta);
        //M2M_Operator**  m2m_operator;
        //Device::custom_alloc(m2m_operator, num_operators*sizeof(M2M_Operator*));
        //assert(d_delta == 1);

        for (size_t ii = 0; ii <= 1; ++ii) {
            for (size_t jj = 0; jj <= 1; ++jj) {
                for (size_t kk = 0; kk <= 1; ++kk) {
                    size_t idx_p = depth_offset_p + make_boxid(0, 0, 0, d_p);
                    size_t idx_c = depth_offset_c + make_boxid(2 * 0 + ii, 2 * 0 + jj, 2 * 0 + kk, d_c);
                    size_t idx_op = make_boxid(ii, jj, kk, d_delta);
                    m2m_operators[d_p][idx_op] = new M2M_Operator(expansion_points[idx_c], expansion_points[idx_p], p);
                    //printf("m2m operator index %d target depth %d\n",idx_op,d_c);
                    //printf("max %e\n", m2m_operator[idx_op]->max());
                    //printf("min %e\n", m2m_operator[idx_op]->min());
                    //dump(*m2m_operator[idx_op],p);
                }
            }
        }
        //m2m_operators[d_p] = m2m_operator;
    }
}

void fmm_algorithm::prepare_m2m_impl(){

    set_m2m_operator();
    size_t box_id = 1;
    for(size_t d = 1; d <= depth; d++)
    {
        size_t d_target = d - 1;
        size_t depth_offset_target = boxes_above_depth(d_target);
        size_t dim = size_t(1) << d;

        for (ssize_t z = 0; z < 2; ++z)
        {
            for (ssize_t y = 0; y < 2; ++y)
            {
                for (ssize_t x = 0; x < 2; ++x)
                {
                    size_t idx_op = make_boxid(z, y, x, 1);
                    for (ssize_t i = z; i < (ssize_t)dim; i+=2)
                    {
                        for (ssize_t j = y; j < (ssize_t)dim; j+=2)
                        {
                            for (ssize_t k = x; k < (ssize_t)dim; k+=2)
                            {
                                 size_t omega_target_id = depth_offset_target + make_boxid(i/2, j/2, k/2, d_target);
                                 box[box_id].set_a_interaction_pairs(0, omega[omega_target_id], m2m_operators[d-1][idx_op], box_id_map[omega_target_id]);
                                 box_id++;
                            }
                        }
                    }
                }
            }
        }
    }
    if(0)
    {
        for(int i = 0; i < num_boxes_tree; ++i)
        {
            printf("box_id = %d\n",&box[i].id);
        }
    }

#ifdef M2M_SOA_OPTIMIZATION
    typedef typename CoeffMatrix::value_type complex;
    dim3 block(1,1,1);
    dim3 grid(1, num_boxes_tree, num_boxes_tree);
    //not performance critical. Only one call at the beginning
    __prepare_a_SoA<CoeffMatrix, CoeffMatrixSoA, Box, complex>
    <<<grid,block, 0, priority_streams[current_priority_stream]>>>
    (box,omegaSoA,box_a_targets_SoA,num_boxes_tree,p1xp2_2,p,1);
    CUDA_CHECK_ERROR();
#endif

}

void fmm_algorithm::m2m_impl(){

    typedef typename CoeffMatrix::value_type complex_type;

    dim3 b(8,1,1);
    dim3 g(1,1,1);
    dim3 grid(p1,p1,8);
    for (size_t d = depth; d > 0; --d)
    {
        size_t depth_offset = boxes_above_depth(d);
        size_t boxes_on_depth_above = boxes_on_depth(d-1);
        size_t dimx = std::max<size_t>(boxes_on_depth_above,p1);
        dim3 block(dimx,1,1);
        g.x = boxes_on_depth_above;

        __activate_boxes_above<<<g,b,0,priority_streams[current_priority_stream]>>>(box, depth_offset);

        if(d<5)
        {
            __M2M_one<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type, INTTYPE, 0>
            <<<grid,block,p1xx2*sizeof(complex_type), priority_streams[current_priority_stream]>>>
            (box,omegaSoA,box_a_targets_SoA,0,boxes_on_depth_above,num_boxes_tree,depth_offset,p,p1,p1xx2,pxp1_2,p1xp2_2);
        }
        else
        {
            grid.z = size_t(std::pow(8,d-4));
            block.x = 512;
            for(int op_id = 0; op_id < 8; op_id++)
            {
                __M2M_one<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type, INTTYPE, 1>
                <<<grid,block,p1xx2*sizeof(complex_type), priority_streams[current_priority_stream]>>>
                (box,omegaSoA,box_a_targets_SoA,op_id,boxes_on_depth_above,num_boxes_tree,depth_offset,p,p1,p1xx2,pxp1_2,p1xp2_2);
            }
        }
        CUDA_CHECK_ERROR();
    }
    //no dipole computed if p < 1
    if(p > 0)
    {
        __get_dipole<<<1,1,0,priority_streams[current_priority_stream]>>>(omegaSoA, dipole);
        cudaEventRecord(priority_events[current_priority_stream], priority_streams[current_priority_stream]);
    }

    if(0)
    {
        __SoA2AoS_omega__(box, omegaSoA, num_boxes_tree, p1xp2_2, priority_streams[current_priority_stream]);
        Device::devSync();
        //for (size_t boxid = 0; boxid < num_boxes_tree - num_boxes_lowest;++boxid)
        for (size_t boxid = 0; boxid < num_boxes_tree;++boxid)
        {
            printf("active = %d\n", box[box_id_map[boxid]].active);
            printf("box index %d  id %d\n", boxid, box_id_map[boxid]);

            dump(*omega[boxid],p);
        }
    }

}

}//namespace end
