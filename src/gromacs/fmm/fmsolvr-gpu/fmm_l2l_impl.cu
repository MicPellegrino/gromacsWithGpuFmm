#include "fmm.hpp"
#include "cuda_L2L.hpp"
#include "cuda_lib.hpp"

namespace gmx_gpu_fmm{

template <typename CoefficientMatrix,typename CoefficientMatrixSoA, typename Box, typename complex_type>
__global__ void
__prepare_c_SoA(Box *box,CoefficientMatrixSoA *mus, complex_type*** targets_SoA, const size_t num_boxes_tree, size_t p1xp2_2,size_t p, size_t num_of_operators)
{
    size_t op_id = blockIdx.x * blockDim.x + threadIdx.x;
    if(op_id >= num_of_operators)
        return;

    size_t target_box_id = blockIdx.y;
    size_t source_box_id = blockIdx.z;

    CoefficientMatrix *ptr = box[source_box_id].c_targets[op_id];

    if( ptr != box[target_box_id].mu)
    {
        return;
    }

    size_t soa_index = op_id*(num_boxes_tree-1) + source_box_id;
    for(size_t index = 0; index < p1xp2_2; ++index)
        targets_SoA[index][soa_index] = mus->get_lin_SoA_ptr(index,target_box_id);
}

void fmm_algorithm::set_l2l_operator(){

    for (size_t d_c = depth; d_c > 0; --d_c){
        size_t d_p = d_c - 1;
        size_t depth_offset_p = boxes_above_depth(d_p);
        size_t depth_offset_c = boxes_above_depth(d_c);
        size_t d_delta = d_c - d_p;
        //size_t num_operators = boxes_on_depth(d_delta);

        //L2L_Operator **  l2l_operator;
        //Device::custom_alloc(l2l_operator, num_operators*sizeof(L2L_Operator*));
        for (size_t ii = 0; ii <= 1; ++ii) {
            for (size_t jj = 0; jj <= 1; ++jj) {
                for (size_t kk = 0; kk <= 1; ++kk) {
                    size_t idx_p = depth_offset_p + make_boxid(0, 0, 0, d_p);
                    size_t idx_c = depth_offset_c + make_boxid(2 * 0 + ii, 2 * 0 + jj, 2 * 0 + kk, d_c);
                    size_t idx_op = make_boxid(ii, jj, kk, d_delta);
                    l2l_operators[d_p][idx_op] = new L2L_Operator(expansion_points[idx_p], expansion_points[idx_c], p);
                    //printf("l2l operator index %d target depth %d\n",idx_op,d_p);
                    //printf("max %e\n", l2l_operator[idx_op]->max());
                    //printf("min %e\n", l2l_operator[idx_op]->min());
                    //dump(*l2l_operator[idx_op],p);
                }
            }
        }
        //l2l_operators[d_p] = l2l_operator;
    }
}

void fmm_algorithm::prepare_l2l_impl(){

    set_l2l_operator();
    size_t box_id = 0;
    for(size_t d = 0; d < depth; d++)
    {
        size_t d_target = d + 1;
        size_t depth_offset_target = boxes_above_depth(d_target);
        size_t dim = size_t(1) << d;

        for (ssize_t z = 0; z < 2; ++z)
        {
            for (ssize_t y = 0; y < 2; ++y)
            {
                for (ssize_t x = 0; x < 2; ++x)
                {
                    for (ssize_t i = z; i < (ssize_t)dim; i+=2)
                    {
                        for (ssize_t j = y; j < (ssize_t)dim; j+=2)
                        {
                            for (ssize_t k = x; k < (ssize_t)dim; k+=2)
                            {
                                size_t local_op_index = 0;
                                for (size_t ii = 0; ii <= 1; ++ii)
                                {
                                    for (size_t jj = 0; jj <= 1; ++jj)
                                    {
                                        for (size_t kk = 0; kk <= 1; ++kk)
                                        {

                                            size_t idx_op = make_boxid(ii, jj, kk, 1);
                                            size_t mu_target_id = depth_offset_target + make_boxid(2*i+ii, 2*j+jj, 2*k+kk, d_target);
                                            box[box_id].set_c_interaction_pairs(local_op_index, mu[mu_target_id], l2l_operators[d][idx_op], box_id_map[mu_target_id]);
                                            ++local_op_index;
                                        }
                                    }
                                }
                                box_id++;
                            }
                        }
                    }
                }
            }
        }
    }
    typedef typename CoeffMatrix::value_type complex;

#ifdef L2L_SOA_OPTIMIZATION
    dim3 block(512,1,1);
    dim3 grid((8-1)/block.x+1,num_boxes_tree, num_boxes_tree - num_boxes_lowest);

    if(grid.z > 0)
    {
        //not performance critical. Only one call at the beginning
        __prepare_c_SoA<CoeffMatrix, CoeffMatrixSoA, Box, complex>
        <<<grid, block, 0, priority_streams[current_priority_stream]>>>
        (box,muSoA,box_c_targets_SoA,num_boxes_tree,p1xp2_2,p,8);
        CUDA_CHECK_ERROR();
    }
#endif
    Device::devSync();

}

void fmm_algorithm::l2l_impl(){

    typedef typename CoeffMatrix::value_type complex_type;

    dim3 grid(p1,p1,8);
    for (size_t d = 1; d <= depth; ++d)
    {
        size_t depth_offset = boxes_above_depth(d-1);
        size_t boxes_on_depth_above = boxes_on_depth(d-1);

        size_t dimx = std::max<size_t>(boxes_on_depth_above,p+1);
        dim3 block(dimx,1,1);

        if(d<5)
        {
            __L2L_one<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type, INTTYPE, 0>
            <<<grid,block,p1xx2*sizeof(complex_type), priority_streams[current_priority_stream]>>>
            (box,muSoA,box_c_targets_SoA,0,boxes_on_depth_above,num_boxes_tree,depth_offset,p,p1,p1xx2,pxp1_2,p1xp2_2);
        }
        else
        {
            block.x = 512;
            for(int op_id = 0; op_id < 8; op_id++)
            {
                __L2L_one<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type, INTTYPE, 1>
                <<<grid,block,p1xx2*sizeof(complex_type), priority_streams[current_priority_stream]>>>
                (box,muSoA,box_c_targets_SoA,op_id,boxes_on_depth_above,num_boxes_tree,depth_offset,p,p1,p1xx2,pxp1_2,p1xp2_2);
            }
            grid.z *=8;

        }
        CUDA_CHECK_ERROR();
    }

    __SoA2AoS_mu__(box, muSoA, num_boxes_tree, p1xp2_2, priority_streams[current_priority_stream]);

    if(0)
    {
        Device::devSync();
        for (size_t boxid = global_offset; boxid < num_boxes_tree;++boxid)
        {
            printf("%d\n", boxid);
            std::cout<<expansion_points[boxid]<<std::endl;
            dump(*mu[boxid],p);
        }

        M2L_Operator* B = new M2L_Operator(expansion_points[52] - expansion_points[49], p);
        dump(*B,2*p);

        M2L_Operator* BB = new M2L_Operator(expansion_points[49] - expansion_points[52], p);
        dump(*BB,2*p);
    }
}

}//namespace end

