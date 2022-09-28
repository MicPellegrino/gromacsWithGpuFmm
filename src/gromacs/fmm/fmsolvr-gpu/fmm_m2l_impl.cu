#include "fmm_m2l_impl.hpp"

namespace gmx_gpu_fmm{

CUDA
size_t make_operator_boxid(size_t x, size_t y, size_t z, unsigned int dim)
{
    size_t size = dim/2;
    return ((z+size)*dim*dim) + ((y+size)*dim) + x + size;
}

template <typename CoefficientMatrix,typename CoefficientMatrixSoA, typename Box, typename complex_type>
__global__ void
__prepare_b_SoA(Box *box,CoefficientMatrixSoA *mus, complex_type*** targets_SoA, const size_t num_boxes_tree, size_t p1xp2_2,size_t p, size_t num_of_operators)
{
    size_t op_id = threadIdx.x;
    size_t target_box_id = blockIdx.y + 1;
    size_t source_box_id = blockIdx.z + 1;

    CoefficientMatrix *ptr = box[source_box_id].b_targets[op_id];

    if( ptr == box[target_box_id].mu )
    {
        size_t soa_index = op_id*(num_boxes_tree-1) + source_box_id - 1;
        //printf("SOA %d\n",soa_index);
        for(size_t index = 0; index < p1xp2_2; ++index)
            targets_SoA[index][soa_index] =  mus->get_lin_SoA_ptr(index,target_box_id);
    }
}


void prepare_lm_map(size_t p, cudaStream_t &stream){
    __prepare_lm_map<<<1,1,0,stream>>>(p);
}

void fmm_algorithm::prepare_m2l_impl()
{
    prepare_lm_map(p, priority_streams[current_priority_stream]);
    CUDA_CHECK_ERROR();
    Device::devSync();
    unsigned int bitset_size = 1;
    if(p>0)
        bitset_size = bitset_offset(2*(p1)*(p1)-1)+1;

    size_t num_of_bitsets = p1xp2_2;

    Device::custom_alloc(bitset_mem, num_of_all_ops * num_of_bitsets * bitset_size * sizeof(unsigned int));

    //printf("num of all ops %d\n", num_of_all_ops);
    bitset = new Bitit<32,5,Device>[num_of_all_ops * num_of_bitsets];
    for (size_t i = 0; i < num_of_all_ops*num_of_bitsets; ++i)
    {
        bitset[i].init(&bitset_mem[i*bitset_size], bitset_size);
    }

    Real3 exp_delta = Real3(1.,1.,1.);

    non_periodic_m2l_dummy = new M2L_Operator(io->abc*exp_delta, p);
    non_periodic_mu_dummy  = new CoeffMatrix(p);

    for (size_t box_id = 1; box_id < num_boxes_tree; ++box_id)
    {
        for (size_t op_id = 0; op_id < num_of_efective_ops; ++op_id)
        {
            box[box_id].set_b_interaction_pairs(op_id, non_periodic_mu_dummy, non_periodic_m2l_dummy, num_boxes_tree);
        }
    }

    for(size_t d = 1; d <= depth; d++)
    {
        size_t d_index = d - 1;
        dummy_m2l_operators[d_index] = new M2L_Operator(Real3(1,1,1), p);
        for (size_t op_id = 0; op_id < num_of_all_ops; ++op_id)
        {
             m2l_operators[d_index][op_id] = dummy_m2l_operators[d_index];
        }
    }

    typedef typename CoeffMatrix::value_type complex;
    size_t box_id = 1;
    for(size_t d = 1; d <= depth; d++)
    {
        //int index = 0;
        ssize_t dim = ssize_t(1) << d;
        Real scale = reciprocal(Real(dim));
        size_t d_index = d - 1;

        for (ssize_t ii = (3+4*ws)/2; ii >= 0; --ii)
        {
            for (ssize_t jj = (3+4*ws)/2; jj >= 0; --jj)
            {
                for (ssize_t kk = (3+4*ws)/2; kk >= 0; --kk)
                {
                    for (ssize_t si = -1; si <= 1; si+=2)
                    {
                        for (ssize_t sj = -1; sj <= 1; sj+=2)
                        {
                            for (ssize_t sk = -1; sk <= 1; sk+=2)
                            {
                                ssize_t i = ii*si;
                                ssize_t j = jj*sj;
                                ssize_t k = kk*sk;

                                if(ii == 0 && si == -1)
                                    continue;
                                if(jj == 0 && sj == -1)
                                    continue;
                                if(kk == 0 && sk == -1)
                                    continue;

                                if(i == 0 && j == 0 && k == 0)
                                    continue;

                                //printf("%d %d %d\n",i,j,k);
                                size_t idx_op = make_operator_boxid(i, j, k, 3+4*ws);
                                Real3 expansionpointdelta = Real3(i, j, k) * scale;
                                m2l_operators[d_index][idx_op] = new M2L_Operator(io->abc*expansionpointdelta, p);
                                //printf("m2l operator index %d depth %d linear_index %d\n",idx_op,d, index++);
                                //printf("max %e\n", m2l_operator[idx_op]->max());
                                //printf("min %e\n", m2l_operator[idx_op]->min());
                                //dump(*m2l_operator[idx_op],2*p);
                                m2l_operators[d_index][idx_op]->init_bitset(&(bitset[idx_op * num_of_bitsets]));
                                m2l_operators[d_index][idx_op]->set_bitset();
                            }
                        }
                    }
                }
            }
        }
    }

    box_id = 1;
    for (size_t d = 1; d <= depth; ++d)
    {
        size_t depth_offset = boxes_above_depth(d);
        ssize_t dim = ssize_t(1) << d;
        //printf("DEPTH %lu DIMENSION %lu\n",d,dim);
        size_t glob_pos = 0;
        for (ssize_t z = 0; z < 2; ++z)
            for (ssize_t y = 0; y < 2; ++y)
                for (ssize_t x = 0; x < 2; ++x)
                {
                    for (ssize_t i = z; i < dim; i+=2)
                    {
                        for (ssize_t j = y; j < dim; j+=2)
                        {
                            for (ssize_t k = x; k < dim; k+=2)
                            {
                                //size_t omega_id = depth_offset + make_boxid(i, j, k, d);
                                std::vector<std::vector<size_t> > positions(4,std::vector<size_t>());

                                // parent box coordinates
                                ssize_t ip = i / 2;
                                ssize_t jp = j / 2;
                                ssize_t kp = k / 2;
                                // outer ws shells
                                ssize_t io_beg = (ip - ws) * 2;
                                ssize_t io_end = (ip + ws + 1) * 2;
                                ssize_t jo_beg = (jp - ws) * 2;
                                ssize_t jo_end = (jp + ws + 1) * 2;
                                ssize_t ko_beg = (kp - ws) * 2;
                                ssize_t ko_end = (kp + ws + 1) * 2;
                                // inner - no far field interactions
                                ssize_t ii_beg = i - ws;
                                ssize_t ii_end = i + ws + 1;
                                ssize_t ji_beg = j - ws;
                                ssize_t ji_end = j + ws + 1;
                                ssize_t ki_beg = k - ws;
                                ssize_t ki_end = k + ws + 1;

                                size_t sum_of_global_indizes = 0;
                                size_t local_op_index = 0;
                                size_t same_op = 0;
                                for (ssize_t ii_ = io_beg-1; ii_ <= i; ++ii_)
                                {
                                    for (ssize_t jj_ = jo_beg-1; jj_ <= j; ++jj_)
                                    {
                                        for (ssize_t kk_ = ko_beg-1; kk_ <= k; ++kk_)
                                        {
                                            if(same_op == 1)
                                            {
                                                positions[0].push_back(local_op_index-1);
                                            }

                                            if(same_op == 2)
                                            {
                                                positions[1].push_back(local_op_index-2);
                                            }

                                            if(same_op == 4)
                                            {
                                                positions[2].push_back(local_op_index-4);
                                            }

                                            if(same_op == 8)
                                            {
                                                positions[3].push_back(local_op_index-8);
                                            }

                                            same_op = 0;
                                            for (ssize_t si = -1; si <= 1; si+=2)
                                            {
                                                for (ssize_t sj = -1; sj <= 1; sj+=2)
                                                {
                                                    for (ssize_t sk = -1; sk <= 1; sk+=2)
                                                    {

                                                        if(ii_ == i && si ==-1)
                                                            continue;
                                                        if(jj_ == j && sj ==-1)
                                                            continue;
                                                        if(kk_ == k && sk ==-1)
                                                            continue;

                                                        ssize_t i_ii = i - ii_;
                                                        ssize_t j_jj = j - jj_;
                                                        ssize_t k_kk = k - kk_;

                                                        ssize_t ii = i + i_ii * si;
                                                        ssize_t jj = j + j_jj * sj;
                                                        ssize_t kk = k + k_kk * sk;

                                                        if(ii >= io_end || jj >= jo_end || kk >= ko_end)
                                                            continue;
                                                        if(ii < io_beg || jj < jo_beg || kk < ko_beg)
                                                            continue;

                                                        if (ii < ii_beg || ii >= ii_end
                                                                || jj < ji_beg || jj >= ji_end
                                                                || kk < ki_beg || kk >= ki_end)
                                                        {
                                                            size_t idx_op;
                                                            if (open_boundary_conditions && (ii < 0 || ii >= dim || jj < 0 || jj >= dim || kk < 0 || kk >= dim))
                                                            {

                                                            }
                                                            else
                                                            {
                                                                ssize_t iif = ii;
                                                                while (iif < 0) {
                                                                    iif += dim;
                                                                }
                                                                while (iif >= dim) {
                                                                    iif -= dim;
                                                                }
                                                                ssize_t jjf = jj;
                                                                while (jjf < 0) {
                                                                    jjf += dim;
                                                                }
                                                                while (jjf >= dim) {
                                                                    jjf -= dim;
                                                                }
                                                                ssize_t kkf = kk;
                                                                while (kkf < 0) {
                                                                    kkf += dim;
                                                                }
                                                                while (kkf >= dim) {
                                                                    kkf -= dim;
                                                                }
                                                                assert(0 <= iif && iif < dim && 0 <= jjf && jjf < dim && 0 <= kkf && kkf < dim);
                                                                size_t mu_id = depth_offset + make_boxid(iif, jjf, kkf, d);
                                                                idx_op = make_operator_boxid(i-ii, j-jj, k-kk, 3+4*ws);
                                                                size_t periodic_mu = box_id_map[mu_id];

                                                                box[box_id].set_b_interaction_pairs(local_op_index, mu[mu_id], m2l_operators[d-1][idx_op], periodic_mu);
                                                            }


                                                            ++local_op_index;
                                                            same_op++;
                                                            sum_of_global_indizes +=idx_op;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                box[box_id].set_permutations(positions);
                                box_id++;
                            }
                        }
                    }
                    glob_pos++;
                }
    }

    for(size_t box_id = 1; box_id < num_boxes_tree; box_id++)
    {
        box[box_id].permute_ops(box_operators_tmp, box_targets_tmp, box_target_ids_tmp, num_of_efective_ops);
    }

    //printf("num_of_efective_ops=%d\n",num_of_efective_ops);
    //dim3 block(512,1,1);
    //dim3 grid((num_of_efective_ops-1)/block.x+1,num_boxes_tree, num_boxes_tree);
    Device::devSync();

#ifdef M2L_SOA_OPTIMIZATION
    dim3 block(189,1,1);
    dim3 grid(1,num_boxes_tree-1,num_boxes_tree-1);
    //not performance critical. Only one call at the beginning
    if(num_boxes_tree > 1)
    {
        __prepare_b_SoA<CoeffMatrix, CoeffMatrixSoA, Box, complex>
        <<<grid, block, 0, priority_streams[current_priority_stream]>>>(box, muSoA, box_b_targets_SoA, num_boxes_tree, p1xp2_2,p, num_of_efective_ops);
        Device::devSync();
        CUDA_CHECK_ERROR();
    }
#endif
}

void fmm_algorithm::m2l_impl(){

    typedef typename CoeffMatrix::value_type complex_type;

    const size_t num_of_streams = STREAMS;

    //wait for m2m in all streams we need for m2l
    for(int i = 0; i < num_of_streams; ++i)
    {
        cudaStreamWaitEvent(priority_streams[i], priority_events[current_priority_stream], 0);
    }

    size_t op_p1xx2 = (2*p1) * (2*p1);
    //m2m will be ready here
    if(open_boundary_conditions)
    {
        __SoA2AoS_omega__(box, omegaSoA, num_boxes_tree, p1xp2_2, priority_streams[current_priority_stream]);
        cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
        int event_to_wait = current_priority_stream;

        for (size_t d = 1; d <= depth; ++d)
        {
            current_priority_stream = (++current_priority_stream)%num_of_streams;
            size_t depth_offset = boxes_above_depth(d);
            size_t boxes_on_depth_above = boxes_on_depth(d-1);
            //printf("%d\n", boxes_on_depth_above);
            dim3 grid(boxes_on_depth_above,189,8);
            dim3 block(p1+1,p1,1);

            cudaStreamWaitEvent(priority_streams[current_priority_stream], priority_events[event_to_wait], 0);

#ifndef GMX_FMM_DOUBLE
            __M2L_one_p2<CoeffMatrix, M2L_Operator, Box, Real, Real3, complex_type>
            <<<grid,block,(p1xx2+op_p1xx2)*sizeof(complex_type),priority_streams[current_priority_stream]>>>
            (box,omega,depth_offset,num_boxes_tree,p,p1,p1xx2,op_p1xx2);

#else
            if(p < 21)
            {
                __M2L_one_p2<CoeffMatrix, M2L_Operator, Box, Real, Real3, complex_type>
                <<<grid,block,(p1xx2+op_p1xx2)*sizeof(complex_type),priority_streams[current_priority_stream]>>>
                (box,omega,depth_offset,num_boxes_tree,p,p1,p1xx2,op_p1xx2);
            }
            else
            {
                CUDA_CHECK_ERROR();
                dim3 griD(boxes_on_depth_above,189,p1);
                dim3 blocK(p1,8,1);
                __M2L_one_p2_no_shared<CoeffMatrix, M2L_Operator, Box, Real, Real3, complex_type>
                <<<griD,blocK,0,priority_streams[current_priority_stream]>>>
                (box,omega,depth_offset,num_boxes_tree,p,p1,p1xx2,op_p1xx2);
                CUDA_CHECK_ERROR();
            }

#endif
            cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
        }
        for(int i = 0; i < num_of_streams; ++i)
        {
            cudaStreamWaitEvent(priority_streams[current_priority_stream], priority_events[i], 0);
        }

        __AoS2SoA_mu__(box, muSoA, num_boxes_tree, p1xp2_2, priority_streams[current_priority_stream]);
        cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
    }
    else
    {

        for (size_t d = 1; d <= depth; ++d)
        {

            size_t depth_offset = boxes_above_depth(d);
            size_t boxes_on_depth_above = boxes_on_depth(d-1);

#ifndef GMX_FMM_DOUBLE
            if(d < 3 || sparse_fmm)
            {

                dim3 grid(boxes_on_depth_above,189,8);
                dim3 block(p1,p1+1,1);

                if(sparse_fmm)
                {
                    __M2L_one_p2_SoA<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type, true>
                    <<<grid,block,(p1xx2+op_p1xx2)*sizeof(complex_type),priority_streams[current_priority_stream]>>>
                    (box,omegaSoA,muSoA,box_b_targets_SoA,depth_offset,num_boxes_tree,p,p1,p1xx2,op_p1xx2);
                }
                else
                {
                    __M2L_one_p2_SoA<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type, false>
                    <<<grid,block,(p1xx2+op_p1xx2)*sizeof(complex_type),priority_streams[current_priority_stream]>>>
                    (box,omegaSoA,muSoA,box_b_targets_SoA,depth_offset,num_boxes_tree,p,p1,p1xx2,op_p1xx2);
                }
                cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
                current_priority_stream = (++current_priority_stream)%num_of_streams;
            }
            else
            {
                size_t dimx = std::max<size_t>(boxes_on_depth_above, p1xx2);
                size_t gridz = 1;

                if(dimx > 512)
                {
                    dimx = 512;
                    gridz = boxes_on_depth_above / dimx;
                }

                dim3 block(dimx,1,1);
                unsigned int size_of_bitset = 1;
                if(p1xx2 > 0)
                    size_of_bitset = bitset_offset(2*p1xx2-1) + 1;

                if(ws==1)
                {
                    dim3 grid_one(p1xp2_2,7,gridz);
                    dim3 grid_two(p1xp2_2,21,gridz);
                    dim3 grid_four(p1xp2_2,21,gridz);
                    dim3 grid_eight(p1xp2_2,7,gridz);

                    for (size_t op_id = 0; op_id < 8; ++op_id)
                    {
                        __M2L_one_shared_v8_all<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type,0,1,INTTYPE>
                                <<<grid_one,block,p1xx2*sizeof(complex_type),priority_streams[current_priority_stream]>>>
                                (box,omegaSoA,muSoA,box_b_targets_SoA,op_id,boxes_on_depth_above,num_boxes_tree,depth_offset,p,p1,p1xx2,pxp1_2,p1xp2_2, size_of_bitset);

                        cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
                        current_priority_stream = (++current_priority_stream)%num_of_streams;

                        __M2L_one_shared_v8_all<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type,7,2,INTTYPE>
                                <<<grid_two,block,p1xx2*sizeof(complex_type)+size_of_bitset*sizeof(unsigned int),priority_streams[current_priority_stream]>>>
                                (box,omegaSoA,muSoA,box_b_targets_SoA,op_id,boxes_on_depth_above,num_boxes_tree,depth_offset,p,p1,p1xx2,pxp1_2,p1xp2_2, size_of_bitset);

                        cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
                        current_priority_stream = (++current_priority_stream)%num_of_streams;

                        __M2L_one_shared_v8_all<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type,49,4,INTTYPE>
                                <<<grid_four,block,p1xx2*sizeof(complex_type)+3*size_of_bitset*sizeof(unsigned int),priority_streams[current_priority_stream]>>>
                                 (box,omegaSoA,muSoA,box_b_targets_SoA,op_id,boxes_on_depth_above,num_boxes_tree,depth_offset,p,p1,p1xx2,pxp1_2,p1xp2_2, size_of_bitset);

                        cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
                        current_priority_stream = (++current_priority_stream)%num_of_streams;

                        __M2L_one_shared_v8_all<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type,133,8,INTTYPE>
                                <<<grid_eight,block,p1xx2*sizeof(complex_type)+7*size_of_bitset*sizeof(unsigned int),priority_streams[current_priority_stream]>>>
                                (box,omegaSoA,muSoA,box_b_targets_SoA,op_id,boxes_on_depth_above,num_boxes_tree,depth_offset,p,p1,p1xx2,pxp1_2,p1xp2_2, size_of_bitset);

                        cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
                        current_priority_stream = (++current_priority_stream)%num_of_streams;
                    }
                }

                if(ws==2)
                {
                    dim3 grid_one(7,p1,p1);
                    dim3 grid_two(42,p1,p1);
                    dim3 grid_four(84,p1,p1);
                    dim3 grid_eight(56,p1,p1);

                    for (size_t op_id = 0; op_id < 8; ++op_id)
                    {
                        __M2L_one_shared_v8_all<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type,0,1,INTTYPE>
                                <<<grid_one,block,p1xx2*sizeof(complex_type),priority_streams[current_priority_stream]>>>
                                (box,omegaSoA,muSoA,box_b_targets_SoA,op_id,boxes_on_depth_above,num_boxes_tree,depth_offset,p,p1,p1xx2,pxp1_2,p1xp2_2, size_of_bitset);

                        cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
                        current_priority_stream = (++current_priority_stream)%num_of_streams;

                        __M2L_one_shared_v8_all<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type,7,2,INTTYPE>
                                <<<grid_two,block,p1xx2*sizeof(complex_type)+size_of_bitset*sizeof(unsigned int),priority_streams[current_priority_stream]>>>
                                (box,omegaSoA,muSoA,box_b_targets_SoA,op_id,boxes_on_depth_above,num_boxes_tree,depth_offset,p,p1,p1xx2,pxp1_2,p1xp2_2, size_of_bitset);

                        cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
                        current_priority_stream = (++current_priority_stream)%num_of_streams;

                        __M2L_one_shared_v8_all<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type,91,4,INTTYPE>
                                <<<grid_four,block,p1xx2*sizeof(complex_type)+3*size_of_bitset*sizeof(unsigned int),priority_streams[current_priority_stream]>>>
                                (box,omegaSoA,muSoA,box_b_targets_SoA,op_id,boxes_on_depth_above,num_boxes_tree,depth_offset,p,p1,p1xx2,pxp1_2,p1xp2_2, size_of_bitset);

                        cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
                        current_priority_stream = (++current_priority_stream)%num_of_streams;

                        __M2L_one_shared_v8_all<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type,427,8,INTTYPE>
                                <<<grid_eight,block,p1xx2*sizeof(complex_type)+7*size_of_bitset*sizeof(unsigned int),priority_streams[current_priority_stream]>>>
                                (box,omegaSoA,muSoA,box_b_targets_SoA,op_id,boxes_on_depth_above,num_boxes_tree,depth_offset,p,p1,p1xx2,pxp1_2,p1xp2_2, size_of_bitset);

                        cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
                        current_priority_stream = (++current_priority_stream)%num_of_streams;
                    }
                }
            }
        }
#else
            if(p < 21)
            {
                if(d < 3 || sparse_fmm)
                {

                    dim3 grid(boxes_on_depth_above,189,8);
                    dim3 block(p1,p1+1,1);

                    if(sparse_fmm)
                    {
                        __M2L_one_p2_SoA<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type, true>
                        <<<grid,block,(p1xx2+op_p1xx2)*sizeof(complex_type),priority_streams[current_priority_stream]>>>
                        (box,omegaSoA,muSoA,box_b_targets_SoA,depth_offset,num_boxes_tree,p,p1,p1xx2,op_p1xx2);
                    }
                    else
                    {
                        __M2L_one_p2_SoA<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type, false>
                        <<<grid,block,(p1xx2+op_p1xx2)*sizeof(complex_type),priority_streams[current_priority_stream]>>>
                        (box,omegaSoA,muSoA,box_b_targets_SoA,depth_offset,num_boxes_tree,p,p1,p1xx2,op_p1xx2);
                    }
                    cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
                    current_priority_stream = (++current_priority_stream)%num_of_streams;
                }
                else
                {
                    size_t dimx = std::max<size_t>(boxes_on_depth_above, p1xx2);
                    size_t gridz = 1;

                    if(dimx > 512)
                    {
                        dimx = 512;
                        gridz = boxes_on_depth_above / dimx;
                    }

                    dim3 block(dimx,1,1);
                    unsigned int size_of_bitset = 1;
                    if(p1xx2 > 0)
                        size_of_bitset = bitset_offset(2*p1xx2-1) + 1;

                    if(ws==1)
                    {
                        dim3 grid_one(p1xp2_2,7,gridz);
                        dim3 grid_two(p1xp2_2,21,gridz);
                        dim3 grid_four(p1xp2_2,21,gridz);
                        dim3 grid_eight(p1xp2_2,7,gridz);

                        for (size_t op_id = 0; op_id < 8; ++op_id)
                        {
                            __M2L_one_shared_v8_all<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type,0,1,INTTYPE>
                                    <<<grid_one,block,p1xx2*sizeof(complex_type),priority_streams[current_priority_stream]>>>
                                    (box,omegaSoA,muSoA,box_b_targets_SoA,op_id,boxes_on_depth_above,num_boxes_tree,depth_offset,p,p1,p1xx2,pxp1_2,p1xp2_2, size_of_bitset);

                            cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
                            current_priority_stream = (++current_priority_stream)%num_of_streams;

                            __M2L_one_shared_v8_all<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type,7,2,INTTYPE>
                                    <<<grid_two,block,p1xx2*sizeof(complex_type)+size_of_bitset*sizeof(unsigned int),priority_streams[current_priority_stream]>>>
                                    (box,omegaSoA,muSoA,box_b_targets_SoA,op_id,boxes_on_depth_above,num_boxes_tree,depth_offset,p,p1,p1xx2,pxp1_2,p1xp2_2, size_of_bitset);

                            cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
                            current_priority_stream = (++current_priority_stream)%num_of_streams;

                            __M2L_one_shared_v8_all<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type,49,4,INTTYPE>
                                    <<<grid_four,block,p1xx2*sizeof(complex_type)+3*size_of_bitset*sizeof(unsigned int),priority_streams[current_priority_stream]>>>
                                     (box,omegaSoA,muSoA,box_b_targets_SoA,op_id,boxes_on_depth_above,num_boxes_tree,depth_offset,p,p1,p1xx2,pxp1_2,p1xp2_2, size_of_bitset);

                            cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
                            current_priority_stream = (++current_priority_stream)%num_of_streams;

                            __M2L_one_shared_v8_all<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type,133,8,INTTYPE>
                                    <<<grid_eight,block,p1xx2*sizeof(complex_type)+7*size_of_bitset*sizeof(unsigned int),priority_streams[current_priority_stream]>>>
                                    (box,omegaSoA,muSoA,box_b_targets_SoA,op_id,boxes_on_depth_above,num_boxes_tree,depth_offset,p,p1,p1xx2,pxp1_2,p1xp2_2, size_of_bitset);

                            cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
                            current_priority_stream = (++current_priority_stream)%num_of_streams;
                        }
                    }

                    if(ws==2)
                    {
                        dim3 grid_one(7,p1,p1);
                        dim3 grid_two(42,p1,p1);
                        dim3 grid_four(84,p1,p1);
                        dim3 grid_eight(56,p1,p1);

                        for (size_t op_id = 0; op_id < 8; ++op_id)
                        {
                            __M2L_one_shared_v8_all<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type,0,1,INTTYPE>
                                    <<<grid_one,block,p1xx2*sizeof(complex_type),priority_streams[current_priority_stream]>>>
                                    (box,omegaSoA,muSoA,box_b_targets_SoA,op_id,boxes_on_depth_above,num_boxes_tree,depth_offset,p,p1,p1xx2,pxp1_2,p1xp2_2, size_of_bitset);

                            cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
                            current_priority_stream = (++current_priority_stream)%num_of_streams;

                            __M2L_one_shared_v8_all<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type,7,2,INTTYPE>
                                    <<<grid_two,block,p1xx2*sizeof(complex_type)+size_of_bitset*sizeof(unsigned int),priority_streams[current_priority_stream]>>>
                                    (box,omegaSoA,muSoA,box_b_targets_SoA,op_id,boxes_on_depth_above,num_boxes_tree,depth_offset,p,p1,p1xx2,pxp1_2,p1xp2_2, size_of_bitset);

                            cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
                            current_priority_stream = (++current_priority_stream)%num_of_streams;

                            __M2L_one_shared_v8_all<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type,91,4,INTTYPE>
                                    <<<grid_four,block,p1xx2*sizeof(complex_type)+3*size_of_bitset*sizeof(unsigned int),priority_streams[current_priority_stream]>>>
                                    (box,omegaSoA,muSoA,box_b_targets_SoA,op_id,boxes_on_depth_above,num_boxes_tree,depth_offset,p,p1,p1xx2,pxp1_2,p1xp2_2, size_of_bitset);

                            cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
                            current_priority_stream = (++current_priority_stream)%num_of_streams;

                            __M2L_one_shared_v8_all<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type,427,8,INTTYPE>
                                    <<<grid_eight,block,p1xx2*sizeof(complex_type)+7*size_of_bitset*sizeof(unsigned int),priority_streams[current_priority_stream]>>>
                                    (box,omegaSoA,muSoA,box_b_targets_SoA,op_id,boxes_on_depth_above,num_boxes_tree,depth_offset,p,p1,p1xx2,pxp1_2,p1xp2_2, size_of_bitset);

                            cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
                            current_priority_stream = (++current_priority_stream)%num_of_streams;
                        }
                    }
                }
            }
            else
            {
                __SoA2AoS_omega__(box, omegaSoA, num_boxes_tree, p1xp2_2, priority_streams[current_priority_stream]);
                CUDA_CHECK_ERROR();
                dim3 griD(boxes_on_depth_above,189,p1);
                dim3 blocK(p1,8,1);
                __M2L_one_p2_no_shared<CoeffMatrix, M2L_Operator, Box, Real, Real3, complex_type>
                <<<griD,blocK,0,priority_streams[current_priority_stream]>>>
                (box,omega,depth_offset,num_boxes_tree,p,p1,p1xx2,op_p1xx2);
                CUDA_CHECK_ERROR();

                __AoS2SoA_mu__(box, muSoA, num_boxes_tree, p1xp2_2, priority_streams[current_priority_stream]);
                cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
                CUDA_CHECK_ERROR();
            }
        }
#endif
    }

    if(0)
    {
        __SoA2AoS_mu__(box, muSoA, num_boxes_tree, p1xp2_2, priority_streams[current_priority_stream]);
        Device::devSync();
        for (size_t boxid = 0; boxid < num_boxes_tree;++boxid)
        {
            printf("MU %d\n",boxid);
            dump(*mu[boxid],p);
        }
    }
}

std::vector<size_t> fmm_algorithm::m2l_impl_gpu_cpu_prep(){

    __SoA2AoS_omega__(box,omegaSoA,num_boxes_tree,p1xp2_2, priority_streams[current_priority_stream]);
    Device::devSync();

    size_t num_boxes = num_boxes_tree;
    size_t threads = 8;
    size_t chunk_size = num_boxes/threads;
    std::vector<size_t> box_chunks(threads+1,0);

    box_chunks[0] = 1;
    box_chunks[threads] = num_boxes;

    for(size_t i = 1; i < threads; i++)
    {
        box_chunks[i] = i*chunk_size;
    }

    return box_chunks;
}

void fmm_algorithm::m2l_impl_gpu_cpu(std::vector<size_t>& box_chunks){


    //std::thread a(&fmm_algorithm::m2l_impl, this);
    /*
    std::thread a(make_m2l_on_chunk<Box>, box_chunks[0], box_chunks[1], num_of_efective_ops, box, p);
    std::thread b(make_m2l_on_chunk<Box>, box_chunks[1], box_chunks[2], num_of_efective_ops, box, p);
    std::thread c(make_m2l_on_chunk<Box>, box_chunks[2], box_chunks[3], num_of_efective_ops, box, p);
    std::thread d(make_m2l_on_chunk<Box>, box_chunks[3], box_chunks[4], num_of_efective_ops, box, p);
    std::thread e(make_m2l_on_chunk<Box>, box_chunks[4], box_chunks[5], num_of_efective_ops, box, p);
    std::thread f(make_m2l_on_chunk<Box>, box_chunks[5], box_chunks[6], num_of_efective_ops, box, p);
    std::thread g(make_m2l_on_chunk<Box>, box_chunks[6], box_chunks[7], num_of_efective_ops, box, p);
    std::thread h(make_m2l_on_chunk<Box>, box_chunks[7], box_chunks[8], num_of_efective_ops, box, p);

    a.join();
    b.join();
    c.join();
    d.join();
    e.join();
    f.join();
    g.join();
    h.join();


    __AoS2SoA_mu__(box,muSoA,num_boxes_tree,p1xp2_2,priority_streams[current_priority_stream]);
    */
}

void fmm_algorithm::m2l_impl_cpu(){

    typedef typename CoeffMatrix::complex_type complex_type;

    __SoA2AoS_omega__(box, omegaSoA, num_boxes_tree, p1xp2_2,priority_streams[current_priority_stream]);
    Device::devSync();

    for (size_t d = 1; d <= depth; ++d)
    {
        size_t depth_offset = boxes_above_depth(d);
        size_t num_boxes_on_this_level = boxes_on_depth(d);

        for (ssize_t l = 0; l <= p; ++l)
        {
            for (ssize_t m = 0; m <= l; ++m)
            {
                for (size_t i = 0; i < num_of_efective_ops; ++i)
                {
                    for (size_t boxid = depth_offset; boxid < num_boxes_on_this_level + depth_offset; ++boxid)
                    {
                        CoeffMatrix* B = box[boxid].b_operators[i];

                        complex_type mu_l_m(0.);

                        for (ssize_t j = 0; j <= p; j+=2 )
                        {
                            for (ssize_t k = -j; k < 0; ++k)
                            {
                                mu_l_m += B->get(j + l, k + m) * box[boxid].omega->get(j, k);
                            }

                            for (ssize_t k = 0; k <= j; ++k)
                            {
                                mu_l_m += B->get(j + l, k + m) * box[boxid].omega->get_upper(j, k);
                            }
                        }

                        for (ssize_t j = 1; j <= p; j+=2)
                        {
                            for (ssize_t k = -j; k < 0; ++k)
                            {
                                mu_l_m -= B->get(j + l, k + m) * box[boxid].omega->get(j, k);
                            }

                            for (ssize_t k = 0; k <= j; ++k)
                            {
                                mu_l_m -= B->get(j + l, k + m) * box[boxid].omega->get_upper(j, k);
                            }
                        }

                        box[boxid].b_targets[i]->operator ()(l, m) += mu_l_m;
                    }
                }
            }
        }
    }

    __AoS2SoA_mu__(box, muSoA, num_boxes_tree, p1xp2_2,priority_streams[current_priority_stream]);
}

}//namespace end
