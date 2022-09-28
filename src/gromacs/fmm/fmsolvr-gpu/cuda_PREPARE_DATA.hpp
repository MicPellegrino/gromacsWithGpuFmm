#ifndef CUDA_PREPARE_DATA_HPP
#define CUDA_PREPARE_DATA_HPP

namespace gmx_gpu_fmm{

template <typename CoefficientMatrixSoA, typename complex_type>
__global__
void __initSoA(CoefficientMatrixSoA* SoA, size_t num_boxes_tree, complex_type val)
{
    size_t box_id = blockIdx.x * blockDim.x + threadIdx.x;
    if(box_id >= num_boxes_tree)
        return;
    size_t index = blockIdx.y;
    (SoA->get_lin_SoA_val(index, box_id)) = val;
}

template <typename outputadapter, typename Real3>
__global__
void __flush_result(outputadapter *result, size_t n)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i<n)
    {
        result->vF_ptr[i] = Real3(0.,0.,0.);
        result->vPhi_ptr[i] = 0.;
    }
}

template <typename CoefficientMatrix>
__global__
void __flush_multipoles(CoefficientMatrix **omega, CoefficientMatrix **mu)
{
    size_t idx = blockIdx.x;
    size_t i = threadIdx.x;
    omega[idx]->flush(i);
    mu[idx]->flush(i);
}

template <typename CoefficientMatrix>
__global__
void __flush_multipoles2(CoefficientMatrix **omega, CoefficientMatrix **mu, size_t p1xp2_2, size_t num_boxes_tree)
{
    size_t idx = blockIdx.y * gridDim.z + blockIdx.z;
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i < p1xp2_2 && idx < num_boxes_tree)
    {
        omega[idx]->flush(i);
        mu[idx]->flush(i);
    }
}

template <typename Real4, typename Box>
__global__
void
__COPY_parts(size_t offset,
             size_t *box_particle_offset,
             Real4* ordered_particles,
             const Real4* input,
             size_t* orig_ids,
             Box* box,
             size_t* id_map,
             size_t *block_map,
             size_t *offset_map)
{
    //maps particles blocks to cuda treadblocks
    const size_t id = block_map[blockIdx.x*blockDim.y + threadIdx.y];
    const size_t i = (blockIdx.x*blockDim.y + threadIdx.y)*blockDim.x + threadIdx.x - offset_map[blockIdx.x*blockDim.y + threadIdx.y];
    size_t idx = offset + id;
    box[id_map[idx]].active = 1;
    if(i < box_particle_offset[id+1])
    {
        size_t orig_ptcl_index =  box[id_map[idx]].orig_ptcl_ids[i - box_particle_offset[id]];
        ordered_particles[i] = input[orig_ptcl_index];
        orig_ids[i] = orig_ptcl_index;
    }
}

}//namespace end

#endif // CUDA_PREPARE_DATA_HPP
