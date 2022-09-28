#ifndef CUDA_LIB_HPP
#define CUDA_LIB_HPP

#include "data_type.hpp"
#include "cuda_keywords.hpp"
#include "xyz.hpp"

namespace gmx_gpu_fmm{

extern void start_dummy_kernel(int grid, int block);

extern void start_dummy_p2p(int grid, int block);

extern void start_async_dummy_kernel(int grid, int block, cudaStream_t &stream);

template <typename Real3>
extern Real3 get_dipole_from_device(Real3* dipole_device, cudaStream_t &stream);

template <typename CoefficientMatrixSoA, typename Box>
extern void __AoS2SoA_omega__(Box *box, CoefficientMatrixSoA* omegaSoA, size_t num_boxes_tree, size_t p1xp2_2, cudaStream_t &stream);

template <typename CoefficientMatrixSoA, typename Box>
extern void __SoA2AoS_omega__(Box *box, CoefficientMatrixSoA* omegaSoA, size_t num_boxes_tree, size_t p1xp2_2, cudaStream_t &stream);

template <typename CoefficientMatrixSoA, typename Box>
extern void __AoS2SoA_mu__(Box *box, CoefficientMatrixSoA* muSoA, size_t num_boxes_tree, size_t p1xp2_2, cudaStream_t &stream);

template <typename CoefficientMatrixSoA, typename Box>
extern void __SoA2AoS_mu__(Box *box, CoefficientMatrixSoA* muSoA, size_t num_boxes_tree, size_t p1xp2_2, cudaStream_t &stream);

template <typename CoefficientMatrixSoA, typename Box>
extern  void __AoS_addto_SoA_omega__(Box *box, CoefficientMatrixSoA* omegaSoA, size_t num_boxes_tree, size_t p1xp2_2, cudaStream_t &stream);

template <typename CoefficientMatrixSoA, typename Box>
extern void __AoS_addto_SoA_mu__(Box *box, CoefficientMatrixSoA* muSoA, size_t num_boxes_tree, size_t p1xp2_2, cudaStream_t &stream);

}//namespace end

/*
template <typename CoefficientMatrix>
__global__
void __check_for_nans(CoefficientMatrix **omega, size_t num_boxes_tree, int *found)
{

    typedef typename CoefficientMatrix::value_type complex_value;
    size_t box_id = blockIdx.x * blockDim.x + threadIdx.x;
    if(box_id >= num_boxes_tree)
        return;
    size_t index = blockIdx.y;
    complex_value val = omega[box_id]->get_vectorized(index);

    if(isnan(val.real()) || isinf(val.real()) || abs(val.real())< 1.0e-32)
    {
        //printf("%e ",val.real());
        omega[box_id]->operator()(index).real(0.0);
        *found = 1;
    }

    if(isnan(val.imag()) || isinf(val.imag()) || abs(val.imag())< 1.0e-32)
    {
        //printf("%e ",val.imag());
        omega[box_id]->operator()(index).imag(0.0);
        *found = 1;
    }
}
*/

#endif



