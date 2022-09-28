#ifndef CUDA_LIB2_HPP
#define CUDA_LIB2_HPP

#include "cuda_lib.hpp"
#include "cuda_keywords.hpp"
#include "box.hpp"
#include "multipole.hpp"
#include "xyz.hpp"

namespace gmx_gpu_fmm{

__global__
void __dummy(){};

__global__
void __dummy_async(){};

__global__
void __dummy_p2p(){};

void start_dummy_kernel(int grid, int block){

    __dummy<<<grid, block, 0, 0>>>();
}

void start_dummy_p2p(int grid, int block){

    __dummy_p2p<<<grid, block, 0, 0>>>();
}

void start_async_dummy_kernel(int grid, int block, cudaStream_t &stream){

    __dummy_async<<<grid, block, 0, stream>>>();
}

template <typename Real3>
extern Real3 get_dipole_from_device(Real3* dipole_device, cudaStream_t &stream)
{
    Real3 dp;
    cudaMemcpyAsync(&dp, dipole_device, sizeof(Real3), cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);
    return dp;
}

template <typename CoefficientMatrixSoA, typename Box>
__global__
void __AoS2SoA_omega(Box *box, CoefficientMatrixSoA* omegaSoA, size_t num_boxes_tree)
{
    size_t box_id = blockIdx.x * blockDim.x + threadIdx.x;
    if(box_id >= num_boxes_tree)
        return;
    size_t index = blockIdx.y;
    (omegaSoA->get_lin_SoA_val(index, box_id)) = box[box_id].omega->get_vectorized(index);
}

template <typename CoefficientMatrixSoA, typename Box>
void __AoS2SoA_omega__(Box *box, CoefficientMatrixSoA *omegaSoA, size_t num_boxes_tree, size_t p1xp2_2, cudaStream_t &stream)
{
    dim3 SoA_block(512,1,1);
    dim3 SoA_grid((num_boxes_tree-1)/SoA_block.x+1,p1xp2_2,1);
    __AoS2SoA_omega<<<SoA_grid, SoA_block,0,stream>>>(box, omegaSoA, num_boxes_tree);
}

template <typename CoefficientMatrixSoA, typename Box>
__global__ void
__SoA2AoS_omega(Box *box, CoefficientMatrixSoA* omegaSoA, size_t num_boxes_tree)
{
    size_t box_id = blockIdx.x * blockDim.x + threadIdx.x;
    if(box_id >= num_boxes_tree)
        return;
    size_t index = blockIdx.y;

    (*box[box_id].omega)(index) =  *(omegaSoA->get_lin_SoA_ptr(index,box_id));
}

template <typename CoefficientMatrixSoA, typename Box>
void __SoA2AoS_omega__(Box *box, CoefficientMatrixSoA* omegaSoA, size_t num_boxes_tree, size_t p1xp2_2, cudaStream_t &stream)
{
    dim3 SoA_block(512,1,1);
    dim3 SoA_grid((num_boxes_tree-1)/SoA_block.x+1,p1xp2_2,1);
    __SoA2AoS_omega<<<SoA_grid,SoA_block,0,stream>>>(box,omegaSoA,num_boxes_tree);
}

template <typename CoefficientMatrixSoA, typename Box>
__global__
void __AoS2SoA_mu(Box *box, CoefficientMatrixSoA* muSoA, size_t num_boxes_tree)
{
    size_t box_id = blockIdx.x * blockDim.x + threadIdx.x;
    if(box_id >= num_boxes_tree)
        return;
    size_t index = blockIdx.y;
    (muSoA->get_lin_SoA_val(index, box_id)) = box[box_id].mu->get_vectorized(index);
}

template <typename CoefficientMatrixSoA, typename Box>
void __AoS2SoA_mu__(Box *box, CoefficientMatrixSoA* muSoA, size_t num_boxes_tree, size_t p1xp2_2, cudaStream_t &stream)
{

    dim3 SoA_block(512,1,1);
    dim3 SoA_grid((num_boxes_tree-1)/SoA_block.x+1,p1xp2_2,1);
    __AoS2SoA_mu<<<SoA_grid,SoA_block,0,stream>>>(box,muSoA,num_boxes_tree);
}

template <typename CoefficientMatrixSoA, typename Box>
__global__ void
__AoS_addto_SoA_omega(Box *box, CoefficientMatrixSoA* omegaSoA, size_t num_boxes_tree)
{
    size_t box_id = blockIdx.x * blockDim.x + threadIdx.x;
    if(box_id >= num_boxes_tree)
        return;
    size_t index = blockIdx.y;
    (omegaSoA->get_lin_SoA_val(index, box_id)) += box[box_id].omega->get_vectorized(index);
}

template <typename CoefficientMatrixSoA, typename Box>
void __AoS_addto_SoA_omega__(Box *box, CoefficientMatrixSoA* omegaSoA, size_t num_boxes_tree, size_t p1xp2_2, cudaStream_t &stream)
{
    dim3 SoA_block(1,1,1);
    dim3 SoA_grid(1,p1xp2_2,1);
    __AoS_addto_SoA_omega<<<SoA_grid,SoA_block,0,stream>>>(box, omegaSoA, 1);
}

template <typename CoefficientMatrixSoA, typename Box>
__global__ void
__AoS_addto_SoA_mu(Box *box, CoefficientMatrixSoA* muSoA, size_t num_boxes_tree)
{
    size_t box_id = blockIdx.x * blockDim.x + threadIdx.x;
    if(box_id >= num_boxes_tree)
        return;
    size_t index = blockIdx.y;
    (muSoA->get_lin_SoA_val(index, box_id)) += box[box_id].mu->get_vectorized(index);
}

template <typename CoefficientMatrixSoA, typename Box>
void __AoS_addto_SoA_mu__(Box *box, CoefficientMatrixSoA* muSoA, size_t num_boxes_tree, size_t p1xp2_2, cudaStream_t &stream)
{
    dim3 SoA_block(1,1,1);
    dim3 SoA_grid(1,p1xp2_2,1);
    __AoS_addto_SoA_mu<<<SoA_grid,SoA_block,0,stream>>>(box, muSoA, 1);
}

template <typename CoefficientMatrixSoA, typename Box>
__global__
void __SoA2AoS_mu(Box *box, CoefficientMatrixSoA* muSoA, size_t num_boxes_tree)
{
    size_t box_id = blockIdx.x * blockDim.x + threadIdx.x;
    if(box_id >= num_boxes_tree)
        return;
    size_t index = blockIdx.y;

    (*box[box_id].mu)(index) =  *(muSoA->get_lin_SoA_ptr(index,box_id));
}

template <typename CoefficientMatrixSoA, typename Box>
void __SoA2AoS_mu__(Box *box, CoefficientMatrixSoA* muSoA, size_t num_boxes_tree, size_t p1xp2_2, cudaStream_t &stream)
{
    dim3 SoA_block(512,1,1);
    dim3 SoA_grid((num_boxes_tree-1)/SoA_block.x+1,p1xp2_2,1);
    __SoA2AoS_mu<<<SoA_grid,SoA_block,0,stream>>>(box,muSoA,num_boxes_tree);
}

typedef MultipoleCoefficientsUpperSoA<REAL, Device<REAL> > CoefficientMatrixSoA;
typedef Box<REAL, Device<REAL>> Boxtype;
typedef XYZ<REAL> Real3;

template void __AoS2SoA_omega__<CoefficientMatrixSoA, Boxtype>(Boxtype *box, CoefficientMatrixSoA* omegaSoA, size_t num_boxes_tree, size_t p1xp2_2, cudaStream_t &stream);
template void __SoA2AoS_omega__<CoefficientMatrixSoA, Boxtype>(Boxtype *box, CoefficientMatrixSoA* omegaSoA, size_t num_boxes_tree, size_t p1xp2_2, cudaStream_t &stream);
template void __AoS2SoA_mu__<CoefficientMatrixSoA, Boxtype>(Boxtype *box, CoefficientMatrixSoA* muSoA, size_t num_boxes_tree, size_t p1xp2_2, cudaStream_t &stream);
template void __SoA2AoS_mu__<CoefficientMatrixSoA, Boxtype>(Boxtype *box, CoefficientMatrixSoA* muSoA, size_t num_boxes_tree, size_t p1xp2_2, cudaStream_t &stream);
template void __AoS_addto_SoA_omega__<CoefficientMatrixSoA, Boxtype>(Boxtype *box, CoefficientMatrixSoA* omegaSoA, size_t num_boxes_tree, size_t p1xp2_2, cudaStream_t &stream);
template void __AoS_addto_SoA_mu__<CoefficientMatrixSoA, Boxtype>(Boxtype *box, CoefficientMatrixSoA* muSoA, size_t num_boxes_tree, size_t p1xp2_2, cudaStream_t &stream);
template Real3 get_dipole_from_device<Real3>(Real3* device_dipole, cudaStream_t &stream);

}//namespace end

#endif // CUDA_LIB_HPP
