#include "fmm.hpp"
#include "cuda_LATTICE.hpp"
#include "cuda_lib.hpp"

namespace gmx_gpu_fmm{

template <typename CoefficientMatrix>
__global__
void __reset_omega(CoefficientMatrix **omega, size_t p1xp2_2)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i < p1xp2_2)
        omega[0]->flush(i);
}

void fmm_algorithm::lattice_impl(){

    const int num_of_streams = STREAMS;
    typedef typename CoeffMatrix::value_type complex_type;

    //wait for m2l
    for(size_t i = 0; i < STREAMS; ++i)
    {
        cudaStreamWaitEvent(priority_streams[current_priority_stream], priority_events[i], 0 );
    }

    size_t op_p1xx2 = ( (2*p+1) * (2*p+1) );

    if(dipole_compensation)
    {
        if(depth == 0)
        {
            __reset_omega<CoeffMatrix><<<(p1xp2_2-1)/512 + 1, 512, 0, priority_streams[current_priority_stream]>>>(omega, p1xp2_2);
        }

        int last_stream = current_priority_stream;
        __P2M_P2L_dipole_corr<<<1, 64, 0, priority_streams[current_priority_stream]>>>(&q0abc[0], omega, &fake_particles[0], mu, p, fake_particle_size);
        cudaEventRecord(priority_events[current_priority_stream],priority_streams[current_priority_stream]);
        current_priority_stream = (++current_priority_stream)%num_of_streams;

        cudaStreamWaitEvent(priority_streams[current_priority_stream], priority_events[last_stream], 0 );
        __AoS_addto_SoA_omega__(box, omegaSoA, 1, p1xp2_2, priority_streams[current_priority_stream]);
        cudaStreamWaitEvent(priority_streams[current_priority_stream], priority_events[last_stream], 0 );
        __AoS_addto_SoA_mu__(box, muSoA, 1, p1xp2_2, priority_streams[current_priority_stream]);
    }

    if (!open_boundary_conditions)
    {

#ifndef GMX_FMM_DOUBLE
        dim3 grid(1,1,1);
        dim3 block(p1+1,p1,1);
        __lattice<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type>
        <<<grid,block,(p1*p1+op_p1xx2)*sizeof(complex_type), priority_streams[current_priority_stream]>>>
        (box, omegaSoA, muSoA, Lattice, 0, num_boxes_tree, p, p1, p1xx2, op_p1xx2, lattice_rescale);
#else
        if(p<21)
        {
            dim3 grid(1,1,1);
            dim3 block(p1+1,p1,1);
            __lattice<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type>
            <<<grid,block,(p1*p1+op_p1xx2)*sizeof(complex_type), priority_streams[current_priority_stream]>>>
            (box, omegaSoA, muSoA, Lattice, 0, num_boxes_tree, p, p1, p1xx2, op_p1xx2, lattice_rescale);
        }
        else
        {
            dim3 griD(p1,1,1);
            dim3 blocK(p1,1,1);
            __lattice_no_shared<CoeffMatrix, CoeffMatrixSoA, Box, Real, Real3, complex_type>
            <<<griD,blocK,0,priority_streams[current_priority_stream]>>>
            (box, omegaSoA, muSoA, Lattice, 0, num_boxes_tree, p, p1, p1xx2, op_p1xx2, lattice_rescale);
        }
#endif
    }
}

}//namespace end
