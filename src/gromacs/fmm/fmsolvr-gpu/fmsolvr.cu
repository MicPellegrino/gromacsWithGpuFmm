//#define CUDADEBUG
//#define THRUST_DEBUG 0

#include <stdio.h>
#include "fmsolvr.h"
#include <omp.h>
#include "cuda_lib.hpp"
#include "timer.hpp"
#include "global_functions.hpp"
#include "input_output.hpp"
#include "data_type.hpp"
#include "testdata.hpp"
#include "fmm.hpp"
#include "ioadapter.hpp"

namespace gmx_gpu_fmm{

#ifdef GMX_FMM_DOUBLE
    __global__
    void
    __cast_coords(real *input, REAL* output, int n_elements)
    {
        int i = blockIdx.x * blockDim.x + threadIdx.x;
        if(i >= n_elements)
            return;
        output[i] = input[i];
    }
#endif


Gpu_fmm::Gpu_fmm(int d, int p_order, const size_t n, std::vector<std::vector<int> > &excl, int bound)
    : depth(get_env_int(d != -1 ? d : 3, "DEPTH")),
      p(get_env_int(p_order != -1 ? p_order : 18, "MULTIPOLEORDER")),
      n_parts(n),
      excl_(excl),
      ws(get_env_int(1, "WS")),
      open_boundary_conditions(bound != 0),
      dipole_compensation(get_env_int(1, "DIPOLE_COMPENSATION") != 0),
      fmm_sparse(get_env_int(0, "FMM_SPARSE") != 0),
      calc(true),
      step(0)

{
#ifdef GMX_FMM_DOUBLE
    printf("Running FMM in double precision !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
#endif
    //cudaSetDevice(0);
    printf("DEPTH = %lu MULTIPOLEORDER = %lu\n",depth, p);
    fm_cuda = new fmm_algorithm(p, ws, depth, open_boundary_conditions, dipole_compensation, fmm_sparse);
}

Gpu_fmm::~Gpu_fmm()
{
    delete fm_cuda;
}

void Gpu_fmm::init(const matrix &box, const rvec* x, const real* q)
{

    if(step==0)
    {
        cuda_settings();
        cudaDeviceSynchronize();
        typedef typename fmm_algorithm::Real3 Real3;
        typedef typename fmm_algorithm::Real33 Real33;

        Real3 a = Real3(box[0][0],box[0][1],box[0][2]);
        Real3 b = Real3(box[1][0],box[1][1],box[1][2]);
        Real3 c = Real3(box[2][0],box[2][1],box[2][2]);
        Real33 abc(a,b,c);

        fm_cuda->init_io(n_parts);
        fm_cuda->update_fmm_boxscale(abc);

        fm_cuda->alloc_exclusions(excl_);
        fm_cuda->alloc_data();

        fm_cuda->set_expansion_points();
        fm_cuda->prepare_m2m();
        fm_cuda->prepare_m2l();
        fm_cuda->prepare_l2l();

        typedef typename fmm_algorithm::Real Real;
        fake_AoS_GMX<Real, real> particle_data(x, q, fm_cuda->io->n);
        fm_cuda->copy_particles(particle_data);

        fm_cuda->alloc_box_particle_memory(5.0);
        fm_cuda->prepare_p2p();

    }
    step++;
}

void Gpu_fmm::run(const gmx_bool            bqChanged,
                  const real               *q,
                  const rvec               *x,
                  nonbonded_verlet_t       *nbv,
                  const gmx_bool            calc_energies,
                  const gmx_bool            gmx_neighbor_search,
                  const matrix              &box,
                  real one4pieps0)

{

    if(step == 1)
    {
        fm_cuda->lj_comb = reinterpret_cast<REAL2*>(Nbnxm::gpu_get_lj_comb(nbv->gpu_nbv));
    }

    typedef typename fmm_algorithm::Real Real;
    typedef typename fmm_algorithm::Real3 Real3;
    typedef typename fmm_algorithm::Real4 Real4;
    typedef typename fmm_algorithm::Real33 Real33;
    fm_cuda->current_stream = 0;
    fm_cuda->current_priority_stream = 0;
    cudaEvent_t* gmx_h2d_ready_event = reinterpret_cast<cudaEvent_t*>(Nbnxm::gpu_get_gmx_h2d_done_event(nbv->gpu_nbv));
    if(gmx_neighbor_search)
    {
        fake_AoS_GMX<Real, real> particle_data(x, q, fm_cuda->io->n);
        fm_cuda->copy_particles(particle_data);
    }
    else
    {
#ifndef GMX_FMM_DOUBLE
        REAL* coordinates = reinterpret_cast<REAL*>(Nbnxm::nbnxn_get_gpu_xrvec(nbv->gpu_nbv));
#else
        typedef typename fmm_algorithm::io_type::potential_vector_type REALVEC;
        REALVEC fmm_coords(fm_cuda->io->n*3,0.0);
        REAL* coordinates = &fmm_coords[0];

        real* gmx_coords = reinterpret_cast<real*>(Nbnxm::nbnxn_get_gpu_xrvec(nbv->gpu_nbv));

        dim3 block(512,1,1);
        dim3 grid((fm_cuda->io->n*3 - 1)/block.x + 1,1,1);
        __cast_coords<<<grid,block>>>(gmx_coords, coordinates, fm_cuda->io->n*3);
        fmm_algorithm::Device::devSync();

#endif
        fm_cuda->gmx_copy_particles_buffer(coordinates, gmx_h2d_ready_event);
    }
#ifndef GMX_FMM_DOUBLE
    const Real3* a = reinterpret_cast<const Real3*>(&box[0]);
    const Real3* b = reinterpret_cast<const Real3*>(&box[1]);
    const Real3* c = reinterpret_cast<const Real3*>(&box[2]);
#else
    const Real3 A(box[0][0], box[0][1], box[0][2]);
    const Real3 B(box[1][0], box[1][1], box[1][2]);
    const Real3 C(box[2][0], box[2][1], box[2][2]);
    const Real3* a = &A;
    const Real3* b = &B;
    const Real3* c = &C;
#endif

    Real33 abc(*a,*b,*c);
    fm_cuda->update_fmm_boxscale(abc);

    fm_cuda->distribute_particles(gmx_neighbor_search);
    fm_cuda->prepare_data();
    fm_cuda->p2p(calc_energies,13,1);
    fm_cuda->p2m();
    fm_cuda->m2m();
    fm_cuda->m2l();
    fm_cuda->compute_dipole_compensating_particles();
    fm_cuda->lattice();
    fm_cuda->l2l();
    fm_cuda->forces(one4pieps0);
}

void Gpu_fmm::getForcesCpu(const Forcesarray   &f_fmm,
                        double              *coulombEnergy,
                        real                 one4pieps0)
{
    *coulombEnergy = fm_cuda->energy(one4pieps0);
    size_t num_particles = fm_cuda->io->n;
    //omp_set_num_threads(6);
#pragma omp parallel for
    for(size_t i = 0; i < num_particles; ++i)
    {

#ifndef GMX_FMM_DOUBLE
        gmx::RVec* force = reinterpret_cast<gmx::RVec*>(&fm_cuda->io->forces_orig_order_host[i]);
#else
        gmx::RVec f(fm_cuda->io->forces_orig_order_host[i].x, fm_cuda->io->forces_orig_order_host[i].y, fm_cuda->io->forces_orig_order_host[i].z);
        gmx::RVec *force = &f;
#endif
        //printf("host gmx %d %e %e %e\n",i, f_fmm[i][0], f_fmm[i][1], f_fmm[i][2]);
        f_fmm[i] += *force;
        //printf("host fmm %d %e %e %e\n",i, fm_cuda->io->forces_orig_order_host[i].x, fm_cuda->io->forces_orig_order_host[i].y, fm_cuda->io->forces_orig_order_host[i].z);
    }
}

void Gpu_fmm::getForcesGpu(nonbonded_verlet_t *nbv)
{
    cudaStream_t* gmx_sync_stream = reinterpret_cast<cudaStream_t*>(Nbnxm::gpu_get_command_stream(nbv->gpu_nbv, Nbnxm::InteractionLocality::Local));
    float3* gmx_forces = reinterpret_cast<float3*>(nbv->getDeviceForces());
    fm_cuda->get_forces_in_orig_order(gmx_forces, gmx_sync_stream);
}

void Gpu_fmm::cuda_settings()
{
    //cuda settings
    //cudaDeviceSynchronize();
    //cudaDeviceReset();
    //cudaSetDevice(0);
    //const char* DEV = "CUDA_VISIBLE_DEVICES=0,1";
    //const char* CUDA_MAX_CON = "CUDA_DEVICE_MAX_CONNECTIONS=32";
    //putenv(const_cast<char*>(DEV));
    //putenv(const_cast<char*>(CUDA_MAX_CON));
    //cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
    //cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
}

}//namespace end






