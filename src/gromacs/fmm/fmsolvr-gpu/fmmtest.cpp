#define CUDADEBUG

#include "data_type.hpp"
#include "fmm_helper_functions.hpp"
#include "fmm_lambda_helper_functions.hpp"
#include "sim_handler.hpp"
#include "fmm_cpu.hpp"
#include "fmm_gpu.hpp"

#include "testdata.hpp"


int main(int argc, char ** argv)
{

    //cuda settings
    cudaDeviceReset();
    const char* DEV = "CUDA_VISIBLE_DEVICES=0,1";
    const char* CUDA_MAX_CON = "CUDA_DEVICE_MAX_CONNECTIONS=32";
    putenv(const_cast<char*>(DEV));
    putenv(const_cast<char*>(CUDA_MAX_CON));
    cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
    //cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

    cudaSetDevice(0);

    std::cout << std::scientific;
    //time to solution measurement
    fmmgpu::papiwrapper tts;
    tts.start();

    //input data -> fmm parameters
    const size_t first = get_env_int(1, "FIRST");
    const bool print_energy = get_env_int(1, "PRINT_ENERGY") != 0;

    const char* sim_type = getenv("SIM"); //COM,GPU,CPU
    size_t p = get_env_int(10, "MULTIPOLEORDER");
    const ssize_t ws = get_env_int(1, "WS");
    const size_t depth = get_env_int(4, "DEPTH");
    const size_t steps = get_env_int(1, "STEPS");
    const bool near_field = get_env_int(0, "P2P");
    const int dataset_id = get_env_int(4, "DATASET");
    const bool open_boundary_conditions = get_env_int(0, "OPENBOUNDARY") != 0;
    const bool dipole_compensation = get_env_int(1, "DIPOLE_COMPENSATION") != 0;

    Sim sim(sim_type);
    if(!sim.init())
        return 0;

    typedef gmx_gpu_fmm::input_output<no_threads> cpu_io_type;
    typedef gmx_gpu_fmm::input_output<cu> gpu_io_type;

    typedef gmx_gpu_fmm::fmm_algorithm<no_threads> fm_cpu;
    typedef gmx_gpu_fmm::fmm_algorithm<cu> fm_gpu;
    typedef typename fm_gpu::Real Real;

    fm_gpu * fm_cuda;
    //fm_cpu * fm;

    Testdata<Real>* testdata = new Testdata<Real>();

    if(sim.sim == SIM::GPU)
    {
        if(!testdata->init_values(dataset_id))
            return 0;

        gpu_io_type* gpu_in_and_out = new gpu_io_type(testdata->n());
        gpu_in_and_out->set_box_info(testdata->abc, testdata->reference_point, open_boundary_conditions);

        fm_cuda  = new fm_gpu(p, ws, depth, open_boundary_conditions, dipole_compensation, gpu_in_and_out->n);
        fm_cuda ->init_io(gpu_in_and_out);

        fake_aos4<Real> particle_data(testdata->x(), testdata->y(), testdata->z(), testdata->q(), testdata->n());
        gpu_in_and_out->set_particles(particle_data);

        fm_cuda->alloc_exclusions();

        fm_cuda->printer.setup(steps);

        fm_cuda->alloc_data();
        fm_cuda->set_expansion_points();
        fm_cuda->prepare_m2m();
        fm_cuda->prepare_m2l();
        fm_cuda->prepare_l2l();
        fm_cuda->prepare_p2p();

        if(first)
        {
            fm_cuda->printer.print_header();
            fm_cuda->printer.print("ALLOC DATA");
            fm_cuda->printer.print("PREPARE M2M");
            fm_cuda->printer.print("PREPARE M2L");
            fm_cuda->printer.print("PREPARE L2L");
            fm_cuda->printer.print("PREPARE P2P");
        }


        for (int i = 0; i < steps; i++)
        {
            fm_cuda->realloc_box_index_mem(i);
            fm_cuda->printer.set_av(i+1);

            fake_aos4<Real> particle_data(testdata->x(), testdata->y(), testdata->z(), testdata->q(), testdata->n());
            gpu_in_and_out->set_particles(particle_data);

            fm_cuda->distribute_particles();
            fm_cuda->printer.print("DISTRIBUTE PARTICLES");
            fm_cuda->prepare_data();
            fm_cuda->printer.print("COPY AND FLUSH");
            if(near_field)
            {
                fm_cuda->p2p();
                fm_cuda->printer.print("P2P");
            }

            fm_cuda->p2m();
            fm_cuda->printer.print("P2M");
            fm_cuda->m2m();
            fm_cuda->printer.print("M2M");
            fm_cuda->m2l();
            fm_cuda->printer.print("M2L");
            fm_cuda->lattice();
            fm_cuda->printer.print("LATTICE");
            fm_cuda->l2l();
            fm_cuda->printer.print("L2L");
            fm_cuda->forces();
            fm_cuda->printer.print("FORCES");

            fm_cuda->printer.print_walltime();

            if(print_energy)
            {
                //fm_cuda->energy_dump();
                printf("Energy host    %.20e\n",fm_cuda->energy());
                //printf("Energy host oo %.20e\n",fm_cuda->energy_orig_order());
            }
        }


        if(print_energy)
        {
            typedef typename fm_gpu::Real3 Real3;
            Real3 force = fm_cuda->force_l2_norm();

            //parameter 0 prints all forces, do not use for large system
            //fm_cuda->force_dump(0);

            //fm_cuda->force_dump(1);
            printf("Forcex host    %.20e\n",force.x);
            //fm_cuda->force_dump(2);
            printf("Forcey host    %.20e\n",force.y);
            //fm_cuda->force_dump(3);
            printf("Forcez host    %.20e\n",force.z);
        }

    }


    return 0;
}
