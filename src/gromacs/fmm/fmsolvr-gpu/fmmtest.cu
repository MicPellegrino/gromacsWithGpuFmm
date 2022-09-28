#include "architecture.hpp"
#include "global_functions.hpp"
#include "input_output.hpp"
#include "data_type.hpp"
#include "testdata.hpp"
#include "fmm.hpp"
#include "ioadapter.hpp"
#include "cuda_lib.hpp"
#include <vector>
#include <omp.h>
#include <unistd.h>

using namespace gmx_gpu_fmm;
/*!
 * \brief The class for executung the FMM and presetting the main simulation parameter.
 */
class Fmm_runner
{

public:

    //! basic data types
    typedef typename fmm_algorithm::Real Real;
    typedef typename fmm_algorithm::Real33 Real33;

    //! choose different p2p versions (0,1) and split between host and device (0,13) - experimental
    int p2p_host_device, p2p_version;
    //! multipole order \p p, depth of the octree \p depth, number of test steps \p,
    const size_t p, depth, steps, dataset_id;
    const bool print_energy, near_field, open_boundary_conditions, dipole_compensation, fmm_sparse, run_benchmark;
    const ssize_t ws;
    fmm_algorithm* fmm;
    Testdata<Real>* testdata;
    int deviceCount, current_device;

    Fmm_runner() :  print_energy(get_env_int(1, "PRINT_ENERGY") != 0),
        p(get_env_int(10, "MULTIPOLEORDER")),
        ws(get_env_int(1, "WS")),
        depth(get_env_int(3, "DEPTH")),
        steps(get_env_int(1, "STEPS")),
        near_field(get_env_int(1, "P2P")),
        dataset_id(get_env_int(4, "DATASET")),
        open_boundary_conditions(get_env_int(0, "OPENBOUNDARY") != 0),
        dipole_compensation(get_env_int(1, "DIPOLE_COMPENSATION") != 0),
        p2p_host_device(get_env_int(13, "P2P_HOST_DEVICE")),
        p2p_version(get_env_int(1, "P2P_VERSION")),
        fmm_sparse(get_env_int(0, "FMM_SPARSE") != 0),
        run_benchmark(get_env_int(0, "RUNBENCHMARK") != 0)
    {
        //cuda settings
        //const char* DEV = "CUDA_VISIBLE_DEVICES=-1";
        //putenv(const_cast<char*>(DEV));

        deviceCount = 0;
        CUDA_SAFE_CALL(cudaGetDeviceCount(&deviceCount));
        if(!run_benchmark)
        {
            printf("Detected %d CUDA Capable device(s)\n", deviceCount);
        }


        for (int dev = 0; dev < deviceCount; ++dev)
        {
            cudaSetDevice(dev);
            cudaDeviceProp deviceProp;
            cudaGetDeviceProperties(&deviceProp, dev);
            if(!run_benchmark)
                printf("Device %d : %s\n", dev, deviceProp.name);
        }

        current_device = 0;
    }

    void init(Testdata<Real>* test_data)
    {
        testdata = test_data;

        CUDA_SAFE_CALL(cudaSetDevice(current_device));
        //const char* CUDA_MAX_CON = "CUDA_DEVICE_MAX_CONNECTIONS=32";
        //putenv(const_cast<char*>(CUDA_MAX_CON));
        //cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
        //cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

        fmm = new fmm_algorithm(p, ws, depth, open_boundary_conditions, dipole_compensation, fmm_sparse);

        fmm->exec_time.init("ALLOC DATA");
        fmm->exec_time.init("PREPARE M2M");
        fmm->exec_time.init("PREPARE M2L");
        fmm->exec_time.init("PREPARE L2L");
        fmm->exec_time.init("PREPARE P2P");
        fmm->exec_time.init("COPY PARTICLES");
        fmm->exec_time.init("DISTRIBUTE PARTICLES");
        fmm->exec_time.init("PREP DATASTRUCTURES");
        fmm->exec_time.init("DIPOLE COMPENSATION");
        fmm->exec_time.init("P2P");
        fmm->exec_time.init("P2P_HOST");
        fmm->exec_time.init("P2M");
        fmm->exec_time.init("M2M");
        fmm->exec_time.init("M2L");
        fmm->exec_time.init("LATTICE");
        fmm->exec_time.init("L2L");
        fmm->exec_time.init("FORCES");
        CUDA_CHECK_ERROR();
        fmm->init_io(testdata->n());
        fmm->update_fmm_boxscale(testdata->abc);
        fmm->alloc_exclusions(testdata->excl);
        fmm->alloc_data();
        CUDA_CHECK_ERROR();
        fmm->set_expansion_points();
        fmm->prepare_m2m();
        fmm->prepare_m2l();
        fmm->prepare_l2l();
        CUDA_CHECK_ERROR();
        fake_AoS<Real> particle_data(testdata->x(), testdata->y(), testdata->z(), testdata->q(), testdata->n());
        fmm->copy_particles(particle_data);
        //param for memory allocated memory size prefactor
        fmm->alloc_box_particle_memory(5.0);
        fmm->prepare_p2p(p2p_version);
        CUDA_CHECK_ERROR();

        if(!run_benchmark)
        {
            printf("number of steps %d\n", steps);
            fmm->printer.print_header();
            fmm->printer.print("ALLOC DATA");
            fmm->printer.print("PREPARE M2M");
            fmm->printer.print("PREPARE M2L");
            fmm->printer.print("PREPARE L2L");
            fmm->printer.print("PREPARE P2P");
        }

        fmm->printer.setup_steps(steps);

    }

    void run()
    {
        //CUDA_SAFE_CALL(cudaSetDevice(current_device));
        double energy_av = 0.0;
        double energy;
        //times of kernels
        fmm->exec_time.init(steps);
        //time to solution
        Timer_ tts;
        int i = 0;
        
        tts.start();
        for (;i < steps; i++)
        {

            CUDA_CHECK_ERROR();
            fake_AoS<Real> particle_data(testdata->x(), testdata->y(), testdata->z(), testdata->q(), testdata->n());
            CUDA_CHECK_ERROR();
            fmm->update_fmm_boxscale(testdata->abc);

            CUDA_CHECK_ERROR();
            fmm->copy_particles(particle_data);
            CUDA_CHECK_ERROR();
            fmm->distribute_particles();
            CUDA_CHECK_ERROR();
            fmm->prepare_data();
            CUDA_CHECK_ERROR();
            if(near_field)
            {
                fmm->p2p(true,13,p2p_version);
            }
            CUDA_CHECK_ERROR();
            fmm->p2m();
            CUDA_CHECK_ERROR();
            fmm->m2m();
            CUDA_CHECK_ERROR();
            fmm->m2l();
            CUDA_CHECK_ERROR();
            fmm->compute_dipole_compensating_particles();
            CUDA_CHECK_ERROR();
            fmm->lattice();
            CUDA_CHECK_ERROR();
            fmm->l2l();
            CUDA_CHECK_ERROR();
            fmm->forces();
            CUDA_CHECK_ERROR();

            //t.join();
            //fmm->gather_results();

            //printf("-----------------------------------------step %d\n",(int)fmm->exec_time.actual_step);
            if(!run_benchmark)
            {
                fmm->printer.print("COPY PARTICLES");
                fmm->printer.print("DISTRIBUTE PARTICLES");
                fmm->printer.print("PREP DATASTRUCTURES");
                fmm->printer.print("DIPOLE COMPENSATION");
                fmm->printer.print("P2P_HOST");
                fmm->printer.print("P2P");
                fmm->printer.print("P2M");
                fmm->printer.print("M2M");
                fmm->printer.print("M2L");
                fmm->printer.print("LATTICE");
                fmm->printer.print("L2L");
                fmm->printer.print("FORCES");
            }


            energy = fmm->energy();
            CUDA_CHECK_ERROR();
            //fmm->dump_result_in_orig_order();
            energy_av += energy;

            if(!run_benchmark)
                fmm->printer.print_walltime(0);
            fmm->exec_time.add_step();
            CUDA_CHECK_ERROR();
            //printf("Energy host   %.20e\n", energy);
            //usleep(100);
        }
        tts.stop();

        if(print_energy)
        {
            Device<Real>::devSync();
            double p2p_reference;
            double fmm_reference;
            energy_av /=((double)fmm->exec_time.steps);


            if(open_boundary_conditions)
            {
                p2p_reference = testdata->reference_energies["open_reference"];
                fmm_reference = testdata->reference_energies["open"];
            }
            else
            {
                if(dipole_compensation)
                {
                    p2p_reference = testdata->reference_energies["dipole_correction_reference"];
                    fmm_reference = testdata->reference_energies["dipole_correction"];
                }
                else
                {
                    p2p_reference = testdata->reference_energies["periodic_reference"];
                    fmm_reference = testdata->reference_energies["periodic"];
                }
            }

            printf("Energy host   %.20e\n", energy);
            printf("Energy error  %.20e\n", (energy - p2p_reference) / p2p_reference);
            printf("FMM    error  %.20e\n", (energy - fmm_reference) / fmm_reference );
            typedef typename fmm_algorithm::Real3 Real3;
            //Real3 force = fmm->force_l2_norm();

            //parameter 0 prints all forces, do not use for large system
            //fm_cuda->force_dump(0);

            //fmm->force_dump(1);
            //printf("Forcex host    %.20e\n",force.x);
            //fmm->force_dump(2);
            //printf("Forcey host    %.20e\n",force.y);
            //fmm->force_dump(3);
            //printf("Forcez host    %.20e\n",force.z);

            //fm_cuda->energy_dump();
            //printf("Energy host oo %.20e\n",fm_cuda->energy_orig_order());
        }

        if(run_benchmark)
        {
            fmm->printer.print_results(fmm->io->n, fmm->depth, fmm->p, "", tts.get_time()/steps);
            printf("\n");
        }
        else
        {
            printf("\n");
            printf("av tts = %f (%f ns/day)\n", tts.get_time()/steps, 24*60*60/((tts.get_time()/steps)*250000));
        }

    }

    void exit()
    {
        fmm->~fmm_algorithm();
    }
};

int main(int argc, char ** argv)
{
    cudaDeviceReset();
    Fmm_runner* fmm_runner = new Fmm_runner();
    typedef typename Fmm_runner::Real Real;
    Testdata<Real>* testdata = new Testdata<Real>();
    if(!testdata->init_values(fmm_runner->dataset_id, fmm_runner->open_boundary_conditions))
    {
        printf("failed to create testdata\n");
        return 1;
    }

    int inits = get_env_int(1, "INITS");
    //printf("INITS %d\n",inits);
    for(int i = 0; i < inits; ++i)
    {
        //printf("INIT!!!!!!!!!!!!!!!! %d\n",i);
        alloc_counter::allocs = 0;
        fmm_runner->init(testdata);
        //printf ("allocated %d \n",alloc_counter::get());
        fmm_runner->run();
        fmm_runner->exit();
        //printf ("freed %d \n",alloc_counter::get());
        //usleep(10000000);
    }

    delete testdata;
    delete fmm_runner;
    return 0;
}











