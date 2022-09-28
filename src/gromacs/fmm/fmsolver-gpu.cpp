#include "gmxpre.h"

#include "fmm.h"

#ifdef GMX_WITH_FMM
#include "fmm-impl-helpers.h"
#include "fmm-impl.h"
#include "fmsolver-gpu.h"
#include "fmsolvr-gpu/fmsolvr.h"

#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"

namespace gmx
{

/*! \brief
 *  Identification string that can be prepended to FMM-related output
 */
constexpr char fmmStr[] = "FMM:";

FmmImpl::FmmImpl(const FmmInputParameters *fmmParameters, const ForceProviderInitOptions *options)
    : n_(options->mtop_->natoms),
      ePBC_(options->ePBC_)
{
    fprintf(stderr, "%s Initializing GPU version of fast multipole solver\n", fmmStr);

    excl_    = makeExclusionList(options->mtop_);
    gpu_fmm_ = std::make_unique<gmx_gpu_fmm::Gpu_fmm>(fmmParameters->override_depth_, fmmParameters->override_order_, options->mtop_->natoms, excl_, ePBC_);

#if GMX_GPU != GMX_GPU_CUDA
    gmx_fatal(FARGS, "GPU FMM requested, but GROMACS was compiled without CUDA.\n");
#endif
}

void FmmImpl::calculateForces(const ForceProviderInput &forceProviderInput,
                              ForceProviderOutput *forceProviderOutput)
{
    if(forceProviderOutput == nullptr)
    {
        if(gpu_fmm_->calc)
        {
            gpu_fmm_->calc = false;
            real* q = forceProviderInput.mdatoms_.chargeA;
            gpu_fmm_->init(forceProviderInput.box_,as_rvec_array(forceProviderInput.x_.data()), q);
            gpu_fmm_->run(1, q, as_rvec_array(forceProviderInput.x_.data()), forceProviderInput.nb_, forceProviderInput.calculateEnergies_, forceProviderInput.doNeighborSearching_, forceProviderInput.box_);
        }
        else
        {
            gpu_fmm_->calc = true;
            if(!forceProviderInput.calculateEnergies_)
            {
                gpu_fmm_->getForcesGpu(forceProviderInput.nb_);
            }
        }
    }
    else
    {
        gpu_fmm_->calc = true;
        if(forceProviderInput.calculateEnergies_)
        {
            const ArrayRef<RVec> force = forceProviderOutput->forceWithVirial_.force_;
            double energy_fmm = 0.0;
            gpu_fmm_->getForcesCpu(force, &energy_fmm);
            // Add FMM Coulomb energy to GROMACS energies, for now to Coulomb-short range
            forceProviderOutput->enerd_.grpp.ener[egCOULSR][0] += energy_fmm;
        }
    }
}

std::unique_ptr<FmmImpl> createFmmImpl(const FmmInputParameters *fmmParameters, const ForceProviderInitOptions *options)
{
    return std::make_unique<FmmImpl>(fmmParameters, options);
}


} // namespace gmx


#endif // GMX_WITH_FMM
