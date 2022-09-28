#ifndef GMX_FMM_FMSOLVER_GPU_H
#define GMX_FMM_FMSOLVER_GPU_H

#include "config.h"

#include "gromacs/fmm/fmm-impl.h"
#include "fmsolvr-gpu/fmsolvr.h"

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"

struct gmx_mtop_t;
struct t_commrec;
struct t_inputrec;

namespace gmx_gpu_fmm
{
class Gpu_fmm;
}

namespace gmx
{
using namespace gmx_gpu_fmm;

class FmmImpl final : public IForceProvider
{
public:
    FmmImpl(const FmmInputParameters *fmmParameters, const ForceProviderInitOptions *options);

    // This is the main FMM routine that wraps calls to lower level routines
    void calculateForces(const ForceProviderInput &forceProviderInput,
                         ForceProviderOutput      *forceProviderOutput) override;

    int                                   n_;                  // Number of atoms
    int                                   ePBC_;               // Type of periodic boundary conditions

private:
    std::vector< std::vector<int> >    excl_;
    std::unique_ptr<gmx_gpu_fmm::Gpu_fmm> gpu_fmm_;

};

}      // namespace gmx

#endif // GMX_FMM_FMSOLVER_GPU_H
