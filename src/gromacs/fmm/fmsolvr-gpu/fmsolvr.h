#ifndef GMX_FMM_FMSOLVR_GPU_H
#define GMX_FMM_FMSOLVR_GPU_H

#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/fmm/fmm-impl.h"
#include "gromacs/fmm/fmm-impl-helpers.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "global_functions.hpp"

namespace gmx_gpu_fmm
{

//forward declaration
class fmm_algorithm;

class Gpu_fmm
{

public:

    const size_t n_parts;
    const size_t p;
    const ssize_t ws;
    const size_t depth;
    const bool open_boundary_conditions;
    const bool dipole_compensation;
    const bool fmm_sparse;
    std::vector<std::vector<int> > excl_;
    size_t step;
    fmm_algorithm* fm_cuda;
    bool calc;

    typedef gmx::ArrayRef<gmx::RVec> Forcesarray;

    Gpu_fmm(int depth, int p_order, const size_t n, std::vector<std::vector<int> > &excl, int bound);

    void init(const matrix &box, const rvec *x, const real *q);

    void run(const gmx_bool            bqChanged,
             const real               *q,
             const rvec               *x,
             nonbonded_verlet_t       *nbv,
             const gmx_bool            calc_energies,
             const gmx_bool            gmx_neighbor_search,
             const matrix              &box,
             real one4pieps0 = ONE_4PI_EPS0);

    void getForcesCpu(const Forcesarray  &f_fmm,
                   double             *coulombEnergy,
                   real one4pieps0 = ONE_4PI_EPS0);

    void getForcesGpu(nonbonded_verlet_t *nbv);

    void cuda_settings();

    ~Gpu_fmm();
};

} // namespace gmx_gpu_fmm

#endif // GMX_FMM_FMSOLVR_GPU_H
