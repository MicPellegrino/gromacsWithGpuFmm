#ifndef _BK_FMM_kernel_pairwise_reducer_hpp_
#define _BK_FMM_kernel_pairwise_reducer_hpp_

#include "../../../xyz.hpp"
#include "../../../cuda_keywords.hpp"

namespace fmm_out {
    namespace pairwise {

        template <typename Real, typename Result>
        struct reducer
        {
            typedef size_t size_type;
            typedef gmx_gpu_fmm::XYZ<Real> Real3;

            Result &result;
            size_type i;
            Real potential;
            Real3 field;

            CUDA
            reducer(Result & result, size_type i)
                : result(result), i(i), potential(0.), field(0., 0., 0.)
            { }

            CUDA
            void operator () (Real p, const Real3 & f)
            {
                potential += p;
                field += f;
            }

            CUDA  ~reducer()
            {
                result.reduce_pf(i, potential, field);
            }
        };

    }  // namespace pairwise
}  // namespace fmm

#endif
// vim: et:ts=4:sw=4
