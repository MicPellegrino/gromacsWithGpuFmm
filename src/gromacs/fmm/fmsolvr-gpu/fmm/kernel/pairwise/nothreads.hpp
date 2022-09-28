#ifndef _BK_FMM_kernel_nothreads_hpp_
#define _BK_FMM_kernel_nothreads_hpp_

#include "../../../xyz.hpp"
#include "reducer.hpp"
#include "../../lib.hpp"

namespace fmm_out {
    namespace pairwise {

        template <typename Real, typename Kernel, typename Result, typename Target, typename Source>
        struct nothreads {
            typedef size_t size_type;
            typedef gmx_gpu_fmm::XYZ<Real> Real3;
            typedef gmx_gpu_fmm::XYZQ<Real> Real4;
            typedef reducer<Real, Result> reducer_type;

            nothreads(Result & result, const Target & target, const Source & source)
            {

            }

            void interaction_sequence(const Real3 & xyz_i, reducer_type & reducer_i, size_type j_begin, size_type j_end, const Source & source) const
            {

                for (size_type j = j_begin; j < j_end; ++j) {
                    Real4 xyzs_j = source.load_xyzs(j);
                    std::pair<Real, Real3> res = single_interaction(xyz_i, xyzs_j);  // returns a pair (potential, field)
                    reducer_i(res.first, res.second);
                }
            }

            // returns a pair (potential, field)
            std::pair<Real, Real3> single_interaction(const Real3 & xyz_i, const Real4 & xyzs_j) const
            {
                return kernel(xyz_i, xyzs_j);
            }

            void overlapping_interactions(size_type i_begin, size_type i_end, size_type j_begin, size_type j_end)
            {
                for (size_type i = i_begin; i < i_end; ++i)
                {
                    Real3 xyz_i = target.load_xyz(i);
                    reducer_type reducer_i(result, i);
                    interaction_sequence(xyz_i, reducer_i, j_begin, std::min(j_end, i), source);
                    interaction_sequence(xyz_i, reducer_i, std::max(j_begin, i + 1), j_end, source);
                }
            }

            void disjoint_interactions(size_type i_begin, size_type i_end, size_type j_begin, size_type j_end, Real3 offset = Real3(0.,0.,0.))
            {
                for (size_type i = i_begin; i < i_end; ++i)
                {
                    Real3 xyz_i = target.load_xyz(i);
                    reducer_type reducer_i(result, i);
                    interaction_sequence(xyz_i, reducer_i, j_begin, j_end, source);
                }
            }
        };

    }  // namespace pairwise
}  // namespace fmm

#endif
// vim: et:ts=4:sw=4
