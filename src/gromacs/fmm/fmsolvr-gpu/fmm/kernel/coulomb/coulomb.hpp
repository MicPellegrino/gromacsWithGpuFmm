#ifndef _BK_FMM_kernel_coulomb_hpp_
#define _BK_FMM_kernel_coulomb_hpp_

#include "../../../xyz.hpp"
#include "../../../cuda_keywords.hpp"


namespace fmm_out {

    // generic Coulomb (or '1/r') kernel for the directed interaction between
    // two distinct particles (target, source)
    template <typename Real>
    struct coulomb_kernel
    {
        typedef gmx_gpu_fmm::XYZ<Real> Real3;

        template <typename Target, typename Source>
        CUDA
        std::pair<Real, Real3> operator () (const Target & target, const Source & source) const
        {
            Real3 x_i(target);
            Real3 x_j(source);
            Real q_j(source.s);

            Real3 diff = x_i - x_j;                 // (x_i - x_j)
            Real rlen = rcplength(diff);            // 1 / |x_i - x_j|
            Real q_j_rlen = q_j * rlen;             // q_j / |x_i - x_j|
            Real rlen2 = rlen * rlen;               // 1 / |x_i - x_j|^2
            Real q_j_rlen3 = q_j_rlen * rlen2;      // q_j / |x_i - x_j|^3
            Real3 f_i = diff * q_j_rlen3;           // (x_i - x_j) * q_j / |x_i - x_j|^3
            Real phi_i = q_j_rlen;
            return std::make_pair(phi_i, f_i);
        }
    };

}  // namespace fmm

#endif
// vim: et:ts=4:sw=4
