#ifndef LATTICE_H
#define LATTICE_H

#include "xyz.hpp"
#include "abc.hpp"
#include "triangular_array.hpp"
#include "fmm_complex.hpp"
#include "multipole.hpp"
#include "multipole2local.hpp"

namespace gmx_gpu_fmm{

template <typename Real33, typename CoefficientMatrix>
extern size_t __attribute__((noinline)) LatticeOperator3D(
        const Real33 & abc,
        const size_t ws,
        CoefficientMatrix & lattice,
        const size_t p,
        const size_t maxlayers,
        const double maxerror = std::numeric_limits<double>::min());

}//namespace end

#endif // LATTICE_H
