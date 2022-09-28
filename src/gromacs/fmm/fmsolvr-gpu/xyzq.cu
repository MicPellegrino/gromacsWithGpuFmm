#include "xyzq.hpp"

namespace gmx_gpu_fmm{

template <typename Real>
CUDA
XYZQ<Real>::XYZQ(): x(value_type()), y(value_type()), z(value_type()), q(value_type())
{}

template <typename Real>
CUDA
XYZQ<Real>::XYZQ(value_type x, value_type y, value_type z, value_type q): x(x), y(y), z(z), q(q)
{}

template class XYZQ<double>;
template class XYZQ<float>;

}//namespace end

