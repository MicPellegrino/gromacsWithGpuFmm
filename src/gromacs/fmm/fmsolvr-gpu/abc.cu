#include "abc.hpp"

namespace gmx_gpu_fmm{

template <typename T>
CUDA
ABC<T>::ABC(const value_type & a, const value_type & b, const value_type & c)
    : a(a), b(b), c(c)
{}

template <typename T>
CUDA
ABC<T>::ABC(){}

template <typename T> template <typename Tother>
CUDA
ABC<T>::ABC(const Tother & a, const Tother & b, const Tother & c)
    : a(value_type(a)), b(value_type(b)), c(value_type(c))
{ }

template <typename T>  template <typename ABCother>
CUDA
ABC<T>::ABC(const ABCother & abc)
{
  a = abc.a;
  b = abc.b;
  c = abc.c;
}

template <typename T>
CUDA
ABC<T>& ABC<T>::operator *=(const Real & scale)
{
  a *= scale;
  b *= scale;
  c *= scale;

  return *this;
}

template <typename T>
CUDA
typename ABC<T>::value_type ABC<T>::operator * (const value_type & xyz) const
{
    return value_type(
            a.x * xyz.x + b.x * xyz.y + c.x * xyz.z,
            a.y * xyz.x + b.y * xyz.y + c.y * xyz.z,
            a.z * xyz.x + b.z * xyz.y + c.z * xyz.z);
}

template <typename T>
CUDA
ABC<T> ABC<T>::half() const
{
    ABC<T> abc(*this);

    abc.a *= 0.5;
    abc.b *= 0.5;
    abc.c *= 0.5;

    return abc;
}

template <typename T>
std::ostream & operator << (std::ostream & os, const ABC<T> & abc)
{
    os << "[" << abc.a << ", " << abc.b << ", " << abc.c << "]";
    return os;
}

template class ABC<XYZ<float> >;
template class ABC<XYZ<double> >;

}//namespace end
