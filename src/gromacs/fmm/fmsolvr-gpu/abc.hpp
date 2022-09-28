#ifndef ABC_HPP
#define ABC_HPP

#include <ostream>
#include "cuda_keywords.hpp"
#include "xyz.hpp"

namespace gmx_gpu_fmm{
/*!
 *  \brief Abstracts the size and the form of the simulation box
 *  \tparam T Datatype
 */
template <typename T>
class ABC {

public:

    typedef T value_type;
    value_type a,b,c;
    typedef typename T::value_type Real;

    CUDA
    ABC(const value_type& a, const value_type& b, const value_type& c);

    CUDA
    ABC();

    template <typename Tother>
    CUDA
    ABC(const Tother& a, const Tother& b, const Tother& c);

    template <typename ABCother>
    CUDA
    ABC(const ABCother& abc);

    CUDA
    ABC &operator *=(const Real& scale);

    CUDA
    value_type operator * (const value_type& xyz) const;

    CUDA
    ABC half() const;
};


}//namespace end
#endif // ABC_HPP
