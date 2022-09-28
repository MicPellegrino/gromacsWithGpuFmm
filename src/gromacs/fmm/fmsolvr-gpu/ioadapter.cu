#include "ioadapter.hpp"

namespace gmx_gpu_fmm{

template <typename T4Vector>
inputadapter<T4Vector>::inputadapter(const vector4_type & v4) : v4(v4), v4_ptr(&(v4[0]))
{}

template <typename T4Vector>
CUDA
inputadapter<T4Vector>::Real3 inputadapter<T4Vector>::load_xyz(size_t i) const
{
#ifdef __CUDA_ARCH__
    return Real3( vectorized_load( v4_ptr, i ) );
    //return Real3(v4_ptr[i]);
#else
    return Real3(v4[i]);
#endif
}

template <typename T4Vector>
CUDA
inputadapter<T4Vector>::Real4 inputadapter<T4Vector>::load_xyzs(size_t i) const
{
#ifdef __CUDA_ARCH__
    return v4_ptr[i];
#else
    return v4[i];
#endif
}

template <typename T4Vector>
CUDA
inputadapter<T4Vector>::Real inputadapter<T4Vector>::load_s(size_t i) const
{
#ifdef __CUDA_ARCH__
    return v4_ptr[i].s;
#else
    return v4[i].s;
#endif
}

template <typename T1Vector, typename T3Vector, typename Arch>
outputadapter<T1Vector, T3Vector, Arch>::outputadapter(vector1_type & vPhi, vector3_type & vF) : vPhi(vPhi), vF(vF),vPhi_ptr(&(vPhi[0])),vF_ptr(&vF[0])
{}

template <typename T1Vector, typename T3Vector, typename Arch>
outputadapter<T1Vector, T3Vector, Arch>::~outputadapter()
{}

template <typename T1Vector, typename T3Vector, typename Arch>
CUDA
void
outputadapter<T1Vector, T3Vector, Arch>::reduce_pf(size_t i, Real p, const Real3 & f)
{

#ifdef __CUDA_ARCH__
    __atomicAdd(&(vPhi_ptr[i]),p);
    vF_ptr[i] += f;
#else
    vPhi[i] += p;
    vF[i] += f;
#endif

}

template <typename T1Vector, typename T3Vector, typename Arch>
DEVICE
void
outputadapter<T1Vector, T3Vector, Arch>::atomic_reduce_pf(size_t i, Real p, const Real3 & f)
{
    __atomicAdd(&(vPhi_ptr[i]), p);
    __atomicAdd(&(vF_ptr[i]).x, f.x);
    __atomicAdd(&(vF_ptr[i]).y, f.y);
    __atomicAdd(&(vF_ptr[i]).z, f.z);
}

template <typename T1Vector, typename T3Vector, typename Arch>
DEVICE
void
outputadapter<T1Vector, T3Vector, Arch>::atomic_reduce_pf(size_t i, Real p, const REAL3 &f)
{
    __atomicAdd(&(vPhi_ptr[i]), p);
    __atomicAdd(&(vF_ptr[i]).x, f.x);
    __atomicAdd(&(vF_ptr[i]).y, f.y);
    __atomicAdd(&(vF_ptr[i]).z, f.z);

    //vF_ptr[i] %= f;
}

template <typename T1Vector, typename T3Vector, typename Arch>
DEVICE
void
outputadapter<T1Vector, T3Vector, Arch>::atomic_reduce_f(size_t i, const Real3 & f)
{
    __atomicAdd(&(vF_ptr[i]).x, f.x);
    __atomicAdd(&(vF_ptr[i]).y, f.y);
    __atomicAdd(&(vF_ptr[i]).z, f.z);
}

template <typename T1Vector, typename T3Vector, typename Arch>
DEVICE
void
outputadapter<T1Vector, T3Vector, Arch>::atomic_reduce_f(size_t i, const REAL3 &f)
{
    __atomicAdd(&(vF_ptr[i]).x, f.x);
    __atomicAdd(&(vF_ptr[i]).y, f.y);
    __atomicAdd(&(vF_ptr[i]).z, f.z);
}

template <typename T1Vector, typename T3Vector, typename Arch>
CUDA
void
outputadapter<T1Vector, T3Vector, Arch>::not_atomic_reduce_pf(size_t i, Real p, const Real3 & f)
{
    vPhi_ptr[i] += p;
    vF_ptr[i] += f;
}

template <typename T1Vector, typename T3Vector, typename Arch>
DEVICE
void
outputadapter<T1Vector, T3Vector, Arch>::not_atomic_reduce_pf(size_t i, Real p, const REAL3 & f)
{
    vPhi_ptr[i] += p;
    vF_ptr[i].x += f.x;
    vF_ptr[i].y += f.y;
    vF_ptr[i].z += f.z;
}


template struct outputadapter<std::vector<REAL, cuda_allocator<REAL> >, std::vector<XYZ<REAL>, cuda_allocator<XYZ<REAL> > >, Device<REAL> >;
template struct outputadapter<std::vector<REAL, cuda_host_allocator<REAL> >, std::vector<XYZ<REAL>, cuda_host_allocator<XYZ<REAL> > >, Host<REAL> >;
template struct outputadapter<cudavector<REAL, cuda_device_allocator<REAL> >, cudavector<XYZ<REAL>, cuda_device_allocator<XYZ<REAL> > >, Device<REAL> >;

}//namespace end
