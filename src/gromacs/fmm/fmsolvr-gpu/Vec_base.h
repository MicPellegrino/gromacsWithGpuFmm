#ifndef _BK_FMM_Vec_base_h_
#define _BK_FMM_Vec_base_h_

namespace gmx_gpu_fmm
{

template <typename Vec, typename Scalar>
struct Vec_base
{
	typedef Vec vec_type;
	typedef Scalar scalar_type;

	enum { VecSize = sizeof(vec_type) / sizeof(scalar_type) };

	vec_type unaligned_load(const scalar_type * p) const
	{
		return *p;
	}

	vec_type from(const scalar_type * p) const
	{
		return *(reinterpret_cast<const vec_type *>(p));
	}

	const vec_type & at(const scalar_type * p) const
	{
		return *(reinterpret_cast<const vec_type *>(p));
	}

	vec_type & at(scalar_type * p)
	{
		return *(reinterpret_cast<vec_type *>(p));
	}

	const vec_type & element(const vec_type & v, size_t i) const
	{
		return reinterpret_cast<const scalar_type *>(&v)[i];
	}

	scalar_type & element(vec_type & v, size_t i)
	{
		return reinterpret_cast<scalar_type *>(&v)[i];
	}

	size_t index_floor(size_t x) const
	{
		size_t skewed = x % VecSize;
		return x - skewed;
	}

	size_t index_ceil(size_t x) const
	{
		size_t skewed = x % VecSize;
		return x + ((skewed == 0) ? (0) : (VecSize - skewed));
	}

};

template <typename Vec>
struct Vec_traits : public Vec_base<Vec, Vec>
{
	typedef Vec_base<Vec, Vec> base;
	typedef typename base::vec_type vec_type;
	typedef typename base::scalar_type scalar_type;

    CUDA
	vec_type zero() const
	{
		return 0;
	}
    CUDA
	vec_type same(scalar_type s) const
	{
		return s;
	}
    CUDA
	scalar_type sum(const vec_type & v) const
	{
		return v;
	}
};

template <typename Vec>
Vec rcplength(const Vec & dx, const Vec & dy, const Vec & dz)
{
	Vec ssq = dx * dx + dy * dy + dz * dz;
	return rsqrt(ssq);
}

}//namespace end

#endif
