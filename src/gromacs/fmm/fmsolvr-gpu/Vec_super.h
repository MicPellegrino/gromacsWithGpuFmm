#ifndef _BK_FMM_Vec_Super_h_
#define _BK_FMM_Vec_Super_h_

#include "Vec_base.h"

namespace gmx_gpu_fmm
{

template <typename SIMD, unsigned k>
struct SuperSIMD
{
	typedef Vec_traits<SIMD> simd_traits;
	typedef typename simd_traits::scalar_type Scalar;
	enum { SIMD_width = simd_traits::VecSize };

	SIMD a[k];

	Scalar & operator[] (size_t i)
	{
		size_t simd_elem = i / SIMD_width;
		size_t scalar_elem = i % SIMD_width;
		return a[simd_elem][scalar_elem];
	}

	const Scalar & operator[] (size_t i) const
	{
		size_t simd_elem = i / SIMD_width;
		size_t scalar_elem = i % SIMD_width;
		return a[simd_elem][scalar_elem];
	}

	SuperSIMD operator - () const
	{
		SuperSIMD c;
		for (size_t i = 0; i < k; ++i)
			c.a[i] = -a[i];
		return c;
	}

	SuperSIMD operator + (const SuperSIMD & b) const
	{
		SuperSIMD c;
		for (size_t i = 0; i < k; ++i)
			c.a[i] = a[i] + b.a[i];
		return c;
	}

	SuperSIMD operator - (const SuperSIMD & b) const
	{
		SuperSIMD c;
		for (size_t i = 0; i < k; ++i)
			c.a[i] = a[i] - b.a[i];
		return c;
	}

	SuperSIMD operator * (const SuperSIMD & b) const
	{
		SuperSIMD c;
		for (size_t i = 0; i < k; ++i)
			c.a[i] = a[i] * b.a[i];
		return c;
	}

	SuperSIMD operator / (const SuperSIMD & b) const
	{
		SuperSIMD c;
		for (size_t i = 0; i < k; ++i)
			c.a[i] = a[i] / b.a[i];
		return c;
	}

	SuperSIMD & operator += (const SuperSIMD & b)
	{
		for (size_t i = 0; i < k; ++i)
			a[i] += b.a[i];
		return *this;
	}

	SuperSIMD & operator -= (const SuperSIMD & b)
	{
		for (size_t i = 0; i < k; ++i)
			a[i] -= b.a[i];
		return *this;
	}

	SuperSIMD & operator *= (const SuperSIMD & b)
	{
		for (size_t i = 0; i < k; ++i)
			a[i] *= b.a[i];
		return *this;
	}
};

template <typename SIMD, unsigned k>
struct Vec_traits<SuperSIMD<SIMD, k> > : public Vec_base<SuperSIMD<SIMD, k>, typename Vec_traits<SIMD>::scalar_type>
{
	typedef Vec_base<SuperSIMD<SIMD, k>, typename Vec_traits<SIMD>::scalar_type> base;
	typedef typename base::vec_type vec_type;
	typedef typename base::scalar_type scalar_type;
	typedef SIMD nested_simd_type;
	typedef Vec_traits<nested_simd_type> nested_simd_traits;

	vec_type unaligned_load(const scalar_type * p) const
	{
		nested_simd_traits nested;
		vec_type z;
		for (size_t i = 0; i < k; ++i, p += nested.VecSize)
			z.a[i] = nested.unaligned_load(p);
		return z;
	}

	vec_type zero() const
	{
		nested_simd_traits nested;
		vec_type z;
		for (size_t i = 0; i < k; ++i)
			z.a[i] = nested.zero();
		return z;
	}

	vec_type same(scalar_type s) const
	{
		nested_simd_traits nested;
		vec_type v;
		for (size_t i = 0; i < k; ++i)
			v.a[i] = nested.same(s);
		return v;
	}

	scalar_type sum(const vec_type & v) const
	{
		nested_simd_traits nested;
		scalar_type s = 0;
		for (size_t i = 0; i < k; ++i)
			s += nested.sum(v.a[i]);
		return s;
	}

	nested_simd_type supersum(const vec_type & v) const
	{
		nested_simd_traits nested;
		nested_simd_type s = v.a[0];
		for (size_t i = 1; i < k; ++i)
			s += v.a[i];
		return s;
	}
};

template <typename SIMD, unsigned k>
SuperSIMD<SIMD, k> rsqrt(const SuperSIMD<SIMD, k> & a)
{
	SuperSIMD<SIMD, k> r;
	for (size_t i = 0; i < k; ++i)
        {
	        r.a[i] = rsqrt(a.a[i]);
	}
	return r;
}

}//namespace end

#endif
