#ifndef _BK_FMM_Vec_Array_h_
#define _BK_FMM_Vec_Array_h_

#include "Vec_base.h"

namespace gmx_gpu_fmm{

template <typename Scalar, unsigned k>
struct TrivialArray
{
	Scalar a[k];

	Scalar & operator[] (size_t i)
	{
		return a[i];
	}

	const Scalar & operator[] (size_t i) const
	{
		return a[i];
	}

	TrivialArray operator + (const TrivialArray & b) const
	{
		TrivialArray c;
		for (size_t i = 0; i < k; ++i)
			c[i] = a[i] + b[i];
		return c;
	}

	TrivialArray operator - (const TrivialArray & b) const
	{
		TrivialArray c;
		for (size_t i = 0; i < k; ++i)
			c[i] = a[i] - b[i];
		return c;
	}

	TrivialArray operator * (const TrivialArray & b) const
	{
		TrivialArray c;
		for (size_t i = 0; i < k; ++i)
			c[i] = a[i] * b[i];
		return c;
	}

	TrivialArray & operator += (const TrivialArray & b)
	{
		for (size_t i = 0; i < k; ++i)
			a[i] += b[i];
		return *this;
	}

	TrivialArray & operator -= (const TrivialArray & b)
	{
		for (size_t i = 0; i < k; ++i)
			a[i] -= b[i];
		return *this;
	}
};

template <typename Scalar, unsigned k>
struct Vec_traits<TrivialArray<Scalar, k> > : public Vec_base<TrivialArray<Scalar, k>, Scalar>
{
	typedef Vec_base<TrivialArray<Scalar, k>, Scalar> base;
	typedef typename base::vec_type vec_type;
	typedef typename base::scalar_type scalar_type;

	vec_type zero() const
	{
		vec_type z;
		for (size_t i = 0; i < k; ++i)
			z[i] = 0;
		return z;
	}

	vec_type same(scalar_type s) const
	{
		vec_type v;
		for (size_t i = 0; i < k; ++i)
			v[i] = s;
		return v;
	}

	scalar_type sum(const vec_type & v) const
	{
		scalar_type s = 0;
		for (size_t i = 0; i < k; ++i)
			s += v[i];
		return s;
	}
};

template <typename Scalar, unsigned k>
TrivialArray<Scalar, k> rsqrt(const TrivialArray<Scalar, k> & a)
{
	TrivialArray<Scalar, k> r;
	for (size_t i = 0; i < k; ++i) {
		r[i] = rsqrt(a[i]);
	}
	return r;
}

}//namespace end

#endif
