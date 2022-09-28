#ifndef _BK_FMM_lib_hpp_
#define _BK_FMM_lib_hpp_


namespace fmm_out {

    namespace lib {

#ifdef __CUDACC__

        template <typename Tp>
	CUDA
	const Tp & min(const Tp & a, const Tp & b)
	{
            return b < a ? b : a;
	}

        template <typename Tp>
	CUDA
	const Tp & max(const Tp & a, const Tp & b)
	{
            return b < a ? a : b;
	}

#else

	using std::min;
	using std::max;
    using std::pair;

#endif

    }

}  // namespace fmm

#endif
// vim: et:ts=4:sw=4
