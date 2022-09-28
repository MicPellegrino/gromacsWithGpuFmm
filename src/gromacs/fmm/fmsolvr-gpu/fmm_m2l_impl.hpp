#include "fmm.hpp"
#include "cuda_M2L.hpp"
#include "cuda_lib.hpp"
#include <thread>
#include <mutex>

namespace gmx_gpu_fmm{

std::mutex atomic_write;

template <typename CoefficientMatrix, typename OperatorMatrix>
void M2L_CPU(const CoefficientMatrix* omega, const OperatorMatrix* B, CoefficientMatrix* mu, size_t p)
{
    typedef typename CoefficientMatrix::complex_type complex_type;

    for (ssize_t l = 0; l <= p; ++l)
    {
        for (ssize_t m = 0; m <= l; ++m)
        {
            complex_type mu_l_m(0.);

            for (ssize_t j = 0; j <= p; j+=2 )
            {
                for (ssize_t k = -j; k < 0; ++k)
                {
                    mu_l_m += B->get(j + l, k + m) * omega->get(j, k);
                }

                for (ssize_t k = 0; k <= j; ++k)
                {
                    mu_l_m += B->get(j + l, k + m) * omega->get_upper(j, k);
                }
            }

            for (ssize_t j = 1; j <= p; j+=2)
            {
                for (ssize_t k = -j; k < 0; ++k)
                {
                    mu_l_m -= B->get(j + l, k + m) * omega->get(j, k);
                }

                for (ssize_t k = 0; k <= j; ++k)
                {
                    mu_l_m -= B->get(j + l, k + m) * omega->get_upper(j, k);
                }
            }

            atomic_write.lock();
            mu->operator ()(l, m) += mu_l_m;
            atomic_write.unlock();
        }
    }
}

template <typename Box>
void make_m2l_on_chunk(size_t boxid_start, size_t boxid_end, size_t num_of_efective_ops, Box* box, size_t p)
{

    for (size_t boxid = boxid_start; boxid < boxid_end; ++boxid)
    {
        for (size_t i = 0; i < num_of_efective_ops; ++i)
        {
            M2L_CPU(box[boxid].omega, box[boxid].b_operators[i], box[boxid].b_targets[i], p);
        }
    }
}

}//namespace end
