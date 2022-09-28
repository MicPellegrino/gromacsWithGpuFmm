#ifndef _BK_FMM_energy_stress_triangular_array_group_hpp_
#define _BK_FMM_energy_stress_triangular_array_group_hpp_

#include <functional>

namespace gmx_gpu_fmm{

template <typename Tp>
class EnergyStressTriangularArrayGroup
{
public:
    typedef Tp value_type;

private:
    struct _ABC
    {
        _ABC(int p)
            : a(p), b(p), c(p)
        { }

        value_type a;
        value_type b;
        value_type c;
    };

    typedef EnergyStressTriangularArrayGroup _Self;
    //EnergyStressTriangularArrayGroup(const _Self &) = delete;
    //_Self & operator = (const _Self &) = delete;

public:
    EnergyStressTriangularArrayGroup(int p)
        : energy(p), stress(p)
    { }

    void zero()
    {
        energy.zero();
        stress.a.zero();
        stress.b.zero();
        stress.c.zero();
    }

    void zero(size_t p)
    {
        energy.zero(p);
        stress.a.zero(p);
        stress.b.zero(p);
        stress.c.zero(p);
    }

    void populate_lower()
    {
        energy.populate_lower();
        stress.a.populate_lower();
        stress.b.populate_lower();
        stress.c.populate_lower();
    }

    value_type energy;
    _ABC stress;
};


template <typename ESGroupStorage, typename Value, template <typename, typename> class Operation>
class EnergyStressStorageOperator : public std::binary_function<ESGroupStorage, Value, void>
{
    typedef ESGroupStorage es_group_storage_type;
    typedef typename es_group_storage_type::value_type storage_type;
    typedef typename storage_type::value_type storage_value_type;
    typedef Value compute_value_type;
    typedef Operation<storage_value_type, compute_value_type> Op;
    // FIXME: use traits to get complex::value_type || type
    typedef typename compute_value_type::value_type scale_type;
public:
    typedef typename storage_type::index_type index_type;

    EnergyStressStorageOperator(scale_type a, scale_type b, scale_type c)

    { scale.a = a;
      scale.b = b;
      scale.c = c;
    }

    void operator () (es_group_storage_type & st, index_type l, index_type m, const compute_value_type & value) const
    {
        const Op op;
        op(st.energy(l, m), value);
        op(st.stress.a(l, m), value * scale.a);
        op(st.stress.b(l, m), value * scale.b);
        op(st.stress.c(l, m), value * scale.c);
    }

private:
    struct { scale_type a, b, c; } scale;
};

}//namespace end

#endif
