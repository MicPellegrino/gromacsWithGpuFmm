#ifndef _BK_FMM_generic_accumulate_hpp_
#define _BK_FMM_generic_accumulate_hpp_

#include <functional>

namespace gmx_gpu_fmm{

template <typename Storage>
class generic_accumulate
{
    typedef Storage storage_type;
public:
    typedef typename storage_type::value_type value_type;
    typedef typename storage_type::index_type index_type;

    generic_accumulate(storage_type & storage) : storage(storage)
    { }

    void operator () (index_type l, index_type m, const value_type & value)
    {
        storage(l, m) += value;
    }

private:
    storage_type & storage;
};

template <typename StorageT, typename ValueT = StorageT>
class just_accumulate : public std::binary_function<StorageT, ValueT, void>
{
    typedef StorageT storage_type;
public:
    typedef ValueT value_type;

    just_accumulate() { }

    void operator () (storage_type & target, const value_type & value) const
    {
        target += value;
    }
};

}//namespace end

#endif
