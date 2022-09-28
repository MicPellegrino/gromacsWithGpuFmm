#ifndef _BK_FMM_iterables_hpp_
#define _BK_FMM_iterables_hpp_

namespace gmx_gpu_fmm{

template <typename Operation>
class bind_op_first
{
    typedef Operation op;
    typedef typename op::first_argument_type target_type;
public:
    typedef typename op::second_argument_type value_type;

    bind_op_first(target_type & target) : target(target)
    { }

    void operator () (const value_type & value)
    {
        op()(target, value);
    }

private:
    target_type & target;
};

template <typename Iterator, typename Operation>
class bound_iterable
{
    typedef bound_iterable self;
    typedef Iterator iterator;
public:

    bound_iterable(iterator iter) : iter(iter)
    { }

    self & operator ++ ()
    {
        ++iter;
        return *this;
    }

    bind_op_first<Operation> operator * ()
    {
        //auto bound = std::bind(op, std::ref(*thing), std::placeholders::_1);

        return bind_op_first<Operation>(*iter);
    }

    void prefetch(ssize_t distance) const
    {
        const char * p = reinterpret_cast<const char *>(&(*iter));
        // (p, r/w, low temporal hint)
        __builtin_prefetch(p + distance, 1, 3);
    }

private:
    iterator iter;
};

template <typename Storage, typename Value, template <typename, typename> class Operation>
class row_iterable
{
    typedef Storage storage_type;
    typedef typename storage_type::value_type storage_value_type;
    typedef typename storage_type::row_iterator storage_row_iterator;
    typedef Value compute_value_type;
    typedef Operation<storage_value_type, compute_value_type> op;
public:
    typedef typename storage_type::index_type index_type;
    typedef compute_value_type value_type;
    typedef bound_iterable<storage_row_iterator, op> row_iterator;

    row_iterable(storage_type & storage) : storage(storage)
    { }

    row_iterator get_row_iterator(index_type l, index_type m)
    {
        return row_iterator(storage.get_row_iterator(l, m));
    }

private:
    storage_type & storage;
};

}//namespace end


#endif
