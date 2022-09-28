#ifndef _BK_FMM_triangular_array_adaptor_hpp_
#define _BK_FMM_triangular_array_adaptor_hpp_

namespace gmx_gpu_fmm{

// or generic?
template <typename Storage, typename Value, template <typename, typename> class Operation>
class simple_triangular_array_adaptor
{
    typedef Storage storage_type;
    typedef typename storage_type::value_type storage_value_type;
    typedef Value compute_value_type;
    typedef Operation<storage_value_type, compute_value_type> Op;
public:
    typedef typename storage_type::index_type index_type;
    typedef compute_value_type value_type;

    simple_triangular_array_adaptor(storage_type & storage)
        : storage(storage)
    {
        //storage = storage_;
    }

    void operator () (index_type l, index_type m, const compute_value_type & value) const
    {
        const Op op;
        op(storage(l, m), value);
    }

private:
    storage_type & storage;
};


template <typename Operation>
class custom_triangular_array_adaptor
{
    typedef typename Operation::first_argument_type storage_type;
    typedef typename Operation::second_argument_type compute_value_type;
    typedef Operation Op;
public:
    typedef typename Operation::index_type index_type;
    typedef compute_value_type value_type;


    custom_triangular_array_adaptor(storage_type & storage_, const Op & op_ = Op())
        : storage(storage_), op(op_)
    {
        //init(storage_,op_);
    }

    template <typename Args>
    custom_triangular_array_adaptor(storage_type & storage, Args args)
        : storage(storage), op(Op(args))
    {
        //init(storage, Op(args...));
        //custom_triangular_array_adaptor(storage, );
    }

    void operator () (index_type l, index_type m, const compute_value_type & value)
    {
        op(storage, l, m, value);
    }

    Op & get_op()
    {
        return op;
    }

private:
    //template <typename storage_type,typename compute_value_type>
    void init(storage_type & st, const compute_value_type & op_ = Op())
    {

        storage = st;
        op = op_;
    }

    storage_type & storage;
    Op op;
};

}//namespace end

#endif
