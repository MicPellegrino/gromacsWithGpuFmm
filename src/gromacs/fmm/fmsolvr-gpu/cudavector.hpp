#ifndef CUDAVECTOR_HPP
#define CUDAVECTOR_HPP

namespace gmx_gpu_fmm{

template<typename T, typename Allocator>
class cudavector
{

public:
    T* p;
    typedef T value_type;

    cudavector(){}

    cudavector(size_t size, T val)
    {
        p = allocator.allocate(size);
        cudaMemset(p, 0, size*sizeof(val));
    }

    ~cudavector()
    {
        allocator.deallocate(p,size());
    }

    void resize(size_t size)
    {
        p = allocator.allocate(size);
    }

    T& operator[](size_t position)
    {
        return p[position];
    }

    size_t size()
    {
        return mysize;
    }

private:
    size_t mysize;
    Allocator allocator;
};

}//namespace end
#endif // CUDAVECTOR_HPP

