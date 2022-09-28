#ifndef _BK_FMM_cuda_alloc_hpp_
#define _BK_FMM_cuda_alloc_hpp_

#include <new>
#include <stdexcept>
#include <cstdlib>
#include <cerrno>
#include <iostream>
#include <cuda.h>
#include <type_traits>
#include "alloc_counter.hpp"

namespace gmx_gpu_fmm{

//c++03 alignof implementation
template <typename T> struct alignof_;

template <int size_diff>
struct helper
{
    template <typename T>
    struct Val {
        enum { value = size_diff };
    };
};

template <>
struct helper<0>
{
    template <typename T>
    struct Val {
        enum { value = alignof_<T>::value };
    };
};

template <typename T>
struct alignof_
{
    struct Big { T x; char c; };

    enum {
        diff = sizeof (Big) - sizeof (T),
        value = helper<diff>::template Val<Big>::value
        };
};

template <class T, unsigned alignment = alignof_<T>::value>
class cuda_allocator;

//specialize for void:
template <unsigned alignment>
class cuda_allocator<void, alignment>
{
    public:
    typedef void* pointer;
    typedef const void* const_pointer;
    // reference-to-void members are impossible.
    typedef void value_type;
    template <class U> struct rebind { typedef cuda_allocator<U, alignment> other; };
};

template <class T, unsigned alignment>
class cuda_allocator : public std::allocator<T>
{
public:
    typedef size_t size_type;
    typedef std::ptrdiff_t difference_type;

    typedef T* pointer;
    typedef const T* const_pointer;

    typedef T& reference;
    typedef const T& const_reference;

    typedef T value_type;

    template <class U>
    struct rebind { typedef cuda_allocator<U, alignment> other; };

    cuda_allocator() throw() { }
    cuda_allocator(const cuda_allocator& a) throw()
        : std::allocator<T>(a) { }
    template <class U>
    cuda_allocator(const cuda_allocator<U, alignment>&) throw() { }
    ~cuda_allocator() throw() { }

    pointer allocate(size_type n, cuda_allocator<void, alignof_<void*>::value>::const_pointer /*hint*/ = 0) // space for n Ts
    {
        void * p;
        unsigned ali = alignment;
        if (ali < alignof_<void*>::value)
            ali = alignof_<void*>::value;
        
        cudaMallocManaged(&p,n*sizeof(T));
        //alloc_counter::up();
        //printf("alloc %d\n",alloc_counter::get());


        return pointer(p);
    }

    void deallocate(pointer p, size_type /*n*/)   // deallocate n Ts, don't destroy
    {
        cudaFree(p);
        //alloc_counter::down();
        //printf("free %d\n",alloc_counter::get());
    }

    void construct(pointer p, const T& val) { new(p) T(val); }          // initialize *p by val
    void destroy(pointer p) { p->~T(); }     // destroy *p but don't deallocate

};

template<class T1, class T2, unsigned alignment>
bool operator==(const cuda_allocator<T1, alignment>&, const cuda_allocator<T2, alignment>&) //noexcept
{
    return true;
}

template<class T1, class T2, unsigned alignment>
bool operator!=(const cuda_allocator<T1, alignment>&, const cuda_allocator<T2, alignment>&) //noexcept
{
    return false;
}


template <class T, unsigned alignment = alignof_<T>::value>
class cuda_host_allocator;

//specialize for void:
template <unsigned alignment>
class cuda_host_allocator<void, alignment>
{
    public:
    typedef void* pointer;
    typedef const void* const_pointer;
    // reference-to-void members are impossible.
    typedef void value_type;
    template <class U> struct rebind { typedef cuda_host_allocator<U, alignment> other; };
};

template <class T, unsigned alignment>
class cuda_host_allocator : public std::allocator<T>
{
public:
    typedef size_t size_type;
    typedef std::ptrdiff_t difference_type;

    typedef T* pointer;
    typedef const T* const_pointer;

    typedef T& reference;
    typedef const T& const_reference;

    typedef T value_type;

    template <class U>
    struct rebind { typedef cuda_host_allocator<U, alignment> other; };

    cuda_host_allocator() throw() { }
    cuda_host_allocator(const cuda_host_allocator& a) throw()
        : std::allocator<T>(a) { }
    template <class U>
    cuda_host_allocator(const cuda_host_allocator<U, alignment>&) throw() { }
    ~cuda_host_allocator() throw() { }

    pointer allocate(size_type n, cuda_host_allocator<void, alignof_<void*>::value>::const_pointer /*hint*/ = 0 ) // space for n Ts
    {
        void * p;
        unsigned ali = alignment;
        if (ali < alignof_<void*>::value)
            ali = alignof_<void*>::value;

        cudaMallocHost(&p,n*sizeof(T));
        //alloc_counter::up();

        return pointer(p);
    }
    void deallocate(pointer p, size_type /*n*/)   // deallocate n Ts, don't destroy
    {
        cudaFreeHost(p);
        //alloc_counter::down();
    }

    void construct(pointer p, const T& val) { new(p) T(val); }          // initialize *p by val
    void destroy(pointer p) { p->~T(); }     // destroy *p but don't deallocate
};

template<class T1, class T2, unsigned alignment>
bool operator==(const cuda_host_allocator<T1, alignment>&, const cuda_host_allocator<T2, alignment>&) //noexcept
{
    return true;
}

template<class T1, class T2, unsigned alignment>
bool operator!=(const cuda_host_allocator<T1, alignment>&, const cuda_host_allocator<T2, alignment>&) //noexcept
{
    return false;
}

template <class T, unsigned alignment = alignof_<T>::value>
class cuda_device_allocator;

//specialize for void:
template <unsigned alignment>
class cuda_device_allocator<void, alignment>
{
    public:
    typedef void* pointer;
    typedef const void* const_pointer;
    // reference-to-void members are impossible.
    typedef void value_type;
    template <class U> struct rebind { typedef cuda_device_allocator<U, alignment> other; };
};

template <class T, unsigned alignment>
class cuda_device_allocator : public std::allocator<T>
{
public:
    typedef size_t size_type;
    typedef std::ptrdiff_t difference_type;

    typedef T* pointer;
    typedef const T* const_pointer;

    typedef T& reference;
    typedef const T& const_reference;

    typedef T value_type;

    template <class U>
    struct rebind { typedef cuda_device_allocator<U, alignment> other; };

    cuda_device_allocator() throw() { }
    cuda_device_allocator(const cuda_device_allocator& a) throw()
        : std::allocator<T>(a) { }
    template <class U>
    cuda_device_allocator(const cuda_device_allocator<U, alignment>&) throw() { }
    ~cuda_device_allocator() throw() { }

    pointer allocate(size_type n, cuda_device_allocator<void, alignof_<void*>::value>::const_pointer /*hint*/ = 0) // space for n Ts
    {
        void * p;
        unsigned ali = alignment;
        if (ali < alignof_<void*>::value)
            ali = alignof_<void*>::value;

        cudaMalloc(&p,n*sizeof(T));
        //alloc_counter::up();

        return pointer(p);
    }

    void deallocate(pointer p, size_type /*n*/)   // deallocate n Ts, don't destroy
    {
        cudaFree(p);
        //alloc_counter::down();
    }

    //void construct(pointer p, const T& val) { new(p) T(val); }          // initialize *p by val
    //void construct(pointer p, const T& val) { }          // initialize *p by val
    void destroy(pointer p) { p->~T(); }     // destroy *p but don't deallocate

};

template<class T1, class T2, unsigned alignment>
bool operator==(const cuda_device_allocator<T1, alignment>&, const cuda_device_allocator<T2, alignment>&) //noexcept
{
    return true;
}

template<class T1, class T2, unsigned alignment>
bool operator!=(const cuda_device_allocator<T1, alignment>&, const cuda_device_allocator<T2, alignment>&) //noexcept
{
    return false;
}

}//namespace end

#endif
