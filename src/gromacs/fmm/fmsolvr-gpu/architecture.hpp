#ifndef ARCHITECTURE_HPP
#define ARCHITECTURE_HPP

#include "data_type.hpp"
#include "cuda_alloc.hpp"
#include "timer.hpp"
#include "alloc_counter.hpp"
#include "cuda_keywords.hpp"

namespace gmx_gpu_fmm{

/*! \brief Abstracts the memory allocation for CUDA Unified Memory used on Host and Device.
 *  \tparam Type Allocation datatype
 */
template<typename Type>
class Device{

public:

    //! Allocation data type
    typedef Type value_type;
    //! Unified memory allocator allocator
    typedef cuda_allocator<value_type> allocator;

    //! Device synchronization only in case of deivce debug modus
    static void devSyncDebug()
    {
#ifdef CUDADEBUG
        cudaDeviceSynchronize();
#endif
    }

    //! Device synchronization
    static void devSync()
    {
        cudaDeviceSynchronize();
    }

    /*! \brief Allocates memory on Device
     *  \tparam T     Data type
     *  \param ptr    Pointer to the memory location
     *  \param size   Allocation size in bytes
     */
    template<typename T>
    static void custom_alloc(T &ptr, size_t size){

        if(size > 0)
        {
            cudaMallocManaged(&ptr, size);
            alloc_counter::up();
        }
    }

    /*! \brief Frees memory on Device
     *  \tparam T  Data type
     *  \param ptr Pointer to the memory location
     */
    template<typename T>
    static void custom_free(T ptr){

        cudaFree(ptr);
        alloc_counter::down();
    }

    /*! \brief Rebind structure to enable different data type allocations
     *  \tparam other_type Datatype to rebind
     */
    template<typename other_type>
    struct rebind{
        typedef cuda_allocator<other_type> other;
    };  
};
/*! \brief Abstracts the memory allocation for on Host
 *  \tparam Type Allocation datatype
 */
template<typename Type>
struct Host{

public:

    //! Allocation data type
    typedef Type value_type;
    //! Host allocator
    typedef std::allocator<value_type> allocator;

    //! Left empty since no need to synchronize on the Host
    static void devSync(){}

    //! Left empty since no need to synchronize on the Host
    static void devSyncDebug(){}

    /*! \brief Allocates memory on Host
     *  \tparam T     Data type
     *  \param ptr    Pointer to the memory location
     *  \param size   Allocation size in bytes
     */
    template<typename T>
    static void custom_alloc(T &ptr, size_t size){

        ptr = (T) malloc(size);
        //alloc_counter::up();
        //printf("malloc\n");
    }

    /*! \brief Frees memory on Host
     *  \tparam T  Data type
     *  \param ptr Pointer to the memory location
     */
    template<typename T>
    static void custom_free(T &ptr){

        free(ptr);
        //alloc_counter::down();
        //printf("free\n");
    }

    /*! \brief Rebind structure to enable different data type allocations
     *  \tparam other_type Datatype to rebind
     */
    template<typename other_type>
    struct rebind{
        typedef std::allocator<other_type> other;
    };
};

/*! \brief Abstracts the memory allocation on Host with pinned memory
 *  \tparam Type Allocation datatype
 */
template<typename Type>
struct CudaHost{

public:

    //! Allocation data type
    typedef Type value_type;

    //! Host allocator
    typedef std::allocator<value_type> allocator;

    //! Left empty since no need to synchronize on the Host
    static void devSync(){}

    //! Left empty since no need to synchronize on the Host
    static void devSyncDebug(){}

    /*! \brief Allocates pinned memory on Host
     *  \tparam T     Data type
     *  \param ptr    Pointer to the memory location
     *  \param size   Allocation size in bytes
     */
    template<typename T>
    static void custom_alloc(T &ptr, size_t size){

        cudaMallocHost(&ptr, size);
        //alloc_counter::up();
        //printf("malloc\n");
    }

    /*! \brief Frees pinned memory on Host
     *  \tparam T  Data type
     *  \param ptr Pointer to the memory location
     */
    template<typename T>
    static void custom_free(T &ptr){

        free(ptr);
        //alloc_counter::down();
        //printf("free\n");
    }

    /*! \brief Rebind structure to enable different data type allocations
     *  \tparam other_type Datatype to rebind
     */
    template<typename other_type>
    struct rebind{
        typedef std::allocator<other_type> other;
    };
};
/*! \brief Abstracts the memory allocation on Host with device memory only
 *  \tparam Type Allocation datatype
 */
template<typename Type>
class DeviceOnly{

public:

    //! Allocation data type
    typedef Type value_type;
    //! device allocator
    typedef cuda_device_allocator<value_type> allocator;

    //! Left empty since no need to synchronize for device only memory
    static void devSyncDebug()
    {
#ifdef CUDADEBUG
        cudaDeviceSynchronize();
#endif
    }

    //! Left empty since no need to synchronize for device only memory
    static void devSync()
    {
        cudaDeviceSynchronize();
    }

    /*! \brief Allocates memory on Device
     *  \tparam T     Data type
     *  \param ptr    Pointer to the memory location
     *  \param size   Allocation size in bytes
     */
    template<typename T>
    static void custom_alloc(T &ptr, size_t size){

        cudaMalloc(&ptr, size );
        //alloc_counter::up();
    }

    /*! \brief Frees pinned memory on Device
     *  \tparam T  Data type
     *  \param ptr Pointer to the memory location
     */
    template<typename T>
    static void custom_free(T ptr){

        cudaFree(ptr);
        //alloc_counter::down();
    }

    /*! \brief Rebind structure to enable different data type allocations
     *  \tparam other_type Datatype to rebind
     */
    template<typename other_type>
    struct rebind{
        typedef cuda_allocator<other_type> other;
    };
};


}//namespace end


#endif // ARCHITECTURE_HPP
