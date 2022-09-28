#ifndef _BK_testdata_hpp_
#define _BK_testdata_hpp_

#include "xyz.hpp"
#include "abc.hpp"
#include <chrono>
#include <vector>
#include <map>
#include <memory>

namespace gmx_gpu_fmm{

template <typename T>
class Testdata {

public:

    const XYZ<T> e1;
    const XYZ<T> e2;
    const XYZ<T> e3;

    T boxlength_;
    XYZ<T> boxsize;
    int n_;
    std::vector<std::vector<int>> excl;
    std::vector<T> x_;
    std::vector<T> y_;
    std::vector<T> z_;
    std::vector<T> q_;

    ABC<XYZ<T> > abc;
    XYZ<T> box_centre;
    std::map<std::string, double> reference_energies;

    Testdata();
    ~Testdata();

    int init_values(int dataset_id, bool open_boundaries);

    template <typename Datatype>
    void set(Datatype* inputdata);

    int n();
    T* x();
    T* y();
    T* z();
    T* q();
};

}//namespace end

#endif
