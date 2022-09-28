#ifndef FMM_PRINTER_HPP
#define FMM_PRINTER_HPP

#include "architecture.hpp"
#include "fmm.hpp"

namespace gmx_gpu_fmm{

template <typename fmm_algorithm>
class fmm_printer{

public:

    fmm_algorithm *fm;
    double all_steps;

    fmm_printer(fmm_algorithm * fm):fm(fm), all_steps(1.0){}

    void print(std::string kernel_name){
        //printf("all steps = %d, fm->exec_time.steps = %f\n", all_steps, fm->exec_time.steps);
#ifdef CUDADEBUG
        if(fm->exec_time.actual_step == all_steps - 1|| fm->exec_time.actual_step == -1)
        {
            print_results(fm->io->n, fm->depth, fm->p, fm->exec_time.name(kernel_name), fm->exec_time.time(kernel_name));
            print_results(fm->exec_time.time(kernel_name+"_"));
        }
#endif
    }

    void print_walltime(int id)
    {
#ifdef CUDADEBUG
        if(fm->exec_time.actual_step == all_steps - 1)
            printf("kernels on device %d = %f\n", id, fm->exec_time.walltime());
#endif
    }

    void print_header(){
#ifdef CUDADEBUG
        //printf("actual_step %d\n",fm->exec_time.actual_step);
        if(fm->exec_time.actual_step < 1.0){
            printf("n        dept multipoleorder   stage                 time\n");
        }
#endif
    }

    void print_com_header(){
        printf("n        dept multipoleorder   type                  ref_time      kernel_time    real_error      imag_error\n");
    }

    void setup_steps(int steps)
    {
         all_steps = (double)steps;
    }

    void print_results(size_t n, size_t depth, int p, std::string kernel_name, double gpu_time){

        std::string output_0("%d");
        std::string output_1("%d");
        std::string output_2("%d");
        std::string output_3("%s");
        std::string output_4("%e");

        int size_of_name = kernel_name.size();
        std::string inter0("");
        std::string inter1("    ");
        std::string inter2("");
        std::string inter3("");
        std::string inter4("");


        size_t n_length = 1;
        size_t nn = 0;
        while(nn<=n){
            nn = std::pow(10,n_length);
            n_length++;
        }

        for(size_t i=0;i<10-n_length;i++){
             inter0+=" ";
        }

        size_t p_length = p/10;
        for(size_t i=0;i<16-p_length;i++){
             inter2+=" ";
        }

        for(int i=0;i<22-size_of_name;i++)
            inter3+=" ";

        std::string output;
        output = output_0+inter0+output_1+inter1+output_2+inter2+output_3+inter3+output_4+inter4;

        const char * char_output = output.c_str();

        printf(char_output,n,  depth, p, kernel_name.c_str(),gpu_time);
    }

    void print_results(size_t n, size_t depth, int p, std::string kernel_name, double cpu_time, double gpu_time, std::pair<double,double> errors){

        std::string output_0("%d");
        std::string output_1("%d");
        std::string output_2("%d");
        std::string output_3("%s");
        std::string output_4("%e  %e   %e    %e");

        int size_of_name = kernel_name.size();
        std::string inter0("");
        std::string inter1("    ");
        std::string inter2("");
        std::string inter3("");
        std::string inter4("");

        size_t n_length = 1;
        size_t nn = 0;
        while(nn<=n){
            nn = std::pow(10,n_length);
            n_length++;
        }
        for(size_t i=0;i<10-n_length;i++){
             inter0+=" ";
        }

        size_t p_length = p/10;
        for(size_t i=0;i<16-p_length;i++){
             inter2+=" ";
        }

        for(int i=0;i<20-size_of_name;i++)
            inter3+=" ";

        std::string output;
        output = output_0+inter0+output_1+inter1+output_2+inter2+output_3+inter3+output_4+inter4;

        const char * char_output = output.c_str();

        printf(char_output,n,  depth, p, kernel_name.c_str(), cpu_time, gpu_time, errors.first, errors.second);

    }

    void print_results(double gpu_time){

        printf(" %e\n",gpu_time);
    }

};

}//namespace end

#endif // FMM_PRINTER_HPP
