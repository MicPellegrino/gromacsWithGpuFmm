#include "global_functions.hpp"
#include "testdata.hpp"
#include <random>
//#include "testdata/qxyz_nacl262144.hpp"
//#include "testdata/qxyz_nacl8.hpp"
//#include "testdata/qxyz_nacl2.hpp"
#include "testdata/qxyz_nacl32768.hpp"
#include "testdata/qxyz_silicamelt.hpp" //
//#include "testdata/qxyz_silicamelt_test.hpp"  //older one
#include "testdata/qxyz_2_particles.hpp"
#include "testdata/qxyz_25_random_data.hpp"
#include "testdata/qxyz_anion_canal_285119.hpp"
#include "testdata/qxyz_water.hpp"
#include "xyz.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <string>

namespace gmx_gpu_fmm{

int get_env_integer(int dflt, const char* vrnm)
{
    if (const char * val = getenv(vrnm))
        return atoi(val);
    else
        return dflt;
}

double get_env_real(double dflt, const char* vrnm)
{
    if (const char * val = getenv(vrnm))
        return atof(val);
    else
        return dflt;
}

template <typename T>
Testdata<T>::Testdata():e1(1,0,0),e2(0,1,0),e3(0,0,1)
{
}

template <typename T>
Testdata<T>::~Testdata()
{
}

template <typename T>
int Testdata<T>::init_values(int dataset_id, bool open_boundaries)
{

    if(dataset_id > 10)
    {
        printf("Invalid Dataset Id, Datasets from 0 to 10 \n");
        return 0;
    }
    if(dataset_id == -1)
    {
        printf("Reading pdb\n");
        const char* c_filename = getenv("PDBFILE");
        std::string filename(c_filename);
        std::ifstream infile(filename);
        if(!infile.is_open()){
           std::cerr<<"Failed to open file"<<std::endl;
           return -1;
         }
        std::string line;
        std::string atom_type;
        std::string x_str;
        std::string y_str;
        std::string z_str;
        std::string q_str;

        size_t size = 0;
        while (std::getline(infile, line))
        {
            if(line.substr(0,4) == "ATOM")
            {
                size++;
            }
        }
        n_ = size;
        x_.resize(n_);
        y_.resize(n_);
        z_.resize(n_);
        q_.resize(n_);

        infile.close();
        infile.open(filename);

        int i = 0;
        while (std::getline(infile, line))
        {
            if(line.substr(0,6) == "CRYST1")
            {
                boxlength_ = (T) stod(line.substr(6,15))/10;
                boxsize.x = boxlength_;
                boxsize.y = (T) stod(line.substr(15,24))/10;
                boxsize.z = (T) stod(line.substr(24,33))/10;

            }
            if(line.substr(0,4) == "ATOM")
            {
                atom_type = line.substr(13,4);
                x_str = line.substr(26,12);
                y_str = line.substr(38,8);
                z_str = line.substr(46,8);
                q_str = line.substr(54,6);
                x_[i] = (T) stod(x_str)/10;
                y_[i] = (T) stod(y_str)/10;
                z_[i] = (T) stod(z_str)/10;
                q_[i] = (T) stod(q_str);
                //std::cout<<atom_type<<std::endl;
                //printf("%s\n", atom_type);
                if(atom_type == "Na  ")
                {
                    q_[i] = 1.0;
                }
                if(atom_type == "Cl  ")
                {
                    q_[i] = -1.0;
                }
                if(atom_type == "OW  ")
                {
                    q_[i] = -.8340000;
                }
                if(atom_type == "HW1 " || atom_type == "HW2 ")
                {
                    q_[i] = .4170000;
                }

                //std::cout<<atom_type<<x_str<<y_str<<z_str<<std::endl;
                i++;
            }
        }
        infile.close();
        excl.resize(n_);
        for(int i = 0; i < n_; i++)
        {
            excl[i].resize(1);
            excl[i][0] = i;
        }

        if(open_boundaries && false)
        {
            //putting atoms in box if needed
            //printf("boxlength %e\n", boxlength_);
            T min_x = boxlength_;
            T min_y = boxlength_;
            T min_z = boxlength_;

            T max_x = 0.0;
            T max_y = 0.0;
            T max_z = 0.0;

            for (int i = 0; i < n_; i++)
            {
                min_x = std::min(min_x, x_[i]);
                min_y = std::min(min_y, y_[i]);
                min_z = std::min(min_z, z_[i]);

                max_x = std::max(max_x, x_[i]);
                max_y = std::max(max_y, y_[i]);
                max_z = std::max(max_z, z_[i]);
            }

            T min = std::min(min_x, min_y);
            min = std::min(min, min_z);

            T max = std::max(max_x, max_y);
            max = std::max(max, max_z);

            for (int i = 0; i < n_; i++)
            {
                if(min_x < 0.0)
                {
                    x_[i] -= min_x;
                }
                if(min_y < 0.0)
                {
                    y_[i] -= min_y;
                }
                if(min_z < 0.0)
                {
                    z_[i] -= min_z;
                }
            }
            if(max > boxlength_)
            {
                boxlength_ = max;
            }
            if(min <  0.0)
            {
                boxlength_ -= min;
            }

            printf("min %f %f %f\n", min_x, min_y, min_z);
            printf("max %f %f %f\n", max_x, max_y, max_z);
            printf("min max %f %f\n", min, max);
        }
        if(!open_boundaries && false)
        {
            for (int i = 0; i < n_; i++)
            {
                if(x_[i] < 0.0)
                {
                    x_[i] += boxsize.x;
                }
                if(y_[i] < 0.0)
                {
                    y_[i] += boxsize.y;
                }
                if(z_[i] < 0.0)
                {
                    z_[i] += boxsize.z;
                }

                if(x_[i] >= boxsize.x)
                {
                    x_[i] -= boxsize.x;
                }
                if(y_[i] >= boxsize.y)
                {
                    y_[i] -= boxsize.y;
                }
                if(z_[i] >= boxsize.z)
                {
                    z_[i] -= boxsize.z;
                }
            }
        }

    }

    if(dataset_id == 0)
    {
        std::default_random_engine generator(get_env_integer(1234, "SEED"));

        boxlength_ = get_env_real(1.0,"RAND_BOX_SIZE");
        int particles_per_box = get_env_integer(256, "RAND_N");
        double squeeze = get_env_real(1.0,"RAND_SQUEEZE");
        int depth = get_env_integer(0, "DEPTH");
        int num_boxes_in_one_dim = std::pow(2,depth);
        //printf("num boxes in one dim %d\n",num_boxes_in_one_dim);
        double box_size = boxlength_/num_boxes_in_one_dim;
        int num_of_all_boxes = std::pow(num_boxes_in_one_dim,3);
        n_ = num_of_all_boxes * particles_per_box;

        T epsilon = 0.0;
        if(squeeze > 1.0)
        {
            epsilon = (box_size / squeeze)/2.0;
        }

        std::uniform_real_distribution<T> distribution(0.0 + epsilon, box_size - epsilon);
        //printf("%f %f\n", 0.0 + epsilon, box_size - epsilon);

        x_.resize(n_);
        y_.resize(n_);
        z_.resize(n_);
        q_.resize(n_);

        int index = 0;
        int cor_z = -1;
        int cor_y = -1;

        int skip = get_env_integer(1, "RAND_SKIP");
        int p_index = 0;
        for(int z = 0; z < num_boxes_in_one_dim; z++)
        {
            cor_z++;
            for(int y = 0; y < num_boxes_in_one_dim; y++)
            {
                cor_y++;
                for(int x = 0; x < num_boxes_in_one_dim; x++)
                {
                    index = num_boxes_in_one_dim*num_boxes_in_one_dim*z + num_boxes_in_one_dim*y + x;
                    if((index + cor_y + cor_z)%skip == 0)
                    {
                        for (int i = 0; i < particles_per_box*skip; ++i)
                        {
                            if(p_index < n_)
                            {
                                x_[p_index] = x * box_size + distribution(generator);
                                y_[p_index] = y * box_size + distribution(generator);
                                z_[p_index] = z * box_size + distribution(generator);
                                //printf("%d  (%f %f %f)\n", p_index, x_[p_index], y_[p_index], z_[p_index]);
                            }

                            ++p_index;
                        }
                    }
                }
            }
        }
        //printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
        if(p_index < n_)
        {
            for (int i = p_index - 1; i < n_; ++i)
            {
                //printf("%d  (%f %f %f)\n", i, x_[p_index], y_[p_index], z_[p_index]);
                x_[i] = (num_boxes_in_one_dim-1) * box_size + distribution(generator);
                y_[i] = (num_boxes_in_one_dim-1) * box_size + distribution(generator);
                z_[i] = (num_boxes_in_one_dim-1) * box_size + distribution(generator);
            }
        }

        std::normal_distribution<T> charge_distribution(0.0, 1.0);
        for (int i = 0; i < n_; ++i)
        {
            q_[i] = 1.0;//charge_distribution(generator);
        }

        reference_energies["open_reference"]               = 1.0;
        reference_energies["open"]                         = 1.0;
        reference_energies["periodic_reference"]           = 1.0;
        reference_energies["periodic"]                     = 1.0;
        reference_energies["dipole_correction_reference"]  = 1.0;
        reference_energies["dipole_correction"]            = 1.0;

        excl.resize(n_);
        for(int i = 0; i < n_; i++)
        {
            excl[i].resize(1);
            excl[i][0] = i;
        }
    }

    if(dataset_id == 1)
    {
        std::default_random_engine generator(get_env_integer(1234, "SEED"));

        boxlength_ = get_env_real(1.0,"RAND_BOX_SIZE");
        n_ = get_env_integer(256, "RAND_N");

        x_.resize(n_);
        y_.resize(n_);
        z_.resize(n_);
        q_.resize(n_);
        std::uniform_real_distribution<T> particle_distribution(0.0, boxlength_);
        //std::normal_distribution<T> charge_distribution(0.0, 1.0);
        std::uniform_real_distribution<T> charge_distribution(-2.0, 2.0);
        T centre = boxlength_/2.0;
        for (int i = 0 ; i < n_; ++i)
        {
            T xi = boxlength_;
            T yi = boxlength_;
            T zi = boxlength_;
            while( std::pow( (xi-centre)*(xi-centre) + (yi-centre)*(yi-centre) + (zi-centre)*(zi-centre) , 0.5) > boxlength_/2.0)
            {
                xi = particle_distribution(generator);
                yi = particle_distribution(generator);
                zi = particle_distribution(generator);
            }
            x_[i] = xi;
            y_[i] = yi;
            z_[i] = zi;
            q_[i] = charge_distribution(generator);
            //q_[i] = 1.0;
        }

        reference_energies["open_reference"]               = 1.0;
        reference_energies["open"]                         = 1.0;
        reference_energies["periodic_reference"]           = 1.0;
        reference_energies["periodic"]                     = 1.0;
        reference_energies["dipole_correction_reference"]  = 1.0;
        reference_energies["dipole_correction"]            = 1.0;

        excl.resize(n_);
        for(int i = 0; i < n_; i++)
        {
            excl[i].resize(1);
            excl[i][0] = i;
        }
    }
    /*
    if(dataset_id == 2)
    {
        random_parts_2<T>* inputdata = new random_parts_2<T>();
        set(inputdata);
        reference_energies["open_reference"]               = -1.23076923076923083755e+00;
        reference_energies["open"]                         = -1.23076923076922462030e+00; // p=20 d=3
        reference_energies["periodic_reference"]           = -4.02820799496271675366e+00; // p2p reference, lattice p=20
        reference_energies["periodic"]                     = -4.02820823584402987194e+00; // p=20 d=3
        reference_energies["dipole_correction_reference"]  = -5.46274348076880045255e+00; // p2p reference, lattice p=20
        reference_energies["dipole_correction"]            = -5.41083645039080707306e+00; // p=20 d=3

        excl.resize(inputdata->n);
        for(int i = 0; i < inputdata->n; i++)
        {
            excl[i].resize(1);
            excl[i][0] = i;
        }
        delete inputdata;
    }
    */
    if(dataset_id == 2)
    {
        parts_25<T>* inputdata = new parts_25<T>();
        set(inputdata);
        reference_energies["open_reference"]               = -9.16537563777904296103e+00;
        reference_energies["open"]                         = -9.16579294232553998256e+00; // p=8 d=3
        reference_energies["periodic_reference"]           =  1.19957236118806633840e+02; // p2p reference, lattice p=50
        reference_energies["periodic"]                     =  1.19961909956159075819e+02; // p=8 d=3
        reference_energies["dipole_correction_reference"]  =  1.09923865145869385174e+02; // p2p reference, lattice p=50
        reference_energies["dipole_correction"]            =  1.10190043180439971593e+02; // p=20 d=3

        excl.resize(inputdata->n);
        for(int i = 0; i < inputdata->n; i++)
        {
            excl[i].resize(1);
            excl[i][0] = i;
        }
        delete inputdata;
    }

    if(dataset_id == 3)
    {
        nacl<T>* inputdata = new nacl<T>();
        set(inputdata);
        reference_energies["open_reference"]               = -2.84273038966285057541e+04;
        reference_energies["open"]                         = -2.84274341017965416540e+04; // p=8 d=3
        reference_energies["periodic_reference"]           = -1.747564594633182190636212e+00 * 32768/2; // madelung constant(-1.747564594633182190636212)
        reference_energies["periodic"]                     = -2.86322883818380651064e+04; // p=8 d=3
        reference_energies["dipole_correction_reference"]  = -2.86320981936067328206e+04; // p2p reference, lattice p=50
        reference_energies["dipole_correction"]            = -2.86322883818380651064e+04; // p=8 d=3

        excl.resize(inputdata->n);
        for(int i = 0; i < inputdata->n; i++)
        {
            excl[i].resize(1);
            excl[i][0] = i;
        }
        delete inputdata;
    }

    if(dataset_id == 4)
    {
        silicamelt<T>* inputdata = new silicamelt<T>();
        set(inputdata);
        reference_energies["open_reference"]               = -1.31499516063500341261e+05;
        reference_energies["open"]                         = -1.31499517738219496096e+05; // p=8 d=3
        reference_energies["periodic_reference"]           = -1.32969898959316225955e+05; // p2p reference, lattice p=50
        reference_energies["periodic"]                     = -1.32969925259017734788e+05; // p=8 d=3
        reference_energies["dipole_correction_reference"]  = -1.32985680484800104750e+05; // p2p reference, lattice p=50
        reference_energies["dipole_correction"]            = -1.32985996586082212161e+05; // p=8 d=3

        excl.resize(inputdata->n);
        for(int i = 0; i < inputdata->n; i++)
        {
            excl[i].resize(1);
            excl[i][0] = i;
        }
        delete inputdata;
    }

    if(dataset_id == 5)
    {
        teha<T>* inputdata = new teha<T>();
        set(inputdata);
        reference_energies["open_reference"]               = -5.99727737601174158044e+05;
        reference_energies["open"]                         = -5.99727736879982054234e+05; // p=8 d=3
        reference_energies["periodic_reference"]           = -6.05178529554840526544e+05; // p2p reference, lattice p=50
        reference_energies["periodic"]                     = -6.05178549943373189308e+05; // p=8 d=3
        reference_energies["dipole_correction_reference"]  = -6.05223191328209126368e+05; // p2p reference, lattice p=50
        reference_energies["dipole_correction"]            = -6.05226432500290917233e+05; // p=8 d=3

        excl.resize(inputdata->n);
        for(int i = 0; i < inputdata->n; i++)
        {
            excl[i].resize(1);
            excl[i][0] = i;
        }
        delete inputdata;
    }

    if(dataset_id == 6)
    {
        parts_2<T>* inputdata = new parts_2<T>();
        set(inputdata);
        reference_energies["open_reference"]               = 1.0;
        reference_energies["open"]                         = 1.0;
        reference_energies["periodic_reference"]           = 1.0;
        reference_energies["periodic"]                     = 1.0;
        reference_energies["dipole_correction_reference"]  = 1.0;
        reference_energies["dipole_correction"]            = 1.0;

        excl.resize(2);
        for(int i = 0; i < inputdata->n; i++)
        {
            excl[i].resize(1);
            excl[i][0] = i;
        }
        excl[0].resize(2);
        excl[1].resize(2);

        excl[0][0] = 0;
        excl[0][1] = 1;

        excl[1][0] = 0;
        excl[1][1] = 1;
        delete inputdata;
    }

    if(dataset_id == 7)
    {
        T shift = 0.125;
        parts_2<T>* inputdata = new parts_2<T>();
        for(int i = 0; i < inputdata->n; i++)
        {
            inputdata->x[i] += shift;
            if(inputdata->x[i] > inputdata->boxlength)
            {
                inputdata->x[i] -= inputdata->boxlength;
            }
            inputdata->y[i] += shift;
            if(inputdata->y[i] > inputdata->boxlength)
            {
                inputdata->y[i] -= inputdata->boxlength;
            }
            inputdata->z[i] += shift;
            if(inputdata->z[i] > inputdata->boxlength)
            {
                inputdata->z[i] -= inputdata->boxlength;
            }
        }
        set(inputdata);
        reference_energies["open_reference"]               = 1.0;
        reference_energies["open"]                         = 1.0;
        reference_energies["periodic_reference"]           = 1.0;
        reference_energies["periodic"]                     = 1.0;
        reference_energies["dipole_correction_reference"]  = 1.0;
        reference_energies["dipole_correction"]            = 1.0;

        excl.resize(2);
        for(int i = 0; i < inputdata->n; i++)
        {
            excl[i].resize(1);
            excl[i][0] = i;
        }
        excl[0].resize(2);
        excl[1].resize(2);

        excl[0][0] = 0;
        excl[0][1] = 1;

        excl[1][0] = 0;
        excl[1][1] = 1;
        delete inputdata;
    }

    if(dataset_id == 8)
    {
        T shift = 0.625;
        parts_2<T>* inputdata = new parts_2<T>();
        for(int i = 0; i < inputdata->n; i++)
        {
            inputdata->x[i] += shift;
            if(inputdata->x[i] > inputdata->boxlength)
            {
                inputdata->x[i] -= inputdata->boxlength;
            }
            inputdata->y[i] += shift;
            if(inputdata->y[i] > inputdata->boxlength)
            {
                inputdata->y[i] -= inputdata->boxlength;
            }
            inputdata->z[i] += shift;
            if(inputdata->z[i] > inputdata->boxlength)
            {
                inputdata->z[i] -= inputdata->boxlength;
            }
        }
        set(inputdata);
        reference_energies["open_reference"]               = 1.0;
        reference_energies["open"]                         = 1.0;
        reference_energies["periodic_reference"]           = 1.0;
        reference_energies["periodic"]                     = 1.0;
        reference_energies["dipole_correction_reference"]  = 1.0;
        reference_energies["dipole_correction"]            = 1.0;

        excl.resize(2);
        for(int i = 0; i < inputdata->n; i++)
        {
            excl[i].resize(1);
            excl[i][0] = i;
        }


        excl[0].resize(2);
        excl[1].resize(2);

        excl[0][0] = 0;
        excl[0][1] = 1;

        excl[1][0] = 0;
        excl[1][1] = 1;
        delete inputdata;
    }


    if(dataset_id == 9)
    {
        water_in<T>* inputdata = new water_in<T>();
        set(inputdata);
        reference_energies["open_reference"]               = 1.0;
        reference_energies["open"]                         = 1.0;
        reference_energies["periodic_reference"]           = 1.0;
        reference_energies["periodic"]                     = 1.0;
        reference_energies["dipole_correction_reference"]  = 1.0;
        reference_energies["dipole_correction"]            = 1.0;

        excl.resize(inputdata->n);
        for(int i = 0; i < inputdata->n; i++)
        {
            excl[i].resize(1);
            excl[i][0] = i;
        }
        for(int i = 0; i < inputdata->n; i+=3)
        {
            excl[i].resize(3);
            for(int j = 0; j < 3; j++)
            {
                excl[i][j] = i+j;
            }
        }
        delete inputdata;
    }

    if(dataset_id == 10)
    {
        T shift = 25.0;
        water_in<T>* inputdata = new water_in<T>();
        int out_of_box_particles = 0;
        for(int i = 0; i < inputdata->n; i++)
        {
            inputdata->x[i] += shift;
            if(inputdata->x[i] > inputdata->boxlength)
            {
                out_of_box_particles++;
                inputdata->x[i] -= inputdata->boxlength;
            }
            inputdata->y[i] += shift;
            if(inputdata->y[i] > inputdata->boxlength)
            {
                out_of_box_particles++;
                inputdata->y[i] -= inputdata->boxlength;
            }
            inputdata->z[i] += shift;
            if(inputdata->z[i] > inputdata->boxlength)
            {
                out_of_box_particles++;
                inputdata->z[i] -= inputdata->boxlength;
            }
        }
        printf("out of box = %d\n", out_of_box_particles);
        set(inputdata);
        reference_energies["open_reference"]               = 1.0;
        reference_energies["open"]                         = 1.0;
        reference_energies["periodic_reference"]           = 1.0;
        reference_energies["periodic"]                     = 1.0;
        reference_energies["dipole_correction_reference"]  = 1.0;
        reference_energies["dipole_correction"]            = 1.0;

        excl.resize(inputdata->n);
        for(int i = 0; i < inputdata->n; i++)
        {
            excl[i].resize(1);
            excl[i][0] = i;
        }
        for(int i = 0; i < inputdata->n; i+=3)
        {
            excl[i].resize(3);
            for(int j = 0; j < 3; j++)
            {
                excl[i][j] = i+j;
            }
        }
        delete inputdata;
    }

    if(boxsize.x == 0.0)
    {
        boxsize.x = boxlength_;
        boxsize.y = boxlength_;
        boxsize.z = boxlength_;
    }

    abc = ABC<XYZ<T> >(e1 * boxsize.x, e2 * boxsize.y, e3 * boxsize.z);
    box_centre = XYZ<T>((abc.a + abc.b + abc.c) * 0.5);

    bool runbenchmark = get_env_int(0, "RUNBENCHMARK") != 0;
    if(!runbenchmark)
    {
        printf("writing output pdb (fmmdata.pdb)\n");
        std::ofstream myfile;
        myfile.open("fmmdata.pdb");
        myfile <<"TITLE     Random t=   0.00000\n";
        myfile <<"REMARK    THIS IS A SIMULATION BOX\n";
        myfile <<"CRYST1";
        myfile<<std::fixed<<std::setprecision(3)<<std::setw(9)<<boxsize.x*10<<std::setw(9)<<boxsize.y*10<<std::setw(9)<<boxsize.z*10;
        myfile<<"  "<<"90.00"<<"  "<<"90.00"<<"  "<<"90.00"<<" P 1           1\n";
        myfile <<"MODEL        1\n";
        for(int i = 1; i < n_+1; i++)
        {
            myfile<<"ATOM"<<std::setw(7)<<i<<"  ";
            if(q_[i] < 0.0)
            {
                myfile<<"Na   Na";
            }
            else
            {
                myfile<<"Cl   Cl";
            }
            myfile<<"     1";
            myfile<<std::fixed<<std::setprecision(3)<<std::setw(12)<<x_[i-1]*10<<std::setw(8)<<y_[i-1]*10<<std::setw(8)<<z_[i-1]*10;
            myfile<<std::fixed<<std::setprecision(2)<<std::setw(6)<<q_[i-1];
            myfile <<"  0.00"<<"            \n";
        }
        myfile <<"TER\n";
        myfile <<"ENDMDL\n";
        myfile.close();
    }

    return 1;
}

template <typename T> template <typename Datatype>
void Testdata<T>::set(Datatype* inputdata )
{
    boxlength_ = inputdata->boxlength;
    n_ = inputdata->n;

    x_.resize(n_);
    y_.resize(n_);
    z_.resize(n_);
    q_.resize(n_);
    x_.assign(inputdata->x,inputdata->x + n_);
    y_.assign(inputdata->y,inputdata->y + n_);
    z_.assign(inputdata->z,inputdata->z + n_);
    q_.assign(inputdata->q,inputdata->q + n_);
}

template <typename T>
int Testdata<T>::n()
{
    return n_;
}

template <typename T>
T* Testdata<T>::x()
{
    return &x_[0];
}

template <typename T>
T* Testdata<T>::y()
{
    return &y_[0];
}

template <typename T>
T* Testdata<T>::z()
{
    return &z_[0];
}

template <typename T>
T* Testdata<T>::q()
{
    return &q_[0];
}

template class Testdata<float>;
template class Testdata<double>;


namespace lambda_testdata {

////#include "testdata/lambda_qxyz_random_data.hpp"
////#include "testdata/lambda_qxyz_random_data_5_sites.hpp"
////#include "testdata/lambda_qxyz_two_sites_two_variable_forms_nacl.hpp"
////#include "testdata/lambda_qxyz_one_site_two_forms_nacl.hpp"
////#include "testdata/lambda_qxyz_two_sites_nacl.hpp"

//const REAL *x,*y,*z,*q;

//ABC<XYZ<REAL> > abc(e1 * boxlength, e2 * boxlength, e3 * boxlength);
//XYZ<REAL> reference_point((abc.a + abc.b + abc.c) * 0.5);

}

}//namespace end
