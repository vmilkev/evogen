#ifndef effects_hpp__
#define effects_hpp__

#include "cs_matrix.hpp"

namespace evo
{
    class Effects
    {

    private:
        matrix<int> i_effect;
        matrix<size_t> l_effect;
        matrix<bool> b_effect;
        matrix<float> f_effect;
        matrix<double> d_effect;

    public:
        int type = 0;

        int set(matrix<int> &in_effect);
        int set(matrix<size_t> &in_effect);
        int set(matrix<bool> &in_effect);
        int set(matrix<float> &in_effect);
        int set(matrix<double> &in_effect);

        int get(matrix<int> &out_effect);
        int get(matrix<size_t> &out_effect);
        int get(matrix<bool> &out_effect);
        int get(matrix<float> &out_effect);
        int get(matrix<double> &out_effect);

        matrix<float> get_float();

        // Interfaces to the matrix class
        void fread();
        void fwrite();
        matrix<float> fget(size_t irow[], size_t icol[]);
        void vect_fget(size_t irow[], size_t icol[], std::vector<std::vector<float>> &vect);
        void vect_fget(size_t irow[], size_t icol[], float **vect);
        int print(std::string msg);
        matrix<size_t> shape();
        size_t size();
        int clear();
    };

}

#endif // effects_hpp__
