/*
    cs_matrix.cpp

    Explicit instantiation declaration of matrix class.

*/

#include "cs_matrix.hpp"

namespace evo
{
    extern template class matrix<float>;
    extern template class matrix<double>;
}
