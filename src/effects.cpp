#include "effects.hpp"

namespace evo
{
    int Effects::set(matrix<int> &in_effect)
    {
        try
        {
            i_effect = in_effect;
            l_effect.resize(1);
            b_effect.resize(1);
            f_effect.resize(1);
            d_effect.resize(1);

            type = 1;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Effects::set(matrix<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Effects::set(matrix<int> &)." << '\n';
            throw;
        }

        return 0;
    }

    int Effects::set(matrix<size_t> &in_effect)
    {
        try
        {
            l_effect = in_effect;
            i_effect.resize(1);
            b_effect.resize(1);
            f_effect.resize(1);
            d_effect.resize(1);

            type = 2;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Effects::set(matrix<size_t> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Effects::set(matrix<size_t> &)." << '\n';
            throw;
        }

        return 0;
    }

    int Effects::set(matrix<bool> &in_effect)
    {
        try
        {
            b_effect = in_effect;
            l_effect.resize(1);
            i_effect.resize(1);
            f_effect.resize(1);
            d_effect.resize(1);

            type = 3;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Effects::set(matrix<bool> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Effects::set(matrix<bool> &)." << '\n';
            throw;
        }

        return 0;
    }

    int Effects::set(matrix<float> &in_effect)
    {
        try
        {
            f_effect = in_effect;
            l_effect.resize(1);
            b_effect.resize(1);
            i_effect.resize(1);
            d_effect.resize(1);

            type = 4;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Effects::set(matrix<float> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Effects::set(matrix<float> &)." << '\n';
            throw;
        }

        return 0;
    }

    int Effects::set(matrix<double> &in_effect)
    {
        try
        {
            d_effect = in_effect;
            l_effect.resize(1);
            b_effect.resize(1);
            f_effect.resize(1);
            i_effect.resize(1);

            type = 5;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Effects::set(matrix<double> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Effects::set(matrix<double> &)." << '\n';
            throw;
        }

        return 0;
    }

    int Effects::get(matrix<int> &out_effect)
    {
        try
        {
            out_effect = i_effect;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Effects::get(matrix<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Effects::get(matrix<int> &)." << '\n';
            throw;
        }

        return 0;
    }

    int Effects::get(matrix<size_t> &out_effect)
    {
        try
        {
            // l_effect.fget();
            out_effect = l_effect;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Effects::get(matrix<size_t> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Effects::get(matrix<size_t> &)." << '\n';
            throw;
        }

        return 0;
    }

    int Effects::get(matrix<bool> &out_effect)
    {
        try
        {
            out_effect = b_effect;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Effects::get(matrix<bool> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Effects::get(matrix<bool> &)." << '\n';
            throw;
        }

        return 0;
    }

    int Effects::get(matrix<float> &out_effect)
    {
        try
        {
            out_effect = f_effect;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Effects::get(matrix<float> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Effects::get(matrix<float> &)." << '\n';
            throw;
        }

        return 0;
    }

    int Effects::get(matrix<double> &out_effect)
    {
        try
        {
            out_effect = d_effect;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Effects::get(matrix<double> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Effects::get(matrix<double> &)." << '\n';
            throw;
        }

        return 0;
    }

    matrix<float> Effects::get_float()
    {
        try
        {
            matrix<float> z;

            if (type == 1)
            {
                matrix<int> i_z;
                get(i_z);
                i_z.fread();
                z = i_z._float();
            }

            if (type == 2)
            {
                matrix<size_t> i_z;
                get(i_z);
                i_z.fread();
                z = i_z._float();
            }

            if (type == 3)
            {
                matrix<bool> i_z;
                get(i_z);
                i_z.fread();
                z = i_z._float();
            }

            if (type == 4)
            {
                get(z);
                z.fread();
            }

            if (type == 5)
            {
                matrix<double> i_z;
                get(i_z);
                i_z.fread();
                z = i_z._float();
            }

            return z;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Effects::get_float()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Effects::get_float()." << '\n';
            throw;
        }
    }

    void Effects::fread()
    {
        try
        {
            if (type == 1)
                i_effect.fread();
            if (type == 2)
                l_effect.fread();
            if (type == 3)
                b_effect.fread();
            if (type == 4)
                f_effect.fread();
            if (type == 5)
                d_effect.fread();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Effects::fread()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Effects::fread()." << '\n';
            throw;
        }
    }

    void Effects::fwrite()
    {
        try
        {
            if (type == 1)
                i_effect.fwrite();
            if (type == 2)
                l_effect.fwrite();
            if (type == 3)
                b_effect.fwrite();
            if (type == 4)
                f_effect.fwrite();
            if (type == 5)
                d_effect.fwrite();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Effects::fwrite()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Effects::fwrite()." << '\n';
            throw;
        }
    }

    matrix<float> Effects::fget(size_t irow[], size_t icol[])
    {
        matrix<float> vect;

        try
        {
            if (type == 1)
                i_effect.cast_fget(irow, icol, vect);
            if (type == 2)
                l_effect.cast_fget(irow, icol, vect);
            if (type == 3)
                b_effect.cast_fget(irow, icol, vect);
            if (type == 4)
                f_effect.cast_fget(irow, icol, vect);
            if (type == 5)
                d_effect.cast_fget(irow, icol, vect);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Effects::fget(size_t, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Effects::fget(size_t, size_t)." << '\n';
            throw;
        }

        return vect;
    }

    void Effects::vect_fget(size_t irow[], size_t icol[], std::vector<std::vector<float>> &vect)
    {
        try
        {
            if (type == 1)
                i_effect.cast_fget(irow, icol, vect);
            if (type == 2)
                l_effect.cast_fget(irow, icol, vect);
            if (type == 3)
                b_effect.cast_fget(irow, icol, vect);
            if (type == 4)
                f_effect.cast_fget(irow, icol, vect);
            if (type == 5)
                d_effect.cast_fget(irow, icol, vect);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Effects::vect_fget(size_t, size_t, std::vector<std::vector<float>> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Effects::vect_fget(size_t, size_t, std::vector<std::vector<float>> &)." << '\n';
            throw;
        }
    }

    void Effects::vect_fget(size_t irow[], size_t icol[], float **vect)
    {
        try
        {
            if (type == 1)
                i_effect.cast_fget(irow, icol, vect);
            if (type == 2)
                l_effect.cast_fget(irow, icol, vect);
            if (type == 3)
                b_effect.cast_fget(irow, icol, vect);
            if (type == 4)
                f_effect.cast_fget(irow, icol, vect);
            if (type == 5)
                d_effect.cast_fget(irow, icol, vect);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Effects::vect_fget(std::vector<std::vector<float>> &, std::vector<std::vector<float>> &, float **)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Effects::vect_fget(std::vector<std::vector<float>> &, std::vector<std::vector<float>> &, float **)." << '\n';
            throw;
        }
    }

    int Effects::print(std::string msg)
    {
        try
        {
            if (type == 1)
                i_effect.print(msg);
            if (type == 2)
                l_effect.print(msg);
            if (type == 3)
                b_effect.print(msg);
            if (type == 4)
                f_effect.print(msg);
            if (type == 5)
                d_effect.print(msg);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Effects::print(std::string)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Effects::print(std::string)." << '\n';
            throw;
        }

        return 0;
    }

    matrix<size_t> Effects::shape()
    {
        try
        {
            matrix<size_t> shape(2, 1);
            if (type == 1)
                shape = i_effect.shape();
            if (type == 2)
                shape = l_effect.shape();
            if (type == 3)
                shape = b_effect.shape();
            if (type == 4)
                shape = f_effect.shape();
            if (type == 5)
                shape = d_effect.shape();

            return shape;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Effects::shape()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Effects::shape()." << '\n';
            throw;
        }
    }

    size_t Effects::size()
    {
        try
        {
            size_t sz = 0;

            if (type == 1)
                sz = i_effect.size();
            if (type == 2)
                sz = l_effect.size();
            if (type == 3)
                sz = b_effect.size();
            if (type == 4)
                sz = f_effect.size();
            if (type == 5)
                sz = d_effect.size();

            return sz;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Effects::size()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Effects::size()." << '\n';
            throw;
        }
    }

    int Effects::clear()
    {
        try
        {
            i_effect.fclear();
            i_effect.clear();
            l_effect.fclear();
            l_effect.clear();
            b_effect.fclear();
            b_effect.clear();
            f_effect.fclear();
            f_effect.clear();
            d_effect.fclear();
            d_effect.clear();

            return 0;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Effects::clear()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Effects::clear()." << '\n';
            throw;
        }
    }

}