#include "model.hpp"

namespace evo
{

#ifdef PYBIND

    int Model::append_residual(pybind11::array_t<float> arr, size_t lda)
    {
        try
        {
            /* lda := is leading diagonal of symmetric matrix.

               NOTE: accepting a full format matrix (not a compact format)!
            */

            matrix<float> residual;

            pybind11::buffer_info buf1 = arr.request();

            if (buf1.ndim != 1)
                throw std::runtime_error("Number of dimensions must be one");

            float *ptr1 = static_cast<float *>(buf1.ptr);

            residual.resize(lda, lda);
            // residual.rectosym();

            for (size_t i = 0; i < buf1.shape[0]; i++)
                residual[i] = ptr1[i];

            residual.fwrite();
            residuals.push_back(residual);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::append_residual(pybind11::array_t<float>, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::append_residual(pybind11::array_t<float>, size_t)." << '\n';
            throw;
        }

        return 0;
    }

    int Model::append_observation(pybind11::array_t<float> arr, size_t lda)
    {
        try
        {
            matrix<float> observation;

            // lda := is size of vector
            pybind11::buffer_info buf1 = arr.request();

            if (buf1.ndim != 1)
                throw std::runtime_error("Number of dimensions must be one");

            float *ptr1 = static_cast<float *>(buf1.ptr);

            observation.resize(lda, 1);

            for (size_t i = 0; i < buf1.shape[0]; i++)
                observation[i] = ptr1[i];

            observation.fwrite();
            observations.push_back(observation);
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Model::append_observation(pybind11::array_t<float>, size_t): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::append_observation(pybind11::array_t<float>, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::append_observation(pybind11::array_t<float>, size_t)." << '\n';
            throw;
        }

        return 0;
    }

    int Model::append_effect(pybind11::array_t<float> arr, size_t lda1, size_t lda2)
    {
        try
        {
            // lda1 := number of rows
            // lda2 := number of cols

            matrix<float> effect;
            Effects e;

            pybind11::buffer_info buf1 = arr.request();

            if (buf1.ndim != 1)
                throw std::runtime_error("Number of dimensions must be one");

            float *ptr1 = static_cast<float *>(buf1.ptr);

            effect.resize(lda1, lda2);

            for (size_t i = 0; i < buf1.shape[0]; i++)
            {
                effect[i] = ptr1[i];
            }

            effect.transpose();

            effect.fwrite();
            e.set(effect);

            // effects.push_back(effect);
            effects.push_back(e);

            return 0;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Model::append_effect(pybind11::array_t<float>, size_t, size_t): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::append_effect(pybind11::array_t<float>, size_t, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::append_effect(pybind11::array_t<float>, size_t, size_t)." << '\n';
            throw;
        }

        return 0;
    }

    int Model::append_corrstruct(pybind11::array_t<float> var, size_t lda1, pybind11::array_t<float> corr, size_t lda2, pybind11::array_t<int> which_effects)
    {
        try
        {
            matrix<float> variance;
            matrix<float> correlation;
            matrix<int> effects;

            pybind11::buffer_info buf1 = var.request();
            pybind11::buffer_info buf2 = corr.request();
            pybind11::buffer_info buf3 = which_effects.request();

            if (buf1.ndim != 1 || buf2.ndim != 1 || buf3.ndim != 1)
                throw std::runtime_error("Number of dimensions must be one");

            float *ptr1 = static_cast<float *>(buf1.ptr);
            float *ptr2 = static_cast<float *>(buf2.ptr);
            int *ptr3 = static_cast<int *>(buf3.ptr);

            effects.resize(buf3.shape[0], 1);
            variance.resize(lda1, lda1);
            // variance.rectosym();
            correlation.resize(lda2, lda2);
            // correlation.rectosym();

            for (size_t i = 0; i < buf3.shape[0]; i++)
                effects[i] = ptr3[i];

            for (size_t i = 0; i < buf1.shape[0]; i++)
                variance[i] = ptr1[i];

            for (size_t i = 0; i < buf2.shape[0]; i++)
                correlation[i] = ptr2[i];

            effects.fwrite();
            correlated_effects.push_back(effects);

            variance.fwrite();
            ;
            variances.push_back(variance);

            correlation.fwrite();
            correlations.push_back(correlation);

            identity_correlations.push_back(false);
            identity_dimension.push_back(0);
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Model::append_corrstruct(pybind11::array_t<float>, size_t, pybind11::array_t<float>, size_t, pybind11::array_t<int>): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::append_corrstruct(pybind11::array_t<float>, size_t, pybind11::array_t<float>, size_t, pybind11::array_t<int>)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::append_corrstruct(pybind11::array_t<float>, size_t, pybind11::array_t<float>, size_t, pybind11::array_t<int>)." << '\n';
            throw;
        }

        return 0;
    }

    int Model::append_corrstruct(pybind11::array_t<float> var, size_t lda1, std::string &identity, size_t lda2, pybind11::array_t<int> which_effects)
    {
        try
        {
            std::string identity_string("I");

            if (identity_string.compare(identity) != 0)
                throw std::runtime_error("Error in the third input argument, the expected is I");

            matrix<float> variance;
            matrix<float> correlation;
            matrix<int> effects;

            pybind11::buffer_info buf1 = var.request();
            // pybind11::buffer_info buf2 = corr.request();
            pybind11::buffer_info buf3 = which_effects.request();

            if (buf1.ndim != 1 || buf3.ndim != 1)
                throw std::runtime_error("Number of dimensions must be one");

            float *ptr1 = static_cast<float *>(buf1.ptr);
            // float *ptr2 = static_cast<float *>(buf2.ptr);
            int *ptr3 = static_cast<int *>(buf3.ptr);

            effects.resize(buf3.shape[0], 1);
            variance.resize(lda1, lda1);
            // variance.rectosym();

            correlation.resize(1, 1);

            for (size_t i = 0; i < buf3.shape[0]; i++)
                effects[i] = ptr3[i];

            for (size_t i = 0; i < buf1.shape[0]; i++)
                variance[i] = ptr1[i];

            correlation[0] = 1.0;

            effects.fwrite();
            correlated_effects.push_back(effects);

            variance.fwrite();
            ;
            variances.push_back(variance);

            correlation.fwrite();
            correlations.push_back(correlation);

            identity_correlations.push_back(true);
            identity_dimension.push_back(lda2);
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Model::append_corrstruct(pybind11::array_t<float>, size_t, std::string &, size_t, pybind11::array_t<int>): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::append_corrstruct(pybind11::array_t<float>, size_t, std::string &, size_t, pybind11::array_t<int>)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::append_corrstruct(pybind11::array_t<float>, size_t, std::string &, size_t, pybind11::array_t<int>)." << '\n';
            throw;
        }

        return 0;
    }

    int Model::append_corrstruct(const std::string &fname_var, std::string &identity, size_t lda2, pybind11::array_t<int> which_effects)
    {
        try
        {
            std::string identity_string("I");

            if (identity_string.compare(identity) != 0)
                throw std::runtime_error("Error in the third input argument, the expected is I");

            IOInterface datstream;

            datstream.set_fname(fname_var);
            std::vector<std::vector<float>> var;
            datstream.fgetdata(var);

            size_t lda1 = var.size();

            if (lda2 == 0)
                throw "The first dimension of the arrey CORR is 0, expected et least 1!";

            if (lda1 == 0)
                throw "The first dimension of the arrey VAR is 0, expected et least 1!";

            if (var[0].size() == 0)
                throw "The second dimension of the arrey VAR is 0, expected et least 1!";

            if (lda1 != var[0].size())
                throw "The arrey VAR is not square, that is expected in this case!";

            matrix<float> variance;
            matrix<float> correlation;
            matrix<int> effects;

            pybind11::buffer_info buf3 = which_effects.request();

            if (buf3.ndim != 1)
                throw std::runtime_error("Number of dimensions must be one");

            int *ptr3 = static_cast<int *>(buf3.ptr);

            effects.resize(buf3.shape[0], 1);

            variance.resize(lda1, lda1);

            correlation.resize(1, 1);

            for (size_t i = 0; i < buf3.shape[0]; i++)
                effects[i] = ptr3[i];

            for (size_t i = 0; i < lda1; i++)
            {
                for (size_t j = 0; j < lda1; j++)
                    variance(i, j) = var[i][j];
            }

            correlation[0] = 1.0;

            effects.fwrite();
            correlated_effects.push_back(effects);

            variance.fwrite();
            ;
            variances.push_back(variance);

            correlation.fwrite();
            correlations.push_back(correlation);

            identity_correlations.push_back(true);
            identity_dimension.push_back(lda2);

            var.clear();
            var.shrink_to_fit();
        }
        catch (std::string err)
        {
            std::cerr << "Exception in Model::append_corrstruct(const std::string &, std::string &, size_t, pybind11::array_t<int>)." << '\n';
            std::cerr << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::append_corrstruct(const std::string &, std::string &, size_t, pybind11::array_t<int>)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::append_corrstruct(const std::string &, std::string &, size_t, pybind11::array_t<int>)." << '\n';
            throw;
        }

        return 0;
    }

    int Model::append_corrstruct(pybind11::array_t<float> var, size_t lda1, const std::string &fname_corr, pybind11::array_t<int> which_effects)
    {
        try
        {
            IOInterface datstream;

            datstream.set_fname(fname_corr);
            std::vector<std::vector<float>> corr;
            datstream.fgetdata(corr);

            size_t lda2 = corr.size();

            if (lda2 == 0)
                throw "The first dimension of the arrey CORR is 0, expected et least 1!";

            if (corr[0].size() == 0)
                throw "The second dimension of the arrey CORR is 0, expected et least 1!";

            if (lda2 != corr[0].size())
                throw "The arrey CORR is not square, that is expected in this case!";

            matrix<float> variance;
            matrix<float> correlation;
            matrix<int> effects;

            pybind11::buffer_info buf1 = var.request();
            // pybind11::buffer_info buf2 = corr.request();
            pybind11::buffer_info buf3 = which_effects.request();

            if (buf1.ndim != 1 || /*buf2.ndim != 1 ||*/ buf3.ndim != 1)
                throw std::runtime_error("Number of dimensions must be one");

            float *ptr1 = static_cast<float *>(buf1.ptr);
            // float *ptr2 = static_cast<float *>(buf2.ptr);
            int *ptr3 = static_cast<int *>(buf3.ptr);

            effects.resize(buf3.shape[0], 1);
            variance.resize(lda1, lda1);
            // variance.rectosym();
            correlation.resize(lda2, lda2);
            // correlation.rectosym();

            for (size_t i = 0; i < buf3.shape[0]; i++)
                effects[i] = ptr3[i];

            for (size_t i = 0; i < buf1.shape[0]; i++)
                variance[i] = ptr1[i];

            for (size_t i = 0; i < lda2; i++)
            {
                for (size_t j = 0; j < lda2; j++)
                    correlation(i, j) = corr[i][j];
            }

            effects.fwrite();
            correlated_effects.push_back(effects);

            variance.fwrite();
            ;
            variances.push_back(variance);

            correlation.fwrite();
            correlations.push_back(correlation);

            identity_correlations.push_back(false);
            identity_dimension.push_back(0);

            corr.clear();
            corr.shrink_to_fit();
        }
        catch (std::string err)
        {
            std::cerr << "Exception in Model::append_corrstruct(pybind11::array_t<float>, size_t, const std::string &, pybind11::array_t<int>)." << '\n';
            std::cerr << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::append_corrstruct(pybind11::array_t<float>, size_t, const std::string &, pybind11::array_t<int>)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::append_corrstruct(pybind11::array_t<float>, size_t, const std::string &, pybind11::array_t<int>)." << '\n';
            throw;
        }

        return 0;
    }

    int Model::append_corrstruct(const std::string &fname_var, const std::string &fname_corr, pybind11::array_t<int> which_effects)
    {
        try
        {
            IOInterface datstream;

            datstream.set_fname(fname_corr);
            std::vector<std::vector<float>> corr;
            datstream.fgetdata(corr);

            datstream.set_fname(fname_var);
            std::vector<std::vector<float>> var;
            datstream.fgetdata(var);

            size_t lda2 = corr.size();
            size_t lda1 = var.size();

            if (lda2 == 0)
                throw "The first dimension of the arrey CORR is 0, expected et least 1!";

            if (corr[0].size() == 0)
                throw "The second dimension of the arrey CORR is 0, expected et least 1!";

            if (lda2 != corr[0].size())
                throw "The arrey CORR is not square, that is expected in this case!";

            if (lda1 == 0)
                throw "The first dimension of the arrey VAR is 0, expected et least 1!";

            if (var[0].size() == 0)
                throw "The second dimension of the arrey VAR is 0, expected et least 1!";

            if (lda1 != var[0].size())
                throw "The arrey VAR is not square, that is expected in this case!";

            matrix<float> variance;
            matrix<float> correlation;
            matrix<int> effects;

            // pybind11::buffer_info buf1 = var.request();
            //  pybind11::buffer_info buf2 = corr.request();
            pybind11::buffer_info buf3 = which_effects.request();

            if (/*buf1.ndim != 1 || buf2.ndim != 1 ||*/ buf3.ndim != 1)
                throw std::runtime_error("Number of dimensions must be one");

            // float *ptr1 = static_cast<float *>(buf1.ptr);
            //  float *ptr2 = static_cast<float *>(buf2.ptr);
            int *ptr3 = static_cast<int *>(buf3.ptr);

            effects.resize(buf3.shape[0], 1);
            variance.resize(lda1, lda1);
            // variance.rectosym();
            correlation.resize(lda2, lda2);
            // correlation.rectosym();

            for (size_t i = 0; i < buf3.shape[0]; i++)
                effects[i] = ptr3[i];

            for (size_t i = 0; i < lda1; i++)
            {
                for (size_t j = 0; j < lda1; j++)
                    variance(i, j) = var[i][j];
            }

            for (size_t i = 0; i < lda2; i++)
            {
                for (size_t j = 0; j < lda2; j++)
                    correlation(i, j) = corr[i][j];
            }

            effects.fwrite();
            correlated_effects.push_back(effects);

            variance.fwrite();
            ;
            variances.push_back(variance);

            correlation.fwrite();
            correlations.push_back(correlation);

            identity_correlations.push_back(false);
            identity_dimension.push_back(0);

            corr.clear();
            corr.shrink_to_fit();

            var.clear();
            var.shrink_to_fit();
        }
        catch (std::string err)
        {
            std::cerr << "Exception in Model::append_corrstruct(const std::string &, const std::string &, pybind11::array_t<int>)." << '\n';
            std::cerr << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::append_corrstruct(const std::string &, const std::string &, pybind11::array_t<int>)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::append_corrstruct(const std::string &, const std::string &, pybind11::array_t<int>)." << '\n';
            throw;
        }

        return 0;
    }

    int Model::append_traitstruct(int obs_id, pybind11::array_t<int> eff_id)
    {
        try
        {
            matrix<int> eff;

            pybind11::buffer_info buf2 = eff_id.request();

            if (buf2.ndim != 1)
                throw std::runtime_error("Number of dimensions must be one");

            int *ptr2 = static_cast<int *>(buf2.ptr);

            eff.resize(buf2.shape[0], 1);

            for (size_t i = 0; i < buf2.shape[0]; i++)
                eff[i] = ptr2[i];

            observation_trait.push_back(obs_id);
            effects_trait.push_back(eff);
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Model::append_traitstruct(int, pybind11::array_t<int>): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::append_traitstruct(int, pybind11::array_t<int>)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::append_traitstruct(int, pybind11::array_t<int>)." << '\n';
            throw;
        }

        return 0;
    }

#else

    int Model::append_residual(const std::vector<float> &arr, size_t lda)
    {
        try
        {
            // lda := is leading diagonal of symmetric matrix, though arr is in a full-store format

            matrix<float> residual;

            size_t n = sizeof(arr) / sizeof(arr[0]);
            residual.resize(lda, lda);
            // residual.rectosym();
            for (size_t i = 0; i < n; i++)
                residual[i] = arr[i];

            residual.fwrite();
            residuals.push_back(residual);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::append_residual(const std::vector<float> &, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::append_residual(const std::vector<float> &, size_t)." << '\n';
            throw;
        }

        return 0;
    }

    int Model::append_observation(const std::vector<float> &arr, size_t lda)
    {
        try
        {
            matrix<float> observation;

            // lda := is size of vector

            observation.resize(lda, 1);

            for (size_t i = 0; i < lda; i++)
                observation[i] = arr[i];

            observation.fwrite();
            observations.push_back(observation);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::append_observation(const std::vector<float> &, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::append_observation(const std::vector<float> &, size_t)." << '\n';
            throw;
        }

        return 0;
    }

    template <typename T>
    int Model::append_effect(const std::vector<T> &arr, size_t lda1, size_t lda2)
    {
        try
        {
            matrix<T> effect;
            Effects e;

            // lda1 := number of rows
            // lda2 := number of cols

            effect.resize(lda1, lda2);

            for (size_t i = 0; i < lda1 * lda2; i++)
                effect[i] = arr[i];

            effect.transpose();

            effect.fwrite();

            e.set(effect);

            effects.push_back(e);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::append_effect(const std::vector<T> &, size_t, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::append_effect(const std::vector<T> &, size_t, size_t)." << '\n';
            throw;
        }

        return 0;
    }

    template int Model::append_effect(const std::vector<int> &arr, size_t lda1, size_t lda2);
    template int Model::append_effect(const std::vector<size_t> &arr, size_t lda1, size_t lda2);
    template int Model::append_effect(const std::vector<bool> &arr, size_t lda1, size_t lda2);
    template int Model::append_effect(const std::vector<float> &arr, size_t lda1, size_t lda2);
    template int Model::append_effect(const std::vector<double> &arr, size_t lda1, size_t lda2);

    int Model::append_corrstruct(const std::vector<float> &var, size_t lda1, const std::vector<float> &corr, size_t lda2, const std::vector<int> &which_effects)
    {
        try
        {
            matrix<float> variance;
            matrix<float> correlation;
            matrix<int> effects;

            effects.resize(which_effects.size(), 1);
            variance.resize(lda1, lda1);
            // variance.rectosym();
            correlation.resize(lda2, lda2);
            // correlation.rectosym();

            for (size_t i = 0; i < which_effects.size(); i++)
                effects[i] = which_effects[i];

            for (size_t i = 0; i < var.size(); i++)
                variance[i] = var[i];

            for (size_t i = 0; i < corr.size(); i++)
                correlation[i] = corr[i];

            variance.fwrite();

            variances.push_back(variance);

            correlation.fwrite();
            correlations.push_back(correlation);

            identity_correlations.push_back(false);
            identity_dimension.push_back(0);

            effects.fwrite();
            correlated_effects.push_back(effects);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::append_corrstruct(const std::vector<float> &, size_t, const std::vector<float> &, size_t, const std::vector<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::append_corrstruct(const std::vector<float> &, size_t, const std::vector<float> &, size_t, const std::vector<int> &)." << '\n';
            throw;
        }

        return 0;
    }

    int Model::append_corrstruct(const std::vector<float> var, size_t lda1, std::string &identity, size_t lda2, const std::vector<int> which_effects)
    {
        try
        {
            std::string identity_string("I");

            if (identity_string.compare(identity) != 0)
                throw std::runtime_error("Error in the third input argument, the expected is I");

            matrix<float> variance;
            matrix<float> correlation;
            matrix<int> effects;

            effects.resize(which_effects.size(), 1);

            variance.resize(lda1, lda1);

            correlation.resize(1, 1);

            for (size_t i = 0; i < which_effects.size(); i++)
                effects[i] = which_effects[i];

            for (size_t i = 0; i < var.size(); i++)
                variance[i] = var[i];

            correlation[0] = 1.0;

            variance.fwrite();

            variances.push_back(variance);

            correlation.fwrite();
            correlations.push_back(correlation);

            identity_correlations.push_back(true);
            identity_dimension.push_back(lda2);

            effects.fwrite();
            correlated_effects.push_back(effects);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::append_corrstruct(const std::vector<float>, size_t, std::string &, size_t, const std::vector<int>)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::append_corrstruct(const std::vector<float>, size_t, std::string &, size_t, const std::vector<int>)." << '\n';
            throw;
        }

        return 0;
    }

    int Model::append_corrstruct(const std::string &fname_var, std::string &identity, size_t lda2, const std::vector<int> which_effects)
    {
        try
        {
            std::string identity_string("I");

            if (identity_string.compare(identity) != 0)
                throw std::runtime_error("Error in the third input argument, the expected is I");

            IOInterface datstream;

            datstream.set_fname(fname_var);
            std::vector<std::vector<float>> var;
            datstream.fgetdata(var);

            size_t lda1 = var.size();

            if (lda2 == 0)
                throw "The first dimension of the arrey CORR is 0, expected et least 1!";

            if (lda1 == 0)
                throw "The first dimension of the arrey VAR is 0, expected et least 1!";

            if (var[0].size() == 0)
                throw "The second dimension of the arrey VAR is 0, expected et least 1!";

            if (lda1 != var[0].size())
                throw "The arrey VAR is not square, that is expected in this case!";

            matrix<float> variance;
            matrix<float> correlation;
            matrix<int> effects;

            effects.resize(which_effects.size(), 1);

            variance.resize(lda1, lda1);

            correlation.resize(1, 1);

            for (size_t i = 0; i < which_effects.size(); i++)
                effects[i] = which_effects[i];

            for (size_t i = 0; i < lda1; i++)
            {
                for (size_t j = 0; j < lda1; j++)
                    variance(i, j) = var[i][j];
            }

            correlation[0] = 1.0;

            variance.fwrite();

            variances.push_back(variance);

            correlation.fwrite();
            correlations.push_back(correlation);

            identity_correlations.push_back(true);
            identity_dimension.push_back(lda2);

            effects.fwrite();
            correlated_effects.push_back(effects);

            var.clear();
            var.shrink_to_fit();
        }
        catch (std::string err)
        {
            std::cerr << "Exception in Model::append_corrstruct(const std::string &, std::string &, size_t, const std::vector<int>)." << '\n';
            std::cerr << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::append_corrstruct(const std::string &, std::string &, size_t, const std::vector<int>)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::append_corrstruct(const std::string &, std::string &, size_t, const std::vector<int>)." << '\n';
            throw;
        }

        return 0;
    }

    int Model::append_corrstruct(const std::vector<float> &var, size_t lda1, const std::string &fname_corr, const std::vector<int> &which_effects)
    {
        try
        {
            IOInterface datstream;

            datstream.set_fname(fname_corr);

            std::vector<std::vector<float>> corr;

            datstream.fgetdata(corr);

            size_t lda2 = corr.size();

            if (lda2 == 0)
                throw "The first dimension of the arrey is 0, expected et least 1!";

            if (corr[0].size() == 0)
                throw "The second dimension of the arrey is 0, expected et least 1!";

            if (lda2 != corr[0].size())
                throw "The arrey is not square, that is expected in this case!";

            matrix<float> variance;
            matrix<float> correlation;
            matrix<int> effects;

            effects.resize(which_effects.size(), 1);
            variance.resize(lda1, lda1);
            // variance.rectosym();
            correlation.resize(lda2, lda2);
            // correlation.rectosym();

            for (size_t i = 0; i < which_effects.size(); i++)
                effects[i] = which_effects[i];

            for (size_t i = 0; i < var.size(); i++)
                variance[i] = var[i];

            for (size_t i = 0; i < lda2; i++)
            {
                for (size_t j = 0; j < lda2; j++)
                    correlation(i, j) = corr[i][j];
            }

            variance.fwrite();

            variances.push_back(variance);

            correlation.fwrite();
            correlations.push_back(correlation);

            identity_correlations.push_back(false);
            identity_dimension.push_back(0);

            effects.fwrite();
            correlated_effects.push_back(effects);

            corr.clear();
            corr.shrink_to_fit();
        }
        catch (std::string err)
        {
            std::cerr << "Exception in Model::append_corrstruct(const std::vector<float> &, size_t, const std::string &, const std::vector<int> &)." << '\n';
            std::cerr << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::append_corrstruct(const std::vector<float> &, size_t, const std::string &, const std::vector<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::append_corrstruct(const std::vector<float> &, size_t, const std::string &, const std::vector<int> &)." << '\n';
            throw;
        }

        return 0;
    }

    int Model::append_corrstruct(const std::string &fname_var, const std::string &fname_corr, const std::vector<int> &which_effects)
    {
        try
        {
            IOInterface datstream;

            datstream.set_fname(fname_corr);
            std::vector<std::vector<float>> corr;
            datstream.fgetdata(corr);

            datstream.set_fname(fname_var);
            std::vector<std::vector<float>> var;
            datstream.fgetdata(var);

            size_t lda2 = corr.size();
            size_t lda1 = var.size();

            if (lda2 == 0)
                throw "The first dimension of the arrey CORR is 0, expected et least 1!";

            if (corr[0].size() == 0)
                throw "The second dimension of the arrey CORR is 0, expected et least 1!";

            if (lda2 != corr[0].size())
                throw "The arrey CORR is not square, that is expected in this case!";

            if (lda1 == 0)
                throw "The first dimension of the arrey VAR is 0, expected et least 1!";

            if (var[0].size() == 0)
                throw "The second dimension of the arrey VAR is 0, expected et least 1!";

            if (lda1 != var[0].size())
                throw "The arrey VAR is not square, that is expected in this case!";

            matrix<float> variance;
            matrix<float> correlation;
            matrix<int> effects;

            effects.resize(which_effects.size(), 1);
            variance.resize(lda1, lda1);
            // variance.rectosym();
            correlation.resize(lda2, lda2);
            // correlation.rectosym();

            for (size_t i = 0; i < which_effects.size(); i++)
                effects[i] = which_effects[i];

            for (size_t i = 0; i < lda1; i++)
            {
                for (size_t j = 0; j < lda1; j++)
                    variance(i, j) = var[i][j];
            }

            for (size_t i = 0; i < lda2; i++)
            {
                for (size_t j = 0; j < lda2; j++)
                    correlation(i, j) = corr[i][j];
            }

            variance.fwrite();
            ;
            variances.push_back(variance);

            correlation.fwrite();
            correlations.push_back(correlation);

            identity_correlations.push_back(false);
            identity_dimension.push_back(0);

            effects.fwrite();
            correlated_effects.push_back(effects);

            corr.clear();
            corr.shrink_to_fit();

            var.clear();
            var.shrink_to_fit();
        }
        catch (std::string err)
        {
            std::cerr << "Exception in Model::append_corrstruct(const std::string &, const std::string &, const std::vector<int> &)." << '\n';
            std::cerr << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::append_corrstruct(const std::string &, const std::string &, const std::vector<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::append_corrstruct(const std::string &, const std::string &, const std::vector<int> &)." << '\n';
            throw;
        }

        return 0;
    }

    int Model::append_traitstruct(int obs_id, const std::vector<int> &eff_id)
    {
        try
        {
            matrix<int> eff;

            eff.resize(eff_id.size(), 1);

            for (size_t i = 0; i < eff_id.size(); i++)
                eff[i] = eff_id[i];

            observation_trait.push_back(obs_id);
            effects_trait.push_back(eff);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::append_traitstruct(int, const std::vector<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::append_traitstruct(int, const std::vector<int> &)." << '\n';
            throw;
        }

        return 0;
    }

#endif

    int Model::append_residual(const std::string &fname)
    {
        try
        {
            IOInterface datstream;

            datstream.set_fname(fname);

            std::vector<std::vector<float>> arr;

            datstream.fgetdata(arr);

            size_t lda = arr.size();

            if (lda == 0)
                throw "The first dimension of the arrey is 0, expected et least 1!";

            if (arr[0].size() == 0)
                throw "The second dimension of the arrey is 0, expected et least 1!";

            if (lda != arr[0].size())
                throw "The arrey is not square, that is expected in this case!";

            matrix<float> residual;

            residual.resize(lda, lda);
            // residual.rectosym();

            for (size_t i = 0; i < lda; i++)
            {
                for (size_t j = 0; j < lda; j++)
                    residual(i, j) = arr[i][j];
            }

            residual.fwrite();
            residuals.push_back(residual);

            arr.clear();
            arr.shrink_to_fit();
        }
        catch (std::string err)
        {
            std::cerr << "Exception in Model::append_residual(const std::string &)." << '\n';
            std::cerr << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::append_residual(const std::string &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::append_residual(const std::string &)." << '\n';
            throw;
        }

        return 0;
    }

    int Model::append_observation(const std::string &fname)
    {
        try
        {
            IOInterface datstream;

            datstream.set_fname(fname);

            std::vector<std::vector<float>> arr;

            datstream.fgetdata(arr);

            size_t lda = arr.size();

            if (lda == 0)
                throw "The first dimension of the arrey is 0, expected et least 1!";

            matrix<float> observation;

            // lda := is size of vector

            observation.resize(lda, 1);

            for (size_t i = 0; i < lda; i++)
                observation[i] = arr[i][0];

            observation.fwrite();
            observations.push_back(observation);

            arr.clear();
            arr.shrink_to_fit();
        }
        catch (std::string err)
        {
            std::cerr << "Exception in Model::append_observation(const std::string &)." << '\n';
            std::cerr << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::append_observation(const std::string &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::append_observation(const std::string &)." << '\n';
            throw;
        }

        return 0;
    }

    template <typename T>
    int Model::append_effect(const std::string &fname, T dummy_var)
    {
        try
        {
            IOInterface datstream;

            datstream.set_fname(fname);

            std::vector<std::vector<T>> arr;

            datstream.fgetdata(arr);

            size_t lda1 = arr.size();
            size_t lda2 = arr[0].size();

            if (lda1 == 0)
                throw "The first dimension of the arrey is 0, expected et least 1!";

            if (lda2 == 0)
                throw "The second dimension of the arrey is 0, expected et least 1!";

            matrix<T> effect;
            Effects e;

            // lda1 := number of rows
            // lda2 := number of cols

            effect.resize(lda1, lda2);

            for (size_t i = 0; i < lda1; i++)
            {
                for (size_t j = 0; j < lda2; j++)
                    effect(i, j) = arr[i][j];
            }

            effect.transpose();

            effect.fwrite();

            e.set(effect);

            effects.push_back(e);

            arr.clear();
            arr.shrink_to_fit();
        }
        catch (std::string err)
        {
            std::cerr << "Exception in Model::append_effect(const std::string &, T)." << '\n';
            std::cerr << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::append_effect(const std::string &, T)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::append_effect(const std::string &, T)." << '\n';
            throw;
        }

        return 0;
    }

    // template int Model::append_effect(const std::string &fname, int dummy_var); The 'int' type is specifically reserved for the overloaded method for reading genotypes
    template int Model::append_effect(const std::string &fname, size_t dummy_var);
    template int Model::append_effect(const std::string &fname, bool dummy_var);
    template int Model::append_effect(const std::string &fname, float dummy_var);
    template int Model::append_effect(const std::string &fname, double dummy_var);

    // ------------------------------------------

    int Model::clear_residuals()
    {
        try
        {
            for (auto &e : residuals)
            {
                e.fclear();
                e.clear();
            }
            std::vector<matrix<float>>().swap(residuals);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::clear_residuals()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::clear_residuals()." << '\n';
            throw;
        }

        return 0;
    }

    int Model::clear_observations()
    {
        try
        {
            for (auto &e : observations)
            {
                e.fclear();
                e.clear();
            }
            std::vector<matrix<float>>().swap(observations);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::clear_observations()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::clear_observations()." << '\n';
            throw;
        }

        return 0;
    }

    int Model::clear_effects()
    {
        try
        {
            for (auto &e : effects)
            {
                e.clear();
                // e.fclear();
                // e.clear();
            }
            // std::vector < matrix<int> >().swap( effects );
            std::vector<Effects>().swap(effects);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::clear_effects()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::clear_effects()." << '\n';
            throw;
        }

        return 0;
    }

    int Model::clear_corrstruct()
    {
        try
        {
            for (auto &e : correlated_effects)
            {
                e.fclear();
                e.clear();
            }
            std::vector<matrix<int>>().swap(correlated_effects);

            for (auto &e : variances)
            {
                e.fclear();
                e.clear();
            }
            std::vector<matrix<float>>().swap(variances);

            for (auto &e : correlations)
            {
                e.fclear();
                e.clear();
            }
            std::vector<matrix<float>>().swap(correlations);

            identity_correlations.clear();
            identity_correlations.shrink_to_fit();
            identity_dimension.clear();
            identity_dimension.shrink_to_fit();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::clear_corrstruct()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::clear_corrstruct()." << '\n';
            throw;
        }

        return 0;
    }

    int Model::clear_traitstruct()
    {
        try
        {
            std::vector<int>().swap(observation_trait);
            std::vector<matrix<int>>().swap(effects_trait);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::clear_traitstruct()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::clear_traitstruct()." << '\n';
            throw;
        }

        return 0;
    }

    int Model::clear()
    {
        try
        {
            clear_residuals();
            clear_observations();
            clear_effects();
            clear_corrstruct();
            clear_traitstruct();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::clear()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::clear()." << '\n';
            throw;
        }

        return 0;
    }

#ifdef UTEST

    int Model::print()
    {
        try
        {
            matrix<size_t> shape(2, 1);

            for (auto i2 = 0; i2 < residuals.size(); i2++)
            {
                residuals[i2].fread();
                residuals[i2].print("residual");
                shape = residuals[i2].shape();
                residuals[i2].fwrite();
                shape.transpose();
                shape.print("Shape:");
                shape.clear();
            }
            for (auto i2 = 0; i2 < observations.size(); i2++)
            {
                observations[i2].fread();
                observations[i2].print("observations");
                shape = observations[i2].shape();
                observations[i2].fwrite();
                shape.transpose();
                shape.print("Shape:");
                shape.clear();
            }
            for (auto i2 = 0; i2 < effects.size(); i2++)
            {
                effects[i2].fread();
                effects[i2].print("effects");
                shape = effects[i2].shape();
                effects[i2].fwrite();
                shape.transpose();
                shape.print("Shape:");
                shape.clear();
            }
            for (auto i2 = 0; i2 < correlated_effects.size(); i2++)
            {
                correlated_effects[i2].fread();
                correlated_effects[i2].print("correlated_effects");
                shape = correlated_effects[i2].shape();
                correlated_effects[i2].fwrite();
                shape.transpose();
                shape.print("Shape:");
                shape.clear();
            }
            for (auto i2 = 0; i2 < variances.size(); i2++)
            {
                variances[i2].fread();
                variances[i2].print("variances");
                shape = variances[i2].shape();
                variances[i2].fwrite();
                shape.transpose();
                shape.print("Shape:");
                shape.clear();
            }
            for (auto i2 = 0; i2 < correlations.size(); i2++)
            {
                correlations[i2].fread();
                correlations[i2].print("correlations");
                shape = correlations[i2].shape();
                correlations[i2].fwrite();
                shape.transpose();
                shape.print("Shape:");
                shape.clear();
            }
            for (auto i2 = 0; i2 < effects_trait.size(); i2++)
            {
                effects_trait[i2].print("effects_trait");
                shape = effects_trait[i2].shape();
                shape.transpose();
                shape.print("Shape:");
                shape.clear();
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::print()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::print()." << '\n';
            throw;
        }

        return 0;
    }

    size_t Model::size_of(const std::string type)
    {
        try
        {
            if (type.compare("res") == 0)
            {
                size_t sz = 0;
                if (residuals.size() > 0)
                {
                    for (auto i = 0; i < residuals.size(); i++)
                        sz = sz + residuals[i].size();
                }
                return sz;
            }

            if (type.compare("obs") == 0)
            {
                size_t sz = 0;
                if (observations.size() > 0)
                {
                    for (auto i = 0; i < observations.size(); i++)
                        sz = sz + observations[i].size();
                }
                return sz;
            }

            if (type.compare("eff") == 0)
            {
                size_t sz = 0;
                if (effects.size() > 0)
                {
                    for (auto i = 0; i < effects.size(); i++)
                        sz = sz + effects[i].size();
                }
                return sz;
            }

            if (type.compare("var") == 0)
            {
                size_t sz = 0;
                if (variances.size() > 0)
                {
                    for (auto i = 0; i < variances.size(); i++)
                    {
                        sz = sz + variances[i].size();
                    }
                }
                return sz;
            }

            if (type.compare("cor") == 0)
            {
                size_t sz = 0;
                if (correlations.size() > 0)
                {
                    for (auto i = 0; i < correlations.size(); i++)
                        sz = sz + correlations[i].size();
                }
                return sz;
            }

            if (type.compare("cor_eff") == 0)
            {
                size_t sz = 0;
                if (correlated_effects.size() > 0)
                {
                    for (auto i = 0; i < correlated_effects.size(); i++)
                        sz = sz + correlated_effects[i].size();
                }
                return sz;
            }

            if (type.compare("obs_trt") == 0)
            {
                size_t sz = 0;
                if (observation_trait.size() > 0)
                {
                    for (auto i = 0; i < observation_trait.size(); i++)
                        sz = sz + 1;
                }
                return sz;
            }

            if (type.compare("eff_trt") == 0)
            {
                size_t sz = 0;
                if (effects_trait.size() > 0)
                {
                    for (auto i = 0; i < effects_trait.size(); i++)
                        sz = sz + effects_trait[i].size();
                }
                return sz;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::size_of(const std::string)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::size_of(const std::string)." << '\n';
            throw;
        }
    }

    std::vector<std::vector<size_t>> Model::shape_of(const std::string type)
    {
        try
        {
            std::vector<size_t> sz(2, 0);
            matrix<size_t> shape(2, 1);
            std::vector<std::vector<size_t>> shapes;

            if (type.compare("res") == 0)
            {
                if (residuals.size() > 0)
                {
                    for (auto i = 0; i < residuals.size(); i++)
                    {
                        shape = residuals[i].shape();
                        sz[0] = shape[0];
                        sz[1] = shape[1];
                        shapes.push_back(sz);
                    }
                    shape.clear();
                }
                return shapes;
            }

            if (type.compare("obs") == 0)
            {
                if (observations.size() > 0)
                {
                    for (auto i = 0; i < observations.size(); i++)
                    {
                        shape = observations[i].shape();
                        sz[0] = shape[0];
                        sz[1] = shape[1];
                        shapes.push_back(sz);
                    }
                    shape.clear();
                }
                return shapes;
            }

            if (type.compare("eff") == 0)
            {
                if (effects.size() > 0)
                {
                    for (auto i = 0; i < effects.size(); i++)
                    {
                        shape = effects[i].shape();
                        // sz[0] = shape[0]; // on not transposed
                        // sz[1] = shape[1]; // on not transposed
                        sz[0] = shape[1]; // on transposed
                        sz[1] = shape[0]; // on transposed

                        shapes.push_back(sz);
                    }
                    shape.clear();
                }
                return shapes;
            }

            if (type.compare("var") == 0)
            {
                if (variances.size() > 0)
                {
                    for (auto i = 0; i < variances.size(); i++)
                    {
                        shape = variances[i].shape();
                        sz[0] = shape[0];
                        sz[1] = shape[1];
                        shapes.push_back(sz);
                    }
                    shape.clear();
                }
                return shapes;
            }

            if (type.compare("cor") == 0)
            {
                if (correlations.size() > 0)
                {
                    for (auto i = 0; i < correlations.size(); i++)
                    {
                        shape = correlations[i].shape();
                        sz[0] = shape[0];
                        sz[1] = shape[1];
                        shapes.push_back(sz);
                    }
                    shape.clear();
                }
                return shapes;
            }

            if (type.compare("cor_eff") == 0)
            {
                if (correlated_effects.size() > 0)
                {
                    for (auto i = 0; i < correlated_effects.size(); i++)
                    {
                        shape = correlated_effects[i].shape();
                        sz[0] = shape[0];
                        sz[1] = shape[1];
                        shapes.push_back(sz);
                    }
                    shape.clear();
                }
                return shapes;
            }

            if (type.compare("obs_trt") == 0)
            {
                if (observation_trait.size() > 0)
                {
                    for (auto i = 0; i < observation_trait.size(); i++)
                    {
                        // shape = observation_trait[i].shape();
                        sz[0] = 1;
                        sz[1] = 1;
                        shapes.push_back(sz);
                    }
                    shape.clear();
                }
                return shapes;
            }

            if (type.compare("eff_trt") == 0)
            {
                if (effects_trait.size() > 0)
                {
                    for (auto i = 0; i < effects_trait.size(); i++)
                    {
                        shape = effects_trait[i].shape();
                        sz[0] = shape[0];
                        sz[1] = shape[1];
                        shapes.push_back(sz);
                    }
                    shape.clear();
                }
                return shapes;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::shape_of(const std::string)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::shape_of(const std::string)." << '\n';
            throw;
        }
    }
    template <typename T>
    std::vector<T> Model::test_effects(size_t which_effect, T dummy_type)
    {
        try
        {
            std::vector<T> out;

            matrix<T> effect;

            Effects e = effects[which_effect];

            e.get(effect);

            effect.fread();

            effect.transpose(); // because we transposing when appending

            for (size_t i = 0; i < effect.size(); i++)
                out.push_back(effect[i]);

            return out;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::test_effects(size_t, T)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::test_effects(size_t, T)." << '\n';
            throw;
        }
    }
    template std::vector<int> Model::test_effects(size_t which_effect, int dummy_type);
    template std::vector<size_t> Model::test_effects(size_t which_effect, size_t dummy_type);
    template std::vector<bool> Model::test_effects(size_t which_effect, bool dummy_type);
    template std::vector<float> Model::test_effects(size_t which_effect, float dummy_type);
    template std::vector<double> Model::test_effects(size_t which_effect, double dummy_type);

    std::vector<float> Model::test_observations(size_t which_observations)
    {
        try
        {
            std::vector<float> out;

            matrix<float> obs;

            obs = observations[which_observations];

            obs.fread();

            for (size_t i = 0; i < obs.size(); i++)
                out.push_back(obs[i]);

            return out;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::test_observations(size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::test_observations(size_t)." << '\n';
            throw;
        }
    }

    std::vector<float> Model::test_residual(size_t which_residual)
    {
        try
        {
            std::vector<float> out;

            matrix<float> resid;

            resid = residuals[which_residual];

            resid.fread();

            for (size_t i = 0; i < resid.size(); i++)
                out.push_back(resid[i]);

            return out;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::test_residual(size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::test_residual(size_t)." << '\n';
            throw;
        }
    }

    std::vector<float> Model::test_variance(size_t which_variance)
    {
        try
        {
            std::vector<float> out;

            matrix<float> var;

            var = variances[which_variance];

            var.fread();

            for (size_t i = 0; i < var.size(); i++)
                out.push_back(var[i]);

            return out;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::test_variance(size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::test_variance(size_t)." << '\n';
            throw;
        }
    }

    std::vector<float> Model::test_correlation(size_t which_correlation)
    {
        try
        {
            std::vector<float> out;

            matrix<float> cor;

            cor = correlations[which_correlation];

            cor.fread();

            for (size_t i = 0; i < cor.size(); i++)
                out.push_back(cor[i]);

            return out;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Model::test_correlation(size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Model::test_correlatino(size_t)." << '\n';
            throw;
        }
    }

#endif

} // end of namespace evo