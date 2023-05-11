#ifndef model_hpp__
#define model_hpp__

#include <string>
#include <vector>
#include <typeinfo>

#include "cs_matrix.hpp"
#include "effects.hpp"
#include "iointerface.hpp"

#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#endif

namespace evo
{

        class Model
        {

        public:
                friend class Solver;

#ifdef PYBIND
                int append_residual(pybind11::array_t<float> arr, size_t lda);
                int append_observation(pybind11::array_t<float> arr, size_t lda);

                int append_effect(pybind11::array_t<float> arr, size_t lda1, size_t lda2);

                int append_traitstruct(int obs_id, pybind11::array_t<int> eff_id);

                int append_corrstruct(pybind11::array_t<float> var, size_t lda1, pybind11::array_t<float> corr, size_t lda2, pybind11::array_t<int> which_effects);
                int append_corrstruct(pybind11::array_t<float> var, size_t lda1, const std::string &fname_corr, pybind11::array_t<int> which_effects);
                int append_corrstruct(const std::string &fname_var, const std::string &fname_corr, pybind11::array_t<int> which_effects);

                int append_corrstruct(pybind11::array_t<float> var, size_t lda1, std::string &identity, size_t lda2, pybind11::array_t<int> which_effects);
                int append_corrstruct(const std::string &fname_var, std::string &identity, size_t lda2, pybind11::array_t<int> which_effects);
#else
                int append_residual(const std::vector<float> &arr, size_t lda);
                int append_observation(const std::vector<float> &arr, size_t lda);

                template <typename T>
                int append_effect(const std::vector<T> &arr, size_t lda1, size_t lda2);

                int append_traitstruct(int obs_id, const std::vector<int> &eff_id);

                int append_corrstruct(const std::vector<float> &var, size_t lda1, const std::vector<float> &corr, size_t lda2, const std::vector<int> &which_effects);
                int append_corrstruct(const std::vector<float> &var, size_t lda1, const std::string &fname_corr, const std::vector<int> &which_effects);
                int append_corrstruct(const std::string &fname_var, const std::string &fname_corr, const std::vector<int> &which_effects);

                int append_corrstruct(const std::vector<float> var, size_t lda1, std::string &identity, size_t lda2, const std::vector<int> which_effects);
                int append_corrstruct(const std::string &fname_var, std::string &identity, size_t lda2, const std::vector<int> which_effects);
#endif

                // overloaded methods
                int append_residual(const std::string &fname);
                int append_observation(const std::string &fname);

                template <typename T>
                int append_effect(const std::string &fname, T dummy_var);
                // ------------------

                int clear_residuals();
                int clear_observations();
                int clear_effects();
                int clear_corrstruct();
                int clear_traitstruct();
                int clear(); // clears all

#ifdef UTEST
                size_t size_of(const std::string type);
                std::vector<std::vector<size_t>> shape_of(const std::string type);
                int print();
                template <typename T>
                std::vector<T> test_effects(size_t which_effect, T dummy_type);
                std::vector<float> test_observations(size_t which_observations);
                std::vector<float> test_residual(size_t which_residual);
                std::vector<float> test_variance(size_t which_variance);
                std::vector<float> test_correlation(size_t which_correlation);
#endif

        private:
                // Data
                std::vector<matrix<float>> residuals;    // assumed to full matrix
                std::vector<matrix<float>> observations; // assumed to be (nx1) vectors
                std::vector<Effects> effects;

                // Correlation structure
                std::vector<matrix<int>> correlated_effects; // assumed to be (nx1) vectors. In matlab: which_correlated_effects
                std::vector<matrix<float>> variances;        // assumed to be full matrix. In matlab: G
                std::vector<matrix<float>> correlations;     // assumed to full matrix. In matlab: covariance_matrix
                std::vector<bool> identity_correlations;     // indicates which provided correlations should be treated as identity matrices
                std::vector<size_t> identity_dimension;

                // Trait structure
                std::vector<int> observation_trait;     // assumed to be (nx1) vectors
                std::vector<matrix<int>> effects_trait; // assumed to be (nx1) vectors
        };

} // end of namespace evo

#endif // model_hpp__