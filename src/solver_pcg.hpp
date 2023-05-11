#ifndef solver_pcg_hpp__
#define solver_pcg_hpp__

#include <fstream>
#include <thread>
#include <ctime>

#include "solver.hpp"

#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#endif

namespace evo
{

        class Pcg : public Solver
        {

        public:
                Pcg()
                {
                        tolerance = 1e-6; // default
                        befault_max_iter = true;
                        amatrix_ondisk = false;
                        amatrix_onmem = false;

                        pipeline_val = 4; // default

                        srand(rdtsc());
                        int iNum = std::rand() % 100000;
                        
                        binFilename = "amatrix_" + std::to_string(iNum);
                }

                ~Pcg()
                {
                        if (fA.is_open())
                                fA.close();
                }

                void solve();
                void solve( int pipeline_index );

                void set_tolerance(double tol);
                void set_maxiter(size_t iter);

#ifdef PYBIND
                pybind11::array_t<float> get_solution();
#else
                std::vector<float> get_solution();
#endif

                int get_solution(const std::string &fname);

#ifdef UTEST
                matrix<float> test_z_dot_y(size_t vect_size, size_t i_trait, size_t j_trait, float r);
                matrix<float> test_rhs();

                size_t test_num_all_levels();
                std::vector<int> test_ordered_levels();
                std::vector<std::vector<size_t>> test_cov_offsets();
                std::vector<float> test_dval();
                std::vector<std::vector<float>> test_A(); // coefficient matrix
                std::vector<std::vector<float>> construct_A(std::vector<std::vector<size_t>> &cov_offsets, size_t num_levels, std::vector<int> &ordered_levels);
#endif

        private:
                void construct_rhs();
                matrix<float> z_dot_y(size_t vect_size, size_t i_trait, size_t j_trait, float r);
                matrix<float> z_dot_z(size_t row, size_t vect_size, size_t i_matr, size_t j_matr, float r);
                float v_dot_v(const matrix<float> &v1, const matrix<float> &v2);

                matrix<float> get_row_cmatr(size_t rhs_size,
                                            size_t i_trate,
                                            size_t i_eff,
                                            std::vector<std::vector<size_t>> &cov_offsets,
                                            size_t num_levels,
                                            std::vector<int> &ordered_levels,
                                            size_t i_row);

                size_t get_levels(size_t which_trait);
                size_t get_all_levels();
                size_t get_all_levels(size_t before_trait);
                size_t num_all_levels();
                std::vector<int> get_ordered_levels();
                std::vector<std::vector<size_t>> get_cov_offsets(const std::vector<int> &ordered_levels);

                void jacobi_pcg(std::vector<std::vector<size_t>> &cov_offsets,
                                size_t num_levels,
                                std::vector<int> &ordered_levels);

                matrix<float> construct_dval(std::vector<std::vector<size_t>> &cov_offsets,
                                             size_t num_levels,
                                             std::vector<int> &ordered_levels);

                void set_amatr(std::vector<std::vector<size_t>> &cov_offsets,
                               size_t num_levels,
                               std::vector<int> &ordered_levels);

                void set_amatr(std::vector<std::vector<size_t>> &cov_offsets,
                               size_t num_levels,
                               std::vector<int> &ordered_levels, bool on_mem);

                void update_vect(std::vector<std::vector<size_t>> &cov_offsets,
                                 size_t num_levels,
                                 std::vector<int> &ordered_levels,
                                 matrix<float> &out_vect,
                                 matrix<float> &in_vect);

                matrix<float> fget_vect(size_t all_tr_levels, size_t row);

                void set_amatrix();
                void memload_amatr();
                void set_pipeline(int which_pipeline);

                unsigned long long rdtsc(); /* Seed for random number generator. */

                matrix<float> rhs;
                matrix<float> sol;
                matrix<float> dval; // vector of inverse of diagonal elements of coefficient maatrix
                matrix<float> amatr;
                double tolerance;
                size_t iterations;
                size_t max_iterations;
                bool befault_max_iter;
                bool amatrix_ondisk;
                bool amatrix_onmem;

                int pipeline_val;

                std::string binFilename; /* Name of binary file to store A on disck. */
                std::fstream fA;
        };

} // end of namespace evo

#endif // solver_pcg_hpp__