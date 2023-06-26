#ifndef solver_hpp__
#define solver_hpp__

#include <vector>
#include <map>
#include "model.hpp"

namespace evo
{

    class Solver
    {

    public:
        virtual void solve() = 0;

        void append_model(const Model &m);

        void remove_model();

        Solver()
        {
            n_trait = 0;
            z_on_memory = false;
            is_model_added = false;
            var_onmem = false;
            cor_onmem = false;
        };

        virtual ~Solver();

#ifdef UTEST
        void print();
        void test_vect_z_uni(const size_t &which_trait, matrix<float> &out);
#endif

    private:
        void process_model();
        void set_complete_z(const matrix<int> &eff_ids); // build full incedense matrix, z_dat
        void construct_union_z(const matrix<int> &eff_ids);
        void construct_union_z(const matrix<int> &eff_ids, bool on_memory);
        void set_y(const int obs_id);
        void z_row_to_uni_col(const size_t &which_trait,
                              const size_t &in_z_row,
                              size_t *out_uni_col,
                              size_t &out_uni_matr);

    protected:
        Model model;                             // instance of the model to be solved
        std::vector<size_t> n_obs;                  // number of observations for each trait
        size_t n_trait;                             // number of traits
        std::vector<std::vector<size_t>> n_lev;     // num of effects' levels for each trait
        std::vector<matrix<float>> z_dat;        // combined incidence matrix for each trait.
        std::vector<matrix<float>> y;            // Observations for each trait.
        std::vector<std::vector<Effects>> z_uni; // Combines a consecutive sets of incidense matrices.
        matrix<float> r;

        bool z_on_memory;
        bool is_model_added;
        bool var_onmem;
        bool cor_onmem;

        std::map<size_t, size_t> adj_effects_order;

        matrix<float> get_vect_z_uni(const size_t &which_trait, const size_t &which_row);
        void get_vect_z_uni2(const size_t &which_trait, const size_t &which_row, std::vector<std::vector<float>> &vect);
        void get_vect_z_uni2(const size_t &which_trait, const size_t &which_row, float **vect);
        size_t corr_size();
        matrix<int> get_corr_effects(size_t which_trait);
        matrix<float> get_variance(size_t which_trait, size_t row, size_t col);
        matrix<float> get_correlation(size_t which_trait, size_t which_row, size_t col_1, size_t col_2);
        bool identity_correlation(size_t which_corr);

        void memload_effects();
        void diskload_effects();
        void memload_var();
        void diskload_var();
        void memload_cor();
        void diskload_cor();
    };

} // end of namespace evo

#endif // solver_hpp__
