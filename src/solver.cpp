#include "solver.hpp"

namespace evo
{

    void Solver::append_model(const Model &m)
    {
        try
        {
            model = m;
            is_model_added = true;
            process_model();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::append_model(const Model &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::append_model(const Model &)." << '\n';
            throw;
        }
    }

    void Solver::remove_model()
    {
        try
        {
            model.clear();

            std::vector<size_t>().swap(n_obs);

            std::vector<std::vector<size_t>>().swap(n_lev);

            for (auto &e : z_dat)
            {
                e.fclear();
                e.clear();
            }
            std::vector<matrix<float>>().swap(z_dat);

            for (auto &e : z_uni)
            {
                for (auto &p : e)
                {
                    p.clear();
                }
            }
            std::vector<std::vector<Effects>>().swap(z_uni);

            for (auto &e : y)
            {
                e.fclear();
                e.clear();
            }
            std::vector<matrix<float>>().swap(y);

            r.clear();

            z_on_memory = false;
            is_model_added = false;

            adj_effects_order.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::remove_model()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::remove_model()." << '\n';
            throw;
        }
    }

    Solver::~Solver()
    {
        try
        {
            remove_model();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::~Solver()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver destructor." << '\n';
        }
    }

    void Solver::process_model()
    {
        try
        {
            n_trait = model.observation_trait.size();

            size_t i_adj_effects = 0;

            for (size_t i = 0; i < n_trait; i++)
            {
                int trait_obs_id = model.observation_trait[i];

                matrix<int> trait_eff_ids = model.effects_trait[i];

                size_t n_trait_effects = trait_eff_ids.size();

                n_obs.push_back(model.observations[trait_obs_id].size());

                std::vector<size_t> trait_levels;

                for (size_t j = 0; j < n_trait_effects; j++)
                {
                    matrix<size_t> shape = model.effects[trait_eff_ids[j]].shape();

                    // trait_levels.push_back(shape[1]); // on not transposed effects
                    trait_levels.push_back(shape[0]); // on transposed effects

                    adj_effects_order[trait_eff_ids[j]] = i_adj_effects;

                    i_adj_effects++;
                }

                n_lev.push_back(trait_levels);

                construct_union_z(trait_eff_ids); // <- everything on disk !

                set_y(trait_obs_id);

                r = model.residuals[0].fget(); // Note! This is only for the very first submitted residual matrix.
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::process_model()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::process_model()." << '\n';
            throw;
        }
    }

    void Solver::memload_effects()
    {
        try
        {
            if (!is_model_added)
                throw std::string("There is no model added to the solver. Use append_model().");

            if (!z_on_memory)
            {
                if (!z_uni.empty())
                {
                    for (auto &e : z_uni)
                    {
                        for (auto &p : e)
                        {
                            p.fread();
                        }
                    }

                    z_on_memory = true;
                }
                else
                {
                    for (size_t i = 0; i < n_trait; i++)
                    {
                        matrix<int> trait_eff_ids = model.effects_trait[i];
                        construct_union_z(trait_eff_ids, true);
                        z_on_memory = true;
                    }
                }
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Solver::memload_effects(): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::memload_effects()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::memload_effects()." << '\n';
            throw;
        }
    }

    void Solver::diskload_effects()
    {
        try
        {
            if (!is_model_added)
                throw std::string("There is no model added to the solver. Use append_model().");

            if (z_on_memory)
            {
                if (!z_uni.empty())
                {
                    for (auto &e : z_uni)
                    {
                        for (auto &p : e)
                        {
                            p.fwrite();
                        }
                    }

                    z_on_memory = false;
                }
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Solver::diskload_effects(): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::diskload_effects()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::diskload_effects()." << '\n';
            throw;
        }
    }

    void Solver::set_complete_z(const matrix<int> &eff_ids)
    {
        /* Constructs complete incidense for a specific trait.

           NOTE: we are not using this type of matrix yet.
        */

        try
        {
            matrix<float> z;

            z = model.effects[eff_ids[0]].get_float(); // very first effect

            z.transpose();

            for (size_t i = 1; i < eff_ids.size(); i++)
            {
                matrix<float> eff = model.effects[eff_ids[i]].get_float();

                eff.transpose();

                z = z >> eff;

                eff.clear();
            }

            z.fwrite();

            z_dat.push_back(z);

            z.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::set_complete_z(const matrix<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::set_complete_z(const matrix<int> &)." << '\n';
            throw;
        }
    }

    void Solver::construct_union_z(const matrix<int> &eff_ids)
    {
        /* Combines a consecutive sets (from left to right) of incidense matrices
           of different (original & unchanged) types for a specific trait.

           Everything is on disk.
        */

        try
        {
            z_on_memory = false;

            std::vector<Effects> eff;

            eff.push_back(model.effects[eff_ids[0]]); // very first effect

            for (size_t i = 1; i < eff_ids.size(); i++)
            {
                eff.push_back(model.effects[eff_ids[i]]);
            }

            z_uni.push_back(eff);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::construct_union_z(const matrix<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::construct_union_z(const matrix<int> &)." << '\n';
            throw;
        }
    }

    void Solver::construct_union_z(const matrix<int> &eff_ids, bool on_memory)
    {
        /* Combines a consecutive sets (from left to right) of incidense matrices
           of different (original & unchanged) types for a specific trait.

           Everything is on memory.
        */

        try
        {
            if (on_memory)
            {
                z_on_memory = true;

                std::vector<Effects> eff;

                Effects e = model.effects[eff_ids[0]];

                e.fread();

                eff.push_back(e); // very first effect

                e.clear();

                for (size_t i = 1; i < eff_ids.size(); i++)
                {
                    e = model.effects[eff_ids[i]];
                    e.fread();

                    eff.push_back(e);

                    e.clear();
                }

                z_uni.push_back(eff);
            }
            else
            {
                construct_union_z(eff_ids);
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::construct_union_z(const matrix<int> &, bool)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::construct_union_z(const matrix<int> &, bool)." << '\n';
            throw;
        }
    }

    matrix<float> Solver::get_vect_z_uni(const size_t &which_trait, const size_t &which_row)
    {
        /* Assumes z_uni consists of not-transposed matrices.
           Therefore, which_row is supposed to address in combined and transposed
           incidence (effects) matrix.

           CHANGED to transposed!
        */

        matrix<float> vect;

        try
        {
            size_t which_eff_matr = 0; // will be changed inside z_row_to_uni_col(...)

            size_t irow[] = {0, n_obs[0] - 1}; // this is constant for any row, due to the same amount of observations

            size_t icol[] = {0, 0}; // will be changed inside z_row_to_uni_col(...)

            z_row_to_uni_col(which_trait, which_row, icol, which_eff_matr);

            // vect = z_uni[which_trait][which_eff_matr].fget(irow, icol); // here is col vector extracted, assumed not-transposed matrices  from the input

            vect = z_uni[which_trait][which_eff_matr].fget(icol, irow); // here is row vector extracted, assumed transposed matrices from the input

            // vect.transpose(); // transforming to a row vector, because assumed of not-transposed matrices  from the input
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::get_vect_z_uni(const size_t &, const size_t &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::get_vect_z_uni(const size_t &, const size_t &)." << '\n';
            throw;
        }

        return vect;
    }

    void Solver::get_vect_z_uni2(const size_t &which_trait, const size_t &which_row, std::vector<std::vector<float>> &vect)
    {

        try
        {
            size_t which_eff_matr = 0; // will be changed inside z_row_to_uni_col(...)

            size_t irow[] = {0, n_obs[0] - 1}; // this is constant for any row, due to the same amount of observations

            size_t icol[] = {0, 0}; // will be changed inside z_row_to_uni_col(...)

            z_row_to_uni_col(which_trait, which_row, icol, which_eff_matr);

            z_uni[which_trait][which_eff_matr].vect_fget(icol, irow, vect); // here is row vector extracted
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::get_vect_z_uni2(const size_t &, const size_t &, std::vector<std::vector<float>> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::get_vect_z_uni2(const size_t &, const size_t &, std::vector<std::vector<float>> &)." << '\n';
            throw;
        }
    }

    void Solver::get_vect_z_uni2(const size_t &which_trait, const size_t &which_row, float **vect)
    {
        try
        {
            size_t which_eff_matr = 0; // will be changed inside z_row_to_uni_col(...)

            size_t irow[] = {0, n_obs[0] - 1}; // this is constant for any row, due to the same amount of observations

            size_t icol[] = {0, 0}; // will be changed inside z_row_to_uni_col(...)

            z_row_to_uni_col(which_trait, which_row, icol, which_eff_matr);

            z_uni[which_trait][which_eff_matr].vect_fget(icol, irow, vect); // here is row vector extracted
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::get_vect_z_uni2(const size_t &, const size_t &, float **)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::get_vect_z_uni2(const size_t &, const size_t &, float **)." << '\n';
            throw;
        }
    }

    void Solver::z_row_to_uni_col(const size_t &which_trait, const size_t &in_z_row, size_t *out_uni_col, size_t &out_uni_matr)
    {
        /* Here it is assumed the set of incidence matrices in z_uni are not transposed.

           NOTE: changed to transposed from the input. Though, this info is not important in here.
        */

        try
        {
            std::vector<size_t> levels;

            size_t sum_levels = 0;
            levels.push_back(sum_levels);

            for (size_t i = 0; i < n_lev[which_trait].size(); i++)
            {
                sum_levels = sum_levels + n_lev[which_trait][i];
                levels.push_back(sum_levels);

                if (sum_levels - 1 >= in_z_row)
                {
                    out_uni_col[0] = out_uni_col[1] = in_z_row - levels[i];
                    out_uni_matr = i;
                    break;
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::z_row_to_uni_col(const size_t &, const size_t &, size_t *, size_t &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::z_row_to_uni_col(const size_t &, const size_t &, size_t *, size_t &)." << '\n';
            throw;
        }
    }

    void Solver::set_y(const int obs_id)
    {
        try
        {
            matrix<float> yi = model.observations[obs_id];
            y.push_back(yi);
            yi.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::set_y(const int)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::set_y(const int)." << '\n';
            throw;
        }
    }

    size_t Solver::corr_size()
    {
        try
        {
            return model.correlations.size();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::corr_size()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::corr_size()." << '\n';
            throw;
        }
    }

    matrix<int> Solver::get_corr_effects(size_t which_trait)
    {
        try
        {
            matrix<int> which_effects = model.correlated_effects[which_trait];

            which_effects.fread();
            matrix<size_t> shape_eff = which_effects.shape();

            for (size_t i = 0; i < shape_eff[0]; i++)
                which_effects(i, 0) = (int)adj_effects_order[which_effects(i, 0)];

            // return model.correlated_effects[which_trait];
            return which_effects;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::get_corr_effects(size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::get_corr_effects(size_t)." << '\n';
            throw;
        }
    }

    matrix<float> Solver::get_variance(size_t which_trait, size_t row, size_t col)
    {
        try
        {
            if (var_onmem)
            {
                // NOTE, variance was appended directly into memory
                size_t irow[2] = {row, row};
                size_t icol[2] = {col, col};

                matrix<float> var; //(1, 1);

                // var[0] = model.variances[which_trait].at(row, col);
                model.variances[which_trait].cast_fget(irow, icol, var);

                return var;
            }
            else if (!var_onmem)
            {
                // NOTE, variance was appended and stored on disk

                size_t irow[2];
                size_t icol[2];

                irow[0] = row;
                irow[1] = row;

                icol[0] = col;
                icol[1] = col;

                return model.variances[which_trait].fget(irow, icol);
            }
            else
                throw std::string("Problem with determining the state of variance. In Solver::get_variance()");
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Solver::get_variance(size_t, size_t, size_t): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::get_variance(size_t, size_t, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::get_variance(size_t, size_t, size_t)." << '\n';
            throw;
        }
    }

    matrix<float> Solver::get_correlation(size_t which_trait, size_t which_row, size_t col_1, size_t col_2)
    {
        matrix<float> cor_out;

        try
        {
            size_t irow[2] = {which_row, which_row};
            size_t icol[2] = {col_1, col_2};

            if (!cor_onmem)
                return model.correlations[which_trait].fget(irow, icol);
            else if (cor_onmem)
            {
                matrix<float> cor;
                model.correlations[which_trait].cast_fget(irow, icol, cor);
                return cor;
            }
            else
                throw std::string("Problem with determining the state of correlation. In Solver::get_correlation()");
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Solver::get_correlation(size_t, size_t, size_t, size_t): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::get_correlation(size_t, size_t, size_t, size_t)." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::get_correlation(size_t, size_t, size_t, size_t)." << '\n';
            throw;
        }

        return cor_out;
    }

    void Solver::memload_var()
    {
        try
        {
            if (!is_model_added)
                throw std::string("There is no model added to the solver. Use append_model().");

            if (model.variances.size() == 0)
                throw std::string("There are no traits in the model. In Solver::memload_var()");

            if (!var_onmem)
            {
                for (size_t i = 0; i < model.variances.size(); i++)
                    model.variances[i].fread();

                var_onmem = true;
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Solver::memload_var(): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::memload_var()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::memload_var()" << '\n';
            throw;
        }
    }

    void Solver::diskload_var()
    {
        try
        {
            if (!is_model_added)
                throw std::string("There is no model added to the solver. Use append_model().");

            if (model.variances.size() == 0)
                throw std::string("There are no variances appended in the model. In Solver::memload_var()");

            if (var_onmem)
            {
                for (size_t i = 0; i < model.variances.size(); i++)
                    model.variances[i].fwrite();

                var_onmem = false;
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Solver::diskload_var(): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::diskload_var()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::diskload_var()" << '\n';
            throw;
        }
    }

    void Solver::memload_cor()
    {
        try
        {
            if (!is_model_added)
                throw std::string("There is no model added to the solver. Use append_model().");

            if (model.correlations.size() == 0)
                throw std::string("There are no traits in the model. In Solver::memload_cor()");

            if (!cor_onmem)
            {
                for (size_t i = 0; i < model.correlations.size(); i++)
                    model.correlations[i].fread();

                cor_onmem = true;
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Solver::memload_cor(): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::memload_cor()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::memload_cor()" << '\n';
            throw;
        }
    }

    void Solver::diskload_cor()
    {
        try
        {
            if (!is_model_added)
                throw std::string("There is no model added to the solver. Use append_model().");

            if (model.correlations.size() == 0)
                throw std::string("There are no traits in the model. In Solver::memload_cor()");

            if (cor_onmem)
            {
                for (size_t i = 0; i < model.correlations.size(); i++)
                    model.correlations[i].fwrite();

                cor_onmem = false;
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Solver::diskload_cor(): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::diskload_cor()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::diskload_cor()" << '\n';
            throw;
        }
    }

    bool Solver::identity_correlation(size_t which_corr)
    {
        try
        {
            return model.identity_correlations[which_corr];
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::identity_correlation(size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::identity_correlation(size_t)." << '\n';
            throw;
        }
    }

#ifdef UTEST

    void Solver::print()
    {
        try
        {
            std::cout << "n_trait: " << n_trait << std::endl;

            for (size_t i = 0; i < n_obs.size(); i++)
                std::cout << "n_obs: " << n_obs[i] << std::endl;

            for (size_t i = 0; i < n_lev.size(); i++)
                for (size_t j = 0; j < n_lev[i].size(); j++)
                    std::cout << "n_lev: " << n_lev[i][j] << std::endl;

            for (size_t i = 0; i < z_dat.size(); i++)
            {
                z_dat[i].fread();
                z_dat[i].print("z:");
            }

            for (size_t i = 0; i < y.size(); i++)
            {
                y[i].fread();
                y[i].print("y:");
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::print()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::print()." << '\n';
            throw;
        }
    }

    void Solver::test_vect_z_uni(const size_t &which_trait, matrix<float> &out)
    {
        try
        {
            size_t levels = 0;
            for (size_t k = 0; k < n_lev[which_trait].size(); k++)
            {
                levels = levels + n_lev[which_trait][k];
            }

            matrix<float> v(levels, n_obs[which_trait]);

            for (size_t j = 0; j < levels; j++)
            {
                matrix<float> v1 = get_vect_z_uni(which_trait, j);
                for (size_t j2 = 0; j2 < (size_t)n_obs[which_trait]; j2++)
                    v(j, j2) = v1(0, j2);
            }

            out = v;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Solver::test_vect_z_uni(const size_t &, matrix<float> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Solver::test_vect_z_uni(const size_t &, matrix<float> &)." << '\n';
            throw;
        }
    }

#endif

}
