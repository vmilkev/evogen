#include "model.hpp"
#include "solver_pcg.hpp"
#include "iointerface.hpp"
#include <string>

void proces_array(float **arr, size_t lda, size_t ldb);

int main(void)
{
    try
    {
        evo::Pcg solver;
        evo::Model model;

        std::vector<float> iR{0.041}; //, 0.041, 0.041, 0.041};

        std::vector<float> iG1{10.1}; //, 0.1, 0.1, 0.1};

        // size_t dim = 1;

        std::cout << "In main."
                  << "\n";

        model.append_residual(iR, 1);

        model.append_observation("tests/data/model_4/obs_1.dat"); // obs := 0

        size_t type1;
        float type2;

        model.append_effect("tests/data/model_4/obs_489_snp_1000.txt", type1); // eff := 0
        model.append_effect("tests/data/model_4/fixed_1.dat", type2);           // eff := 1

        std::vector<int> corr_eff{0};

        std::vector<int> eff_trate{1, 0};
        int obs_trate = 0;

        std::string identity("I");

        model.append_corrstruct(iG1, 1, identity, 1000, corr_eff);

        model.append_traitstruct(obs_trate, eff_trate);
        // model.append_traitstruct(obs_trate, eff_trate);

        solver.append_model(model);

        solver.solve(2);

        /*

                        std::vector<float> sol = solver.get_solution();

                        solver.get_solution("cpp_solution_model_4.dat");

                        solver.remove_model();

                        model.clear();

                        evo::matrix<float> a(1,3);
                        evo::matrix<float> c(1,3);

                        for(auto i = 0; i < 3; i++)
                            a[i] = i+1;

                        float *b;

                        b = a.return_array();

                        c.insert_array(b);

                        for(auto i = 0; i < 3; i++)
                        {
                            std::cout<<"b = "<<b[i]<<"\n";
                            std::cout<<"c = "<<c[i]<<"\n";
                        }
                */
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
    }
}

void proces_array(float **arr, size_t lda, size_t ldb)
{
    float counter = 1.0;

    for (size_t i = 0; i < lda; i++)
    {
        for (size_t j = 0; j < ldb; j++)
        {
            arr[i][j] = counter;
            counter++;
        }
    }
}