/* test_matrix_mytype.cpp

 Tests of the matrix class from a static library.
 Class instance of type 'mytype'.
*/

#include "catch.hpp"
#include "model.hpp"
#include "solver_pcg.hpp"
#include "iointerface.hpp"
#include <string>

TEST_CASE("Small data test: model 1")
{

    // ---------------------------
    // model:
    // effect: 1       2
    // y1 = b1*X1 + a1*Z1 + e1;
    // effect: 3       4
    // y2 = b2*X2 + a2*Z2 + e2;
    // ---------------------------

    // DATA

    std::vector<float> iR{0.0278, -0.0102,
                          -0.0102, 0.0371}; // full matrix

    std::vector<float> y1{4.5, 2.9, 3.9, 3.5, 5.0};

    std::vector<float> y2{6.8, 5.0, 6.8, 6.0, 7.5};

    std::vector<int> x1{1, 0,
                        0, 1,
                        0, 1,
                        1, 0,
                        1, 0};

    std::vector<int> x2{1, 0,
                        0, 1,
                        0, 1,
                        1, 0,
                        1, 0};

    std::vector<int> z1{0, 0, 0, 1, 0, 0, 0, 0,
                        0, 0, 0, 0, 1, 0, 0, 0,
                        0, 0, 0, 0, 0, 1, 0, 0,
                        0, 0, 0, 0, 0, 0, 1, 0,
                        0, 0, 0, 0, 0, 0, 0, 1};

    std::vector<int> z2{0, 0, 0, 1, 0, 0, 0, 0,
                        0, 0, 0, 0, 1, 0, 0, 0,
                        0, 0, 0, 0, 0, 1, 0, 0,
                        0, 0, 0, 0, 0, 0, 1, 0,
                        0, 0, 0, 0, 0, 0, 0, 1};

    std::vector<float> iA{1.833, 0.5, 0, -0.667, 0, -1, 0, 0,
                          0.5, 2, 0.5, 0, -1, -1, 0, 0,
                          0, 0.5, 2, 0, -1, 0.5, 0, -1,
                          -0.667, 0, 0, 1.833, 0.5, 0, -1, 0,
                          0, -1, -1, 0.5, 2.5, 0, -1, 0,
                          -1, -1, 0.5, 0, 0, 2.5, 0, -1,
                          0, 0, 0, -1, -1, 0, 2, 0,
                          0, 0, -1, 0, 0, -1, 0, 2}; // full matrix

    std::vector<float> iG1{0.084, -0.0378,
                           -0.0378, 0.042}; // full matrix

    std::vector<std::vector<float>> true_z{
        {1, 0, 0, 1, 1},
        {0, 1, 1, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0},
        {0, 1, 0, 0, 0},
        {0, 0, 1, 0, 0},
        {0, 0, 0, 1, 0},
        {0, 0, 0, 0, 1}};

    std::vector<float> _rhs{
        0.1543399,
        0.068767377,
        0,
        0,
        0,
        0.0557924,
        0.0296570,
        0.03911028,
        0.036144578,
        0.062557,
        0.620018,
        0.36811862,
        0,
        0,
        0,
        0.20620945,
        0.1557924,
        0.212326227,
        0.18674698,
        0.22706209};

    std::vector<float> _sol{
        4.360827,
        3.397239,
        0.150977,
        -0.015388,
        -0.078379,
        -0.010185,
        -0.270317,
        0.275839,
        -0.316081,
        0.243784,
        6.799822,
        5.880252,
        0.279713,
        -0.007601,
        -0.170316,
        -0.01257,
        -0.477801,
        0.517297,
        -0.478915,
        0.392017};

    std::vector<float> dval_true{
        0.08341,
        0.0556,
        0.154,
        0.168,
        0.168,
        0.1818,
        0.2378,
        0.2378,
        0.1958,
        0.1958,
        0.1112,
        0.0741,
        0.077,
        0.084,
        0.084,
        0.114,
        0.1421,
        0.1421,
        0.1211,
        0.1211};

    size_t n_all_levels = 4;

    std::vector<int> ordered_random_levels{2, 8, 2, 8};

    std::vector<std::vector<size_t>> rcov_offsets{{0, 0}, {0, 2}, {0, 10}, {0, 12}, {2, 0}, {2, 2}, {2, 10}, {2, 12}, {10, 0}, {10, 2}, {10, 10}, {10, 12}, {12, 0}, {12, 2}, {12, 10}, {12, 12}};

    std::vector<std::vector<float>> CoeffMatrix{{0.08341, 0, 0, 0, 0, 0.0278, 0, 0, 0.0278, 0.0278, -0.030583, 0, 0, 0, 0, -0.0101, 0, 0, -0.0101, -0.0101},
                                                {0, 0.0556, 0, 0, 0, 0, 0.0278, 0.0278, 0, 0, 0, -0.0203, 0, 0, 0, 0, -0.0101, -0.0101, 0, 0},
                                                {0, 0, 0.154, 0.042, 0, -0.056, 0, -0.084, 0, 0, 0, 0, -0.06931, -0.0189, 0, 0.02522, 0, 0.03781, 0, 0},
                                                {0, 0, 0.042, 0.168, 0.042, 0, -0.084, -0.084, 0, 0, 0, 0, -0.0189, -0.07563, -0.0189, 0, 0.03781, 0.03781, 0, 0},
                                                {0, 0, 0, 0.042, 0.168, 0, -0.084, 0.042, 0, -0.084, 0, 0, 0, -0.0189, -0.07563, 0, 0.03781, -0.0189, 0, 0.0378},
                                                {0.0278, 0, -0.056, 0, 0, 0.18183, 0.042, 0, -0.084, 0, -0.0101, 0, 0.02522, 0, 0, -0.0795, -0.0189, 0, 0.03781, 0},
                                                {0, 0.0278, 0, -0.084, -0.084, 0.042, 0.237887, 0, -0.084, 0, 0, -0.0101, 0, 0.03781, 0.03781, -0.0189, -0.10473, 0, 0.03781, 0},
                                                {0, 0.0278, -0.084, -0.084, 0.042, 0, 0, 0.2378, 0, -0.084, 0, -0.0101, 0.0378, 0.0378, -0.0189, 0, 0, -0.10473, 0, 0.0378},
                                                {0.0278, 0, 0, 0, 0, -0.084, -0.084, 0, 0.19587, 0, -0.0101, 0, 0, 0, 0, 0.0378, 0.0378, 0, -0.0858, 0},
                                                {0.0278, 0, 0, 0, -0.08403, 0, 0, -0.084, 0, 0.19587, -0.0101, 0, 0, 0, 0.0378, 0, 0, 0.0378, 0, -0.0858},
                                                {-0.030583, 0, 0, 0, 0, -0.0101, 0, 0, -0.0101, -0.0101, 0.1112, 0, 0, 0, 0, 0.03707, 0, 0, 0.03707, 0.03707},
                                                {0, -0.020389, 0, 0, 0, 0, -0.0101, -0.0101, 0, 0, 0, 0.07414, 0, 0, 0, 0, 0.037, 0.037, 0, 0},
                                                {0, 0, -0.0693, -0.0189, 0, 0.0252, 0, 0.0378, 0, 0, 0, 0, 0.077, 0.021, 0, -0.028, 0, -0.04201, 0, 0},
                                                {0, 0, -0.0189, -0.07563, -0.0189, 0, 0.0378, 0.0378, 0, 0, 0, 0, 0.021, 0.084, 0.021, 0, -0.04201, -0.042, 0, 0},
                                                {0, 0, 0, -0.0189, -0.07563, 0, 0.0378, -0.0189, 0, 0.0378, 0, 0, 0, 0.021, 0.084, 0, -0.042, 0.021, 0, -0.042},
                                                {-0.0101, 0, 0.0252, 0, 0, -0.0795, -0.0189, 0, 0.0378, 0, 0.037, 0, -0.028, 0, 0, 0.114, 0.021, 0, -0.042, 0},
                                                {0, -0.0101, 0, 0.0378, 0.0378, -0.0189, -0.10473, 0, 0.0378, 0, 0, 0.037, 0, -0.042, -0.042, 0.021, 0.1421, 0, -0.042, 0},
                                                {0, -0.0101, 0.0378, 0.0378, -0.0189, 0, 0, -0.10473, 0, 0.0378, 0, 0.037, -0.042, -0.042, 0.021, 0, 0, 0.1421, 0, -0.042},
                                                {-0.0101, 0, 0, 0, 0, 0.0378, 0.0378, 0, -0.0858, 0, 0.037, 0, 0, 0, 0, -0.042, -0.042, 0, 0.121, 0},
                                                {-0.0101, 0, 0, 0, 0.0378, 0, 0, 0.0378, 0, -0.0858, 0.037, 0, 0, 0, -0.042, 0, 0, -0.042, 0, 0.121}};
    // DATA for testing IO interface:

    std::vector<std::vector<int>> allele{
        {2, 0, 1, 1, 0, 0, 0, 2, 1, 2},
        {5, 0, 0, 0, 0, 2, 0, 2, 1, 0},
        {1, 5, 2, 1, 1, 0, 0, 2, 1, 2},
        {0, 0, 2, 1, 0, 1, 0, 2, 2, 1},
        {0, 1, 1, 2, 0, 0, 0, 2, 1, 2},
        {1, 1, 5, 1, 0, 2, 0, 2, 2, 1}};

    std::vector<std::vector<float>> matr{
        {12.2, 20, 51.1},
        {15.5, 30, 10},
        {21.0, 45, 562},
        {30.5, 50, 452},
        {40, 61, 231},
        {51.3, 71, 125},
        {60.6, 80, 121},
        {70.001, 91, 121},
        {82.012, 10, 110.0}};

    std::vector<std::vector<int>> allelebin{
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};

    std::vector<std::vector<int>> alleledat{
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0},
        {0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0},
        {0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0},
        {0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0}};

    // ========================================================================================

    SECTION("The Model class: test size")
    {
        try
        {
            evo::Model model;

            CHECK(model.size_of("res") == 0);
            CHECK(model.size_of("obs") == 0);
            CHECK(model.size_of("eff") == 0);
            CHECK(model.size_of("var") == 0);
            CHECK(model.size_of("cor") == 0);
            CHECK(model.size_of("cor_eff") == 0);

            std::vector<int> corr_eff{2, 3};

            std::vector<int> eff_trate_1{0, 2};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{1, 3};
            int obs_trate_2 = 1;

            model.append_residual(iR, 2);
            model.append_observation(y1, 5); // obs := 0

            CHECK(model.size_of("res") == 4);
            CHECK(model.size_of("obs") == 5);

            model.append_residual(iR, 2);
            model.append_observation(y2, 5); // obs := 1

            CHECK(model.size_of("res") == 8);
            CHECK(model.size_of("obs") == 10);

            model.append_effect(x1, 5, 2); // eff:= 0
            CHECK(model.size_of("eff") == 10);

            model.append_effect(x2, 5, 2); // eff:= 1
            CHECK(model.size_of("eff") == 20);

            model.append_effect(z1, 5, 8); // eff:= 2
            CHECK(model.size_of("eff") == 60);

            model.append_effect(z2, 5, 8); // eff:= 3
            CHECK(model.size_of("eff") == 100);

            model.append_corrstruct(iG1, 2, iA, 8, corr_eff);
            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            CHECK(model.size_of("var") == 4);
            CHECK(model.size_of("cor") == 64);
            CHECK(model.size_of("cor_eff") == 2);
            CHECK(model.size_of("obs_trt") == 2);
            CHECK(model.size_of("eff_trt") == 4);

            model.clear_residuals();
            model.clear_observations();
            model.clear_effects();
            model.clear_corrstruct();
            model.clear_traitstruct();

            CHECK(model.size_of("res") == 0);
            CHECK(model.size_of("obs") == 0);
            CHECK(model.size_of("eff") == 0);
            CHECK(model.size_of("var") == 0);
            CHECK(model.size_of("cor") == 0);
            CHECK(model.size_of("cor_eff") == 0);
            CHECK(model.size_of("obs_trt") == 0);
            CHECK(model.size_of("eff_trt") == 0);
        }
        catch (const std::exception &e)
        {
            std::cerr << " Model 1. Exit from the Model tests size due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << e.what() << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (const std::string &err)
        {
            std::cerr << " Model 1. Exit from the Model tests size due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << err << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (int ierr)
        {
            std::cerr << " Model 1. Exit from the Model tests size due to the exception with code " << ierr << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (...)
        {
            std::cerr << " Model 1. Exit from the Model tests size due to an unknown exception." << '\n';
            //exit(EXIT_FAILURE);
        }
    }

    // ========================================================================================

    SECTION("The Model class: test shapes")
    {
        try
        {
            evo::Model model;

            std::vector<std::vector<size_t>> shapes;

            std::vector<int> corr_eff{2, 3};

            std::vector<int> eff_trate_1{0, 2};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{1, 3};
            int obs_trate_2 = 1;

            model.append_residual(iR, 2);
            model.append_observation(y1, 5); // obs := 0

            model.append_residual(iR, 2);
            model.append_observation(y2, 5); // obs := 1

            model.append_effect(x1, 5, 2); // eff:= 0

            model.append_effect(x2, 5, 2); // eff:= 1

            model.append_effect(z1, 5, 8); // eff:= 2

            model.append_effect(z2, 5, 8); // eff:= 3

            model.append_corrstruct(iG1, 2, iA, 8, corr_eff);
            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            shapes = model.shape_of("eff");

            CHECK(shapes[0].at(0) == 5);
            CHECK(shapes[0].at(1) == 2);
            CHECK(shapes[1].at(0) == 5);
            CHECK(shapes[1].at(1) == 2);
            CHECK(shapes[2].at(0) == 5);
            CHECK(shapes[2].at(1) == 8);
            CHECK(shapes[3].at(0) == 5);
            CHECK(shapes[3].at(1) == 8);

            shapes.clear();

            shapes = model.shape_of("obs");
            for (auto &e : shapes)
            {
                CHECK(e[0] == 5);
                CHECK(e[1] == 1);
            }
            shapes.clear();

            shapes = model.shape_of("res");
            for (auto &e : shapes)
            {
                CHECK(e[0] == 2);
                CHECK(e[1] == 2);
            }
            shapes.clear();

            shapes = model.shape_of("var");
            for (auto &e : shapes)
            {
                CHECK(e[0] == 2);
                CHECK(e[1] == 2);
            }
            shapes.clear();

            shapes = model.shape_of("cor");
            for (auto &e : shapes)
            {
                CHECK(e[0] == 8);
                CHECK(e[1] == 8);
            }
            shapes.clear();

            shapes = model.shape_of("cor_eff");
            for (auto &e : shapes)
            {
                CHECK(e[0] == 2);
                CHECK(e[1] == 1);
            }
            shapes.clear();

            shapes = model.shape_of("eff_trt");
            for (auto &e : shapes)
            {
                CHECK(e[0] == 2);
                CHECK(e[1] == 1);
            }
            shapes.clear();

            shapes = model.shape_of("obs_trt");
            for (auto &e : shapes)
            {
                CHECK(e[0] == 1);
                CHECK(e[1] == 1);
            }
            shapes.clear();

            model.clear_residuals();
            model.clear_observations();
            model.clear_effects();
            model.clear_corrstruct();
            model.clear_traitstruct();
        }
        catch (const std::exception &e)
        {
            std::cerr << " Model 1. Exit from the Model tests shapes due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << e.what() << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (const std::string &err)
        {
            std::cerr << " Model 1. Exit from the Model tests shapes due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << err << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (int ierr)
        {
            std::cerr << " Model 1. Exit from the Model tests shapes due to the exception with code " << ierr << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (...)
        {
            std::cerr << " Model 1. Exit from the Model tests shapes due to an unknown exception." << '\n';
            //exit(EXIT_FAILURE);
        }
    }

    // ========================================================================================

    SECTION("The Pcg class: test vect_z_uni")
    {
        try
        {
            evo::Pcg solver;
            evo::Model model;

            std::vector<int> eff_trate_1{0, 1};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{2, 3};
            int obs_trate_2 = 1;

            // define the model
            model.append_residual(iR, 2);

            model.append_observation(y1, 5); // obs := 0
            model.append_observation(y2, 5); // obs := 1

            model.append_effect(x1, 5, 2); // eff := 0
            model.append_effect(z1, 5, 8); // eff := 1
            model.append_effect(x2, 5, 2); // eff := 2
            model.append_effect(z2, 5, 8); // eff := 3

            std::vector<int> corr_eff{1, 3};

            model.append_corrstruct(iG1, 2, iA, 8, corr_eff);
            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            solver.append_model(model);

            evo::matrix<float> z;

            solver.test_vect_z_uni(0, z);

            for (auto i = 0; i < true_z.size(); i++)
            {
                for (auto j = 0; j < true_z[i].size(); j++)
                {
                    CHECK(z(i, j) == true_z[i][j]);
                }
            }

            z.clear();

            solver.test_vect_z_uni(1, z);

            for (auto i = 0; i < true_z.size(); i++)
            {
                for (auto j = 0; j < true_z[i].size(); j++)
                {
                    CHECK(z(i, j) == true_z[i][j]);
                }
            }

            z.clear();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test vect_z_uni] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << e.what() << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (const std::string &err)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test vect_z_uni] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << err << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (int ierr)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test vect_z_uni] due to the exception with code " << ierr << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (...)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test vect_z_uni] due to an unknown exception." << '\n';
            //exit(EXIT_FAILURE);
        }
    }

    // ========================================================================================

    SECTION("The Pcg class: test RHS")
    {
        try
        {
            evo::Pcg solver;
            evo::Model model;

            std::vector<int> eff_trate_1{0, 1};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{2, 3};
            int obs_trate_2 = 1;

            // define the model
            model.append_residual(iR, 2);

            model.append_observation(y1, 5); // obs := 0
            model.append_observation(y2, 5); // obs := 1

            model.append_effect(x1, 5, 2); // eff := 0
            model.append_effect(z1, 5, 8); // eff := 1
            model.append_effect(x2, 5, 2); // eff := 2
            model.append_effect(z2, 5, 8); // eff := 3

            std::vector<int> corr_eff{1, 3};

            model.append_corrstruct(iG1, 2, iA, 8, corr_eff);
            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            solver.append_model(model);

            evo::matrix<float> model_rhs = solver.test_rhs();

            for (auto i = 0; i < _rhs.size(); i++)
                CHECK((_rhs[i]) == Approx(model_rhs[i]).margin(0.0001).epsilon(1e-3));

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test RHS] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << e.what() << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (const std::string &err)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test RHS] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << err << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (int ierr)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test RHS] due to the exception with code " << ierr << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (...)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test RHS] due to an unknown exception." << '\n';
            //exit(EXIT_FAILURE);
        }
    }

    // ========================================================================================

    SECTION("The Pcg class: test dval")
    {
        try
        {
            evo::Pcg solver;
            evo::Model model;

            std::vector<int> eff_trate_1{0, 1};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{2, 3};
            int obs_trate_2 = 1;

            // define the model
            model.append_residual(iR, 2);

            model.append_observation(y1, 5); // obs := 0
            model.append_observation(y2, 5); // obs := 1

            model.append_effect(x1, 5, 2); // eff := 0
            model.append_effect(z1, 5, 8); // eff := 1
            model.append_effect(x2, 5, 2); // eff := 2
            model.append_effect(z2, 5, 8); // eff := 3

            std::vector<int> corr_eff{1, 3};

            model.append_corrstruct(iG1, 2, iA, 8, corr_eff);
            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            solver.append_model(model);

            std::vector<float> dval = solver.test_dval();

            CHECK(dval.size() == dval_true.size());

            for (size_t i = 0; i < dval_true.size(); i++)
                CHECK(dval[i] == Approx(dval_true[i]).margin(0.0001).epsilon(1e-4));

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test dval] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << e.what() << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (const std::string &err)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test dval] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << err << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (int ierr)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test dval] due to the exception with code " << ierr << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (...)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test dval] due to an unknown exception." << '\n';
            //exit(EXIT_FAILURE);
        }
    }

    // ========================================================================================

    SECTION("The Pcg class: test CoeffMatrix")
    {
        try
        {
            evo::Pcg solver;
            evo::Model model;

            std::vector<int> eff_trate_1{0, 1};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{2, 3};
            int obs_trate_2 = 1;

            // define the model
            model.append_residual(iR, 2);

            model.append_observation(y1, 5); // obs := 0
            model.append_observation(y2, 5); // obs := 1

            model.append_effect(x1, 5, 2); // eff := 0
            model.append_effect(z1, 5, 8); // eff := 1
            model.append_effect(x2, 5, 2); // eff := 2
            model.append_effect(z2, 5, 8); // eff := 3

            std::vector<int> corr_eff{1, 3};

            model.append_corrstruct(iG1, 2, iA, 8, corr_eff);
            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            solver.append_model(model);

            std::vector<std::vector<float>> A = solver.test_A();

            CHECK(A.size() == CoeffMatrix.size());
            CHECK(A[0].size() == CoeffMatrix[0].size());

            for (size_t i = 0; i < A.size(); i++)
            {
                for (size_t j = 0; j < A[0].size(); j++)
                {
                    CHECK(A[i][j] == Approx(CoeffMatrix[i][j]).margin(0.001).epsilon(1e-3));
                }
            }

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test CoeffMatrix] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << e.what() << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (const std::string &err)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test CoeffMatrix] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << err << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (int ierr)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test CoeffMatrix] due to the exception with code " << ierr << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (...)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test CoeffMatrix] due to an unknown exception." << '\n';
            //exit(EXIT_FAILURE);
        }
    }

    // ========================================================================================

    SECTION("The Pcg class: test intermediate data structures")
    {
        try
        {
            evo::Pcg solver;
            evo::Model model;

            std::vector<int> eff_trate_1{0, 1};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{2, 3};
            int obs_trate_2 = 1;

            // define the model
            model.append_residual(iR, 2);

            model.append_observation(y1, 5); // obs := 0
            model.append_observation(y2, 5); // obs := 1

            model.append_effect(x1, 5, 2); // eff := 0
            model.append_effect(z1, 5, 8); // eff := 1
            model.append_effect(x2, 5, 2); // eff := 2
            model.append_effect(z2, 5, 8); // eff := 3

            std::vector<int> corr_eff{1, 3};

            model.append_corrstruct(iG1, 2, iA, 8, corr_eff);
            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            solver.append_model(model);

            size_t n_all_levels_checking = solver.test_num_all_levels();

            CHECK(n_all_levels_checking == n_all_levels);

            std::vector<int> ordered_random_levels_checking = solver.test_ordered_levels();

            CHECK(ordered_random_levels_checking.size() == ordered_random_levels.size());

            for (size_t i = 0; i < ordered_random_levels_checking.size(); i++)
                CHECK(ordered_random_levels_checking[i] == ordered_random_levels[i]);

            std::vector<std::vector<size_t>> rcov_offsets_checking = solver.test_cov_offsets();

            CHECK(rcov_offsets_checking.size() == rcov_offsets.size());
            CHECK(rcov_offsets_checking[0].size() == rcov_offsets[0].size());

            for (size_t i = 0; i < rcov_offsets.size(); i++)
            {
                for (size_t j = 0; j < rcov_offsets[0].size(); j++)
                {
                    CHECK(rcov_offsets_checking[i][j] == rcov_offsets[i][j]);
                }
            }

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test intermediate data structures] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << e.what() << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (const std::string &err)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test intermediate data structures] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << err << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (int ierr)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test intermediate data structures] due to the exception with code " << ierr << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (...)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test intermediate data structures] due to an unknown exception." << '\n';
            //exit(EXIT_FAILURE);
        }
    }

    // ========================================================================================

    SECTION("The Pcg class: test solution")
    {
        try
        {
            evo::Pcg solver;
            evo::Model model;

            std::vector<int> eff_trate_1{0, 1};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{2, 3};
            int obs_trate_2 = 1;

            // define the model
            model.append_residual(iR, 2);

            model.append_observation(y1, 5); // obs := 0
            model.append_observation(y2, 5); // obs := 1

            model.append_effect(x1, 5, 2); // eff := 0
            model.append_effect(z1, 5, 8); // eff := 1
            model.append_effect(x2, 5, 2); // eff := 2
            model.append_effect(z2, 5, 8); // eff := 3

            std::vector<int> corr_eff{1, 3};

            model.append_corrstruct(iG1, 2, iA, 8, corr_eff);
            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            solver.append_model(model);

            solver.solve();

            std::vector<float> sol = solver.get_solution();

            solver.get_solution("cpp_solution_model_1.dat");

            for (auto i = 0; i < sol.size(); i++)
                CHECK((_sol[i]) == Approx(sol[i]).margin(0.0001).epsilon(1e-3));

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test solution] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << e.what() << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (const std::string &err)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test solution] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << err << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (int ierr)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test solution] due to the exception with code " << ierr << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (...)
        {
            std::cerr << " Model 1. Exit from [The Pcg class: test solution] due to an unknown exception." << '\n';
            //exit(EXIT_FAILURE);
        }
    }

    // ========================================================================================

    SECTION("The IO interface class: test base class")
    {
        try
        {
            evo::IOInterface datstream;

            datstream.set_fname("tests/data/allele.dat");

            std::vector<std::vector<int>> in;

            datstream.fgetdata(true, in);

            REQUIRE(in.empty() == false);
            CHECK(in.size() == 6);
            CHECK(in[0].size() == 10);

            for (auto i = 0; i < in.size(); i++)
            {
                for (auto j = 0; j < in[0].size(); j++)
                {
                    CHECK(in[i][j] == allele[i][j]);
                }
            }

            in.clear();
            in.shrink_to_fit();

            datstream.set_fname("tests/data/allele2.dat");

            datstream.fgetdata(in);

            REQUIRE(in.empty() == false);
            CHECK(in.size() == 6);
            CHECK(in[0].size() == 10);

            for (auto i = 0; i < in.size(); i++)
            {
                for (auto j = 0; j < in[0].size(); j++)
                {
                    CHECK(in[i][j] == allele[i][j]);
                }
            }

            in.clear();
            in.shrink_to_fit();

            //-------------------------------------------------

            datstream.set_fname("tests/data/data_matr.dat");

            std::vector<std::vector<float>> in2;

            datstream.fgetdata(in2);

            REQUIRE(in2.empty() == false);
            CHECK(in2.size() == 9);
            CHECK(in2[0].size() == 3);

            for (auto i = 0; i < in2.size(); i++)
            {
                for (auto j = 0; j < in2[0].size(); j++)
                {
                    CHECK(in2[i][j] == matr[i][j]);
                }
            }

            in2.clear();
            in2.shrink_to_fit();

            //-------------------------------------------------

            datstream.set_fname("tests/data/1000G.EUR.QC.22.bed"); // Expected rows: variants (SNPs); columns: samples (observations).

            datstream.fgetdata(489, 141123, in);

            REQUIRE(in.empty() == false);
            CHECK(in[0].size() == 489);
            CHECK(in.size() == 141123);

            for (auto i = 0; i < 10; i++)
            {
                for (auto j = 0; j < 20; j++)
                {
                    CHECK(in[i][j] == allelebin[i][j]);
                }
            }

            in.clear();
            in.shrink_to_fit();

            //-------------------------------------------------

            datstream.set_fname("tests/data/model_4/obs_489_snp_141123.txt"); // Expected columns: variants (SNPs); rows: samples (observations).

            datstream.fgetdata(in);

            REQUIRE(in.empty() == false);
            CHECK(in.size() == 489);
            CHECK(in[0].size() == 141123);

            for (auto i = 0; i < 10; i++)
            {
                for (auto j = 0; j < 20; j++)
                {
                    CHECK(in[i][j] == alleledat[i][j]);
                }
            }

            in.clear();
            in.shrink_to_fit();
        }
        catch (const std::exception &e)
        {
            std::cerr << " Model 1. Exit from [The IO interface class: test base class] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << e.what() << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (const std::string &err)
        {
            std::cerr << " Model 1. Exit from [The IO interface class: test base class] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << err << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (int ierr)
        {
            std::cerr << " Model 1. Exit from [The IO interface class: test base class] due to the exception with code " << ierr << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (...)
        {
            std::cerr << " Model 1. Exit from [The IO interface class: test base class] due to an unknown exception." << '\n';
            //exit(EXIT_FAILURE);
        }
    }

    // ========================================================================================

    SECTION("The overloaded methods (Model & IO interface): testing IN -> data")
    {
        try
        {
            evo::Model model;

            size_t dummy_type1;
            float dummy_type2;

            model.append_effect("tests/data/z1.dat", dummy_type1); // eff_0
            model.append_effect("tests/data/z2.dat", dummy_type1); // eff_1
            model.append_effect("tests/data/x1.dat", dummy_type2); // eff_2
            model.append_effect("tests/data/x2.dat", dummy_type2); // eff_3

            std::vector<size_t> eff_0 = model.test_effects(0, dummy_type1);
            std::vector<size_t> eff_1 = model.test_effects(1, dummy_type1);
            std::vector<float> eff_2 = model.test_effects(2, dummy_type2);
            std::vector<float> eff_3 = model.test_effects(3, dummy_type2);

            //std::cout<<"eff_0.size(): "<<eff_0.size()<<"\n";
            for (size_t i = 0; i < eff_0.size(); i++)
            {
                CHECK(z1[i] == eff_0[i]);
                //std::cout<<"eff_0[i]: i "<<i<<" -> "<<eff_0[i]<<"\n";
            }

            for (size_t i = 0; i < eff_1.size(); i++)
                CHECK(z2[i] == eff_1[i]);

            for (size_t i = 0; i < eff_2.size(); i++)
                CHECK(x1[i] == eff_2[i]);

            for (size_t i = 0; i < eff_3.size(); i++)
                CHECK(x2[i] == eff_3[i]);

            model.append_observation("tests/data/y1.dat"); // obs_0
            model.append_observation("tests/data/y2.dat"); // obs_1

            std::vector<float> obs_0 = model.test_observations(0);
            std::vector<float> obs_1 = model.test_observations(1);

            for (size_t i = 0; i < obs_0.size(); i++)
                CHECK(y1[i] == obs_0[i]);

            for (size_t i = 0; i < obs_1.size(); i++)
                CHECK(y2[i] == obs_1[i]);

            CHECK(model.size_of("eff") == 100);
            CHECK(model.size_of("obs") == 10);

            model.append_residual("tests/data/iR.dat");

            std::vector<int> corr_eff{0, 1};

            model.append_corrstruct("tests/data/iG.dat", "tests/data/ainv.dat", corr_eff);

            std::vector<float> res_test = model.test_residual(0);
            std::vector<float> var_test = model.test_variance(0);
            std::vector<float> cor_test = model.test_correlation(0);

            CHECK(res_test.size() == iR.size());
            CHECK(var_test.size() == iG1.size());
            CHECK(cor_test.size() == iA.size());

            for (size_t i = 0; i < res_test.size(); i++)
                CHECK(iR[i] == res_test[i]);

            for (size_t i = 0; i < var_test.size(); i++)
                CHECK(iG1[i] == var_test[i]);

            for (size_t i = 0; i < cor_test.size(); i++)
                CHECK(iA[i] == cor_test[i]);

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing IN -> data] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << e.what() << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (const std::string &err)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing IN -> data] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << err << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (int ierr)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing IN -> data] due to the exception with code " << ierr << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (...)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing IN -> data] due to an unknown exception." << '\n';
            //exit(EXIT_FAILURE);
        }
    }

    // ========================================================================================

    SECTION("The overloaded methods (Model & IO interface): testing vect_z_uni")
    {
        try
        {
            evo::Model model;
            evo::Pcg solver;

            size_t dummy_type1;
            float dummy_type2;

            model.append_effect("tests/data/z1.dat", dummy_type1); // eff_0
            model.append_effect("tests/data/z2.dat", dummy_type1); // eff_1
            model.append_effect("tests/data/x1.dat", dummy_type2); // eff_2
            model.append_effect("tests/data/x2.dat", dummy_type2); // eff_3

            model.append_observation("tests/data/y1.dat"); // obs_0
            model.append_observation("tests/data/y2.dat"); // obs_1

            model.append_residual("tests/data/iR.dat");

            std::vector<int> corr_eff{0, 1};

            model.append_corrstruct("tests/data/iG.dat", "tests/data/ainv.dat", corr_eff);

            std::vector<int> eff_trate_1{2, 0};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{3, 1};
            int obs_trate_2 = 1;

            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            solver.append_model(model);

            evo::matrix<float> z;

            solver.test_vect_z_uni(0, z);

            for (auto i = 0; i < true_z.size(); i++)
            {
                for (auto j = 0; j < true_z[i].size(); j++)
                {
                    CHECK(z(i, j) == true_z[i][j]);
                }
            }

            z.clear();

            solver.test_vect_z_uni(1, z);

            for (auto i = 0; i < true_z.size(); i++)
            {
                for (auto j = 0; j < true_z[i].size(); j++)
                {
                    CHECK(z(i, j) == true_z[i][j]);
                }
            }

            z.clear();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing vect_z_uni] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << e.what() << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (const std::string &err)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing vect_z_uni] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << err << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (int ierr)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing vect_z_uni] due to the exception with code " << ierr << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (...)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing vect_z_uni] due to an unknown exception." << '\n';
            //exit(EXIT_FAILURE);
        }
    }

    // ========================================================================================

    SECTION("The overloaded methods (Model & IO interface): testing dval")
    {
        try
        {
            evo::Model model;
            evo::Pcg solver;

            size_t dummy_type1;
            float dummy_type2;

            model.append_effect("tests/data/x1.dat", dummy_type2); // eff_0
            model.append_effect("tests/data/z1.dat", dummy_type1); // eff_1
            model.append_effect("tests/data/x2.dat", dummy_type2); // eff_2
            model.append_effect("tests/data/z2.dat", dummy_type1); // eff_3

            model.append_observation("tests/data/y1.dat"); // obs_0
            model.append_observation("tests/data/y2.dat"); // obs_1

            model.append_residual("tests/data/iR.dat");

            std::vector<int> corr_eff{1, 3};

            model.append_corrstruct("tests/data/iG.dat", "tests/data/ainv.dat", corr_eff);

            std::vector<int> eff_trate_1{0, 1};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{2, 3};
            int obs_trate_2 = 1;

            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            solver.append_model(model);

            evo::matrix<float> model_rhs = solver.test_rhs();

            for (auto i = 0; i < _rhs.size(); i++)
                CHECK((_rhs[i]) == Approx(model_rhs[i]).margin(0.0001).epsilon(1e-3));

            std::vector<float> dval = solver.test_dval();

            CHECK(dval.size() == dval_true.size());

            for (size_t i = 0; i < dval_true.size(); i++)
                CHECK(dval[i] == Approx(dval_true[i]).margin(0.0001).epsilon(1e-4));

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing dval] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << e.what() << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (const std::string &err)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing dval] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << err << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (int ierr)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing dval] due to the exception with code " << ierr << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (...)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing dval] due to an unknown exception." << '\n';
            //exit(EXIT_FAILURE);
        }
    }

    // ========================================================================================

    SECTION("The overloaded methods (Model & IO interface): testing CoeffMatrix")
    {
        try
        {
            evo::Model model;
            evo::Pcg solver;

            size_t dummy_type1;
            float dummy_type2;

            model.append_effect("tests/data/x1.dat", dummy_type1); // eff_0
            model.append_effect("tests/data/z1.dat", dummy_type1); // eff_1
            model.append_effect("tests/data/x2.dat", dummy_type2); // eff_2
            model.append_effect("tests/data/z2.dat", dummy_type2); // eff_3

            model.append_observation("tests/data/y1.dat"); // obs_0
            model.append_observation("tests/data/y2.dat"); // obs_1

            model.append_residual("tests/data/iR.dat");

            std::vector<int> corr_eff{1, 3};

            model.append_corrstruct("tests/data/iG.dat", "tests/data/ainv.dat", corr_eff);

            std::vector<int> eff_trate_1{0, 1};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{2, 3};
            int obs_trate_2 = 1;

            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            solver.append_model(model);

            std::vector<std::vector<float>> A = solver.test_A();

            CHECK(A.size() == CoeffMatrix.size());
            CHECK(A[0].size() == CoeffMatrix[0].size());

            for (size_t i = 0; i < A.size(); i++)
            {
                for (size_t j = 0; j < A[0].size(); j++)
                {
                    CHECK(A[i][j] == Approx(CoeffMatrix[i][j]).margin(0.001).epsilon(1e-3));
                }
            }

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing CoeffMatrix] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << e.what() << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (const std::string &err)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing CoeffMatrix] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << err << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (int ierr)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing CoeffMatrix] due to the exception with code " << ierr << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (...)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing CoeffMatrix] due to an unknown exception." << '\n';
            //exit(EXIT_FAILURE);
        }
    }

    // ========================================================================================

    SECTION("The overloaded methods (Model & IO interface): testing intermediate data structures")
    {
        try
        {
            evo::Model model;
            evo::Pcg solver;

            size_t dummy_type1;
            float dummy_type2;

            model.append_effect("tests/data/z1.dat", dummy_type1); // eff_0
            model.append_effect("tests/data/z2.dat", dummy_type1); // eff_1
            model.append_effect("tests/data/x1.dat", dummy_type2); // eff_2
            model.append_effect("tests/data/x2.dat", dummy_type2); // eff_3

            model.append_observation("tests/data/y1.dat"); // obs_0
            model.append_observation("tests/data/y2.dat"); // obs_1

            model.append_residual("tests/data/iR.dat");

            std::vector<int> corr_eff{0, 1};

            model.append_corrstruct("tests/data/iG.dat", "tests/data/ainv.dat", corr_eff);

            std::vector<int> eff_trate_1{2, 0};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{3, 1};
            int obs_trate_2 = 1;

            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            solver.append_model(model);

            size_t n_all_levels_checking = solver.test_num_all_levels();

            CHECK(n_all_levels_checking == n_all_levels);

            std::vector<int> ordered_random_levels_checking = solver.test_ordered_levels();

            CHECK(ordered_random_levels_checking.size() == ordered_random_levels.size());

            for (size_t i = 0; i < ordered_random_levels_checking.size(); i++)
                CHECK(ordered_random_levels_checking[i] == ordered_random_levels[i]);

            std::vector<std::vector<size_t>> rcov_offsets_checking = solver.test_cov_offsets();

            CHECK(rcov_offsets_checking.size() == rcov_offsets.size());
            CHECK(rcov_offsets_checking[0].size() == rcov_offsets[0].size());

            for (size_t i = 0; i < rcov_offsets.size(); i++)
            {
                for (size_t j = 0; j < rcov_offsets[0].size(); j++)
                {
                    CHECK(rcov_offsets_checking[i][j] == rcov_offsets[i][j]);
                }
            }

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing intermediate data structures] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << e.what() << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (const std::string &err)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing intermediate data structures] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << err << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (int ierr)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing intermediate data structures] due to the exception with code " << ierr << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (...)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing intermediate data structures] due to an unknown exception." << '\n';
            //exit(EXIT_FAILURE);
        }
    }

    // ========================================================================================

    SECTION("The overloaded methods (Model & IO interface): testing RHS")
    {
        try
        {
            evo::Model model;
            evo::Pcg solver;

            size_t dummy_type1;
            float dummy_type2;

            model.append_effect("tests/data/z1.dat", dummy_type1); // eff_0
            model.append_effect("tests/data/z2.dat", dummy_type1); // eff_1
            model.append_effect("tests/data/x1.dat", dummy_type2); // eff_2
            model.append_effect("tests/data/x2.dat", dummy_type2); // eff_3

            model.append_observation("tests/data/y1.dat"); // obs_0
            model.append_observation("tests/data/y2.dat"); // obs_1

            model.append_residual("tests/data/iR.dat");

            std::vector<int> corr_eff{0, 1};

            model.append_corrstruct("tests/data/iG.dat", "tests/data/ainv.dat", corr_eff);

            std::vector<int> eff_trate_1{2, 0};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{3, 1};
            int obs_trate_2 = 1;

            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            solver.append_model(model);

            evo::matrix<float> model_rhs = solver.test_rhs();

            for (auto i = 0; i < _rhs.size(); i++)
                CHECK((_rhs[i]) == Approx(model_rhs[i]).margin(0.0001).epsilon(1e-3));

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing RHS] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << e.what() << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (const std::string &err)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing RHS] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << err << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (int ierr)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing RHS] due to the exception with code " << ierr << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (...)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing RHS] due to an unknown exception." << '\n';
            //exit(EXIT_FAILURE);
        }
    }

    // ========================================================================================

    SECTION("The overloaded methods (Model & IO interface): testing solution")
    {
        try
        {
            evo::Model model;
            evo::Pcg solver;

            size_t dummy_type1;
            float dummy_type2;

            model.append_effect("tests/data/z2.dat", dummy_type1); // eff_0
            model.append_effect("tests/data/x1.dat", dummy_type1); // eff_1
            model.append_effect("tests/data/z1.dat", dummy_type2); // eff_2
            model.append_effect("tests/data/x2.dat", dummy_type2); // eff_3

            model.append_observation("tests/data/y1.dat"); // obs_0
            model.append_observation("tests/data/y2.dat"); // obs_1

            model.append_residual("tests/data/iR.dat");

            std::vector<int> eff_trate_1{1, 2};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{3, 0};
            int obs_trate_2 = 1;

            std::vector<int> corr_eff{1, 3};
            std::vector<int> corr_eff2{0, 2}; // same matrix but diff order of effects

            // model.append_corrstruct("tests/data/iG.dat", "tests/data/ainv.dat", corr_eff);
            model.append_corrstruct("tests/data/iG_diff_order.dat", "tests/data/ainv.dat", corr_eff2);

            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            solver.append_model(model);

            solver.solve();

            std::vector<float> sol = solver.get_solution();

            solver.get_solution("cpp_solution_model_1.dat");

            for (auto i = 0; i < sol.size(); i++)
            {
                CHECK((_sol[i]) == Approx(sol[i]).margin(0.0001).epsilon(1e-3));
                // std::cout<<sol[i]<<"\n";
            }

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing solution] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << e.what() << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (const std::string &err)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing solution] due to the exception:" << '\n';
            std::cerr << " Model 1.  => " << err << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (int ierr)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing solution] due to the exception with code " << ierr << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (...)
        {
            std::cerr << " Model 1. Exit from [The overloaded methods (Model & IO interface): testing solution] due to an unknown exception." << '\n';
            //exit(EXIT_FAILURE);
        }
    }

    // ========================================================================================
    
}

TEST_CASE("Small data test: model 2")
{
    // one trait setup
    // EXAMPLE: 7.1. p.110.
    // y = Xb + Ua + Wm + Sp + e
    // corr1: [a m]
    // corr2: [p]

    // changes to model interpretation:
    // X moves to random effects but without covariance structure
    // corr1: [a m]
    // corr2: [p]
    // corr0 == 0: [b], because it is a fixed effect

    // ----------- DATA ----------------------------
    std::vector<int> x{
        1, 0, 0, 1, 0,
        1, 0, 0, 0, 1,
        1, 0, 0, 0, 1,
        1, 0, 0, 1, 0,
        0, 1, 0, 1, 0,
        0, 1, 0, 0, 1,
        0, 1, 0, 0, 1,
        0, 0, 1, 0, 1,
        0, 0, 1, 1, 0,
        0, 0, 1, 0, 1};

    std::vector<int> s{
        1, 0, 0, 0,
        1, 0, 0, 0,
        0, 0, 1, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        1, 0, 0, 0,
        0, 0, 0, 1,
        0, 0, 0, 1,
        1, 0, 0, 0,
        0, 0, 1, 0};

    std::vector<int> w{
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0};

    std::vector<int> u{
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

    std::vector<float> iA{
        2, 0.5, 0, 0, -1, 0.5, 0, 0, -1, 0, 0, 0, 0, 0,
        0.5, 3, 1, 0, -1, -1, 0, 0, 0.5, -1, 0, 0, -1, 0,
        0, 1, 3.5, 0, 0.5, -0.5, 0.5, -1, 0, -1, -1, 0, 0, -1,
        0, 0, 0, 1.5, 0, 0.5, -1, 0, 0, 0, 0, 0, 0, 0,
        -1, -1, 0.5, 0, 2.5, 0, 0, -1, 0, 0, 0, 0, 0, 0,
        0.5, -1, -0.5, 0.5, 0, 3.5, -1, 0, -1, 0, 0, 0, 0, -1,
        0, 0, 0.5, -1, 0, -1, 3, 0.5, 0, 0, -1, -1, 0, 0,
        0, 0, -1, 0, -1, 0, 0.5, 2.5, 0, 0, 0, -1, 0, 0,
        -1, 0.5, 0, 0, 0, -1, 0, 0, 2.5, 0, 0, 0, -1, 0,
        0, -1, -1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0,
        0, 0, -1, 0, 0, 0, -1, 0, 0, 0, 2, 0, 0, 0,
        0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 2, 0, 0,
        0, -1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 2, 0,
        0, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 2};

    std::vector<float> corS{
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1};

    std::vector<float> y{35, 20, 25, 40, 42, 22, 35, 34, 20, 40};

    std::vector<float> iR{0.0029};

    std::vector<float> iG1{1};

    std::vector<float> iG2{0.0076, 0.0034, 0.0034, 0.0126};

    std::vector<float> iG3{0.025};

    std::vector<float> _sol{
        14.492235,
        17.878691,
        15.926014,
        20.047707,
        13.198704,
        0.56389,
        -1.244389,
        1.164936,
        -0.484436,
        0.629533,
        -0.858629,
        -1.155968,
        1.917404,
        -0.553263,
        -1.055081,
        0.385358,
        0.863317,
        -2.979573,
        1.751301,
        0.261565,
        -1.583161,
        0.735732,
        0.585864,
        -0.507471,
        0.841041,
        1.299317,
        -0.157915,
        0.659541,
        -0.152954,
        0.915958,
        0.442008,
        0.093056,
        0.362213,
        -1.70083,
        0.415397,
        0.824915,
        0.460519};

    // ---------------------------------------------

    evo::Pcg solver;
    evo::Model model;

    // ----- define the model -------

    model.append_residual(iR, 2);

    model.append_observation(y, 10); // obs := 0

    model.append_effect(x, 10, 5);  // eff := 0
    model.append_effect(s, 10, 4);  // eff := 1
    model.append_effect(w, 10, 14); // eff := 2
    model.append_effect(u, 10, 14); // eff := 3

    std::vector<int> corr_eff_1{3, 2};
    std::vector<int> corr_eff_2{1};

    std::vector<int> eff_trate{0, 3, 2, 1};
    int obs_trate = 0;

    model.append_corrstruct(iG2, 2, iA, 14, corr_eff_1);

    std::string identity("I");
    // model.append_corrstruct(iG3, 1, corS, 4, corr_eff_2);
    model.append_corrstruct(iG3, 1, identity, 4, corr_eff_2);

    model.append_traitstruct(obs_trate, eff_trate);

    solver.append_model(model);

    // ========================================================================================

    SECTION("The Pcg class: test solution")
    {
        try
        {
            solver.solve();

            std::vector<float> sol = solver.get_solution();

            solver.get_solution("cpp_solution_model_2.dat");

            for (auto i = 0; i < sol.size(); i++)
                CHECK((_sol[i]) == Approx(sol[i]).margin(0.03).epsilon(1e-3));

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << " Model 2. Exit from [The Pcg class: test solution] due to the exception:" << '\n';
            std::cerr << " Model 2.  => " << e.what() << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (const std::string &err)
        {
            std::cerr << " Model 2. Exit from [The Pcg class: test solution] due to the exception:" << '\n';
            std::cerr << " Model 2.  => " << err << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (int ierr)
        {
            std::cerr << " Model 2. Exit from [The Pcg class: test solution] due to the exception with code " << ierr << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (...)
        {
            std::cerr << " Model 2. Exit from [The Pcg class: test solution] due to an unknown exception." << '\n';
            //exit(EXIT_FAILURE);
        }
    }

    // ========================================================================================
}

TEST_CASE("Small data test: model 3")
{
    // one trait setup
    // EXAMPLE: 11.2. p.183.
    // y = Xb + Za + e

    // corr: [a], identity matrix

    // ----------- DATA ----------------------------
    std::vector<int> x{
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1};

    std::vector<float> z{
        1.357, -0.357, 0.286, 0.286, -0.286, -1.214, -0.143, 0.071, -0.143, 1.214,
        0.357, -0.357, -0.714, -0.714, -0.286, 0.786, -0.143, 0.071, -0.143, -0.786,
        0.357, 0.643, 1.286, 0.286, 0.714, -1.214, -0.143, 0.071, -0.143, 1.214,
        -0.643, -0.357, 1.286, 0.286, -0.286, -0.214, -0.143, 0.071, 0.857, 0.214,
        -0.643, 0.643, 0.286, 1.286, -0.286, -1.214, -0.143, 0.071, -0.143, 1.214,
        0.357, 0.643, -0.714, 0.286, -0.286, 0.786, -0.143, 0.071, 0.857, 0.214,
        -0.643, -0.357, 0.286, 0.286, -0.286, 0.786, -0.143, 0.071, 0.857, -0.786,
        -0.643, 0.643, 0.286, -0.714, -0.286, -0.214, -0.143, 0.071, 0.857, -0.786};

    std::vector<float> corZ{
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

    std::vector<float> y{9.0, 13.4, 12.7, 15.4, 5.9, 7.7, 10.2, 4.8};

    std::vector<float> iR{0.0041};

    std::vector<float> iG1{0.1004};

    std::vector<float> _sol{
        9.94415,
        0.0873718,
        -0.312064,
        0.263549,
        -0.0805778,
        0.110685,
        0.139673,
        -2.3004e-07,
        -1.62874e-07,
        -0.0609667,
        -0.0158181};

    // ---------------------------------------------

    evo::Pcg solver;
    evo::Model model;

    // ----- define the model -------

    model.append_residual(iR, 1);

    model.append_observation(y, 8); // obs := 0

    model.append_effect(x, 8, 1);  // eff := 0
    model.append_effect(z, 8, 10); // eff := 1

    std::vector<int> corr_eff{1};

    std::vector<int> eff_trate{0, 1};
    int obs_trate = 0;

    model.append_corrstruct(iG1, 1, corZ, 10, corr_eff);

    model.append_traitstruct(obs_trate, eff_trate);

    solver.append_model(model);

    // ========================================================================================

    SECTION("The Pcg class: test solution")
    {
        try
        {
            solver.solve();

            std::vector<float> sol = solver.get_solution();

            solver.get_solution("cpp_solution_model_3.dat");

            for (auto i = 0; i < sol.size(); i++)
                CHECK((_sol[i]) == Approx(sol[i]).margin(0.003).epsilon(1e-3));

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << " Model 3. Exit from [The Pcg class: test solution] due to the exception:" << '\n';
            std::cerr << " Model 3.  => " << e.what() << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (const std::string &err)
        {
            std::cerr << " Model 3. Exit from [The Pcg class: test solution] due to the exception:" << '\n';
            std::cerr << " Model 3.  => " << err << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (int ierr)
        {
            std::cerr << " Model 3. Exit from [The Pcg class: test solution] due to the exception with code " << ierr << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (...)
        {
            std::cerr << " Model 3. Exit from [The Pcg class: test solution] due to an unknown exception." << '\n';
            //exit(EXIT_FAILURE);
        }
    }

    // ========================================================================================

    // ========================================================================================

    SECTION("The Pcg class: test identity corr structure")
    {
        try
        {
            model.clear();
            solver.remove_model();

            model.append_residual(iR, 1);

            model.append_observation(y, 8); // obs := 0

            model.append_effect(x, 8, 1);  // eff := 0
            model.append_effect(z, 8, 10); // eff := 1

            std::vector<int> corr_eff{1};

            std::vector<int> eff_trate{0, 1};
            int obs_trate = 0;

            std::string identity("I");

            model.append_corrstruct(iG1, 1, identity, 10, corr_eff);
            // model.append_corrstruct(iG1, 1, corZ, 10, corr_eff);
            // model.append_corrstruct("tests/data/model_3/iG.dat", identity, 10, corr_eff);

            model.append_traitstruct(obs_trate, eff_trate);

            // model.print();

            solver.append_model(model);

            solver.solve();

            std::vector<float> sol = solver.get_solution();

            solver.get_solution("cpp_solution_model_3(I).dat");

            for (auto i = 0; i < sol.size(); i++)
                CHECK((_sol[i]) == Approx(sol[i]).margin(0.003).epsilon(1e-3));

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << " Model 3. Exit from [The Pcg class: test identity corr structure] due to the exception:" << '\n';
            std::cerr << " Model 3.  => " << e.what() << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (const std::string &err)
        {
            std::cerr << " Model 3. Exit from [The Pcg class: test identity corr structure] due to the exception:" << '\n';
            std::cerr << " Model 3.  => " << err << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (int ierr)
        {
            std::cerr << " Model 3. Exit from [The Pcg class: test identity corr structure] due to the exception with code " << ierr << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (...)
        {
            std::cerr << " Model 3. Exit from [The Pcg class: test identity corr structure] due to an unknown exception." << '\n';
            //exit(EXIT_FAILURE);
        }
    }

    // ========================================================================================
}

TEST_CASE("Small data test: model 4")
{
    SECTION("The Pcg class: test full SNP blup")
    {
        try
        {
            evo::Pcg solver;
            evo::Model model;

            std::vector<float> iR{0.041};

            std::vector<float> iG1{0.1};

            model.append_residual(iR, 1);

            model.append_observation("tests/data/model_4/obs_1.dat"); // obs := 0

            size_t type1;
            float type2;

            model.append_effect("tests/data/model_4/obs_489_snp_100.txt", type1); // eff := 0
            model.append_effect("tests/data/model_4/fixed_1.dat", type2);          // eff := 1

            std::vector<int> corr_eff{0};

            std::vector<int> eff_trate{1, 0};
            int obs_trate = 0;

            std::string identity("I");

            model.append_corrstruct(iG1, 1, identity, 100, corr_eff);

            model.append_traitstruct(obs_trate, eff_trate);

            solver.append_model(model);

            solver.solve();

            //std::vector<float> sol = solver.get_solution();

            solver.get_solution("cpp_solution_model_4.dat");

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << " Model 4. Exit from [The Pcg class: test full SNP blup] due to the exception:" << '\n';
            std::cerr << " Model 4.  => " << e.what() << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (const std::string &err)
        {
            std::cerr << " Model 4. Exit from [The Pcg class: test full SNP blup] due to the exception:" << '\n';
            std::cerr << " Model 4.  => " << err << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (int ierr)
        {
            std::cerr << " Model 4. Exit from [The Pcg class: test full SNP blup] due to the exception with code " << ierr << '\n';
            //exit(EXIT_FAILURE);
        }
        catch (...)
        {
            std::cerr << " Model 4. Exit from [The Pcg class: test full SNP blup] due to an unknown exception." << '\n';
            //exit(EXIT_FAILURE);
        }
    }
}
