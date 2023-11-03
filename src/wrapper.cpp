#include "solver_pcg.hpp"
#include "model.hpp"

PYBIND11_MODULE(evolm, m)
{
    m.doc() = "pybind11 example plugin"; // Optional module docstring

    pybind11::class_<evo::Model>(m, "Model")

        .def(pybind11::init<>())
#ifdef UTEST
        .def("size_of", &evo::Model::size_of)
        .def("shape_of", &evo::Model::shape_of)
        .def("print", &evo::Model::print)
#endif
        .def("clear_residuals", &evo::Model::clear_residuals)
        .def("clear_observations", &evo::Model::clear_observations)
        .def("clear_effects", &evo::Model::clear_effects)
        .def("clear_corrstruct", &evo::Model::clear_corrstruct)
        .def("clear_traitstruct", &evo::Model::clear_traitstruct)
        .def("clear", &evo::Model::clear)        

        .def( "append_residual", static_cast< int (evo::Model::*)(pybind11::array_t<float>, size_t) >(&evo::Model::append_residual) )
        .def( "append_residual", static_cast< int (evo::Model::*)(const std::string &) >(&evo::Model::append_residual) )        
        .def( "append_observation", static_cast< int (evo::Model::*)(pybind11::array_t<float>, size_t) >(&evo::Model::append_observation) )
        .def( "append_observation", static_cast< int (evo::Model::*)(const std::string &) >(&evo::Model::append_observation) )
        
        .def( "append_effect", static_cast< int (evo::Model::*)(pybind11::array_t<float>, size_t, size_t) >(&evo::Model::append_effect) )

        .def( "append_effect", static_cast< int (evo::Model::*)(const std::string &, size_t) >(&evo::Model::append_effect) )
        .def( "append_effect", static_cast< int (evo::Model::*)(const std::string &, bool) >(&evo::Model::append_effect) )
        .def( "append_effect", static_cast< int (evo::Model::*)(const std::string &, float) >(&evo::Model::append_effect) )
        .def( "append_effect", static_cast< int (evo::Model::*)(const std::string &, double) >(&evo::Model::append_effect) )

        .def( "append_corrstruct", static_cast< int (evo::Model::*)(pybind11::array_t<float>, size_t, pybind11::array_t<float>, size_t, pybind11::array_t<int>) >(&evo::Model::append_corrstruct) )
        .def( "append_corrstruct", static_cast< int (evo::Model::*)(pybind11::array_t<float>, size_t, const std::string &, pybind11::array_t<int>) >(&evo::Model::append_corrstruct) )
        .def( "append_corrstruct", static_cast< int (evo::Model::*)(const std::string &, const std::string &, pybind11::array_t<int>) >(&evo::Model::append_corrstruct) )

        .def( "append_corrstruct", static_cast< int (evo::Model::*)(pybind11::array_t<float>, size_t, std::string &, size_t, pybind11::array_t<int>) >(&evo::Model::append_corrstruct) )
        .def( "append_corrstruct", static_cast< int (evo::Model::*)(const std::string &, std::string &, size_t, pybind11::array_t<int>) >(&evo::Model::append_corrstruct) )
        
        .def( "append_traitstruct", static_cast< int (evo::Model::*)(int, pybind11::array_t<int>) >(&evo::Model::append_traitstruct) );
    
    pybind11::class_<evo::Pcg>(m, "Pcg")
        .def(pybind11::init<>())
        .def("append_model", &evo::Pcg::append_model)
        .def("solve", pybind11::overload_cast<>(&evo::Pcg::solve))
        .def("solve", pybind11::overload_cast<int>(&evo::Pcg::solve))
        .def( "get_solution", static_cast< int (evo::Pcg::*)(const std::string &) >(&evo::Pcg::get_solution) )
        .def( "get_solution", static_cast< pybind11::array_t<float> (evo::Pcg::*)() >(&evo::Pcg::get_solution) );

}
