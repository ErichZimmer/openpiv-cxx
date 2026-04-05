#include "piv/piv_common.h"

// pybind
#include <pybind11/pybind11.h>
#include <pybind11/native_enum.h>

namespace py = pybind11;

using namespace openpiv;


void add_piv_method(py::module& m)
{
    py::enum_<openpiv::piv::deform_method>(m, "deform_method")
        .value("SINC", piv::deform_method::SINC)
        .value("LAGRANGE", piv::deform_method::LAGRANGE)
        .export_values();
}

void add_piv_order(py::module& m)
{
    py::enum_<openpiv::piv::deform_order>(m, "deform_order")
        .value("FORWARD", piv::deform_order::FORWARD)
        .value("SYMMETRIC", piv::deform_order::SYMMETRIC)
        .export_values();
}