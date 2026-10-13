#include "core/image.h"
#include "core/vector_field.h"

#include "validate/impute_nans.h"
#include "validate/substitute_nans.h"

// pybind
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

namespace py = pybind11;

using namespace openpiv;


void add_validate_replacement(py::module& m)
{

     m.def("substitute_nans",
          [](core::grid_data field_data)
        {
            uint32_t failed_substitutions = validate::substitute_nans<core::grid_data_t>(
                field_data
            );

            return failed_substitutions;
        },

        py::arg("field_data")
    );

    m.def("impute_nans",
          [](core::grid_data field_data,
             bool include_corners,
             double tolerance,
             uint32_t max_iterations
            )
        {

            core::grid_data_t tolerance_c = core::grid_data_t(tolerance);

            bool converged = validate::impute_nans<core::grid_data_t>(
                field_data,
                include_corners,
                tolerance_c,
                max_iterations
            );

            return converged;
        },

        py::arg("field_data"),
        py::arg("include_corners") = false,
        py::arg("tolerance") = 0.0001,
        py::arg("max_iterations") = 10
    );
}