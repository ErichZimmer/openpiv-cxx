#include "core/image.h"
#include "core/vector_field.h"

#include "validate/validate_median.h"

// pybind
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

namespace py = pybind11;

using namespace openpiv;


void add_validate_median(py::module& m)
{

     m.def("validate_median",
          [](const core::grid_data field_data,
             float threshold_u,
             float threshold_v)
          {
              auto flag = validate::validate_median<core::grid_data_t>(
                field_data,
                threshold_u,
                threshold_v
              );

            return flag;
        },

        py::arg("field_data"),
        py::arg("threshold_u"),
        py::arg("threshold_v")
    );
}