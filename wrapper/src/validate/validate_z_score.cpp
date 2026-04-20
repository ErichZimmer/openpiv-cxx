#include "core/image.h"
#include "core/vector_field.h"

#include "validate/validate_z_score.h"

// pybind
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

namespace py = pybind11;

using namespace openpiv;


void add_validate_z_score(py::module& m)
{

     m.def("validate_z_score",
          [](const core::grid_data field_data,
             float threshold)
          {
              auto flag = validate::validate_z_score<core::grid_data_t>(
                field_data,
                threshold
              );

            return flag;
        },

        py::arg("field_data"),
        py::arg("threshold")
    );
}