#include "core/image.h"
#include "core/vector_field.h"

#include "interp/interp2d.h"

// pybind
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

namespace py = pybind11;

using namespace openpiv;


void add_interp2d(py::module& m)
{

     m.def("interp2d",
          [](const core::grid_coords& field_coords,
            const core::image<core::grid_data_t>& data,
            const core::grid_coords& mappings,
            uint32_t k,
            int32_t threads)
        {
            auto out = interp::interp2d<core::image, core::grid_data_t>(
                field_coords,
                data,
                mappings,
                k,
                threads
            );

            return out;
        },

        py::arg("field_coords"),
        py::arg("data"),
        py::arg("mappings"),
        py::arg("k"),
        py::arg("threads")
    );
}