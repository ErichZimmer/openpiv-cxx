#include <cstdint>
#include <array>
#include <tuple>

#include "core/image.h"
#include "core/point.h"
#include "core/vector_field.h"

#include "filters/gaussian_lowpass.h"

#include "piv/correlation_utils.h"

// pybind
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>


namespace py = pybind11;

using namespace openpiv;


void add_gaussian_lowpass(py::module& m)
{
     m.def("gaussian_lowpass",
          [](const piv::ImageT& src, 
             piv::ImageT& out,
             piv::FloatT sigma,
             piv::FloatT truncate)
          {
              filter::gaussian_lowpass(
                src,
                out,
                sigma,
                truncate
              );
        },

        py::arg("src"),
        py::arg("out"),
        py::arg("sigma"),
        py::arg("truncate") = 4.0
    );

    m.def("gaussian_lowpass",
          [](const piv::ImageT& src, 
             piv::FloatT sigma,
             piv::FloatT truncate) -> piv::ImageT
          {
              auto out = filter::gaussian_lowpass(
                src,
                sigma,
                truncate
              );

              return out;
        },

        py::arg("src"),
        py::arg("sigma"),
        py::arg("truncate") = 4.0
    );
}