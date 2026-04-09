#include <cstdint>
#include <array>
#include <tuple>

#include "core/image.h"
#include "core/point.h"
#include "core/vector_field.h"

#include "filters/gaussian_variance.h"

#include "piv/correlation_utils.h"

// pybind
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>


namespace py = pybind11;

using namespace openpiv;


void add_gaussian_variance(py::module& m)
{
     m.def("gaussian_variance",
          [](const piv::ImageT& src, 
             piv::ImageT& out,
             piv::FloatT sigma1,
             piv::FloatT sigma2,
             piv::FloatT truncate,
            bool clip)
          {
              filter::gaussian_variance(
                src,
                out,
                sigma1,
                sigma2,
                truncate,
                clip
              );
        },

        py::arg("src"),
        py::arg("out"),
        py::arg("sigma1"),
        py::arg("sigma2"),
        py::arg("truncate") = 4.0,
        py::arg("clip") = true
    );

    m.def("gaussian_variance",
          [](const piv::ImageT& src, 
             piv::FloatT sigma1,
             piv::FloatT sigma2,
             piv::FloatT truncate,
            bool clip) -> piv::ImageT
          {
              auto out = filter::gaussian_variance(
                src,
                sigma1,
                sigma2,
                truncate,
                clip
              );

              return out;
        },

        py::arg("src"),
        py::arg("sigma1"),
        py::arg("sigma2"),
        py::arg("truncate") = 4.0,
        py::arg("clip") = true
    );
}