#include <cstdint>
#include <array>
#include <tuple>

#include "core/image.h"
#include "core/point.h"
#include "core/vector_field.h"

#include "piv/firstpass.h"
#include "piv/multipass.h"
#include "piv/correlation_utils.h"
#include "piv/piv_common.h"

// pybind
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

namespace py = pybind11;

using namespace openpiv;


void add_piv_firstpass(py::module& m)
{
     m.def("process_images_standard",
          [](const piv::ImageT& image_a, 
             const piv::ImageT& image_b,
             std::array<uint32_t,2> window_size,
             std::array<uint32_t,2> overlap_size,
             bool step,
             bool zero_pad, 
             bool centered, 
             bool limit_search, 
             int32_t threads) -> py::tuple
          {
              auto [coords, data] = piv::process_images_standard(
                image_a, 
                image_b,
                window_size, 
                overlap_size,
                step,
                zero_pad, 
                centered, 
                limit_search, 
                threads
            );

            return py::make_tuple(std::move(coords), std::move(data));
        },

        py::arg("image_a"),
        py::arg("image_b"),
        py::arg("window_size"),
        py::arg("overlap_size"),
        py::arg("step") = false,
        py::arg("zero_pad") = false,
        py::arg("centered") = false,
        py::arg("limit_search") = false,
        py::arg("threads") = 1
    );
}


void add_piv_multipass(py::module& m)
{
    using ImageT = core::image_gf32;

     m.def("process_images_multipass",
          [](const piv::ImageT& image_a, 
             const piv::ImageT& image_b,
             const core::grid_coords& old_coords,
             const core::grid_data&   old_data,
             std::array<uint32_t, 2> window_size,
             std::array<uint32_t, 2> overlap_size,
             piv::deform_method method,
             piv::deform_order order,
             int32_t k,
             bool step,
             bool zero_pad, 
             bool centered, 
             bool limit_search, 
             int32_t threads) -> py::tuple
          {
              auto [coords, data] = piv::process_images_standard_multi(
                image_a, 
                image_b,
                old_coords,
                old_data,
                window_size, 
                overlap_size,
                method,
                order,
                k,
                step,
                zero_pad, 
                centered, 
                limit_search, 
                threads
            );

            return py::make_tuple(std::move(coords), std::move(data));
        },

        py::arg("image_a"),
        py::arg("image_b"),
        py::arg("old_coords"),
        py::arg("old_data"),
        py::arg("window_size"),
        py::arg("overlap_size"),
        py::arg("method") = piv::deform_method::LAGRANGE,
        py::arg("order") = piv::deform_order::FORWARD,
        py::arg("k") = 3,
        py::arg("step") = false,
        py::arg("zero_pad") = false,
        py::arg("centered") = false,
        py::arg("limit_search") = false,
        py::arg("threads") = 1
    );


    /*
    m.def("process_images_nsqe",
          [](ImageT image_a, 
             ImageT image_b,
             std::array<uint32_t,2> window_size,
             std::array<uint32_t,2> overlap_size,
             bool zero_pad, 
             bool centered, 
             int32_t threads) -> py::array_t<double>
          {
              auto [coords, data] = piv::process_images_nsqe(
                image_a, 
                image_b,
                window_size, 
                overlap_size,
                zero_pad, 
                centered, 
                limit_search, 
                threads
            );

            return openpiv_to_numpy(coords, data);
        },

        py::arg("image_a"),
        py::arg("image_b"),
        py::arg("window_size"),
        py::arg("overlap_size"),
        py::arg("zero_pad") = false,
        py::arg("centered") = false,
        py::arg("limit_search") = false,
        py::arg("threads") = 1
    */
};