#include <cstdint>
#include <array>
#include <tuple>

#include "core/image.h"
#include "core/point.h"
#include "core/vector_field.h"

#include "piv/evaluation.h"

// pybind
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

namespace py = pybind11;

using namespace openpiv;


py::array_t<double> openpiv_to_numpy(
        const core::grid_coords& field_coords,
        const core::grid_data& field_data
    ) {
    const uint32_t h = field_coords.height();
    const uint32_t w = field_coords.width();
    const size_t H = static_cast<size_t>(h);
    const size_t W = static_cast<size_t>(w);
    const size_t C = 6;
    const size_t total = C * H * W;

    // create array with shape (6, H, W)
    py::array_t<double> arr({C, H, W});

    // fill buffer directly
    double *buf = arr.mutable_data();

    // precompute strides in elements: stride_C = H*W, stride_H = W, stride_W = 1
    const size_t strideC = H * W;
    const size_t strideH = W;

    for (size_t y = 0; y < H; ++y) {
        for (size_t x = 0; x < W; ++x) {
            const size_t base = y * W + x;
            
            auto coords = core::point2<uint32_t>(static_cast<uint32_t>(x), static_cast<uint32_t>(y));

            auto p = field_coords[coords];

            buf[0 * strideC + base] = p[0];
            buf[1 * strideC + base] = p[1];
            buf[2 * strideC + base] = field_data.u[coords];
            buf[3 * strideC + base] = field_data.v[coords];
            buf[4 * strideC + base] = field_data.s2n[coords];
            buf[5 * strideC + base] = field_data.peak[coords];
        }
    }

    return arr;
}


void add_piv_firstpass(py::module& m)
{
    using ImageT = core::image_gf64;

     m.def("process_images_standard",
          [](ImageT image_a, 
             ImageT image_b,
             std::array<uint32_t,2> window_size,
             std::array<uint32_t,2> overlap_size,
             bool zero_pad, 
             bool centered, 
             bool limit_search, 
             int32_t threads) -> py::array_t<double>
          {
              auto [coords, data] = piv::process_images_standard(
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