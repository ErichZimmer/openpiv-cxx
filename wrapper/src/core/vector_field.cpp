// openpiv
#include "core/image.h"
#include "core/util.h"
#include "core/vector_field.h"

using namespace openpiv;
using namespace openpiv::core;

// local
#include "pyutils.h"

// std
#include <sstream>
#include <vector>

// pybind
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

namespace py = pybind11;

void add_grid_coords(py::module& m)
{
    using U = double;
    using PU = core::point2<U>;

    py::class_<core::grid_coords>(m, "grid_coords", py::buffer_protocol())
        .def_buffer([](core::grid_coords& coords) -> py::buffer_info {
            const auto dimensions = image_type_trait<core::image, core::point2, U>::dimensions(coords);
            const auto strides = image_type_trait<core::image, core::point2, U>::strides(coords);
            return py::buffer_info(
                coords.data(),                      // pointer to raw data
                sizeof( U ),                        // size of a single underlying data element
                py::format_descriptor<U>::format(), // mapping to python char code for underlying
                dimensions.size(),                  // number of dimensions
                dimensions,                         // dimensions in matrix convention
                strides
            ); 
        })
        .def(py::init())
        .def(py::init<core::grid_coords>())
        .def(py::init<size>())
        .def(py::init<uint32_t, uint32_t>())
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def("width", &core::grid_coords::width)
        .def("height", &core::grid_coords::height)
        .def("size", &core::grid_coords::size)
        .def("pixel_count", &core::grid_coords::pixel_count)
        .def("__getitem__",
             [](const core::grid_coords& p, size_t i) {
                 return p[i];
             })
        .def("__setitem__",
             [](core::grid_coords& p, size_t i, PU v) {
                 p[i] = v;
             })
        .def("__getitem__",
             [](const core::grid_coords& p, const point2<uint32_t>& i) {
                 return p[i];
             })
        .def("__setitem__",
             [](core::grid_coords& p, const point2<uint32_t>& i, PU v) {
                 p[i] = v;
             });             
}

void add_grid_data(py::module& m)
{
    py::enum_<core::FLAG>(m, "FLAG")
        .value("VALID", core::FLAG::VALID)
        .value("INVALID", core::FLAG::INVALID)
        .value("SUBSTITUTED", core::FLAG::SUBSTITUTED)
        .value("INTERPOLATED", core::FLAG::INTERPOLATED)
        .value("MASKED", core::FLAG::MASKED);

    py::class_<grid_data>(m, "grid_data")
        .def(py::init<const grid_data&>())
        .def(py::init<uint32_t, uint32_t>())
        .def(py::init<core::size>())
        // expose public members as read/write attributes
        .def_readwrite("u",    &grid_data::u)
        .def_readwrite("v",    &grid_data::v)
        .def_readwrite("s2n",  &grid_data::s2n)
        .def_readwrite("p2p",  &grid_data::p2p)
        .def_readwrite("peak", &grid_data::peak)
        .def_readwrite("u1",    &grid_data::u1)
        .def_readwrite("v1",    &grid_data::v1)
        .def_readwrite("s2n1",  &grid_data::s2n1)
        .def_readwrite("peak1", &grid_data::peak1)
        .def_readwrite("u2",    &grid_data::u2)
        .def_readwrite("v2",    &grid_data::v2)
        .def_readwrite("s2n2",  &grid_data::s2n2)
        .def_readwrite("peak2", &grid_data::peak2)
        .def_readwrite("u3",    &grid_data::u3)
        .def_readwrite("v3",    &grid_data::v3)
        .def_readwrite("s2n3",  &grid_data::s2n3)
        .def_readwrite("peak3", &grid_data::peak3)
        .def_readwrite("peak_ind", &grid_data::peak_ind)
        .def_readwrite("flag",  &grid_data::flag)
        .def("width", [](const grid_data& data) {
            return data.u.width();
        })
        .def("height", [](const grid_data& data) {
            return data.u.height();
        })
        .def("size", [](const grid_data& data) {
            return data.u.size();
        })
        .def("pixel_count", [](const grid_data& data) {
            return data.u.pixel_count();
        });
}
    


