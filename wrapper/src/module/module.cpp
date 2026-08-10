
// pybind
#include <pybind11/pybind11.h>

namespace py = pybind11;

void add_grid(py::module &);
void add_log(py::module &);
void add_point(py::module &);
void add_rect(py::module &);
void add_size(py::module &);
void add_vector(py::module &);
void add_pixel_types(py::module &);
void add_image(py::module &);
void add_gaussian_lowpass(py::module &);
void add_gaussian_highpass(py::module &);
void add_gaussian_variance(py::module &);
void add_grid_coords(py::module &);
void add_grid_data(py::module &);
void add_interp2d(py::module &);
void add_ducc_simd_backend(py::module &);
void add_piv_firstpass(py::module &);
void add_piv_multipass(py::module &);
void add_piv_method(py::module &);
void add_piv_order(py::module &);
void add_validate_difference(py::module &);
void add_validate_median(py::module &);
void add_validate_normalized_median(py::module &);
void add_validate_peak_ratio(py::module &);
void add_validate_range(py::module &);
void add_validate_s2n(py::module &);
void add_validate_z_score(py::module &);

PYBIND11_MODULE(pyopenpivcore, m) {
    m.doc() = R"pbdoc(
        openpivcore python plugin
        -------------------------

        .. currentmodule:: pyopenpivcore

        .. autosummary::
           :toctree: _generate

        size
    )pbdoc";

    // add each binding chunk here
    add_grid(m);
    add_log(m);
    add_point(m);
    add_rect(m);
    add_size(m);
    add_vector(m);
    add_pixel_types(m);
    add_image(m);
    add_gaussian_lowpass(m);
    add_gaussian_highpass(m);
    add_gaussian_variance(m);
    add_grid_coords(m);
    add_grid_data(m);
    add_interp2d(m);
    add_ducc_simd_backend(m);
    add_piv_method(m);
    add_piv_order(m);
    add_piv_firstpass(m);
    add_piv_multipass(m);
    add_validate_difference(m);
    add_validate_median(m);
    add_validate_normalized_median(m);
    add_validate_peak_ratio(m);
    add_validate_range(m);
    add_validate_s2n(m);
    add_validate_z_score(m);


#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
