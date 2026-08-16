sources = [
    'core/grid.cpp',
    'core/log.cpp',
    'core/point.cpp',
    'core/rect.cpp',
    'core/size.cpp',
    'core/vector.cpp',
    'core/pixel_types.cpp',
    'core/image.cpp',
    'core/vector_field.cpp',
    'filter/gaussian_lowpass.cpp',
    'filter/gaussian_highpass.cpp',
    'filter/gaussian_variance.cpp',
    'interp/interp2d.cpp',
    'piv/evaluation_standard.cpp',
    'piv/piv_enums.cpp',
    'validate/validate_difference.cpp',
    'validate/validate_median.cpp',
    'validate/validate_normalized_median.cpp',
    'validate/validate_peak_ratio.cpp',
    'validate/validate_range.cpp',
    'validate/validate_s2n.cpp',
    'validate/validate_z_score.cpp',
    'module/module.cpp'
]

py.extension_module('pyopenpivcore',
    sources,
    include_directories: [
        openpiv_include, 
        threadlib_include,
        pybind_include
    ],
    cpp_args: [
        '-D_USE_MATH_DEFINES'
    ],
    dependencies: [
        fmt_dep,
        thread_dep,
        duccfft_dep
    ],
    link_with: openpiv_lib,
    install: true,
    subdir: 'openpiv_cxx'
)
