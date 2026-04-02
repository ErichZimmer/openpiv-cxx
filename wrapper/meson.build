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
    'piv/firstpass.cpp',
    'module/module.cpp'
]

py.extension_module('pyopenpivcore',
    sources,
    include_directories: [
        openpiv_include, 
        pybind_include
    ],
    cpp_args: [
        '-D_USE_MATH_DEFINES'
    ],
    dependencies: [
        fmt_dep,
        thread_dep
    ],
    link_with: openpiv_lib,
    install: true
)
