openpiv_src = [
    'core/size.cpp',
    'core/rect.cpp',
    'core/util.cpp',
    'core/detail/grid.impl.cpp',
    'loaders/image_loader.cpp',
    'filters/detail/filter_common.impl.cpp',
    'interp/detail/interp_common.impl.cpp',
    'interp/detail/map_polynomial.impl.cpp',
    'interp/detail/map_sinc.impl.cpp',
    'interp/detail/map_lanczos.impl.cpp',
    'piv/detail/deformation.impl.cpp',
    'piv/detail/firstpass.impl.cpp',
    'piv/detail/multipass.impl.cpp',
]

if not disable_image_loaders
  openpiv_src += [
    'loaders/pnm_image_loader.cpp',
    'loaders/tiff_image_loader.cpp',
    libtiff_cpp_include
  ]
endif

cpp_args = [
  '-D_USE_MATH_DEFINES'
]

if cxx_compiler == 'gcc'
  cpp_args += cpp.get_supported_arguments(
    '-export-symbols'
  )
elif cxx_compiler == 'msvc'
  cpp_args += ['-DEXPORT_DLL_SYMBOLS']
endif

openpiv_lib = shared_library(
  'openpivcore',
  openpiv_src,
  include_directories: [openpiv_include, pocketfft_include, threadlib_include],
  dependencies: libopenpivcore_deps,
  cpp_args: cpp_args,
  install: true
)
