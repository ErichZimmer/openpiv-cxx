openpiv_src = [
    'core/size.cpp',
    'core/rect.cpp',
    'core/util.cpp',
    'core/detail/grid.impl.cpp',
    'loaders/image_loader.cpp',
    'interp/detail/map_polynomial.impl.cpp',
    'interp/detail/map_sinc.impl.cpp',
    'interp/detail/map_lanczos.impl.cpp',
    'algos/duccfft_backend.cpp',
    'piv/detail/deformation.impl.cpp',
    'piv/detail/firstpass.impl.cpp',
    'piv/detail/multipass.impl.cpp',
    'piv/detail/extract_with_padding.impl.cpp'
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

if get_option('export_dll_symbols')
  cpp_args += ['-DEXPORT_DLL_SYMBOLS=1']
endif

openpiv_lib = shared_library(
  'openpivcore',
  openpiv_src,
  include_directories: [openpiv_include, pocketfft_include, threadlib_include],
  dependencies: libopenpivcore_deps,
  cpp_args: cpp_args,
  install: true,
#  install_dir: py.get_install_dir() / 'openpiv_cxx'

)
