duccfft_cpp_syntax = cpp.get_argument_syntax()
duccfft_define_prefix = duccfft_cpp_syntax == 'msvc' ? '/D' : '-D'

duccfft_cpp_args = [
  duccfft_define_prefix + 'DUCC_FFT_BUILDING_LIBRARY=1',
  duccfft_define_prefix + 'DUCC0_NO_FFT_CACHE=1',
  duccfft_define_prefix + 'DUCC0_NO_LOWLEVEL_THREADING=1',
  duccfft_define_prefix + 'HWY_HEADER_ONLY=1',
  duccfft_define_prefix + 'HWY_BROKEN_EMU128=1',
  duccfft_define_prefix + 'HWY_COMPILE_ALL_ATTAINABLE=1',
]

if duccfft_cpp_syntax == 'msvc'
  duccfft_cpp_args += cpp.get_supported_arguments([
    '/bigobj',
    '/EHsc',
  ])
else
  duccfft_cpp_args += cpp.get_supported_arguments([
    '-fno-strict-aliasing',
  ])
endif

duccfft_public_inc = include_directories('.')
duccfft_ducc_inc = include_directories('../external/ducc/src')
duccfft_highway_inc = include_directories('../external/highway')

duccfft_lib = shared_library(
  'duccfft',
  files(
    'ducc_fft.cc',
    '../external/ducc/src/ducc0/infra/threading.cc',
    '../external/highway/hwy/targets.cc',
  ),
  name_prefix: 'lib',
  include_directories: [
    duccfft_public_inc,
    duccfft_ducc_inc,
    duccfft_highway_inc,
  ],
  cpp_args: duccfft_cpp_args,
  gnu_symbol_visibility: 'hidden',
  install: true,
)

duccfft_dep = declare_dependency(
  link_with: duccfft_lib,
  include_directories: duccfft_public_inc,
)

install_headers('ducc_fft.h')

message('libduccfft: compiling every target Highway considers attainable')
message('libduccfft: Highway uses runtime dispatch where supported and the best static or scalar fallback elsewhere')
