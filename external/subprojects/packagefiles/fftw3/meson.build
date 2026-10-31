# Native Meson build for the unmodified FFTW 3.3.11 release sources.
# Only single precision is built. OpenPIV supplies the worker threads.
project(
  'fftw3',
  'c',
  version: '3.3.11',
  license: 'GPL-2.0-or-later',
  meson_version: '>= 1.1.0',
  default_options: ['warning_level=0'],
)

cc = meson.get_compiler('c')
fs = import('fs')
msvc_syntax = cc.get_argument_syntax() == 'msvc'
define_prefix = msvc_syntax ? '/D' : '-D'
math_dep = cc.find_library('m', required: false)
fftw_deps = math_dep.found() ? [math_dep] : []
config = configuration_data()
if math_dep.found()
  config.set('HAVE_LIBM', 1)
endif
config.set('FFTW_SINGLE', 1)
config.set('DISABLE_FORTRAN', 1)
config.set('STDC_HEADERS', 1)
config.set_quoted('PACKAGE', 'fftw')
config.set_quoted('PACKAGE_NAME', 'fftw')
config.set_quoted('PACKAGE_VERSION', meson.project_version())
config.set_quoted('PACKAGE_STRING', 'fftw ' + meson.project_version())
config.set_quoted('VERSION', meson.project_version())
config.set_quoted('FFTW_CC', cc.get_id() + ' ' + cc.version())

foreach header : [
  'alloca.h', 'c_asm.h', 'dlfcn.h', 'fcntl.h', 'fenv.h', 'intrinsics.h',
  'inttypes.h', 'limits.h', 'malloc.h', 'mach/mach_time.h', 'stddef.h',
  'stdint.h', 'stdio.h', 'stdlib.h', 'string.h', 'strings.h',
  'sys/stat.h', 'sys/time.h', 'sys/types.h', 'unistd.h',
]
  if cc.has_header(header)
    config.set('HAVE_' + header.underscorify().to_upper(), 1)
  endif
endforeach

if cc.has_header('sys/time.h') and cc.has_header('time.h')
  config.set('TIME_WITH_SYS_TIME', 1)
endif

foreach probe : [
  ['abort', '#include <stdlib.h>'],
  ['clock_gettime', '#include <time.h>'],
  ['cosl', '#include <math.h>'],
  ['drand48', '#include <stdlib.h>'],
  ['getpagesize', '#include <unistd.h>'],
  ['gettimeofday', '#include <sys/time.h>'],
  ['isnan', '#include <math.h>'],
  ['mach_absolute_time', '#include <mach/mach_time.h>'],
  ['memalign', '#include <malloc.h>'],
  ['memmove', '#include <string.h>'],
  ['memset', '#include <string.h>'],
  ['posix_memalign', '#include <stdlib.h>'],
  ['sinl', '#include <math.h>'],
  ['snprintf', '#include <stdio.h>'],
  ['sqrt', '#include <math.h>'],
  ['strchr', '#include <string.h>'],
  ['vprintf', '#include <stdio.h>'],
]
  if cc.has_function(probe[0], prefix: probe[1], dependencies: fftw_deps)
    config.set('HAVE_' + probe[0].to_upper(), 1)
  endif
endforeach

foreach probe : [
  ['sinl', '#include <math.h>'],
  ['cosl', '#include <math.h>'],
  ['memalign', '#include <malloc.h>'],
  ['drand48', '#include <stdlib.h>'],
  ['srand48', '#include <stdlib.h>'],
  ['posix_memalign', '#include <stdlib.h>'],
]
  declared = cc.compiles(
    probe[1] + '\nvoid f(void) { (void) ' + probe[0] + '; }',
    name: 'declaration of ' + probe[0],
  )
  config.set('HAVE_DECL_' + probe[0].to_upper(), declared ? 1 : 0)
endforeach

foreach type : [
  'int', 'unsigned int', 'long', 'unsigned long', 'long long',
  'unsigned long long', 'size_t', 'ptrdiff_t', 'float', 'double',
]
  config.set('SIZEOF_' + type.underscorify().to_upper(),
    cc.sizeof(type, prefix: '#include <stddef.h>'))
endforeach
config.set('SIZEOF_VOID_P', cc.sizeof('void *'))
config.set('SIZEOF_FFTW_R2R_KIND', cc.sizeof('int'))
if cc.has_type('uintptr_t', prefix: '#include <stdint.h>')
  config.set('HAVE_UINTPTR_T', 1)
endif
if cc.has_type('ptrdiff_t', prefix: '#include <stddef.h>')
  config.set('HAVE_PTRDIFF_T', 1)
endif

alloca_probe = '''
#include <stddef.h>
#if defined(__GNUC__)
# define alloca __builtin_alloca
#elif defined(_MSC_VER)
# include <malloc.h>
# define alloca _alloca
#else
# include <alloca.h>
#endif
int main(void) { void *p = alloca(16); return p == 0; }
'''
if cc.links(alloca_probe, name: 'alloca')
  config.set('HAVE_ALLOCA', 1)
  config.set('FFTW_ENABLE_ALLOCA', 1)
endif

# FFTW provides its own aligned allocator for MinGW and other Windows C
# toolchains. MSVC uses its existing _aligned_malloc/_aligned_free path.
if host_machine.system() == 'windows' and not msvc_syntax
  config.set('WITH_OUR_MALLOC', 1)
endif

# Compile only codelets with extended ISA flags; dispatch stays at baseline.
isas = []
if get_option('simd')
  cpu = host_machine.cpu_family()
  if cpu == 'x86' or cpu == 'x86_64'
    sse2_args = ['-msse2']
    if msvc_syntax
      sse2_args = cpu == 'x86' ? ['/arch:SSE2'] : []
    endif
    avx_args = msvc_syntax ? ['/arch:AVX'] : ['-mavx']
    avx2_args = msvc_syntax ? ['/arch:AVX2'] : ['-mavx2', '-mfma']
    isas = [
      ['sse2', 'HAVE_SSE2', sse2_args, ['sse2'],
       '#include <emmintrin.h>\nint main(void) { __m128 a = _mm_set1_ps(1.0f); return (int)_mm_cvtss_f32(a); }'],
      ['avx', 'HAVE_AVX', avx_args, ['avx'],
       '#include <immintrin.h>\nint main(void) { __m256 a = _mm256_set1_ps(1.0f); return (int)_mm256_cvtss_f32(a); }'],
      ['avx2', 'HAVE_AVX2', avx2_args, ['avx2', 'avx2-128'],
       '#include <immintrin.h>\nint main(void) { __m256 a = _mm256_set1_ps(1.0f); a = _mm256_fmadd_ps(a, a, a); return (int)_mm256_cvtss_f32(a); }'],
    ]
  elif cpu == 'aarch64'
    isas = [
      ['neon', 'HAVE_NEON', [], ['neon'],
       '#include <arm_neon.h>\nint main(void) { float32x4_t a = vdupq_n_f32(1.0f); return (int)vgetq_lane_f32(a, 0); }'],
    ]
  endif
endif

enabled_isas = []
isa_names = []
foreach isa : isas
  if cc.compiles(isa[4], args: isa[2], name: 'FFTW ' + isa[0] + ' codelets')
    config.set(isa[1], 1)
    enabled_isas += [isa]
    isa_names += isa[0]
  endif
endforeach

config_h = configure_file(output: 'config.h', configuration: config)
fftw_inc = include_directories('.')
fftw_api_inc = include_directories('api')
fftw_args = cc.get_supported_arguments(['-fomit-frame-pointer', '-fstrict-aliasing'])
public_args = []
if host_machine.system() == 'windows' and get_option('default_library') != 'static'
  fftw_args += define_prefix + 'FFTW_DLL=1'
  public_args += define_prefix + 'FFTW_DLL=1'
endif

core_sources = fs.read('sources.txt').strip().split('\n')
dft_codelets = fs.read('dft-codelets.txt').strip().split('\n')
rdft_codelets = fs.read('rdft-codelets.txt').strip().split('\n')
simd_objects = []
foreach isa : enabled_isas
  isa_sources = []
  foreach directory : isa[3]
    foreach source : dft_codelets
      isa_sources += 'dft/simd/' + directory + '/' + source
    endforeach
    foreach source : rdft_codelets
      isa_sources += 'rdft/simd/' + directory + '/' + source
    endforeach
  endforeach
  codelets = static_library(
    'fftw3f_' + isa[0] + '_codelets',
    isa_sources,
    config_h,
    include_directories: fftw_inc,
    c_args: fftw_args + isa[2],
    override_options: ['b_lto=false'],
    pic: true,
    install: false,
  )
  simd_objects += codelets.extract_all_objects(recursive: false)
endforeach

fftw_lib = library(
  'fftw3f',
  core_sources,
  config_h,
  objects: simd_objects,
  include_directories: fftw_inc,
  dependencies: fftw_deps,
  c_args: fftw_args,
  override_options: ['b_lto=false'],
  pic: true,
  version: host_machine.system() == 'linux' ? '3' : '3.7.11',
  soversion: '3',
  install: true,
)

fftw_dep = declare_dependency(
  link_with: fftw_lib,
  include_directories: fftw_api_inc,
  dependencies: fftw_deps,
  compile_args: public_args,
  version: meson.project_version(),
)
meson.override_dependency('fftw3f', fftw_dep)
install_headers('api/fftw3.h')
install_data('COPYING', install_dir: get_option('datadir') / 'licenses' / 'fftw3f')
summary({'Precision': 'float', 'SIMD codelets': isa_names, 'Internal threads': false}, section: 'FFTW')
