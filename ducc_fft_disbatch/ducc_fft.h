#ifndef DUCC_FFT_PUBLIC_H_
#define DUCC_FFT_PUBLIC_H_

#include <complex>
#include <cstddef>
#include <vector>

#if defined(_WIN32) && defined(DUCC_FFT_BUILDING_LIBRARY)
#define DUCC_FFT_API __declspec(dllexport)
#elif defined(__GNUC__) && defined(DUCC_FFT_BUILDING_LIBRARY)
#define DUCC_FFT_API __attribute__((visibility("default")))
#else
#define DUCC_FFT_API
#endif

namespace ducc_fft {

    using shape_t = std::vector<std::size_t>;

    // Name of the Highway target selected for the current CPU.
    DUCC_FFT_API const char* backend_name() noexcept;

    // All arrays are contiguous and row-major. Every transform is serial and
    // constructs its DUCC plan for the call; this library has no plan cache.
    DUCC_FFT_API void c2c(
        const std::complex<float>* input,
        std::complex<float>* output, 
        const shape_t& shape,
        const shape_t& axes, 
        bool forward, 
        float scale = 1.0f
    );

    DUCC_FFT_API void c2c(
        const std::complex<double>* input,
        std::complex<double>* output, 
        const shape_t& shape,
        const shape_t& axes, 
        bool forward, 
        double scale = 1.0
    );

    DUCC_FFT_API void r2c(
        const float* input, 
        std::complex<float>* output,
        const shape_t& shape, 
        const shape_t& axes,
        bool forward = true, 
        float scale = 1.0f
    );

    DUCC_FFT_API void r2c(
        const double* input, 
        std::complex<double>* output,
        const shape_t& shape, 
        const shape_t& axes,
        bool forward = true, 
        double scale = 1.0
    );

    DUCC_FFT_API void c2r(
        const std::complex<float>* input, 
        float* output,
        const shape_t& real_shape, 
        const shape_t& axes,
        bool forward = false, 
        float scale = 1.0f
    );

    DUCC_FFT_API void c2r(
        const std::complex<double>* input, 
        double* output,
        const shape_t& real_shape, 
        const shape_t& axes,
        bool forward = false, 
        double scale = 1.0
    );

}  // namespace ducc_fft

#endif  // DUCC_FFT_PUBLIC_H_
