// This is intentionally a Highway toggle guard, not an ordinary include guard.
#if defined(DUCC_FFT_TARGET_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef DUCC_FFT_TARGET_INL_H_
#undef DUCC_FFT_TARGET_INL_H_
#else
#define DUCC_FFT_TARGET_INL_H_
#endif

#include "hwy/highway.h"

// Reparse the complete SIMD-dependent DUCC template stack.  In particular,
// Cmplx must be compiled under the same target attribute as vtp; otherwise a
// butterfly crosses the Highway target boundary for every complex operation
// and the compiler cannot inline the SIMD operators.
#undef DUCC0_CMPLX_H
#undef DUCC0_UNITY_ROOTS_H
#undef DUCC0_FFT_H
#undef DUCC0_FFT1D_IMPL_H
#undef DUCC0_FFTND_IMPL_H

HWY_BEFORE_NAMESPACE();

namespace ducc_fft_detail
{
namespace HWY_NAMESPACE
{
    // DUCC opens namespace ducc0 itself. Seed each target's namespace with the
    // target-independent MAV/threading/error infrastructure from ::ducc0. The
    // include below then creates target-local Cmplx, roots, plans and kernels.
    namespace ducc0
    {
        using namespace ::ducc0;
    }  // namespace ducc0

    #include "ducc0/fft/fftnd_impl.h"

    namespace {

        template <class T>
        void C2CImpl(
            const std::complex<T>* input, std::complex<T>* output,
            const std::vector<std::size_t>& shape,
            const std::vector<std::size_t>& axes, 
            bool forward, 
            T scale
        ) {
            const ducc0::cfmav<std::complex<T>> in(input, shape);
            const ducc0::vfmav<std::complex<T>> out(output, shape);
            ducc0::detail_fft::c2c(in, out, axes, forward, scale, 1);
        }

        template <class T>
        void R2CImpl(
            const T* input, std::complex<T>* output,
            const std::vector<std::size_t>& shape,
            const std::vector<std::size_t>& axes,
            bool forward, 
            T scale
        ) {
            std::vector<std::size_t> complex_shape(shape);
            complex_shape[axes.back()] = shape[axes.back()] / 2 + 1;
            const ducc0::cfmav<T> in(input, shape);
            const ducc0::vfmav<std::complex<T>> out(output, complex_shape);
            ducc0::detail_fft::r2c(in, out, axes, forward, scale, 1);
        }

        template <class T>
        void C2RImpl(
            const std::complex<T>* input, T* output,
            const std::vector<std::size_t>& real_shape,
            const std::vector<std::size_t>& axes,
            bool forward,
            T scale
        ) {
            std::vector<std::size_t> complex_shape(real_shape);
            complex_shape[axes.back()] = real_shape[axes.back()] / 2 + 1;
            const ducc0::cfmav<std::complex<T>> in(input, complex_shape);
            const ducc0::vfmav<T> out(output, real_shape);
            ducc0::detail_fft::c2r(in, out, axes, forward, scale, 1);
        }

    }  // namespace

    const char* BackendName() { return hwy::TargetName(HWY_TARGET); }

    void C2CF32(
        const std::complex<float>* input,
        std::complex<float>* output,
        const std::vector<std::size_t>& shape,
        const std::vector<std::size_t>& axes,
        bool forward,
        float scale
    ) {
        C2CImpl(input, output, shape, axes, forward, scale);
    }

    void C2CF64(
        const std::complex<double>* input,
        std::complex<double>* output,
        const std::vector<std::size_t>& shape,
        const std::vector<std::size_t>& axes,
        bool forward,
        double scale
    ) {
        C2CImpl(input, output, shape, axes, forward, scale);
    }

    void R2CF32(
        const float* input, 
        std::complex<float>* output,
        const std::vector<std::size_t>& shape,
        const std::vector<std::size_t>& axes, 
        bool forward, 
        float scale
    ) {
        R2CImpl(input, output, shape, axes, forward, scale);
    }

    void R2CF64(
        const double* input, 
        std::complex<double>* output,
        const std::vector<std::size_t>& shape,
        const std::vector<std::size_t>& axes, 
        bool forward, 
        double scale
    ) {
        R2CImpl(input, output, shape, axes, forward, scale);
    }

    void C2RF32(
        const std::complex<float>* input, 
        float* output,
        const std::vector<std::size_t>& real_shape,
        const std::vector<std::size_t>& axes, 
        bool forward, 
        float scale
    ) {
        C2RImpl(input, output, real_shape, axes, forward, scale);
    }

    void C2RF64(
        const std::complex<double>* input, 
        double* output,
        const std::vector<std::size_t>& real_shape,
        const std::vector<std::size_t>& axes, 
        bool forward, 
        double scale
    ) {
        C2RImpl(input, output, real_shape, axes, forward, scale);
    }

}  // namespace HWY_NAMESPACE
}  // namespace ducc_fft_detail

HWY_AFTER_NAMESPACE();

#endif  // toggle guard
