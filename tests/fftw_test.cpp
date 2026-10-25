#include "algos/fft.h"
#include "algos/fftw.h"
#include "core/image_view.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>

namespace
{
    using namespace openpiv;
    using real_image = core::image<core::g_f32>;
    using complex_image = core::image<core::c_f32>;
    using algos::direction;

    void require( bool condition, const std::string& message )
    {
        if ( !condition )
            throw std::runtime_error(message);
    }

    template <typename ImageT>
    void fill_real( ImageT& image, float phase )
    {
        for ( std::size_t i = 0; i < image.pixel_count(); ++i )
        {
            const float x = static_cast<float>(i % image.width());
            const float y = static_cast<float>(i / image.width());
            image[i] = std::sin(0.17f * x + phase) * std::cos(0.23f * y)
                     + 0.27f * std::sin(0.071f * x * y + phase);
        }
    }

    template <typename ImageT>
    void fill_complex( ImageT& image )
    {
        for ( std::size_t i = 0; i < image.pixel_count(); ++i )
        {
            const float x = static_cast<float>(i % image.width());
            const float y = static_cast<float>(i / image.width());
            image[i] = core::c_f32{
                std::sin(0.19f * x + 0.31f * y),
                std::cos(0.37f * x - 0.11f * y)
            };
        }
    }

    void close( const real_image& actual, const real_image& expected,
                const std::string& label )
    {
        require(actual.size() == expected.size(), label + ": shape");
        double maximum_error = 0.0;
        double maximum_value = 1.0;
        for ( std::size_t i = 0; i < actual.pixel_count(); ++i )
        {
            const double a = actual[i].v;
            const double b = expected[i].v;
            require(std::isfinite(a) && std::isfinite(b), label + ": nonfinite");
            maximum_error = std::max(maximum_error, std::abs(a - b));
            maximum_value = std::max(maximum_value, std::abs(b));
        }
        require(maximum_error <= 5e-5 * maximum_value,
                label + ": maximum error " + std::to_string(maximum_error));
    }

    void close( const complex_image& actual, const complex_image& expected,
                const std::string& label, std::size_t count = 0 )
    {
        require(actual.size() == expected.size(), label + ": shape");
        if ( count == 0 )
            count = actual.pixel_count();
        double maximum_error = 0.0;
        double maximum_value = 1.0;
        for ( std::size_t i = 0; i < count; ++i )
        {
            const auto a = actual[i];
            const auto b = expected[i];
            require(std::isfinite(a.real) && std::isfinite(a.imag) &&
                    std::isfinite(b.real) && std::isfinite(b.imag),
                    label + ": nonfinite");
            maximum_error = std::max(maximum_error,
                std::hypot(double(a.real) - b.real, double(a.imag) - b.imag));
            maximum_value = std::max(maximum_value,
                std::hypot(double(b.real), double(b.imag)));
        }
        require(maximum_error <= 5e-5 * maximum_value,
                label + ": maximum error " + std::to_string(maximum_error));
    }

    void check_transforms()
    {
        const core::size shape{ 32, 32 };
        algos::FFTW<float> fftw{ shape, true };
        algos::FFT<double> reference{ shape };
        real_image a{ shape }, b{ shape };
        complex_image input{ shape };
        fill_real(a, 0.2f);
        fill_real(b, 0.8f);
        fill_complex(input);

        // Copy results because each FFT object reuses its internal output.
        const complex_image actual{ fftw.transform(input, direction::FORWARD) };
        const complex_image expected{ reference.transform(input, direction::FORWARD) };
        close(actual, expected, "forward transform");

        const complex_image inverse{ fftw.transform(actual, direction::REVERSE) };
        const complex_image expected_inverse{ reference.transform(expected, direction::REVERSE) };
        close(inverse, expected_inverse, "inverse transform");

        const complex_image actual_a{ fftw.transform(a) };
        const complex_image expected_a{ reference.transform(a) };
        close(actual_a, expected_a, "transform of real input");

        auto pair = fftw.transform_real(a, b, direction::FORWARD);
        const complex_image spectrum_a{ std::get<0>(pair) };
        const complex_image spectrum_b{ std::get<1>(pair) };
        const complex_image expected_b{ reference.transform(b) };

        // FFTW stores only the nonredundant part of each real spectrum.
        const std::size_t spectrum_count = a.width() * (a.height() / 2 + 1);
        close(spectrum_a, expected_a, "two real transforms: A", spectrum_count);
        close(spectrum_b, expected_b, "two real transforms: B", spectrum_count);

        const real_image restored = fftw.transform_real(spectrum_a, direction::REVERSE);
        const real_image expected_restored{ core::real(reference.transform(expected_a, direction::REVERSE)) };
        close(restored, expected_restored, "inverse real transform");
    }

    void check_correlations()
    {
        real_image image{ core::size{64, 64} };
        fill_real(image, 0.2f);
        const core::size shape{ 32, 32 };
        const auto a = core::create_image_view(image, core::rect{ {20, 20}, shape });
        const auto b = core::create_image_view(image, core::rect{ {20, 25}, shape });
        algos::FFTW<float> fftw{ shape, false };
        algos::FFTW<float> fftw_complex{ shape, true };
        algos::FFT<double> reference{ shape };

        const real_image expected = reference.cross_correlate(a, b);
        close(fftw_complex.cross_correlate(a, b), expected, "cross correlation of views");
        close(fftw.cross_correlate_real(a, b), expected, "real cross correlation of views");

        bool threw = false;
        try { (void)fftw.cross_correlate(a, b); }
        catch ( const std::runtime_error& error )
        {
            threw = std::string(error.what()).find("complex plans have not been generated") != std::string::npos;
        }
        require(threw, "complex correlation without C2C plans must throw");

        real_image output{ shape };
        fftw_complex.cross_correlate_inplace(a, b, output);
        close(output, expected, "in-place cross correlation");

        fftw.cross_correlate_real_inplace(a, b, output);
        close(output, expected, "in-place real cross correlation");

        const complex_image actual_auto{ fftw_complex.auto_correlate(a) };
        const complex_image expected_auto{ reference.auto_correlate(a) };
        close(actual_auto, expected_auto, "autocorrelation");

        // PIV workers share a correlator and supply their own output images.
        real_image worker_output{ shape };
        std::thread worker([&] { fftw.cross_correlate_real_inplace(a, b, worker_output); });
        fftw.cross_correlate_real_inplace(a, b, output);
        worker.join();
        close(worker_output, expected, "shared correlator: worker output");
        close(output, expected, "shared correlator: main output");
    }

    void check_input_size()
    {
        algos::FFTW<float> fftw{ core::size{32, 32}, true };
        real_image wrong{ core::size{16, 16} };
        bool threw = false;
        try { (void)fftw.transform(wrong); }
        catch ( const std::runtime_error& ) { threw = true; }
        require(threw, "different input size must throw");
    }
}

int main()
{
    try
    {
        check_transforms();
        check_correlations();
        check_input_size();
        std::cout << "FFTW transform and correlation checks passed.\n" << std::flush;
    }
    catch ( const std::exception& error )
    {
        std::cerr << "FFTW class check failed: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
