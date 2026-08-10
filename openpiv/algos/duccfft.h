#pragma once

// std
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

// dynamically dispatched DUCC FFT shared-library API
#include <ducc_fft.h>

// local
#include "algos/fft_common.h"
#include "core/enum_helper.h"
#include "core/exception_builder.h"
#include "core/image.h"
#include "core/image_utils.h"
#include "core/pixel_types.h"
#include "core/util.h"

namespace openpiv::algos {

    using namespace core;

    /// Wrapper for the dynamically dispatched DUCC FFT shared library
    ///
    /// This class is thread-safe
    template <typename T>
    class DuccFFT
    {
        static_assert(
            std::is_same_v<T, float> || std::is_same_v<T, double>,
            "DuccFFT only supports float or double (image_gf32, image_gf64)"
        );

        using value_t = T;
        using std_complex_t = std::complex<value_t>;
        using image_cf_t = core::image< core::complex< value_t > >;  // complex image of T
        using image_gf_t = core::image< core::g< value_t > >;  // real image of T

        const size size_;
        const ducc_fft::shape_t shape_;
        const ducc_fft::shape_t axes_{ 1, 0 };

        /// storage for intermediate data
        struct data_t
        {
            image_cf_t output;
            std::vector<std_complex_t> fft_buffer_a;
            image_cf_t temp;
            std::vector<std_complex_t> fft_buffer_b;
            std::vector<value_t> real_buffer;
        };

        /// helpers to allow TLS for intermediate storage
        using storage_t = std::vector< std::tuple<DuccFFT*, data_t> >;
        storage_t& storage() const
        {
            thread_local static storage_t static_data;
            return static_data;
        }

        /// \fn cache contains a per-thread, per-instance copy of data
        /// that is lazily initialized; this allows a single instance
        /// of FFT to be called from multiple threads without locking
        data_t& cache() const
        {
            DuccFFT* self = const_cast<DuccFFT*>(this);
            for ( auto& [fft, data] : storage() )
            {
                if ( fft == self && data.output.size() == size_ )
                    return data;
            }

            data_t data;
            const std::size_t N = static_cast<std::size_t>( size_.area() );
            data.output.resize( size_ );
            data.temp.resize( size_ );
            data.fft_buffer_a.resize( N );
            data.fft_buffer_b.resize( N );
            data.real_buffer.resize( N );
            auto& [fft, result] = storage().emplace_back(self, std::move(data));

            return result;
        }

        std::size_t spectrum_size() const
        {
            return static_cast<std::size_t>( size_.width() )
                 * (static_cast<std::size_t>( size_.height() ) / 2 + 1);
        }

        static void copy_complex_to_image(
            const std::vector<std_complex_t>& input,
            image_cf_t& output,
            std::size_t count )
        {
            for ( std::size_t i = 0; i < count; ++i )
                output[i] = core::complex<value_t>{ input[i].real(), input[i].imag() };
        }

    public:
        DuccFFT( const core::size& size )
            : size_(size)
            , shape_{
                static_cast<std::size_t>( size.height() ),
                static_cast<std::size_t>( size.width() )
              }
        { }

        /// Perform a 2-D FFT; will always produce a complex floating point image output
        template < template <typename> class ImageT,
                   typename ContainedT,
                   typename = typename std::enable_if_t< is_imagetype_v<ImageT<ContainedT>> >
                   >
        const image_cf_t& transform( const ImageT<ContainedT>& input, direction d = direction::FORWARD ) const
        {
            DECLARE_ENTRY_EXIT
            if ( input.size() != size_ )
            {
                exception_builder< std::runtime_error >()
                    << "image size is different from expected: " << input.size() << ", " << size_;
            }

            // copy data, converting to complex
            cache().temp = input;
            cache().output.resize( input.size() );

            auto& data = cache();
            const std::size_t N = static_cast<std::size_t>( size_.area() );
            for ( std::size_t i = 0; i < N; ++i )
            {
                data.fft_buffer_a[i] = std_complex_t{
                    data.temp[i].real,
                    data.temp[i].imag
                };
            }

            ducc_fft::c2c(
                data.fft_buffer_a.data(),
                data.fft_buffer_b.data(),
                shape_,
                axes_,
                d == direction::FORWARD,
                value_t{1} );

            copy_complex_to_image( data.fft_buffer_b, data.output, N );
            return data.output;
        }

        /// Perform a 2-D FFT of two real images; will produce two
        /// output images
        template < template <typename> class ImageT,
                   typename ContainedT,
                   typename = typename std::enable_if_t<
                       is_imagetype_v<ImageT<ContainedT>> &&
                       is_real_mono_pixeltype_v<ContainedT>
                       >
                   >
        std::tuple<image_cf_t&, image_cf_t&>
        transform_real( const ImageT<ContainedT>& a,
                        const ImageT<ContainedT>& b,
                        direction d = direction::FORWARD ) const
        {
            DECLARE_ENTRY_EXIT
            if ( a.size() != size_ || b.size() != size_ )
            {
                exception_builder< std::runtime_error >()
                    << "image size is different from expected: " << a.size()
                    << ", " << size_;
            }

            cache().output.resize( a.size() );
            cache().temp.resize( b.size() );

            auto& data = cache();
            const std::size_t N = static_cast<std::size_t>( size_.area() );
            const std::size_t spectrum_N = spectrum_size();

            for ( std::size_t i = 0; i < N; ++i )
                data.real_buffer[i] = static_cast<value_t>( a[i].v );

            ducc_fft::r2c(
                data.real_buffer.data(),
                data.fft_buffer_a.data(),
                shape_,
                axes_,
                d == direction::FORWARD,
                value_t{1} );

            for ( std::size_t i = 0; i < N; ++i )
                data.real_buffer[i] = static_cast<value_t>( b[i].v );

            ducc_fft::r2c(
                data.real_buffer.data(),
                data.fft_buffer_b.data(),
                shape_,
                axes_,
                d == direction::FORWARD,
                value_t{1} );

            // DUCC reads and writes only the leading Hermitian-spectrum prefix.
            // Both backing vectors and both OpenPIV images remain full-sized,
            // matching the existing PocketFFT wrapper's storage behavior.
            copy_complex_to_image( data.fft_buffer_a, data.output, spectrum_N );
            copy_complex_to_image( data.fft_buffer_b, data.temp, spectrum_N );

            return { data.output, data.temp };
        }

        template < template <typename> class ImageT,
                   typename ContainedT,
                   typename ValueT = typename ContainedT::value_t,
                   typename OutT = image<g<ValueT>>,
                   typename = typename std::enable_if_t<
                       is_imagetype_v<ImageT<ContainedT>> &&
                       is_complex_mono_pixeltype_v<ContainedT>
                       >
                   >
        OutT
        transform_real( const ImageT<ContainedT>& in,
                        direction d = direction::FORWARD ) const
        {
            DECLARE_ENTRY_EXIT
            if ( in.size() != size_ )
            {
                exception_builder< std::runtime_error >()
                    << "image size is different from expected: " << in.size()
                    << ", " << size_;
            }

            auto& data = cache();
            const std::size_t spectrum_N = spectrum_size();
            for ( std::size_t i = 0; i < spectrum_N; ++i )
            {
                data.fft_buffer_a[i] = std_complex_t{
                    static_cast<value_t>( in[i].real ),
                    static_cast<value_t>( in[i].imag )
                };
            }

            ducc_fft::c2r(
                data.fft_buffer_a.data(),
                data.real_buffer.data(),
                shape_,
                axes_,
                d == direction::FORWARD,
                value_t{1} );

            OutT out{ in.size() };
            const std::size_t N = static_cast<std::size_t>( size_.area() );
            for ( std::size_t i = 0; i < N; ++i )
                out[i] = static_cast<ValueT>( data.real_buffer[i] );

            return out;
        }

        template < template <typename> class ImageT,
                   typename ContainedT,
                   typename ValueT = typename ContainedT::value_t,
                   typename OutT = image<g<ValueT>>,
                   typename = typename std::enable_if_t< is_imagetype_v<ImageT<ContainedT>> >
                   >
        OutT
        cross_correlate( const ImageT<ContainedT>& a,
                         const ImageT<ContainedT>& b ) const
        {
            image_cf_t a_fft{ transform( a, direction::FORWARD ) };
            image_cf_t b_fft{ transform( b, direction::FORWARD ) };

            a_fft = b_fft * conj( a_fft );
            OutT output{ real( transform( a_fft, direction::REVERSE ) ) };
            swap_quadrants( output );

            return output;
        }

        template < template <typename> class ImageT,
                   typename ContainedT,
                   typename ValueT = typename ContainedT::value_t,
                   typename OutT = image<g<ValueT>>,
                   typename = typename std::enable_if_t<
                       is_imagetype_v<ImageT<ContainedT>> &&
                       is_real_mono_pixeltype_v<ContainedT>
                       >
                   >
        OutT
        cross_correlate_real( const ImageT<ContainedT>& a,
                              const ImageT<ContainedT>& b ) const
        {
            auto [a_fft, b_fft] = transform_real( a, b, direction::FORWARD );
            a_fft = b_fft * conj( a_fft );
            OutT output = transform_real( a_fft, direction::REVERSE );
            swap_quadrants( output );

            return output;
        }

        template < template <typename> class ImageT,
                   typename ContainedT,
                   typename = typename std::enable_if_t< is_imagetype_v<ImageT<ContainedT>> >
                   >
        const image_cf_t& auto_correlate( const ImageT<ContainedT>& a ) const
        {
            image_cf_t a_fft{ transform( a, direction::FORWARD ) };

            a_fft = abs_sqr( a_fft );
            cache().output = real( transform( a_fft, direction::REVERSE ) );
            swap_quadrants( cache().output );

            return cache().output;
        }
    };

}
