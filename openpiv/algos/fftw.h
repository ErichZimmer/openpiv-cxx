#pragma once

// std
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <new>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

// FFTW single-precision API
#include <fftw3.h>

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

    /// Wrapper for FFTW single-precision transforms
    ///
    /// This class is thread-safe
    template <typename T>
    class FFTW
    {
        static_assert(
            std::is_same_v<T, float>,
            "FFTW only supports float (image_gf32)"
        );

        using value_t = T;
        using image_cf_t = core::image< core::complex< value_t > >;
        using image_gf_t = core::image< core::g< value_t > >;
        using real_buffer_t = std::unique_ptr<float[], decltype(&fftwf_free)>;
        using complex_buffer_t = std::unique_ptr<fftwf_complex[], decltype(&fftwf_free)>;

        static std::mutex& planner_mutex()
        {
            static std::mutex mutex;
            return mutex;
        }

        static void destroy_plan( fftwf_plan plan ) noexcept
        {
            std::lock_guard<std::mutex> lock( planner_mutex() );
            fftwf_destroy_plan( plan );
        }

        using plan_t = std::unique_ptr<
            std::remove_pointer_t<fftwf_plan>, decltype(&destroy_plan)
        >;

        const size size_;
        plan_t c2c_forward_plan_{ nullptr, &destroy_plan };
        plan_t c2c_reverse_plan_{ nullptr, &destroy_plan };
        plan_t r2c_plan_{ nullptr, &destroy_plan };
        plan_t c2r_plan_{ nullptr, &destroy_plan };

        /// storage for intermediate data
        struct data_t
        {
            image_cf_t output;
            complex_buffer_t fft_buffer_a{ nullptr, &fftwf_free };
            image_cf_t temp;
            complex_buffer_t fft_buffer_b{ nullptr, &fftwf_free };
            real_buffer_t real_buffer{ nullptr, &fftwf_free };
        };

        /// helpers to allow TLS for intermediate storage
        using storage_t = std::vector< std::tuple<FFTW*, data_t> >;
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
            FFTW* self = const_cast<FFTW*>(this);
            for ( auto& [fft, data] : storage() )
            {
                if ( fft == self && data.output.size() == size_ )
                    return data;
            }

            data_t data;
            const std::size_t N = static_cast<std::size_t>( size_.area() );
            data.output.resize( size_ );
            data.temp.resize( size_ );
            data.fft_buffer_a.reset( fftwf_alloc_complex(N) );
            data.fft_buffer_b.reset( fftwf_alloc_complex(N) );
            data.real_buffer.reset( fftwf_alloc_real(N) );

            if ( !data.fft_buffer_a || !data.fft_buffer_b || !data.real_buffer )
                exception_builder< std::bad_alloc >()
                    << "FFTW buffers failed to allocate";

            auto& [fft, result] = storage().emplace_back(self, std::move(data));

            return result;
        }

        static core::size validate_size( const core::size& size )
        {
            const uint64_t width = size.width();
            const uint64_t height = size.height();
            const uint64_t N = width * height;

            if ( width == 0 ||
                 height == 0 ||
                 width > 100000 ||
                 height > 100000 ||
                 N > std::numeric_limits<uint32_t>::max() ||
                 N > std::numeric_limits<std::size_t>::max() / sizeof(fftwf_complex) )
            {
                exception_builder< std::invalid_argument >()
                    << "invalid FFTW image size: " << size;
            }

            return size;
        }

        std::size_t spectrum_size() const
        {
            return static_cast<std::size_t>( size_.width() )
                 * (static_cast<std::size_t>( size_.height() ) / 2 + 1);
        }

        static void copy_complex_to_image(
            const complex_buffer_t& input,
            image_cf_t& output,
            std::size_t count )
        {
            for ( std::size_t i = 0; i < count; ++i )
                output[i] = core::complex<value_t>{ input[i][0], input[i][1] };
        }

        template < template <typename> class ImageT, typename ContainedT >
        static void copy_image_to_complex(
            const ImageT<ContainedT>& input,
            complex_buffer_t& output )
        {
            for ( std::size_t i = 0; i < input.pixel_count(); ++i )
            {
                core::complex<value_t> value;
                core::convert( input[i], value );
                output[i][0] = value.real;
                output[i][1] = value.imag;
            }
        }

    public:
        FFTW( const core::size& size, bool complex = false, unsigned flags = FFTW_ESTIMATE )
            : size_(validate_size(size))
        {
            const std::size_t N = static_cast<std::size_t>( size_.area() );
            real_buffer_t real_buffer{ fftwf_alloc_real(N), &fftwf_free };
            complex_buffer_t fft_buffer_a{ fftwf_alloc_complex(N), &fftwf_free };

            if ( !real_buffer || !fft_buffer_a )
                exception_builder< std::bad_alloc >()
                    << "FFTW buffers failed to allocate";

            const int width = static_cast<int>( size_.width() );
            const int height = static_cast<int>( size_.height() );

            // Reduce axis 0, preserving the existing Hermitian-spectrum prefix.
            fftwf_iodim dims[2]{ {width, 1, 1}, {height, width, width} };
            {
                std::lock_guard<std::mutex> lock( planner_mutex() );
                if ( complex )
                {
                    c2c_forward_plan_.reset( fftwf_plan_dft_2d(
                        height, width, fft_buffer_a.get(), fft_buffer_a.get(),
                        FFTW_FORWARD, flags ) );
                    c2c_reverse_plan_.reset( fftwf_plan_dft_2d(
                        height, width, fft_buffer_a.get(), fft_buffer_a.get(),
                        FFTW_BACKWARD, flags ) );
                }
                r2c_plan_.reset( fftwf_plan_guru_dft_r2c(
                    2, dims, 0,
                    nullptr,
                    real_buffer.get(),
                    fft_buffer_a.get(),
                    flags 
                ));
                c2r_plan_.reset( fftwf_plan_guru_dft_c2r(
                    2, dims, 0, 
                    nullptr, 
                    fft_buffer_a.get(), 
                    real_buffer.get(), 
                    flags 
                ));
            }

            if ( (complex && (!c2c_forward_plan_ || !c2c_reverse_plan_)) ||
                 !r2c_plan_ || !c2r_plan_ )
                exception_builder< std::runtime_error >()
                    << "FFTW plans failed to generate";
        }

        FFTW( const FFTW& ) = delete;
        FFTW& operator=( const FFTW& ) = delete;
        FFTW( FFTW&& ) = delete;
        FFTW& operator=( FFTW&& ) = delete;

        /// Perform a 2-D FFT; will always produce a complex floating point image output
        template < template <typename> class ImageT,
                   typename ContainedT,
                   typename = typename std::enable_if_t< is_imagetype_v<ImageT<ContainedT>> >
                   >
        const image_cf_t& transform( const ImageT<ContainedT>& input, direction d = direction::FORWARD ) const
        {
            DECLARE_ENTRY_EXIT
            transform_inplace( input, d );
            auto& data = cache();
            const std::size_t N = static_cast<std::size_t>( size_.area() );
            copy_complex_to_image( data.fft_buffer_a, data.output, N );
            return data.output;
        }

        template < template <typename> class ImageT,
                   typename ContainedT,
                   typename = typename std::enable_if_t< is_imagetype_v<ImageT<ContainedT>> >
                   >
        void transform_inplace( const ImageT<ContainedT>& input, direction d = direction::FORWARD ) const
        {
            DECLARE_ENTRY_EXIT
            if ( !c2c_forward_plan_ )
                exception_builder< std::runtime_error >()
                    << "FFTW complex plans have not been generated; construct with complex=true";

            if ( input.size() != size_ )
            {
                exception_builder< std::runtime_error >()
                    << "image size is different from expected: " << input.size() << ", " << size_;
            }

            auto& data = cache();
            copy_image_to_complex( input, data.fft_buffer_a );
            transform_inplace( d );
        }

        template < template <typename> class ImageT,
                   typename ContainedT,
                   typename = typename std::enable_if_t< is_imagetype_v<ImageT<ContainedT>> >
                   >
        void transform_inplace( const ImageT<ContainedT>& a,
                                 const ImageT<ContainedT>& b,
                                 direction d = direction::FORWARD ) const
        {
            DECLARE_ENTRY_EXIT
            if ( !c2c_forward_plan_ )
                exception_builder< std::runtime_error >()
                    << "FFTW complex plans have not been generated; construct with complex=true";

            if ( a.size() != size_ || b.size() != size_ )
            {
                exception_builder< std::runtime_error >()
                    << "image size is different from expected: " << a.size()
                    << ", " << size_;
            }

            auto& data = cache();
            copy_image_to_complex( a, data.fft_buffer_a );
            copy_image_to_complex( b, data.fft_buffer_b );
            const auto plan = d == direction::FORWARD ? c2c_forward_plan_.get() : c2c_reverse_plan_.get();
            fftwf_execute_dft( plan, data.fft_buffer_a.get(), data.fft_buffer_a.get() );
            fftwf_execute_dft( plan, data.fft_buffer_b.get(), data.fft_buffer_b.get() );
        }

        void transform_inplace( direction d ) const
        {
            DECLARE_ENTRY_EXIT
            if ( !c2c_forward_plan_ )
                exception_builder< std::runtime_error >()
                    << "FFTW complex plans have not been generated; construct with complex=true";

            auto& data = cache();
            fftwf_execute_dft(
                d == direction::FORWARD ? c2c_forward_plan_.get() : c2c_reverse_plan_.get(),
                data.fft_buffer_a.get(), data.fft_buffer_a.get() );
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
            transform_real_inplace( a, b, d );
            auto& data = cache();
            const std::size_t spectrum_N = spectrum_size();
            copy_complex_to_image( data.fft_buffer_a, data.output, spectrum_N );
            copy_complex_to_image( data.fft_buffer_b, data.temp, spectrum_N );

            return { data.output, data.temp };
        }

        template < template <typename> class ImageT,
                   typename ContainedT,
                   typename = typename std::enable_if_t<
                       is_imagetype_v<ImageT<ContainedT>> &&
                       is_real_mono_pixeltype_v<ContainedT>
                       >
                   >
        void transform_real_inplace( const ImageT<ContainedT>& a,
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

            auto& data = cache();
            const std::size_t N = static_cast<std::size_t>( size_.area() );
            const std::size_t spectrum_N = spectrum_size();

            for ( std::size_t i = 0; i < N; ++i )
                data.real_buffer[i] = static_cast<value_t>( a[i].v );

            fftwf_execute_dft_r2c(
                r2c_plan_.get(), data.real_buffer.get(), data.fft_buffer_a.get() );

            for ( std::size_t i = 0; i < N; ++i )
                data.real_buffer[i] = static_cast<value_t>( b[i].v );

            fftwf_execute_dft_r2c(
                r2c_plan_.get(), data.real_buffer.get(), data.fft_buffer_b.get() );

            if ( d == direction::REVERSE )
            {
                for ( std::size_t i = 0; i < spectrum_N; ++i )
                {
                    data.fft_buffer_a[i][1] = -data.fft_buffer_a[i][1];
                    data.fft_buffer_b[i][1] = -data.fft_buffer_b[i][1];
                }
            }
        }

        void transform_real_inplace( direction d ) const
        {
            DECLARE_ENTRY_EXIT
            auto& data = cache();
            if ( d == direction::FORWARD )
            {
                const std::size_t spectrum_N = spectrum_size();
                for ( std::size_t i = 0; i < spectrum_N; ++i )
                    data.fft_buffer_a[i][1] = -data.fft_buffer_a[i][1];
            }

            fftwf_execute_dft_c2r(
                c2r_plan_.get(), data.fft_buffer_a.get(), data.real_buffer.get() );
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
                data.fft_buffer_a[i][0] = static_cast<value_t>( in[i].real );
                data.fft_buffer_a[i][1] = static_cast<value_t>( in[i].imag );
            }

            transform_real_inplace( d );

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
            OutT output{ a.size() };
            cross_correlate_inplace( a, b, output );
            return output;
        }

        /// Write the correlation into caller-owned output, reusing its storage.
        template < template <typename> class ImageT,
                   typename ContainedT,
                   typename ValueT = typename ContainedT::value_t,
                   typename OutT = image<g<ValueT>>,
                   typename = typename std::enable_if_t< is_imagetype_v<ImageT<ContainedT>> >
                   >
        void
        cross_correlate_inplace( const ImageT<ContainedT>& a,
                                 const ImageT<ContainedT>& b,
                                 OutT& output ) const
        {
            transform_inplace( a, b, direction::FORWARD );

            auto& data = cache();
            const std::size_t N = static_cast<std::size_t>( size_.area() );
            for ( std::size_t i = 0; i < N; ++i )
            {
                const value_t ar = data.fft_buffer_a[i][0];
                const value_t ai = data.fft_buffer_a[i][1];
                const value_t br = data.fft_buffer_b[i][0];
                const value_t bi = data.fft_buffer_b[i][1];
                data.fft_buffer_a[i][0] = br * ar + bi * ai;
                data.fft_buffer_a[i][1] = bi * ar - br * ai;
            }

            transform_inplace( direction::REVERSE );

            output.resize( size_ );
            for ( std::size_t i = 0; i < N; ++i )
                output[i] = static_cast<ValueT>( data.fft_buffer_a[i][0] );
            swap_quadrants( output );
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
            OutT output{ a.size() };
            cross_correlate_real_inplace( a, b, output );
            return output;
        }

        /// Write the real correlation into caller-owned output, reusing its storage.
        template < template <typename> class ImageT,
                   typename ContainedT,
                   typename ValueT = typename ContainedT::value_t,
                   typename OutT = image<g<ValueT>>,
                   typename = typename std::enable_if_t<
                       is_imagetype_v<ImageT<ContainedT>> &&
                       is_real_mono_pixeltype_v<ContainedT>
                       >
                   >
        void
        cross_correlate_real_inplace( const ImageT<ContainedT>& a,
                                      const ImageT<ContainedT>& b,
                                      OutT& output ) const
        {
            transform_real_inplace( a, b, direction::FORWARD );

            auto& data = cache();
            const std::size_t spectrum_N = spectrum_size();
            for ( std::size_t i = 0; i < spectrum_N; ++i )
            {
                const value_t ar = data.fft_buffer_a[i][0];
                const value_t ai = data.fft_buffer_a[i][1];
                const value_t br = data.fft_buffer_b[i][0];
                const value_t bi = data.fft_buffer_b[i][1];
                data.fft_buffer_a[i][0] = br * ar + bi * ai;
                data.fft_buffer_a[i][1] = bi * ar - br * ai;
            }

            transform_real_inplace( direction::REVERSE );

            output.resize( size_ );
            const std::size_t N = static_cast<std::size_t>( size_.area() );
            for ( std::size_t i = 0; i < N; ++i )
                output[i] = static_cast<ValueT>( data.real_buffer[i] );
            swap_quadrants( output );
        }

        template < template <typename> class ImageT,
                   typename ContainedT,
                   typename = typename std::enable_if_t< is_imagetype_v<ImageT<ContainedT>> >
                   >
        const image_cf_t& auto_correlate( const ImageT<ContainedT>& a ) const
        {
            transform_inplace( a, direction::FORWARD );

            auto& data = cache();
            const std::size_t N = static_cast<std::size_t>( size_.area() );
            for ( std::size_t i = 0; i < N; ++i )
            {
                const value_t ar = data.fft_buffer_a[i][0];
                const value_t ai = data.fft_buffer_a[i][1];
                data.fft_buffer_a[i][0] = ar * ar + ai * ai;
                data.fft_buffer_a[i][1] = value_t{0};
            }

            transform_inplace( direction::REVERSE );

            for ( std::size_t i = 0; i < N; ++i )
                data.output[i] = core::complex<value_t>{ data.fft_buffer_a[i][0], value_t{0} };
            swap_quadrants( data.output );

            return data.output;
        }
    };

}
