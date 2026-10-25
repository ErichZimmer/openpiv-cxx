#include "piv/firstpass.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

#include "algos/fftw.h"
#include "algos/pocket_fft.h"
#include "algos/stats.h"

#include "core/exception_builder.h"
#include "core/grid.h"
#include "core/image.h"
#include "core/image_expression.h"
#include "core/image_utils.h"
#include "core/pixel_types.h"
#include "core/vector_field.h"

#include "threadpool.hpp"

#include "piv/correlation_utils.h"
#include "piv/detail/extract_with_padding.impl.h"
#include "piv/detail/robust.impl.h"


namespace openpiv::piv
{

    using namespace openpiv::core;

    // Normalized squared-error cross-correlation
    std::tuple<core::grid_coords, core::grid_data> process_images_robust(
        const ImageT& image_a,
        const ImageT& image_b,
        std::array<uint32_t, 2> window_size,
        std::array<uint32_t, 2> overlap_size,
        bool step,
        bool zero_pad,
        bool centered,
        bool limit_search,
        bool simd,
        int32_t threads
    ) {
        // Force zero_pad to always be true due to nsqe implementation
        zero_pad = true;

        const uint32_t min_window_size = 8;
        const uint32_t max_window_size = 1024;
        const uint32_t min_search_size = 16; // if window_size is below this, disable limit_search

        // Make sure frame a and frame b are the same size
        if (image_a.size() != image_b.size())
        {
            core::exception_builder<std::runtime_error>()
                << "image_a and image_b sizes must be same dimensions";
        }
        
        if ((window_size[0] > max_window_size) ||
            (window_size[1] > max_window_size))
        {
            core::exception_builder<std::runtime_error>()
                << "window size must be less than "
                << max_window_size
                << " pixels";
        }

        if ((window_size[0] < min_window_size) ||
            (window_size[1] < min_window_size))
        {
            core::exception_builder<std::runtime_error>()
                << "window size must be greater than "
                << min_window_size
                << " pixels";
        }

        if ((window_size[0] % 2) || (window_size[1] % 2))
        {
            core::exception_builder<std::runtime_error>()
                << "window size must be even";
        }

        uint32_t thread_count = std::thread::hardware_concurrency() - 1;
        if ((threads > 0) && (static_cast<uint32_t>(threads) < thread_count))
            thread_count = static_cast<uint32_t>(threads);

        if (!step)
        {
            overlap_size[0] = window_size[0] - overlap_size[0];
            overlap_size[1] = window_size[1] - overlap_size[1];
        }

        // create a grid for processing
        const auto ia_size = core::size{window_size[0], window_size[1]};
        const auto grid = core::generate_cartesian_grid(
            image_b.size(),
            ia_size,
            overlap_size,
            centered
        );

        const auto field_shape = core::generate_grid_shape(
            image_b.size(),
            ia_size,
            overlap_size
        );

        // Zero pad by 2N,if requested
        auto corr_window_size = core::size{window_size[0], window_size[1]};
        if (zero_pad)
            corr_window_size = core::size{window_size[0] * 2, window_size[1] * 2};


        auto pocketfft_algo = algos::PocketFFT<FloatT>(corr_window_size);
        auto pocketfft_correlator =
            &algos::PocketFFT<FloatT>::cross_correlate_real<
                core::image,
                ContainerT
            >;

        auto fftw_algo = algos::FFTW<FloatT>(corr_window_size, false);
        auto fftw_correlator =
            &algos::FFTW<FloatT>::cross_correlate_real<
                core::image,
                ContainerT
            >;

        auto field_coords = core::grid_coords(field_shape);
        auto field_data = core::grid_data(field_shape);

        struct scratch_memory
        {
            ImageT a;
            ImageT b;
            std::vector<FloatT> int_arr;

            explicit scratch_memory(const core::size& size)
                : a(size), b(size)
            {}
        };

        auto processor = [
            &image_a,
            &image_b,
            &ia_size,
            &corr_window_size,
            &pocketfft_algo,
            &pocketfft_correlator,
            &fftw_algo,
            &fftw_correlator,
            simd,
            zero_pad,
            limit_search,
            min_search_size,
            &field_coords,
            &field_data
        ](std::size_t i, const core::rect& ia, scratch_memory& scratch_local)
        {
            auto& iw_a = scratch_local.a;
            auto& iw_b = scratch_local.b;

            detail::extract_interrogation_window(image_a, ia, iw_a);
            detail::extract_interrogation_window(image_b, ia, iw_b);

            detail::squared_integral_array<core::image, ContainerT>(
                iw_b, 
                scratch_local.int_arr
            );

            ImageT output = simd
                ? (fftw_algo.*fftw_correlator)(iw_a, iw_b)
                : (pocketfft_algo.*pocketfft_correlator)(iw_a, iw_b);

            detail::normalize_nsqecc<core::image, ContainerT>(
                output,
                iw_a,
                scratch_local.int_arr,
                ia_size
            );

            double dilation_ratio = 1.0;

            if (zero_pad)
                dilation_ratio *= 0.5;

            if (
                limit_search && 
                (ia.width() >= min_search_size || ia.height() >= min_search_size)
            )
                dilation_ratio *= 0.5;

            auto valid_corr = core::create_image_view(
                output,
                output.rect().dilate(dilation_ratio)
            );
            const auto corr_mean = algos::find_mean(valid_corr);

            constexpr uint16_t num_peaks = 2;
            constexpr uint16_t min_peak_count = 1;
            constexpr uint16_t radius = 1;

            const auto peaks = core::find_peaks_brute(
                valid_corr,
                num_peaks,
                radius
            );

            auto bl = ia.bottomLeft();
            const auto midpoint = ia.midpoint();
            field_coords[i] = midpoint;

            // Correct zero padding offsets
            if ( zero_pad )
            {
                bl[0] -= static_cast<int32_t>( (iw_a.width() - ia.width()) / 2 );
                bl[1] -= static_cast<int32_t>( (iw_a.height() - ia.height()) / 2 );
            }

            if (peaks.size() < min_peak_count)
            {
                if ( peaks.size() < min_peak_count )
            {
                field_data.u[i]     = std::numeric_limits<grid_data_t::value_t>::quiet_NaN();
                field_data.v[i]     = std::numeric_limits<grid_data_t::value_t>::quiet_NaN();
                field_data.s2n[i]   = grid_data_t(0);
                field_data.p2p[i]   = grid_data_t(0);
                field_data.peak[i]  = grid_data_t(0);
                field_data.u2[i]    = std::numeric_limits<core::grid_data_t::value_t>::quiet_NaN();
                field_data.v2[i]    = std::numeric_limits<core::grid_data_t::value_t>::quiet_NaN();
                field_data.s2n2[i]  = grid_data_t(0);
                field_data.peak2[i] = grid_data_t(0);
                field_data.flag[i]  = static_cast<uint8_t>(FLAG::INVALID);

                return;
            }
            }

            // Get subpixel information and add it to vector field data
            openpiv::core::point2<openpiv::core::grid_coords_t> subpixel_interp[num_peaks] = {};

            for (uint32_t peak_ind=0; peak_ind<peaks.size(); peak_ind++)
            {
                auto peak = peaks[peak_ind];
                core::point2<core::grid_coords_t> peak_location{};

                // Since robust correlation is always (0..1], there is no need to check if the plane is log safe
                peak_location = core::fit_simple_gaussian( peak );

                subpixel_interp[peak_ind] = peak_location;
            }

            // u and v signs are swapped to match openpiv
            field_data.u[i] = -(midpoint[0] - (bl[0] + subpixel_interp[0][0]));
            field_data.v[i] = -(midpoint[1] - (bl[1] + subpixel_interp[0][1]));
            field_data.s2n[i]  = peaks[0][{1, 1}] / corr_mean;
            field_data.p2p[i]  = (peaks.size() > min_peak_count) ? peaks[0][{1, 1}] / peaks[1][{1, 1}] : 1.0;
            field_data.peak[i] = peaks[0][{1, 1}];

            // Now populate secondary peak info
            if (peaks.size() > 1)
            {
                field_data.u2[i] = -(midpoint[0] - (bl[0] + subpixel_interp[1][0]));
                field_data.v2[i] = -(midpoint[1] - (bl[1] + subpixel_interp[1][1]));
                field_data.s2n2[i]  = peaks[1][{1, 1}] / corr_mean;
                field_data.peak2[i] = peaks[1][{1, 1}];
            }
        };

        if (thread_count > 1)
        {
            ThreadPool pool( thread_count );

            // - split the grid into thread_count chunks
            // - wrap each chunk into a processing for loop and push to thread

            // ensure we don't miss grid locations due to rounding
            size_t chunk_size = grid.size() / thread_count;
            std::vector<size_t> chunk_sizes( thread_count, chunk_size );
            chunk_sizes.back() = grid.size() - (thread_count-1)*chunk_size;

            size_t i = 0;
            for ( const auto& chunk_size_ : chunk_sizes )
            {
                pool.enqueue(
                    [i, chunk_size_, &grid, &processor, &corr_window_size]() {
                        scratch_memory scratch_local_storage(corr_window_size);

                        for ( size_t j=i; j<i + chunk_size_; ++j )
                            processor(j, grid[j], scratch_local_storage);
                    } );
                i += chunk_size_;
            }
        }
        else
        {
            scratch_memory scratch_local_storage(corr_window_size);
            
            for (size_t i = 0; i < grid.size(); ++i)
                processor(i, grid[i], scratch_local_storage);
        }

        return {std::move(field_coords), std::move(field_data)};
    }

} // namespace openpiv::piv