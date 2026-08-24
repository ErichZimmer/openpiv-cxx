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

#include "algos/duccfft.h"
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

        auto duccfft_algo = algos::DuccFFT<FloatT>(corr_window_size);
        auto duccfft_correlator =
            &algos::DuccFFT<FloatT>::cross_correlate_real<
                core::image,
                ContainerT
            >;

        auto corr_weights = ImageT(corr_window_size);

        if (zero_pad)
        {
            core::apply(
                corr_weights,
                [
                    w = corr_weights.width(),
                    h = corr_weights.height(),
                    wx = static_cast<std::size_t>(window_size[0]),
                    wy = static_cast<std::size_t>(window_size[1])
                ](auto i, auto) -> ContainerT
                {
                    const std::size_t index = static_cast<std::size_t>(i);
                    const std::size_t x = index % w;
                    const std::size_t y = index / w;
                    const std::size_t x0 = (w - wx) / 2;
                    const std::size_t y0 = (h - wy) / 2;
                    const bool inside =
                        x >= x0 && x < x0 + wx &&
                        y >= y0 && y < y0 + wy;

                    return inside ? ContainerT(1.0) : ContainerT(0.0);
                }
            );
        }
        else
        {
            core::fill(corr_weights, ContainerT(1.0));
        }

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
            &corr_weights,
            &pocketfft_algo,
            &pocketfft_correlator,
            &duccfft_algo,
            &duccfft_correlator,
            simd,
            zero_pad,
            limit_search,
            min_search_size,
            &field_coords,
            &field_data
        ](std::size_t i, const core::rect& ia, scratch_memory& scratch_local)
        {
            const core::rect extract_window = zero_pad ? ia.dilate(2.0) : ia;

            auto& iw_a = scratch_local.a;
            auto& iw_b = scratch_local.b;

            detail::extract_with_padding(image_a, extract_window, iw_a);
            detail::extract_with_padding(image_b, extract_window, iw_b);

            // Mean subtraction
            //auto view_a_mean = algos::find_mean(iw_a);
            //auto view_b_mean = algos::find_mean(iw_b);

            //iw_a = iw_a - ContainerT(view_a_mean);
            //iw_b = iw_b - ContainerT(view_b_mean);

            iw_a = iw_a * corr_weights;
            iw_b = iw_b * corr_weights;

            detail::squared_integral_array<core::image, ContainerT>(
                iw_b, 
                scratch_local.int_arr
            );

            ImageT output = simd
                ? (duccfft_algo.*duccfft_correlator)(iw_a, iw_b)
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

            const auto bl = extract_window.bottomLeft();
            const auto midpoint = ia.midpoint();
            field_coords[i] = midpoint;

            const bool invalid_correlation =
                !std::isfinite(static_cast<double>(corr_mean)) ||
                static_cast<double>(corr_mean) <=
                    std::numeric_limits<double>::epsilon();

            if (peaks.size() < min_peak_count || invalid_correlation)
            {
                field_data.u[i] = std::numeric_limits<grid_data_t::value_t>::quiet_NaN();
                field_data.v[i] = std::numeric_limits<grid_data_t::value_t>::quiet_NaN();
                field_data.s2n[i] = grid_data_t(0);
                field_data.p2p[i] = grid_data_t(0);
                field_data.peak[i] = grid_data_t(0);
                field_data.flag[i] = static_cast<uint8_t>(FLAG::INVALID);
                return;
            }

            // Correlation intensities are alwats (0, 1], so no need to check if it is safe to use logs
            const auto peak_location = core::fit_simple_gaussian(peaks[0]);

            field_data.u[i] = -(midpoint[0] - (bl[0] + peak_location[0]));
            field_data.v[i] = -(midpoint[1] - (bl[1] + peak_location[1]));
            field_data.s2n[i] = peaks[0][{1, 1}] / corr_mean;
            field_data.p2p[i] = peaks.size() > min_peak_count
                ? peaks[0][{1, 1}] / peaks[1][{1, 1}]
                : 1.0;
            field_data.peak[i] = peaks[0][{1, 1}];
        };

        if (thread_count > 1)
        {
            ThreadPool pool(thread_count);

            const std::size_t chunk_size = grid.size() / thread_count;
            std::vector<std::size_t> chunk_sizes(thread_count, chunk_size);
            chunk_sizes.back() =
                grid.size() - (thread_count - 1) * chunk_size;

            std::size_t i = 0;
            for (const auto chunk_size_ : chunk_sizes)
            {
                pool.enqueue(
                    [
                        i,
                        chunk_size_,
                        &grid,
                        &processor,
                        &corr_window_size
                    ]()
                    {
                        scratch_memory scratch_local_storage(corr_window_size);

                        for (std::size_t j = i; j < i + chunk_size_; ++j)
                            processor(j, grid[j], scratch_local_storage);
                    }
                );
                i += chunk_size_;
            }
        }
        else
        {
            scratch_memory scratch_local_storage(corr_window_size);

            for (std::size_t i = 0; i < grid.size(); ++i)
                processor(i, grid[i], scratch_local_storage);
        }

        return {std::move(field_coords), std::move(field_data)};
    }

} // namespace openpiv::piv