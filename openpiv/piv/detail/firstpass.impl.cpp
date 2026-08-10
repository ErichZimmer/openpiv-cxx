#include "piv/firstpass.h"

#include <atomic>
#include <thread>
#include <exception>
#include <cmath>
#include <vector>
#include <array>
#include <limits>
#include <tuple>

#include "algos/pocket_fft.h"
#include "algos/duccfft.h"
#include "algos/stats.h"

#include "core/enumerate.h"
#include "core/exception_builder.h"
#include "core/grid.h"
#include "core/image.h"
#include "core/pixel_types.h"
#include "core/image_utils.h"
#include "core/image_expression.h"
#include "core/stream_utils.h"
#include "core/vector.h"
#include "core/vector_field.h"

#include "threadpool.hpp"

#include "piv/correlation_utils.h"
#include "piv/detail/extract_with_padding.impl.h"


namespace openpiv::piv
{

    using namespace openpiv::core;
    
    // Basic cross-correlation
    std::tuple<core::grid_coords, core::grid_data> process_images_standard(
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
    ){
        const uint32_t min_window_size = 8;
        const uint32_t max_window_size = 1024;

        // Make sure that window sizes are some sane value
        if ((window_size[0] > max_window_size) || (window_size[1] > max_window_size))
            core::exception_builder<std::runtime_error>() << "window size must be less than " << max_window_size << " pixels";

        if ((window_size[0] < min_window_size) || (window_size[1] < min_window_size))
            core::exception_builder<std::runtime_error>() << "window size must be greater than " << min_window_size << " pixels";

        // assert that the window size is even. Odd stuff throughs off the offsets
        if ((window_size[0] % 2) || (window_size[1] % 2))
            core::exception_builder<std::runtime_error>() << "window size must be even";

        // Setup thread counts - 1 =  no threading; 0 = auto-select thread count; >1 = manually select thread count
        uint32_t thread_count = std::thread::hardware_concurrency()-1;
        if ((threads > 0) && (static_cast<uint32_t>(threads) < thread_count))
            thread_count = static_cast<uint32_t>(threads);

        if (!step)
        {
            overlap_size[0] = window_size[0] - overlap_size[0];
            overlap_size[1] = window_size[1] - overlap_size[1];
        }

        // create a grid for processing
        auto ia_size = core::size{window_size[0], window_size[1]};
        auto grid = core::generate_cartesian_grid(
            image_b.size(), 
            ia_size, 
            overlap_size,
            centered
        );

        auto field_shape = core::generate_grid_shape(
            image_b.size(), 
            ia_size, 
            overlap_size
        );

        // Zero pad by 2N,if requested
        auto corr_window_size = core::size{window_size[0], window_size[1]};
        if (zero_pad)
            corr_window_size = core::size{window_size[0] * 2, window_size[1] * 2};

        // Get FFT correlator (this is somewhat ugly due to pointer to function, but is the most concise?)
        auto pocketfft_algo = algos::PocketFFT<FloatT>( corr_window_size );
        auto pocketfft_correlator = &algos::PocketFFT<FloatT>::cross_correlate_real<core::image, ContainerT>;

        auto duccfft_algo = algos::DuccFFT<FloatT>( corr_window_size );
        auto duccfft_correlator = &algos::DuccFFT<FloatT>::cross_correlate_real<core::image, ContainerT>;


        // Unitform weights for FFT correlation
        // TODO: Add Gaussian weights
        auto corr_weights = ImageT(corr_window_size);
        
        if (zero_pad)
        {
            // Gaussian weights for linear correlation to attenuate the periodic signals
            // While not really zero padding, it is done this way to avoid spectral leakage            
/*          
            core::apply(
                corr_weights,
                [
                    w = corr_weights.width(),
                    h = corr_weights.height()
                ]
                (auto i, auto) -> ContainerT
                {
                    const size_t x = i % w;
                    const size_t y = i / w;

                    const double cx = (static_cast<double>(w) - 1.0) / 2.0;
                    const double cy = (static_cast<double>(h) - 1.0) / 2.0;

                    const double dx = (static_cast<double>(x) - cx) / static_cast<double>(w);
                    const double dy = (static_cast<double>(y) - cy) / static_cast<double>(h);

                    return ContainerT( std::exp(-16.0 * (dx * dx + dy * dy)) );
                }
            );
*/

            // Currently, little variation was seen in gaussian vs linear. May go back to a simpler lienar correlation format

            core::apply(
                corr_weights,
                [
                    w  = corr_weights.width(),
                    h  = corr_weights.height(),
                    wx = static_cast<size_t>(window_size[0]),
                    wy = static_cast<size_t>(window_size[1])
                ]
                (auto i, auto) -> ContainerT
                {
                    const size_t index = static_cast<size_t>(i);
                    const size_t x = index % w;
                    const size_t y = index / w;

                    const size_t x0 = (w - wx) / 2;
                    const size_t y0 = (h - wy) / 2;

                    const bool inside =
                        x >= x0 && x < x0 + wx &&
                        y >= y0 && y < y0 + wy;

                    return inside ? ContainerT(1.0) : ContainerT(0.0);
                }
            );
        }
        else
        {
            // Uniform weights for circular correlation
            core::fill(corr_weights, ContainerT(1.0));
        }

        // Container for vector field
        auto field_coords = core::grid_coords(field_shape);
        auto field_data = core::grid_data(field_shape);

        // Container for reusud memory
        struct scratch_memory
        {
            ImageT a;
            ImageT b;

            explicit scratch_memory(const core::size& size)
                : a(size), b(size)
            {}
        };

        // Lamba func to process PIV image pairs
        auto processor = [
            &image_a,
            &image_b,
            &corr_window_size,
            &corr_weights,
            &pocketfft_algo,
            &pocketfft_correlator,
            &duccfft_algo,
            &duccfft_correlator,
            simd,
            zero_pad,
            limit_search,
            &field_coords,
            &field_data
        ]( size_t i, const core::rect& ia, scratch_memory& scratch_local)
        {
            const core::rect extract_window = zero_pad ? ia.dilate(2.0) : ia;
            
            auto& iw_a = scratch_local.a;
            auto& iw_b = scratch_local.b;

            // Get relavant data from the images
            extract_with_padding(image_a, extract_window, iw_a);
            extract_with_padding(image_b, extract_window, iw_b);

            // Standardize the image
            auto view_a_mean = algos::find_mean(iw_a);
            auto view_b_mean = algos::find_mean(iw_b);

            iw_a = iw_a - ContainerT(view_a_mean);
            iw_b = iw_b - ContainerT(view_b_mean);

            iw_a = iw_a * corr_weights;
            iw_b = iw_b * corr_weights;

            // On gaussian weights, re center around mean
/*
            if (zero_pad)
            {
                auto view_a_mean2 = algos::find_mean(iw_a);
                auto view_b_mean2 = algos::find_mean(iw_b);

                iw_a = iw_a - ContainerT(view_a_mean2);
                iw_b = iw_b - ContainerT(view_b_mean2);
            }
*/

            // Correlate the image extracts
            ImageT output = simd
                ? (duccfft_algo.*duccfft_correlator)(iw_a, iw_b)
                : (pocketfft_algo.*pocketfft_correlator)(iw_a, iw_b);

            // Reduce output correlation matrix size to only contain valid values
            double dilation_ratio = 1.0;

            if (zero_pad)
                dilation_ratio *= 0.5;
            
            if (limit_search)
                dilation_ratio *= 0.5;

            auto valid_corr = core::create_image_view( output, output.rect().dilate(dilation_ratio) );

            // Get mean of valid_corr to calculate s2n ratio
            auto corr_mean = algos::find_mean(valid_corr);
            
            // find peaks
            // core::peaks_t<core::g_f64> peaks;
            constexpr uint16_t num_peaks = 2;
            constexpr uint16_t min_peak_count = 1;
            constexpr uint16_t radius = 1;

            auto peaks = core::find_peaks_brute( valid_corr, num_peaks, radius );

            // Add grid to data
            auto bl = extract_window.bottomLeft();
            auto midpoint = ia.midpoint();

            field_coords[i] = midpoint;
            //field_coords[i][1] = image_a.height() - midpoint[1];
            

            // Early escape if not enough peaks were found
            if ( peaks.size() < min_peak_count )
            {
                field_data.u[i]    = std::numeric_limits<grid_data_t::value_t>::quiet_NaN();
                field_data.v[i]    = std::numeric_limits<grid_data_t::value_t>::quiet_NaN();
                field_data.s2n[i]  = grid_data_t(0);
                field_data.p2p[i]  = grid_data_t(0);
                field_data.peak[i] = grid_data_t(0);
                field_data.flag[i] = static_cast<uint8_t>(FLAG::INVALID);

                return;
            }
            
            // Get subpixel information and add it to vector field data
            auto peak = peaks[0];
            auto peak_location = core::fit_simple_gaussian( peak );

            // u and v signs are swapped to match openpiv
            field_data.u[i] = -(midpoint[0] - (bl[0] + peak_location[0]));
            field_data.v[i] = -(midpoint[1] - (bl[1] + peak_location[1]));
            field_data.s2n[i]  = peaks[0][{1, 1}] / corr_mean;
            field_data.p2p[i]  = (peaks.size() > min_peak_count) ? peaks[0][{1, 1}] / peaks[1][{1, 1}] : 1.0;
            field_data.peak[i] = peaks[0][{1, 1}];
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

        return {field_coords, field_data};
    }

} // end of namespace