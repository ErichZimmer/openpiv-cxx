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
#include "algos/fftw.h"
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
        const uint32_t min_search_size = 16; // if window_size is below this, disable limit_search

        // Make sure frame a and frame b are the same size
        if (image_a.size() != image_b.size())
        {
            core::exception_builder<std::runtime_error>()
                << "image_a and image_b sizes must be same dimensions";
        }

        // Make sure that window sizes are some sane value
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

        // Setup thread counts; 1 =  no threading; 0 = auto-select thread count; >1 = manually select thread count
        uint32_t thread_count = std::thread::hardware_concurrency()-1;
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
        auto corr_window_size = ia_size;
        if (zero_pad)
            corr_window_size = core::size{window_size[0] * 2, window_size[1] * 2};

        // Get FFT correlator (this is somewhat ugly due to pointer to function, but is the most concise?)
        auto pocketfft_algo = algos::PocketFFT<FloatT>( corr_window_size );
        auto pocketfft_correlator = &algos::PocketFFT<FloatT>::cross_correlate_real_inplace<core::image, ContainerT>;

        auto fftw_algo = algos::FFTW<FloatT>( corr_window_size, false );
        auto fftw_correlator = &algos::FFTW<FloatT>::cross_correlate_real_inplace<core::image, ContainerT>;


        // Uniform weights for FFT correlation
        // TODO: Add Gaussian weights
        auto corr_weights = ImageT(corr_window_size);

        // Container for vector field
        auto field_coords = core::grid_coords(field_shape);
        auto field_data = core::grid_data(field_shape);

        // Container for reusud memory
        struct scratch_memory
        {
            ImageT a;
            ImageT b;
            ImageT output;

            explicit scratch_memory(const core::size& size)
                : a(size), b(size), output(size)
            {}
        };

        // Lamba func to process PIV image pairs
        auto processor = [
            &image_a,
            &image_b,
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
        ]( size_t i, const core::rect& ia, scratch_memory& scratch_local)
        {            
            auto& iw_a = scratch_local.a;
            auto& iw_b = scratch_local.b;
            auto& output = scratch_local.output;

            // Get relavant data from the images
            detail::extract_interrogation_window(image_a, ia, iw_a);
            detail::extract_interrogation_window(image_b, ia, iw_b);

            // Mean subtraction to improve cross-correlation quality
            auto valid_iw_rect = iw_a.rect();

            if (zero_pad)
                valid_iw_rect = iw_a.rect().dilate(0.5);

            auto iw_a_view = core::create_image_view(iw_a, valid_iw_rect);
            auto iw_b_view = core::create_image_view(iw_b, valid_iw_rect);

            const auto view_a_mean = algos::find_mean(iw_a_view);
            const auto view_b_mean = algos::find_mean(iw_b_view);

            // TODO: Put this in its own function
            for (uint32_t y = 0; y < iw_a_view.height(); ++y)
            {
                auto* line_a = iw_a_view.line(y);
                auto* line_b = iw_b_view.line(y);

                for (uint32_t x = 0; x < iw_a_view.width(); ++x)
                {
                    line_a[x].v -= view_a_mean.v;
                    line_b[x].v -= view_b_mean.v;
                }
            }

            // Correlate the image extracts
            (simd)
                ? (fftw_algo.*fftw_correlator)(iw_a, iw_b, output)
                : (pocketfft_algo.*pocketfft_correlator)(iw_a, iw_b, output);

            // Reduce output correlation matrix size to only contain valid values
            double dilation_ratio = 1.0;

            if (zero_pad)
                dilation_ratio *= 0.5;
            
            if (
                limit_search && 
                (ia.width() >= min_search_size || ia.height() >= min_search_size)
            )
                dilation_ratio *= 0.5;

            auto valid_corr = core::create_image_view( output, output.rect().dilate(dilation_ratio) );

            // Get mean of valid_corr to calculate s2n ratio
            auto corr_mean = algos::find_mean(valid_corr);

            // Compute denominator to normalize correlatioon plane
            const ContainerT energy_a = detail::find_squared_energy(iw_a_view);
            const ContainerT energy_b = detail::find_squared_energy(iw_b_view);

            ContainerT denominator =
                static_cast<ContainerT>(iw_a.pixel_count()) *
                std::sqrt(energy_a) *
                std::sqrt(energy_b);
            
            // Guard against zero denominator in case of blank windows
            if (denominator == ContainerT(0.0))
                denominator = ContainerT(1.0);

            // find peaks
            // core::peaks_t<core::g_f64> peaks;
            constexpr uint16_t num_peaks = 3;
            constexpr uint16_t min_peak_count = 1;
            constexpr uint16_t radius = 1;

            auto peaks = core::find_peaks_brute(
                valid_corr, 
                num_peaks,
                radius 
            );

            // Add grid to data
            auto bl = ia.bottomLeft();
            const auto midpoint = ia.midpoint();

            field_coords[i] = midpoint;
            //field_coords[i][1] = image_a.height() - midpoint[1];

            // Correct zero padding offsets
            if ( zero_pad )
            {
                bl[0] -= static_cast<int32_t>( (iw_a.width() - ia.width()) / 2 );
                bl[1] -= static_cast<int32_t>( (iw_a.height() - ia.height()) / 2 );
            }
            

            // Early escape if not enough peaks were found
            if ( peaks.size() < min_peak_count )
            {
                field_data.u[i]     = std::numeric_limits<grid_data_t::value_t>::quiet_NaN();
                field_data.v[i]     = std::numeric_limits<grid_data_t::value_t>::quiet_NaN();
                field_data.s2n[i]   = std::numeric_limits<grid_data_t::value_t>::quiet_NaN();
                field_data.p2p[i]   = std::numeric_limits<grid_data_t::value_t>::quiet_NaN();
                field_data.peak[i]  = std::numeric_limits<grid_data_t::value_t>::quiet_NaN();
                field_data.peak_ind[i] = 0;
                field_data.flag[i]  = static_cast<uint8_t>(FLAG::INVALID);

                return;
            }
            
            // Get subpixel information and add it to vector field data
            core::point2<core::grid_coords_t> subpixel_interp[num_peaks] = {};

            for (uint32_t peak_ind=0; peak_ind<peaks.size(); peak_ind++)
            {
                auto peak = peaks[peak_ind];
                core::point2<core::grid_coords_t> peak_location{};

                if ( core::fit_log_safe(peak, radius) )
                    peak_location = core::fit_simple_gaussian( peak );
                else
                    peak_location = core::fit_simple_parabolic( peak );

                subpixel_interp[peak_ind] = peak_location;
            }

            // u and v signs are swapped to match openpiv
            field_data.u1[i] = -(midpoint[0] - (bl[0] + subpixel_interp[0][0]));
            field_data.v1[i] = -(midpoint[1] - (bl[1] + subpixel_interp[0][1]));
            field_data.s2n1[i]  = peaks[0][{1, 1}] / corr_mean;
            field_data.peak1[i] = peaks[0][{1, 1}] / denominator;

            // Now populate secondary peak info
            if (peaks.size() > 1)
            {
                field_data.u2[i] = -(midpoint[0] - (bl[0] + subpixel_interp[1][0]));
                field_data.v2[i] = -(midpoint[1] - (bl[1] + subpixel_interp[1][1]));
                field_data.s2n2[i]  = peaks[1][{1, 1}] / corr_mean;
                field_data.peak2[i] = peaks[1][{1, 1}] / denominator;
            }

            // Now populate tertiary peak info
            if (peaks.size() > 2)
            {
                field_data.u3[i] = -(midpoint[0] - (bl[0] + subpixel_interp[2][0]));
                field_data.v3[i] = -(midpoint[1] - (bl[1] + subpixel_interp[2][1]));
                field_data.s2n3[i]  = peaks[2][{1, 1}] / corr_mean;
                field_data.peak3[i] = peaks[2][{1, 1}] / denominator;
            }

            // Initialize the working fields from the primary peak
            field_data.u[i]    = field_data.u1[i];
            field_data.v[i]    = field_data.v1[i];
            field_data.s2n[i]  = field_data.s2n1[i];
            field_data.peak[i] = field_data.peak1[i];
            field_data.p2p[i]  = (peaks.size() > min_peak_count) ? peaks[0][{1, 1}] / peaks[1][{1, 1}] : 1.0;
            field_data.peak_ind[i] = 0;
            field_data.flag[i] = static_cast<uint8_t>(FLAG::VALID);
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

} // end of namespace