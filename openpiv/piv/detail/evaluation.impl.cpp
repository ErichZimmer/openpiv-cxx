#include "evaluation.impl.h"

#include <atomic>
#include <thread>
#include <cmath>
#include <vector>
#include <array>
#include <tuple>

#include "algos/pocket_fft.h"
#include "algos/stats.h"
#include "core/enumerate.h"
#include "core/grid.h"
#include "core/image.h"
#include "core/pixel_types.h"
#include "core/image_utils.h"
#include "core/image_expression.h"
#include "core/stream_utils.h"
#include "core/vector.h"
#include "Core/vector_field.h"

#include "threadpool.hpp"


namespace openpiv::piv
{

    using namespace openpiv::core;


    // Basic cross-correlation
    std::tuple<core::grid_coords, core::grid_data> process_images_standard(
        ImageT image_a,
        ImageT image_b,
        std::array<uint32_t, 2> window_size,
        std::array<uint32_t, 2> overlap_size,
        bool zero_pad = true,
        bool centered = true,
        bool limit_search = false,
        int32_t threads = 0
    ){
        // Setup thread counts - 1 =  no threading; 0 = auto-select thread count; >1 = manually select thread count
        uint32_t thread_count = std::thread::hardware_concurrency()-1;
        if (threads > 0)
            thread_count = static_cast<uint32_t>(threads);

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

        // Zero pad by 2N if requested
        auto corr_window_size = core::size{window_size[0], window_size[1]};
        if (zero_pad)
            corr_window_size = core::size{window_size[0] * 2, window_size[1] * 2};
        
        // Get FFT correlator (this is somewhat ugly due to pointer to function, but is the most concise?)
        auto fft_algo = algos::PocketFFT( corr_window_size );
        auto correlator = &algos::PocketFFT::cross_correlate_real<core::image, ContainerT>;

        // Now get the correlation normalization matrix (note, this isn't zero-normalized CC)
        auto divisor = ImageT(corr_window_size);
        divisor = divisor + ContainerT(1);

        const auto divisor_corr = (fft_algo.*correlator)( divisor, divisor );

        // Container for vector field
        auto field_coords = core::grid_coords(field_shape);
        auto field_data = core::grid_data(field_shape);

        // Lamba func to process PIV image pairs
        auto processor = [
            &image_a,
            &image_b,
            &corr_window_size,
            &fft_algo,
            &correlator,
            &divisor_corr,
            &limit_search,
            &field_coords,
            &field_data
        ]( size_t i, const core::rect& ia)
        {
            // Get relavant data from the images
            const ImageT view_a{ core::extract( image_a, ia ) };
            const ImageT view_b{ core::extract( image_b, ia ) };

            // Standardize the image
            auto [view_a_mean, view_a_std] = algos::find_mean_std(view_a);
            auto [view_b_mean, view_b_std] = algos::find_mean_std(view_b);

            // Inset image extract into correlation window
            ImageT iw_a{corr_window_size};
            ImageT iw_b{corr_window_size};

            // Make sure iw_a and iw_b are zeroed out
            //core::fill(iw_a, ContainerT(0));
            //core::fill(iw_b, ContainerT(0));

            // if the corr window size is different than ia size, adjust
            std::array<size_t, 2> offset{0};
            if (corr_window_size != ia.size())
            {
                offset[0] = ia.width() / 2;
                offset[1] = ia.height() / 2;
            }

            for (size_t j = 0; j < view_a.height(); ++j)
            {
                for (size_t i = 0; i < view_a.width(); ++i)
                {
                    iw_a[{i + offset[0], j + offset[1]}] = (view_a[{i,j}] - view_a_mean); // / view_a_std;
                    iw_b[{i + offset[0], j + offset[1]}] = (view_b[{i,j}] - view_b_mean); // / view_b_std;
                }
            }

            // Correlate the image extracts
            ImageT output{ (fft_algo.*correlator)( iw_a, iw_b ) };

            // Normalize to 0..1
            output = output / divisor_corr;

            // find peaks
            constexpr uint16_t num_peaks = 2;
            constexpr uint16_t radius = 1;

            core::peaks_t<core::g_f64> peaks;

            // Get default dilation
            double dilation_ratio = 1.0;
            if (corr_window_size != ia.size())
                dilation_ratio *= 0.5;

            if (limit_search)
            {
                // reduce search radius
                auto centre = core::create_image_view( output, output.rect().dilate(0.5) );
                peaks = core::find_peaks( centre, num_peaks, radius );
            } else {
                peaks = core::find_peaks( output, num_peaks, radius );
            }

            // Add grid to data
            auto bl = ia.bottomLeft();
            auto midpoint = ia.midpoint();

            field_coords[i] = midpoint;
            //field_coords[i][1] = image_a.height() - midpoint[1];

            // Early escape if not enough peaks were found
            if ( peaks.size() != num_peaks )
            {
                return;
            }

            // Get subpixel information and add it to vector field data
            auto peak = peaks[0];
            auto peak_location = core::fit_simple_gaussian( peak );
            // u and v signs are swapped to match openpiv
            field_data.u[i] = -(midpoint[0] - (bl[0] + peak_location[0] - offset[0]));
            field_data.v[i] = -(midpoint[1] - (bl[1] + peak_location[1] - offset[1]));
            field_data.s2n[i]  = peaks[0][{1, 1}] / peaks[1][{1, 1}];
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
                        for ( size_t j=i; j<i + chunk_size_; ++j )
                            processor(j, grid[j]);
                    } );
                i += chunk_size_;
            }
        }
        else
        {
            for (size_t i = 0; i < grid.size(); ++i)
                processor(i, grid[i]);
        }

        return {field_coords, field_data};
    }

} // end of namespace