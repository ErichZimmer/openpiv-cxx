#pragma once

#include <cstdint>
#include <cmath>
#include <type_traits>
#include <thread>

#include "threadpool.hpp"

#include "core/vector_field.h"
#include "core/image.h"
#include "core/image_type_traits.h"
#include "core/pixel_types.h"

#include "interp_common.impl.h"


namespace openpiv::interp
{
    using namespace openpiv::core;
    
    double lanczos_weighting(double x, int32_t a);


    template <
        template <typename> class ImageT,
        typename ContainedT,
        typename,
        typename
    >
    void lanczos_interp2d(
        const core::image<ContainedT>& src,
        const core::grid_coords& mappings,
        core::image<ContainedT>& out,
        int32_t kernel_half_size,
        int32_t threads
    )
    {
        // Make sure kernel size is supported
        if ( !((kernel_half_size == 3) || (kernel_half_size == 5)) )
            core::exception_builder<std::runtime_error>() << "kernel_half_size must be either 3 or 5";

        // Setup thread counts - 1 =  no threading; 0 = auto-select thread count; >1 = manually select thread count
        uint32_t thread_count = std::thread::hardware_concurrency()-1;
        if ((threads > 0) && (static_cast<uint32_t>(threads) < thread_count))
            thread_count = static_cast<uint32_t>(threads);

        const uint32_t src_width  = src.width();
        const uint32_t src_height = src.height();

        const int32_t kernel_full_size = 2*kernel_half_size + 1;

        // Precompute Lanczos weights
        auto weighting_func = [kernel_half_size](double r) -> std::vector<double>
        {
            const int32_t k = 2*kernel_half_size + 1;

            std::vector<double> w;
            w.reserve(k);

            // integer-centered nodes: -kernel_half_size to kernel_half_size-1
            for (int32_t n = -kernel_half_size; n <= kernel_half_size; ++n) {
                double x = static_cast<double>(n) - r;
                w.push_back(lanczos_weighting(x, kernel_half_size));
            }

            return w;
        };

        double r_min = static_cast<double>(-kernel_half_size);
        double r_max = static_cast<double>(kernel_half_size);
        size_t steps = 512u;

        const lut_v weights_table = make_lut(
            weighting_func,
            steps,
            r_min,
            r_max
        );

        out.resize(mappings.size());

        auto processor = [&src, &mappings, &out, &weights_table, src_width, src_height, kernel_half_size, kernel_full_size]( uint32_t ind )
        {
            const uint32_t x = (ind % mappings.width());
            const uint32_t y = (ind / mappings.width());

            const auto& point_xy = mappings[{x,y}];

            const double grid_coord_x = point_xy[0];
            const double grid_coord_y = point_xy[1];

            const int32_t cell_ix = static_cast<int32_t>(std::floor(grid_coord_x));
            const int32_t cell_iy = static_cast<int32_t>(std::floor(grid_coord_y));

            // fractional offsets in [0,1)
            const double offset_x = grid_coord_x - static_cast<double>(cell_ix);
            const double offset_y = grid_coord_y - static_cast<double>(cell_iy);

            const auto wx = get_weights(offset_x, weights_table);
            const auto wy = get_weights(offset_y, weights_table);

            double value = 0.0;

            for (int32_t j = 0; j < kernel_full_size; ++j)
            {
                const size_t jj = mirror_index(cell_iy - kernel_half_size + j, src_height);

                const ContainedT* row = src.line(jj);

                for (int32_t i = 0; i < kernel_full_size; ++i)
                {
                    const size_t ii = mirror_index(cell_ix - kernel_half_size + i, src_width);

                    value += static_cast<double>(row[ii]) * wx[i] * wy[j];
                }
            }

            out[{x,y}] = static_cast<ContainedT>(value);
        };

        // check execution
        if (thread_count <= 1)
        {
            for (uint32_t i=0; i < mappings.pixel_count(); i++)
            {
                processor(i);
            }
        }
        else
        {
            ThreadPool pool( thread_count );

            // - split the pixels into thread_count chunks
            // - wrap each chunk into a processing for loop and push to thread

            // ensure we don't miss pixel locations due to rounding
            uint32_t chunk_size = mappings.pixel_count()/thread_count;
            std::vector<uint32_t> chunk_sizes( thread_count, chunk_size );
            chunk_sizes.back() = mappings.pixel_count() - (thread_count-1)*chunk_size;

            uint32_t i = 0;
            for ( const auto& chunk_size_ : chunk_sizes )
            {
                pool.enqueue(
                    [i, chunk_size_, &processor]() {
                        for ( uint32_t j=i; j<i + chunk_size_; ++j )
                            processor(j);
                    } );
                i += chunk_size_;
            }
        }
    }

} // end of namespace