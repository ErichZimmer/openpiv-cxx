#pragma once

#include <cstdint>
#include <vector>
#include <cmath>
#include <type_traits>
#include <thread>

#include "threadpool.hpp"

#include "core/image.h"
#include "core/image_type_traits.h"
#include "core/pixel_types.h"
#include "core/vector_field.h"

#include "interp_common.impl.h"


namespace openpiv::interp
{
    using namespace openpiv::core;

    std::vector<double> poly_weights(double r, int32_t k);


    template <
        template <typename> class ImageT,
        typename ContainedT,
        typename,
        typename
    >
    void lagrange_interp2d(
        const core::image<ContainedT>& src,
        const core::grid_coords& mappings,
        core::image<ContainedT>& out,
        int32_t kernel_half_size,
        int32_t threads
    )
    {
        // Setup thread counts - 1 =  no threading; 0 = auto-select thread count; >1 = manually select thread count
        uint32_t thread_count = std::thread::hardware_concurrency()-1;
        if ((threads > 0) && (static_cast<uint32_t>(threads) < thread_count))
            thread_count = static_cast<uint32_t>(threads);

        const uint32_t src_width  = src.width();
        const uint32_t src_height = src.height();
        
        const int32_t kernel_full_size = 2 * kernel_half_size;

        out.resize(mappings.size());

        auto processor = [&src, &mappings, &out, src_width, src_height, kernel_half_size, kernel_full_size]( uint32_t ind )
        {
            const uint32_t x = (ind % mappings.width());
            const uint32_t y = (ind / mappings.width());

            const auto& point_xy = mappings[{x,y}];

            double grid_coord_x = point_xy[0];
            double grid_coord_y = point_xy[1];

            double cell_ix = std::floor(grid_coord_x);
            double cell_iy = std::floor(grid_coord_y);

            double centered_offset_x = grid_coord_x - (cell_ix + 0.5);
            double centered_offset_y = grid_coord_y - (cell_iy + 0.5);

            std::vector<double> wx = poly_weights(centered_offset_x, kernel_half_size);
            std::vector<double> wy = poly_weights(centered_offset_y, kernel_half_size);

            double value = 0.0;

            for (int32_t j = 0; j < kernel_full_size; ++j)
            {
                const size_t jj = mirror_index(cell_iy - kernel_half_size + 1 + j, src_height);

                const ContainedT* row = src.line(jj);

                for (int32_t i = 0; i < kernel_full_size; ++i)
                {
                    const size_t ii = mirror_index(cell_ix - kernel_half_size + 1 + i, src_width);

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