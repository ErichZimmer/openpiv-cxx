#pragma once

#include <cstdint>
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
    
    // TODO: optimize using lookup tables (LUT) descritized to 256 with linear interp
    double sinc_weighting(double x);


    template <
        template <typename> class ImageT,
        typename ContainedT,
        typename,
        typename
    >
    void sinc_interp2d(
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

        out.resize(mappings.size());

        auto processor = [&src, &mappings, &out, src_width, src_height, kernel_half_size]( uint32_t ind )
        {
            const uint32_t x = (ind % mappings.width());
            const uint32_t y = (ind / mappings.width());

            const auto& point_xy = mappings[{x,y}];

            const double bx = point_xy[0];
            const double by = point_xy[1];

            const int32_t xn = static_cast<int>(bx);
            const int32_t yn = static_cast<int>(by);

            double value = 0.0;

            for (int32_t j = yn - kernel_half_size; j <= yn + kernel_half_size; ++j)
            {
                const size_t jj = mirror_index(j, src_height);

                const double dy = j - by;
                const double sy = sinc_weighting(dy);

                const ContainedT* row = src.line(jj);

                for (int32_t i = xn - kernel_half_size; i <= xn + kernel_half_size; ++i)
                {
                    const size_t ii = mirror_index(i, src_width);

                    const double dx = i - bx;
                    const double sx = sinc_weighting(dx);

                    value += static_cast<double>(row[ii]) * sx * sy;
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