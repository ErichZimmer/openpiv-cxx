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
#include "core/exception_builder.h"
#include "core/image_utils.h"
#include "core/dll_export.h"

namespace openpiv::interp
{
    using namespace openpiv::core;

    DLL_EXPORT void poly_weights(
        double r, 
        int32_t k, 
        std::vector<double>& w
    );


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
        if ((kernel_half_size < 1) || (kernel_half_size > 5))
            core::exception_builder<std::runtime_error>()
                << "kernel_half_size must be between 1 and 5";

        if (src.pixel_count() == 0)
            core::exception_builder<std::runtime_error>()
                << "source image must not be empty";

        if (&src == &out)
            core::exception_builder<std::runtime_error>()
                << "source and output must be different images";

        // Setup thread counts - 1 =  no threading; 0 = auto-select thread count; >1 = manually select thread count
        uint32_t thread_count = std::thread::hardware_concurrency()-1;
        if ((threads > 0) && (static_cast<uint32_t>(threads) < thread_count))
            thread_count = static_cast<uint32_t>(threads);

        const uint32_t src_width  = src.width();
        const uint32_t src_height = src.height();
        
        const int32_t kernel_full_size = 2 * kernel_half_size;

        out.resize(mappings.size());

        auto processor = [&src, &mappings, &out, src_width, src_height, kernel_half_size, kernel_full_size]( uint32_t ind, std::vector<double>& wx, std::vector<double>& wy )
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

            poly_weights(centered_offset_x, kernel_half_size, wx);
            poly_weights(centered_offset_y, kernel_half_size, wy);

            const int32_t start_x = static_cast<int32_t>(cell_ix) - kernel_half_size + 1;
            const int32_t start_y = static_cast<int32_t>(cell_iy) - kernel_half_size + 1;

            const bool interior =
                start_x >= 0 &&
                start_y >= 0 &&
                start_x + kernel_full_size <= static_cast<int32_t>(src_width) &&
                start_y + kernel_full_size <= static_cast<int32_t>(src_height);

            double value = 0.0;
            size_t jj = 0;
            size_t ii = 0;

            if (interior)
            {
                for (int32_t j = 0; j < kernel_full_size; ++j)
                {
                    jj = static_cast<size_t>(start_y + j);

                    const ContainedT* row = src.line(jj);
                    double row_sum = 0.0;

                    for (int32_t i = 0; i < kernel_full_size; ++i)
                    {
                        ii = static_cast<size_t>(start_x + i);

                        row_sum += static_cast<double>(row[ii]) * wx[static_cast<size_t>(i)];
                    }

                    value += row_sum * wy[static_cast<size_t>(j)];
                }

                out[{x,y}] = static_cast<ContainedT>(value);
            }
            else
            {
                for (int32_t j = 0; j < kernel_full_size; ++j)
                {
                    jj = static_cast<size_t>(core::mirror_index<int32_t>(start_y + j, static_cast<int32_t>(src_height)));

                    const ContainedT* row = src.line(jj);
                    double row_sum = 0.0;

                    for (int32_t i = 0; i < kernel_full_size; ++i)
                    {
                        ii = static_cast<size_t>(core::mirror_index<int32_t>(start_x + i, static_cast<int32_t>(src_width)));

                        row_sum += static_cast<double>(row[ii]) * wx[static_cast<size_t>(i)];
                    }

                    value += row_sum * wy[static_cast<size_t>(j)];
                }

                out[{x,y}] = static_cast<ContainedT>(value);
            }
        };
    
        // check execution
        if (thread_count <= 1)
        {
            std::vector<double> wx(kernel_full_size);
            std::vector<double> wy(kernel_full_size);

            for (uint32_t i=0; i < mappings.pixel_count(); i++)
            {
                processor(i, wx, wy);
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
                    [i, chunk_size_, &processor, kernel_full_size]() {
                        std::vector<double> wx(kernel_full_size);
                        std::vector<double> wy(kernel_full_size);

                        for ( uint32_t j=i; j<i + chunk_size_; ++j )
                            processor(j, wx, wy);
                    } );
                i += chunk_size_;
            }
        }
    }


    template <
        template <typename> class ImageT,
        typename ContainedT,
        typename,
        typename
    >
    void lagrange_interp2n(
        const core::image<ContainedT>& src,
        const core::grid_coords& mappings,
        core::image<ContainedT>& out,
        int32_t kernel_half_size,
        int32_t threads
    ) {
        if ((kernel_half_size < 1) || (kernel_half_size > 5))
            core::exception_builder<std::runtime_error>()
                << "kernel_half_size must be between 1 and 5";

        if (src.pixel_count() == 0)
            core::exception_builder<std::runtime_error>()
                << "source image must not be empty";

        if (&src == &out)
            core::exception_builder<std::runtime_error>()
                << "source and output must be different images";

        // Early return in case if mappings is not initialized
        out.resize(mappings.size());
        if (mappings.pixel_count() == 0)
            return;

        const int32_t kernel_full_size = 2 * kernel_half_size;
        const int32_t src_width = static_cast<int32_t>(src.width());
        const int32_t src_height = static_cast<int32_t>(src.height());

        core::image<core::g_f64> horizontal(mappings.width(), src.height());

         auto horizontal_processor = [&](
            uint32_t y, 
            std::vector<double>& wx, 
            uint32_t begin_x, 
            uint32_t end_x
        ) {
            const auto* input = src.line(static_cast<size_t>(y));
            auto* output = horizontal.line(static_cast<size_t>(y));

            for (uint32_t x = begin_x; x < end_x; ++x)
            {
                const double coordinate = mappings[{x,0}][0];
                const int32_t cell_x = static_cast<int32_t>(std::floor(coordinate));
                const int32_t start_x = cell_x - kernel_half_size + 1;

                poly_weights(coordinate - (static_cast<double>(cell_x) + 0.5), kernel_half_size, wx);

                double value = 0.0;
                if (start_x >= 0 && start_x + kernel_full_size <= src_width)
                {
                    for (int32_t i = 0; i < kernel_full_size; ++i)
                        value += static_cast<double>(input[static_cast<size_t>(start_x + i)]) * wx[static_cast<size_t>(i)];
                }
                else
                {
                    for (int32_t i = 0; i < kernel_full_size; ++i)
                    {
                        const int32_t xx = core::mirror_index<int32_t>(start_x + i, src_width);
                        value += static_cast<double>(input[static_cast<size_t>(xx)]) * wx[static_cast<size_t>(i)];
                    }
                }

                output[x] = value;
            }
        };

        auto vertical_processor = [&](
            uint32_t y, 
            std::vector<double>& wy, 
            uint32_t begin_x, 
            uint32_t end_x
        ) {
            const double coordinate = mappings[{0,y}][1];
            const int32_t cell_y = static_cast<int32_t>(std::floor(coordinate));
            const int32_t start_y = cell_y - kernel_half_size + 1;

            poly_weights(coordinate - (static_cast<double>(cell_y) + 0.5), kernel_half_size, wy);

            auto* output = out.line(static_cast<size_t>(y));

            if (start_y >= 0 && start_y + kernel_full_size <= src_height)
            {
                for (uint32_t x = begin_x; x < end_x; ++x)
                {
                    double value = 0.0;
                    for (int32_t j = 0; j < kernel_full_size; ++j)
                        value += static_cast<double>(horizontal.line(static_cast<size_t>(start_y + j))[x]) * wy[static_cast<size_t>(j)];
                    
                    output[x] = static_cast<ContainedT>(value);
                }
            }
            else
            {
                for (uint32_t x = begin_x; x < end_x; ++x)
                {
                    double value = 0.0;
                    for (int32_t j = 0; j < kernel_full_size; ++j)
                    {
                        const int32_t yy = core::mirror_index<int32_t>(start_y + j, src_height);
                        value += static_cast<double>(horizontal.line(static_cast<size_t>(yy))[x]) * wy[static_cast<size_t>(j)];
                    }

                    output[x] = static_cast<ContainedT>(value);
                }
            }
        };

        uint32_t thread_count = std::thread::hardware_concurrency();
        thread_count = thread_count > 1 ? thread_count - 1 : 1;

        if ((threads > 0) && (static_cast<uint32_t>(threads) < thread_count))
            thread_count = static_cast<uint32_t>(threads);

        if (thread_count <= 1)
        {
            std::vector<double> wx(kernel_full_size);
            std::vector<double> wy(kernel_full_size);

            for (uint32_t y = 0; y < src.height(); ++y)
                horizontal_processor(y, wx, 0, out.width());

            for (uint32_t y = 0; y < out.height(); ++y)
                vertical_processor(y, wy, 0, out.width());

            return;
        }

        {
            ThreadPool pool(thread_count);

            for (uint32_t thread = 0; thread < thread_count; ++thread)
            {
                const uint32_t begin_x =
                    out.width() * thread / thread_count;
                const uint32_t end_x =
                    out.width() * (thread + 1) / thread_count;

                pool.enqueue([&, begin_x, end_x]()
                {
                    std::vector<double> wx(kernel_full_size);
                    std::vector<double> wy(kernel_full_size);

                    for (uint32_t y = 0; y < src.height(); ++y)
                        horizontal_processor(y, wx, begin_x, end_x);

                    for (uint32_t y = 0; y < out.height(); ++y)
                        vertical_processor(y, wy, begin_x, end_x);
                });
            }
        }
    }

} // end of namespace
