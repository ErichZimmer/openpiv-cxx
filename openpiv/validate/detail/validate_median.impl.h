#pragma once

#include <type_traits>
#include <cmath>

#include "core/image.h"
#include "core/image_type_traits.h"
#include "core/image_utils.h"
#include "core/vector_field.h"

#include "algos/stats.h"


namespace openpiv::validate
{
    using namespace openpiv::core;
    using namespace openpiv::algos;

    template < 
        typename ContainedT,
        typename ValueT,
        typename ResultT,
        typename
    >
    ResultT validate_median(
        const core::vector_point_data<ContainedT>& field_data,
        ValueT threshold_u,
        ValueT threshold_v
    ) {
        // kernel_half_size=1 --> 3x3 kernel; kernel_half_size=2 --> 5x5 kernel
        const int32_t kernel_half_size = 1;
        const int32_t kernel_full_size = 2*kernel_half_size + 1;

        // Minimum amount of elements in kernel
        const uint32_t min_kernel_size = kernel_full_size;

        ResultT invalid{ field_data.u.size() };

        // Lambda function to process each vector point
        auto processor = [&]( uint32_t ind )
        {
            const int32_t x = static_cast<int32_t>(ind % invalid.width());
            const int32_t y = static_cast<int32_t>(ind / invalid.width());

            const auto vector_u = field_data.u[{x,y}];
            const auto vector_v = field_data.v[{x,y}];

            // If u or v is nan, skip
            if ( !std::isfinite(static_cast<ValueT>(vector_u)) ||
                 !std::isfinite(static_cast<ValueT>(vector_v))
            ) {
                invalid[{x,y}] = 1;
                return;
            }

            // Kernel for medians
            std::vector<ValueT> kernel_u;
            std::vector<ValueT> kernel_v;
            kernel_u.reserve(kernel_full_size*kernel_full_size);
            kernel_v.reserve(kernel_full_size*kernel_full_size);
            
            // Iterate over a 3x3 neighborhood of the vector
            for (int32_t j = 0; j < kernel_full_size; ++j)
            {
                const size_t jj = core::mirror_index<int32_t>(y- kernel_half_size + j, invalid.height());
                const ContainedT* row_u = field_data.u.line(jj);
                const ContainedT* row_v = field_data.v.line(jj);

                for (int32_t i = 0; i < kernel_full_size; ++i)
                {
                    // If i,j is in the center of kernel, skip
                    if ( (i==kernel_half_size) && (j==kernel_half_size) )
                        continue;
                    
                    const size_t ii = core::mirror_index<int32_t>(x - kernel_half_size + i, invalid.width());

                    const ValueT neighbor_u = static_cast<ValueT>(row_u[ii]);
                    const ValueT neighbor_v = static_cast<ValueT>(row_v[ii]);

                    // If any value is nan, skip that value
                    if ( std::isfinite(neighbor_u) )
                        kernel_u.push_back(neighbor_u);

                    if ( std::isfinite(neighbor_v) )
                        kernel_v.push_back(neighbor_v);
                }
            }

            // Calculate median
            auto u_med = algos::median(kernel_u);
            auto v_med = algos::median(kernel_v);
            
            // If not enough values are found within kernels, flag as invalid
            if ( (kernel_u.size() < min_kernel_size) || (kernel_v.size() < min_kernel_size) )
            {
                // invalid[{x,y}] = 1;
                return;
            }

            // Check if anything is above threshold
            if ( (std::abs(vector_u - u_med) > threshold_u) || 
                 (std::abs(vector_v - v_med) > threshold_v) )
            {
                invalid[{x,y}] = 1;
            }

        };

        // Now iterate each vector point
        for (uint32_t i=0; i < invalid.pixel_count(); i++)
        {
            processor(i);
        }

        return invalid;
    }

} // end of namespace