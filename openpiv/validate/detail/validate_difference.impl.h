#pragma once

#include <type_traits>
#include <cmath>
#include <limits>

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
        typename
    >
    void validate_difference(
        core::vector_point_data<ContainedT>& field_data,
        ValueT threshold_u,
        ValueT threshold_v
    ) {
        // kernel_half_size=1 --> 3x3 kernel; kernel_half_size=2 --> 5x5 kernel
        const int32_t kernel_half_size = 1;
        const int32_t kernel_full_size = 2*kernel_half_size + 1;

        // Minimum amount of elements in kernel after skipping nans
        const uint32_t min_kernel_size = kernel_full_size;

        const core::flag_t invalid_flag = static_cast<core::flag_t>(core::FLAG::INVALID);
        const core::flag_t masked_flag = static_cast<core::flag_t>(core::FLAG::MASKED);
        const ValueT nan = std::numeric_limits<ValueT>::quiet_NaN();
        core::image<core::g<core::flag_t>> invalid{ field_data.u.size() };

        // Lambda function to process each vector point
        auto processor = [&]( uint32_t ind )
        {
            const int32_t x = static_cast<int32_t>(ind % invalid.width());
            const int32_t y = static_cast<int32_t>(ind / invalid.width());

            const core::flag_t center_flag = field_data.flag[{x,y}];
            if (center_flag == masked_flag || center_flag == invalid_flag)
                return;

            const auto vector_u = field_data.u[{x,y}];
            const auto vector_v = field_data.v[{x,y}];

            // If u or v is nan, skip
            if ( !std::isfinite(static_cast<ValueT>(vector_u)) ||
                 !std::isfinite(static_cast<ValueT>(vector_v))
            ) {
                invalid[{x,y}] = invalid_flag;
                return;
            }

            // flag to see if vector is invalid
            uint32_t flag = 0;
            uint32_t count = 0;

            // Iterate over a 3x3 neighborhood of the vector
            for (int32_t j = 0; j < kernel_full_size; ++j)
            {
                const size_t jj = core::mirror_index<int32_t>(y - kernel_half_size + j, invalid.height());
                const ContainedT* row_u = field_data.u.line(jj);
                const ContainedT* row_v = field_data.v.line(jj);
                const auto* row_flag = field_data.flag.line(jj);

                for (int32_t i = 0; i < kernel_full_size; ++i)
                {
                    // if i,j is in the center of kernel, skip
                    if ( (i==kernel_half_size) && (j==kernel_half_size) )
                        continue;
                                        
                    const size_t ii = core::mirror_index<int32_t>(x - kernel_half_size + i, invalid.width());

                    const core::flag_t neighbor_flag = row_flag[ii];
                    if (neighbor_flag == masked_flag || neighbor_flag == invalid_flag)
                        continue;

                    const ValueT neighbor_u = static_cast<ValueT>(row_u[ii]);
                    const ValueT neighbor_v = static_cast<ValueT>(row_v[ii]);

                    ContainedT difference_u = std::abs(vector_u - neighbor_u);
                    ContainedT difference_v = std::abs(vector_v - neighbor_v);

                    // If any value is nan, skip that value
                    if ( !std::isfinite(neighbor_u) || !std::isfinite(neighbor_v) )
                        continue;
                    
                    count += 1;
                    
                    if ( (difference_u > threshold_u) || (difference_v > threshold_v) )
                        flag += 1;
                }
            }

            // If not enough values are found within kernels, flag as invalid
            if (count < min_kernel_size)
            {
                // invalid[{x,y}] = invalid_flag;
                return;
            }
            // If the flag is greater than 50% of kernel elements, set invalid to true
            if ( flag > ((kernel_full_size * kernel_full_size) / 2) )
                invalid[{x,y}] = invalid_flag;
        };

        // Now iterate each vector point
        for (uint32_t i=0; i < invalid.pixel_count(); i++)
        {
            processor(i);
        }

        // Apply rejections after all neighborhoods have been evaluated
        for (uint32_t i=0; i < invalid.pixel_count(); i++)
        {
            if (invalid[i])
            {
                field_data.flag[i] = invalid_flag;
                field_data.u[i] = nan;
                field_data.v[i] = nan;
            }
        }
    }

} // end of namespace