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
    ResultT validate_difference(
        const core::vector_point_data<ContainedT>& field_data,
        ValueT threshold_u,
        ValueT threshold_v
    ) {
        // 3x3 kernel size
        const uint32_t kernel_half_size = 1;
        const uint32_t kernel_full_size = 2*kernel_half_size + 1;

        ResultT invalid{ field_data.u.size() };

        // Lambda function to process each vector point
        auto processor = [&]( uint32_t ind )
        {
            const uint32_t x = (ind % invalid.width());
            const uint32_t y = (ind / invalid.width());

            const auto vector_u = field_data.u[{x,y}];
            const auto vector_v = field_data.v[{x,y}];

            // If u or v is nan, skip
            if ( !std::isfinite(static_cast<ValueT>(vector_u)) ||
                 !std::isfinite(static_cast<ValueT>(vector_v))
            ) {
                invalid[{x,y}] = 1;
                return;
            }

            // flag to see if vector is invalid
            uint32_t flag = 0;

            // Iterate over a 3x3 neighborhood of the vector
            for (int32_t j = 0; j < kernel_full_size; ++j)
            {
                const size_t jj = core::mirror_index<int32_t>(y- kernel_half_size + j, invalid.height());
                const ContainedT* row_u = field_data.u.line(jj);
                const ContainedT* row_v = field_data.v.line(jj);

                for (int32_t i = 0; i < kernel_full_size; ++i)
                {
                    // if i,j is in the center of kernel, skip
                    if ( (i==kernel_half_size) && (j==kernel_half_size) )
                        continue;
                    
                    const size_t ii = core::mirror_index<int32_t>(x - kernel_half_size + i, invalid.width());

                    ContainedT difference_u = std::abs(vector_u - row_u[ii]);
                    ContainedT difference_v = std::abs(vector_v - row_v[ii]);
                    
                    if ( (difference_u > threshold_u) || (difference_v > threshold_v) )
                        flag += 1;
                }
            }
            
            // If the flag is greater than 50% of kernel elements, set invalid to true
            if ( flag > ((kernel_full_size * kernel_full_size) / 2) )
                invalid[{x,y}] = 1;
        };

        // Now iterate each vector point
        for (uint32_t i=0; i < invalid.pixel_count(); i++)
        {
            processor(i);
        }

        return invalid;
    }

} // end of namespace