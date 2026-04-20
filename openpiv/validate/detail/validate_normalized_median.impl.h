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
    ResultT validate_normalized_median(
        const core::vector_point_data<ContainedT>& field_data,
        ValueT threshold
    ) {
        // Minimum kernel size after skipping NaNs
        const uint32_t min_kernel_size = 2;

        // EPS for normalized kernel
        const ValueT EPS = 0.1;


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

                    // If any value is nan, skip that value
                    if ( std::isfinite(row_u[ii]) )
                        kernel_u.push_back(row_u[ii]);

                    if ( std::isfinite(row_v[ii]) )
                        kernel_v.push_back(row_v[ii]);
                }
            }

            // Calculate median
            auto med_u = algos::median(kernel_u);
            auto med_v = algos::median(kernel_v);
            
            // If not enough values are found within kernels, flag as invalid
            if ( (kernel_u.size() < min_kernel_size) || (kernel_v.size() < min_kernel_size) )
            {
                invalid[{x,y}] = 1;
                return;
            }

            // Obtain residual
            for (size_t i=0; i < kernel_u.size(); ++i)
                kernel_u[i] = std::abs(kernel_u[i] - med_u);
                
            for (size_t i=0; i < kernel_v.size(); ++i)
                kernel_v[i] = std::abs(kernel_v[i] - med_v);

            auto res_u = algos::median(kernel_u);
            auto res_v = algos::median(kernel_v); 

            // Calculate normalized median
            auto norm_u = std::abs(vector_u - med_u) / (res_u + EPS);
            auto norm_v = std::abs(vector_v - med_v) / (res_v + EPS);

            auto norm_disp = std::sqrt(
                (norm_u * norm_u) + (norm_v * norm_v)
            );

            // Check if anything is above threshold
            if (norm_disp > threshold)
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