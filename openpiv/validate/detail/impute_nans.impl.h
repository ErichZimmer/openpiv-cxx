#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include "core/exception_builder.h"
#include "core/image.h"
#include "core/pixel_types.h"
#include "core/vector_field.h"


namespace openpiv::validate
{
    using namespace openpiv::core;
    
    // Impute nans through iterative imputation.
    // This results in nearly identical imputations compared to spring analogies
    // as long as the weights are the same (4 node equal weights, or an additional 4 corners with 1/sqrt(2) weights)
    // This way avoids having to solve a systen of (sparse) linear equations or using conjugate gradient method.
    // Difference between spring systen implementation and this implementation is usually around 1e-6.
    template < 
        typename ContainedT,
        typename ValueT,
        typename
    >
    bool impute_nans(
        core::vector_point_data<ContainedT>& field_data,
        bool include_corners,
        ContainedT tolerance,
        uint32_t max_iterations
    ) {
        if ( field_data.u.size() != field_data.v.size() ||
             field_data.u.size() != field_data.s2n.size() ||
             field_data.u.size() != field_data.peak.size() ||
             field_data.u.size() != field_data.flag.size() )
        {
            core::exception_builder<std::invalid_argument>()
                << "u, v, s2n, peak, and flag must have the same size";
        }

        if ( !std::isfinite(static_cast<ValueT>(tolerance)) || tolerance < 0.0 )
        {
            core::exception_builder<std::invalid_argument>()
                << "tolerance must be finite and non-negative";
        }

        if ( max_iterations == 0 )
        {
            core::exception_builder<std::invalid_argument>()
                << "max_iterations must be greater than zero";
        }

        // This is mostly to guard against integer overflow
        if ( max_iterations > 500 )
        {
            core::exception_builder<std::invalid_argument>()
                << "max_iterations must be less than 500";
        }

        const core::flag_t invalid_flag = static_cast<core::flag_t>(core::FLAG::INVALID);
        const core::flag_t masked_flag = static_cast<core::flag_t>(core::FLAG::MASKED);
        const core::flag_t interpolated_flag = static_cast<core::flag_t>(core::FLAG::INTERPOLATED);
        const ValueT nan = std::numeric_limits<ValueT>::quiet_NaN();

        core::image<core::g<core::flag_t>> invalid{ field_data.u.size() };
        uint32_t invalid_count = 0;
        uint32_t valid_count = 0;
        ContainedT fill_value = ContainedT{0};

        // Find invalid vectors
        for ( uint32_t i = 0; i < field_data.u.pixel_count(); ++i )
        {
            const core::flag_t flag = static_cast<core::flag_t>(field_data.flag[i]);

            if ( flag == invalid_flag )
            {
                invalid[i] = invalid_flag;
                ++invalid_count;

                continue;
            }

            if ( flag == masked_flag )
                continue;

            const ValueT u = static_cast<ValueT>(field_data.u[i]);
            const ValueT v = static_cast<ValueT>(field_data.v[i]);

            if ( std::isfinite(u) && std::isfinite(v) )
            {
                ++valid_count;
            }
        }

        if ( invalid_count == 0 )
            return true;

        if ( valid_count == 0 )
        {
            core::exception_builder<std::runtime_error>()
                << "cannot replace invalid vectors without a finite valid vector";
        }

        // Set invalid values to an initial guess
        for ( uint32_t i = 0; i < invalid.pixel_count(); ++i )
        {
            if ( static_cast<core::flag_t>(invalid[i]) == 0 )
                continue;

            field_data.u[i] = fill_value;
            field_data.v[i] = fill_value;
        }

        const int32_t width = static_cast<int32_t>(field_data.u.width());
        const int32_t height = static_cast<int32_t>(field_data.u.height());
        const ValueT diagonal_weight = ValueT{1} / std::sqrt(ValueT{2});
        bool converged = false;

        // Iterate over the invalid vectors to propogate valid values to initial estimates
        for ( uint32_t iteration = 0; iteration < max_iterations; ++iteration )
        {
            ValueT maximum_change = ValueT{0};

            for ( int32_t y = 0; y < height; ++y )
            {
                for ( int32_t x = 0; x < width; ++x )
                {
                    if ( static_cast<core::flag_t>(invalid[{x,y}]) == 0 )
                        continue;

                    ValueT sum_u = 0.0;
                    ValueT sum_v = 0.0;
                    ValueT weight_sum = 0.0;

                    for ( int32_t yy = -1; yy <= 1; ++yy )
                    {
                        const int32_t neighbor_y = y + yy;

                        if ( neighbor_y < 0 || neighbor_y >= height )
                            continue;

                        for ( int32_t xx = -1; xx <= 1; ++xx )
                        {
                            if ( xx == 0 && yy == 0 )
                                continue;

                            if ( !include_corners && (xx != 0) && (yy != 0) )
                            {
                                continue;
                            }

                            const int32_t neighbor_x = x + xx;

                            if ( neighbor_x < 0 || neighbor_x >= width )
                                continue;

                            if ( static_cast<core::flag_t>(field_data.flag[{neighbor_x, neighbor_y}]) == masked_flag )
                            {
                                continue;
                            }

                            // Use ValueT raw type to suppress conversion warnings im clang on Windows
                            const ValueT neighbor_u = static_cast<ValueT>(
                                field_data.u[{neighbor_x, neighbor_y}]
                            );
                            const ValueT neighbor_v = static_cast<ValueT>(
                                field_data.v[{neighbor_x, neighbor_y}]
                            );

                            if ( !std::isfinite(neighbor_u) ||
                                 !std::isfinite(neighbor_v) )
                            {
                                continue;
                            }

                            const ValueT weight =
                                ((xx != 0) && (yy != 0))
                                    ? diagonal_weight
                                    : ValueT{1};

                            sum_u += weight * neighbor_u;
                            sum_v += weight * neighbor_v;
                            weight_sum += weight;
                        }
                    }

                    if ( weight_sum == 0.0 )
                        continue;

                    const ValueT old_u = static_cast<ValueT>(field_data.u[{x,y}]);
                    const ValueT old_v = static_cast<ValueT>(field_data.v[{x,y}]);
                    const ValueT new_u = sum_u / weight_sum;
                    const ValueT new_v = sum_v / weight_sum;
                    const ValueT change = std::hypot(new_u - old_u, new_v - old_v);

                    field_data.u[{x,y}] = ContainedT{new_u};
                    field_data.v[{x,y}] = ContainedT{new_v};

                    maximum_change = std::max(maximum_change, change);
                }
            }

            if ( maximum_change <= tolerance )
            {
                converged = true;
                break;
            }
        }

//        if ( converged )
//        {
            for ( uint32_t i = 0; i < invalid.pixel_count(); ++i )
            {
                if ( static_cast<core::flag_t>(invalid[i]) != 0 )
                {
                    field_data.s2n[i] = nan;
                    field_data.peak[i] = nan;
                    field_data.flag[i] = interpolated_flag;
                }
            }
 //       }

        return converged;
    }

} // end of namespace