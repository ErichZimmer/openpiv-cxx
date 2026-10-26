#pragma once

#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>

#include "core/exception_builder.h"
#include "core/vector_field.h"


namespace openpiv::validate
{
    using namespace openpiv::core;

    template < 
        typename ContainedT,
        typename ValueT,
        typename
    >
    uint32_t substitute_nans(
        core::vector_point_data<ContainedT>& field_data
    ) {
        const core::flag_t invalid_flag = static_cast<core::flag_t>(core::FLAG::INVALID);
        const core::flag_t masked_flag = static_cast<core::flag_t>(core::FLAG::MASKED);
        const core::flag_t substituted_flag = static_cast<core::flag_t>(core::FLAG::SUBSTITUTED);
        const ValueT nan = std::numeric_limits<ValueT>::quiet_NaN();
        uint32_t unsucessfull_substitutions = 0;
        
        for ( uint32_t i = 0; i < field_data.u.pixel_count(); ++i )
        {
            const core::flag_t flag = static_cast<core::flag_t>(field_data.flag[i]);

            if ( flag == masked_flag )
                continue;

            if ( flag != invalid_flag )
                continue;

            const auto peak_ind = field_data.peak_ind[i];
            if ( peak_ind >= 2 )
            {
                field_data.u[i] = nan;
                field_data.v[i] = nan;
                field_data.flag[i] = invalid_flag;

                unsucessfull_substitutions += 1;

                continue;
            }

            // Advance one candidate per call, including missing candidates
            const auto next_peak_ind = peak_ind + 1;
            field_data.peak_ind[i] = next_peak_ind;

            const ContainedT next_u = (next_peak_ind == 1) ? field_data.u2[i] : field_data.u3[i];
            const ContainedT next_v = (next_peak_ind == 1) ? field_data.v2[i] : field_data.v3[i];

            if ( !std::isfinite(next_u.v) ||
                 !std::isfinite(next_v.v) )
            {
                field_data.u[i] = nan;
                field_data.v[i] = nan;
                field_data.flag[i] = invalid_flag;

                unsucessfull_substitutions += 1;

                continue;
            }

            field_data.u[i] = next_u;
            field_data.v[i] = next_v;
            field_data.s2n[i]  = (next_peak_ind == 1) ? field_data.s2n2[i] : field_data.s2n3[i];
            field_data.peak[i] = (next_peak_ind == 1) ? field_data.peak2[i] : field_data.peak3[i];
            field_data.flag[i] = substituted_flag;
        }

        return unsucessfull_substitutions;
    }

} // end of namespace