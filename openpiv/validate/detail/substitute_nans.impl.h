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
        uint32_t unsucessfull_substitutions = 0;
        
        for ( uint32_t i = 0; i < field_data.u.pixel_count(); ++i )
        {
            const core::flag_t flag = static_cast<core::flag_t>(field_data.flag[i]);

            if ( flag == masked_flag )
                continue;

            const ContainedT u = field_data.u[i];
            const ContainedT v = field_data.v[i];

            // Basic sanity check, if vector is not flagged as invalid and is finite value, skip
            if ( (flag != invalid_flag) && std::isfinite(u.v) && std::isfinite(v.v) )
            {
                continue;
            }

            const ContainedT secondary_u = field_data.u2[i];
            const ContainedT secondary_v = field_data.v2[i];

            // If secondary peak is not valid, then no need to do further processing
            if ( !std::isfinite(secondary_u.v) ||
                 !std::isfinite(secondary_v.v) )
            {
                field_data.u[i] = std::numeric_limits<ValueT>::quiet_NaN();
                field_data.v[i] = std::numeric_limits<ValueT>::quiet_NaN();
                field_data.flag[i] = invalid_flag;

                unsucessfull_substitutions += 1;

                continue;
            }

            field_data.u[i] = secondary_u;
            field_data.v[i] = secondary_v;
            field_data.s2n[i]  = field_data.s2n2[i];
            field_data.peak[i] = field_data.peak2[i];
            field_data.flag[i] = substituted_flag;
        }

        return unsucessfull_substitutions;
    }

} // end of namespace