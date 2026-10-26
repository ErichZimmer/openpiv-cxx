#pragma once

#include <type_traits>
#include <cmath>
#include <limits>

#include "core/image.h"
#include "core/vector_field.h"


namespace openpiv::validate
{
    using namespace openpiv::core;

    template < 
        typename ContainedT,
        typename ValueT,
        typename
    >
    void validate_s2n(
        core::vector_point_data<ContainedT>& field_data,
        ValueT threshold
    ) {
        const core::flag_t invalid_flag = static_cast<core::flag_t>(core::FLAG::INVALID);
        const core::flag_t masked_flag = static_cast<core::flag_t>(core::FLAG::MASKED);
        const ValueT nan = std::numeric_limits<ValueT>::quiet_NaN();

        for (uint32_t i=0; i<field_data.u.pixel_count(); i++)
        {
            const core::flag_t flag = field_data.flag[i];
            if (flag == masked_flag || flag == invalid_flag)
                continue;

            const auto s2n = field_data.s2n[i];
            if (!std::isfinite(static_cast<ValueT>(field_data.u[i])) ||
                !std::isfinite(static_cast<ValueT>(field_data.v[i])) ||
                !std::isfinite(static_cast<ValueT>(s2n)) || s2n < threshold)
            {
                field_data.flag[i] = invalid_flag;
                field_data.u[i] = nan;
                field_data.v[i] = nan;
            }
        }
    }

} // end of namespace