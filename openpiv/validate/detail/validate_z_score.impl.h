#pragma once

#include <cmath>
#include <limits>

#include "core/image.h"
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
    void validate_z_score(
        core::vector_point_data<ContainedT>& field_data,
        ValueT threshold
    ) {
        const core::flag_t invalid_flag = static_cast<core::flag_t>(core::FLAG::INVALID);
        const core::flag_t masked_flag = static_cast<core::flag_t>(core::FLAG::MASKED);
        const ValueT nan = std::numeric_limits<ValueT>::quiet_NaN();

        uint32_t valid_count = 0;
        for (uint32_t i=0; i<field_data.u.pixel_count(); i++)
        {
            const core::flag_t flag = field_data.flag[i];
            if (flag == masked_flag || flag == invalid_flag)
                continue;

            if (std::isfinite(static_cast<ValueT>(field_data.u[i])) &&
                std::isfinite(static_cast<ValueT>(field_data.v[i])))
                ++valid_count;
        }

        // Compact samples so excluded vectors do not contribute to the divisor
        core::image<ContainedT> valid_u(valid_count, 1);
        core::image<ContainedT> valid_v(valid_count, 1);
        uint32_t sample_ind = 0;
        for (uint32_t i=0; i<field_data.u.pixel_count(); i++)
        {
            const core::flag_t flag = field_data.flag[i];
            if (flag == masked_flag || flag == invalid_flag)
                continue;

            if (std::isfinite(static_cast<ValueT>(field_data.u[i])) &&
                std::isfinite(static_cast<ValueT>(field_data.v[i])))
            {
                valid_u[sample_ind] = field_data.u[i];
                valid_v[sample_ind] = field_data.v[i];
                ++sample_ind;
            }
        }

        auto [mean_u, std_u] = algos::find_nanmean_nanstd(valid_u);
        auto [mean_v, std_v] = algos::find_nanmean_nanstd(valid_v);

        for (uint32_t i=0; i<field_data.u.pixel_count(); i++)
        {
            const core::flag_t flag = field_data.flag[i];
            if (flag == masked_flag || flag == invalid_flag)
                continue;

            const auto u = field_data.u[i];
            const auto v = field_data.v[i];
            if (!std::isfinite(static_cast<ValueT>(u)) ||
                !std::isfinite(static_cast<ValueT>(v)) ||
                std::abs(u - mean_u) > std_u * threshold ||
                std::abs(v - mean_v) > std_v * threshold)
            {
                field_data.flag[i] = invalid_flag;
                field_data.u[i] = nan;
                field_data.v[i] = nan;
            }
        }
    }

} // end of namespace