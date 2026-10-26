#pragma once

#include <type_traits>
#include <array>
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
    void validate_range(
        core::vector_point_data<ContainedT>& field_data,
        std::array<ValueT, 2> threshold_u,
        std::array<ValueT, 2> threshold_v
    ) {
        const core::flag_t invalid_flag = static_cast<core::flag_t>(core::FLAG::INVALID);
        const core::flag_t masked_flag = static_cast<core::flag_t>(core::FLAG::MASKED);
        const ValueT nan = std::numeric_limits<ValueT>::quiet_NaN();

        auto [min_u, max_u] = threshold_u;
        auto [min_v, max_v] = threshold_v;

        for (uint32_t i=0; i<field_data.u.pixel_count(); i++)
        {
            const core::flag_t flag = field_data.flag[i];
            if (flag == masked_flag || flag == invalid_flag)
                continue;

            const auto u = field_data.u[i];
            const auto v = field_data.v[i];
            if (!std::isfinite(static_cast<ValueT>(u)) ||
                !std::isfinite(static_cast<ValueT>(v)) ||
                u > max_u || u < min_u || v > max_v || v < min_v)
            {
                field_data.flag[i] = invalid_flag;
                field_data.u[i] = nan;
                field_data.v[i] = nan;
            }
        }
    }

} // end of namespace