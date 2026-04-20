#pragma once

#include <type_traits>
#include <array>

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
        typename ResultT,
        typename
    >
    ResultT validate_range(
        const core::vector_point_data<ContainedT>& field_data,
        std::array<ValueT, 2> threshold_u,
        std::array<ValueT, 2> threshold_v
    ) {
        ResultT invalid{ field_data.u.size() };

        auto [min_u, max_u] = threshold_u;
        auto [min_v, max_v] = threshold_v;

        for (uint32_t i=0; i<invalid.pixel_count(); i++)
        {
            invalid[i] = invalid[i] || (field_data.u[i] > max_u) || (field_data.u[i] < min_u);
            invalid[i] = invalid[i] || (field_data.v[i] > max_v) || (field_data.v[i] < min_v);
        }

        return invalid;
    }

} // end of namespace