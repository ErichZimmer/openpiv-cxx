#pragma once

#include <cmath>

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
    ResultT validate_z_score(
        const core::vector_point_data<ContainedT>& field_data,
        ValueT threshold
    ) {
        ResultT invalid{ field_data.u.size() };

        auto [mean_u, std_u] = algos::find_nanmean_nanstd(field_data.u);
        auto [mean_v, std_v] = algos::find_nanmean_nanstd(field_data.v);

        for (uint32_t i=0; i<invalid.pixel_count(); i++)
        {
            // Make sure nothing is NaN
            if ( !std::isfinite(static_cast<ValueT>(field_data.u[i])) ||
                 !std::isfinite(static_cast<ValueT>(field_data.v[i]))
            ) {
                invalid[i] = 1;
                continue;
            }

            invalid[i] = invalid[i] || ( std::abs(field_data.u[i] - mean_u) > (std_u * threshold) );
            invalid[i] = invalid[i] || ( std::abs(field_data.v[i] - mean_v) > (std_v * threshold) );
        }

        return invalid;
    }

} // end of namespace