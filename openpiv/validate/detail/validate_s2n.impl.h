#pragma once

#include <type_traits>
#include <cmath>

#include "core/image.h"
#include "core/vector_field.h"


namespace openpiv::validate
{
    using namespace openpiv::core;

    template < 
        typename ContainedT,
        typename ValueT,
        typename ResultT,
        typename
    >
    ResultT validate_s2n(
        const core::vector_point_data<ContainedT>& field_data,
        ValueT threshold
    ) {
        ResultT invalid{ field_data.u.size() };

        for (uint32_t i=0; i<invalid.pixel_count(); i++)
        {
            invalid[i] = field_data.s2n[i] < threshold;
        }

        return invalid;
    }

} // end of namespace