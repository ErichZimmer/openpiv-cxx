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
    ResultT validate_p2p(
        const core::vector_point_data<ContainedT>& field_data,
        ValueT threshold
    ) {
        ResultT invalid{ field_data.p2p.size() };

        for (uint32_t i=0; i<invalid.pixel_count(); i++)
        {
            invalid[i] = field_data.p2p[i] < threshold;
        }

        return invalid;
    }

} // end of namespace