#pragma once

#include <type_traits>
#include <array>

#include "core/image.h"
#include "core/image_type_traits.h"
#include "core/vector_field.h"


namespace openpiv::validate
{
    using namespace openpiv::core;

    template < 
        typename ContainedT,
        typename ValueT = typename ContainedT::value_t,
        typename ResultT = core::image<g_8>,
        typename = typename std::enable_if_t<
            is_real_mono_pixeltype_v<ContainedT> &&
            std::is_floating_point_v<ValueT>
        >
    >
    ResultT validate_range(
        const core::vector_point_data<ContainedT>& field_data,
        std::array<ValueT, 2> threshold_u,
        std::array<ValueT, 2> threshold_v
    );

} // end of namespace

#include "validate/detail/validate_range.impl.h"