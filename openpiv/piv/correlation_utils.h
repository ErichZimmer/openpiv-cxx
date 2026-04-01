#pragma once

#include <type_traits>
#include <limits>
#include <cstdint>

#include "core/image.h"
#include "core/image_type_traits.h"
#include "core/pixel_types.h"

namespace openpiv::piv
{

    using namespace openpiv::core;

    using FloatT = double;
    using ContainerT = core::g<FloatT>;
    using ImageT = core::image<ContainerT>;

} // end of namespace