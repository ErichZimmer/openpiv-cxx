#pragma once

#include <cstdint>
#include <array>
#include <tuple>

#include "core/image.h"
#include "core/pixel_types.h"
#include "core/vector_field.h"
#include "core/dll_export.h"

#include "piv/correlation_utils.h"


namespace openpiv::piv
{
    using namespace openpiv::core;

    // TODO: Add templating to support floats and doubles
    DLL_EXPORT std::tuple<core::grid_coords, core::grid_data> process_images_standard_multi(
        const ImageT& image_a,
        const ImageT& image_b,
        const core::grid_coords& old_coords,
        const core::grid_data&   old_data,
        std::array<uint32_t, 2> window_size,
        std::array<uint32_t, 2> overlap_size,
        int method,
        int order,
        int k,
        bool step,
        bool zero_pad,
        bool centered,
        bool limit_search,
        int32_t threads
    );


} // end of namespace