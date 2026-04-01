#pragma once

#include <cstdint>
#include <array>
#include <tuple>

#include "core/image.h"
#include "core/pixel_types.h"
#include "core/vector_field.h"
#include "piv/correlation_utils.h"
#include "core/dll_export.h"

namespace openpiv::piv
{
    using namespace openpiv::core;

    // TODO: Add future templating to allow float or double precision
    DLL_EXPORT std::tuple<core::grid_coords, core::grid_data> process_images_standard(
        const ImageT& image_a,
        const ImageT& image_b,
        std::array<uint32_t, 2> window_size,
        std::array<uint32_t, 2> overlap_size,
        bool step,
        bool zero_pad,
        bool centered,
        bool limit_search,
        int32_t threads
    );

    /*
    std::tuple<core::grid_coords, core::grid_data> process_images_nsqe(
        ImageT image_a,
        ImageT image_b,
        std::array<uint32_t, 2> window_size,
        std::array<uint32_t, 2> overlap_size,
        bool step,
        bool zero_pad,
        bool centered,
        bool limit_search,
        int32_t threads
    );
    */

} // end of namespace