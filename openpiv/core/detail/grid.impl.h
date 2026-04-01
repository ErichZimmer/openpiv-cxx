
#pragma once

#include <vector>

// openpiv
#include "core/size.h"
#include "core/rect.h"
#include "core/util.h"
#include "core/dll_export.h"


namespace openpiv::core {

    DLL_EXPORT core::size
    generate_grid_shape( const core::size& image_size,
                          const core::size& interrogation_size,
                          double percentage_offset );

    DLL_EXPORT core::size
    generate_grid_shape( const core::size& image_size,
                          const core::size& interrogation_size,
                          std::array<uint32_t, 2> offsets );

    DLL_EXPORT std::vector<core::rect>
    generate_cartesian_grid( const core::size& image_size,
                             const core::size& interrogation_size,
                             double percentage_offset,
                             bool centered );

    DLL_EXPORT std::vector<core::rect>
    generate_cartesian_grid( const core::size& image_size,
                             const core::size& interrogation_size,
                             std::array<uint32_t, 2> offsets,
                             bool centered );

}
