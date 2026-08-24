#pragma once

#include "core/image.h"
#include "core/rect.h"
#include "core/dll_export.h"

#include "piv/correlation_utils.h"


namespace openpiv::piv::detail
{
    using namespace openpiv::core;

    DLL_EXPORT void extract_with_padding(
        const ImageT& source,
        const core::rect& extract_window,
        ImageT& destination
    );

} // namespace piv