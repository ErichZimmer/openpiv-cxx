#pragma once

#include "core/image.h"
#include "core/rect.h"

#include "piv/correlation_utils.h"


namespace openpiv::piv
{
    using namespace openpiv::core;

    void extract_with_padding(
        const ImageT& source,
        const core::rect& extract_window,
        ImageT& destination
    );

} // namespace piv