#include "extract_with_padding.impl.h"

#include <cstdint>

#include "core/image.h"
#include "core/rect.h"
#include "core/size.h"
#include "core/image_utils.h"
#include "core/util.h"

#include "piv/correlation_utils.h"


namespace openpiv::piv::detail
{
    using namespace openpiv::core;

    void
    extract_interrogation_window(
        const ImageT& source,
        const core::rect& ia,
        ImageT& iw
    ) {
        const std::size_t destination_x = (iw.width() - ia.width()) / 2;
        const std::size_t destination_y =  (iw.height() - ia.height()) / 2;

        const std::size_t source_x = static_cast<std::size_t>(ia.left());
        const std::size_t source_y = static_cast<std::size_t>(ia.bottom());

        for (std::size_t y = 0; y < ia.height(); ++y)
        {
            core::typed_memcpy(
                iw.line(destination_y + y) + destination_x,
                source.line(source_y + y) + source_x,
                ia.width()
            );
        }
    }

} // namespace piv