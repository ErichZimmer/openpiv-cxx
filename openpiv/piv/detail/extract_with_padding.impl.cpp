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

    void extract_with_padding(
        const ImageT& source,
        const core::rect& extract_window,
        ImageT& destination
    ) {
        const core::rect image_bounds = core::rect::from_size(source.size());

        if (image_bounds.contains(extract_window))
        {
            const size_t src_x = static_cast<size_t>(extract_window.left());
            const size_t src_y = static_cast<size_t>(extract_window.bottom());

            const size_t copy_width = extract_window.width();

            for (size_t y = 0; y < extract_window.height(); ++y)
            {
                core::typed_memcpy(
                    destination.line(y),
                    source.line(src_y + y) + src_x,
                    copy_width
                );
            }

            return;
        }

        // Clear values left from the previous window.
        core::fill(destination, ContainerT(0));

        const int32_t src_left = std::max(
            extract_window.left(),
            image_bounds.left()
        );

        const int32_t src_bottom = std::max(
            extract_window.bottom(),
            image_bounds.bottom()
        );

        const int32_t src_right = std::min(
            extract_window.right(),
            image_bounds.right()
        );

        const int32_t src_top = std::min(
            extract_window.top(),
            image_bounds.top()
        );

        // The requested rectangle does not intersect the source.
        if (src_left >= src_right || src_bottom >= src_top)
            return;

        const size_t copy_width = static_cast<size_t>(src_right - src_left);
        const size_t copy_height = static_cast<size_t>(src_top - src_bottom);

        const size_t dst_x = static_cast<size_t>(src_left - extract_window.left());
        const size_t dst_y = static_cast<size_t>(src_bottom - extract_window.bottom());

        const size_t source_x = static_cast<size_t>(src_left);
        const size_t source_y = static_cast<size_t>(src_bottom);

        for (size_t y = 0; y < copy_height; ++y)
        {
            core::typed_memcpy(
                destination.line(dst_y + y) + dst_x,
                source.line(source_y + y) + source_x,
                copy_width
            );
        }
    }

} // namespace piv