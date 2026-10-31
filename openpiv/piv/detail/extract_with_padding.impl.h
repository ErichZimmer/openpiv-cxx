#pragma once

#include "core/image.h"
#include "core/image_type_traits.h"
#include "core/rect.h"
#include "core/dll_export.h"

#include "piv/correlation_utils.h"


namespace openpiv::piv::detail
{
    using namespace openpiv::core;

    DLL_EXPORT void extract_interrogation_window(
        const ImageT& source,
        const core::rect& ia,
        ImageT& iw
    );


    // Beware that squared_energy is also in detail:: namespace of robust.impl.h
    // Find squared energy, optimized for core::image and core::image_view
    template < template<typename> class ImageT,
            typename ContainedT,
            typename ValueT = typename ContainedT::value_t,
            typename = typename std::enable_if_t<
                is_imagetype_v<ImageT<ContainedT>> &&
                is_real_mono_pixeltype_v<ContainedT> &&
                std::is_floating_point_v<ValueT>
                >
            >
    ContainedT
    find_squared_energy( const ImageT<ContainedT>& im )
    {
        ValueT energy = ValueT{0};

        for (uint32_t y = 0; y < im.height(); ++y)
        {
            const auto* line = im.line(y);

            for (uint32_t x = 0; x < im.width(); ++x)
            {
                const ValueT value = static_cast<ValueT>(line[x]);

                energy += value * value;
            }
        }

        return ContainedT(energy);
    }

} // namespace piv