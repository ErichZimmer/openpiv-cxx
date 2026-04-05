#pragma once

#include <cstdint>
#include <array>
#include <tuple>

#include "core/image.h"
#include "core/pixel_types.h"
#include "core/vector_field.h"

#include "piv/piv_common.h"


namespace openpiv::piv
{
    using namespace openpiv::core;

    // TODO: Add future templating to allow float or double precision
    using FloatT = double;
    using ContainerT = core::g<FloatT>;
    using ImageT = core::image<ContainerT>;

    template<
        template <typename> class ImageT,
        typename ContainedT,
        typename ValueT = typename ContainedT::value_t,
        typename = typename std::enable_if_t<
            is_imagetype_v<ImageT<ContainedT>> &&
            is_real_mono_pixeltype_v<ContainedT> &&
            std::is_same_v<ValueT,double>
        >
    >
    std::tuple<core::image<ContainedT>, core::image<ContainedT>> deform_images(
        const core::image<ContainedT>& frame_a,
        const core::image<ContainedT>& frame_b,
        const core::grid_coords& coarse_grid,
        const core::grid_data& coarse_data,
        deform_method method,
        deform_order order,
        int32_t k,
        int32_t threads
    );

} // end of namespace


#include "piv/detail/deformation.impl.h"