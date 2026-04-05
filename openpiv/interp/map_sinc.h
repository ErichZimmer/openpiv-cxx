#pragma once

#include <type_traits>

#include "core/image.h"
#include "core/image_type_traits.h"
#include "core/pixel_types.h"
#include "core/vector_field.h"


namespace openpiv::interp {
    using namespace openpiv::core;

    /**
     * @brief Map new coordinates of an image using sinc interpolation.
     * 
     * @param image The input image to interpolate.
     * @param mappings The grid coordinates at which the image is evaluated.
     * @param k The half-size (radius) of the kernel (e.g., k=3 --> 7x7 kernel).
     * @param threads The amount of threads to use during interpolation
     * @return output The interpolated image.
     *
     * @note Only kernel half-sizes of 3 and 5 are allowed due to stability reasons.
     * @note Image borders are handled using a mirror scheme (dcb|abcd|cba)
     * @note The interpolation function evaluates everything using g_f64 (doubles) under-the-hood. This may or may not require type conversions depending on the image pixel types.
     */
     template <
        template<typename> class ImageT,
        typename ContainedT,
        typename ValueT = typename ContainedT::value_t,
        typename = std::enable_if_t<
            is_imagetype_v<ImageT<ContainedT>> &&
            is_real_mono_pixeltype_v<ContainedT> &&
            std::is_same_v<ValueT, double>
        >
    >
    void sinc_interp2d(
        const core::image<ContainedT>& src,
        const core::grid_coords& mappings,
        core::image<ContainedT>& out,
        int32_t k,
        int32_t threads
    );   

} // end of namespace


#include "interp/detail/map_sinc.impl.h"