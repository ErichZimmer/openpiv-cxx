#pragma once

#include <type_traits>

#include "core/image.h"
#include "core/image_type_traits.h"
#include "core/pixel_types.h"
#include "core/vector_field.h"


namespace openpiv::filter {
    using namespace openpiv::core;

    /**
     * @brief Generate a 1d Gaussian kernel.
     * 
     * @param sigma The standard deviation of the Gaussian kernel.
     * @param truncate The cutoff for the gaussian kernel (radius = sigma * truncate).
     * @return void.
     *
     * @note Only floating point data types are allowed.
     * @note Image borders are handled using a mirror scheme (dcb|abcd|cba).
     */
     template <
        typename ContainedT,
        typename ValueT = typename ContainedT::value_t,
        typename = std::enable_if_t<
            is_real_mono_pixeltype_v<ContainedT> &&
            std::is_floating_point_v<ValueT>
        >
    >
    std::vector<ValueT> generate_gaussian_kernel1d(
         ValueT sigma,
         ValueT truncate
    );

} // end of namespace


#include "filters/detail/filter_kernels.impl.h"