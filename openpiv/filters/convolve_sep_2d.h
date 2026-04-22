#pragma once

#include <type_traits>

#include "core/image.h"
#include "core/image_type_traits.h"
#include "core/pixel_types.h"
#include "core/vector_field.h"


namespace openpiv::filter
{
    using namespace openpiv::core;

    /**
     * @brief 1D convolution filters over the x and y axes of a grayscale image.
     * 
     * @param src The input image to filter.
     * @param out The output of the filtered image.
     * @param kernel_x_orig The 1D vector of weight for the x-axis.
     * @param kernel_y_orig The 1D vector of weight for the y-axis.
     * @return void
     *
     * @note Only floating point data types are allowed.
     * @note The 1D kernels are assumed to form a square kernel.
     * @note Image borders are handled using a mirror scheme (dcb|abcd|cba).
     * @note Implementation based from https://github.com/chaowang15/fast-image-convolution-cpp
     */
     template <
        template<typename> class ImageT,
        typename ContainedT,
        typename ValueT = typename ContainedT::value_t,
        typename = std::enable_if_t<
            is_imagetype_v<ImageT<ContainedT>> &&
            is_real_mono_pixeltype_v<ContainedT> &&
            std::is_floating_point_v<ValueT>
        >
    >
    void convolve_2d_sep(
        const ImageT<ContainedT>& src,
        ImageT<ContainedT>& out,
        const std::vector<ValueT>& kernel_x_orig,
        const std::vector<ValueT>& kernel_y_orig
    );

} // end of namespace

#include "filters/detail/convolve_2d_sep.impl.h"