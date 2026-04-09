#pragma once

#include <algorithm>
#include <type_traits>

#include "core/image.h"
#include "core/image_type_traits.h"
#include "core/image_utils.h"
#include "core/pixel_types.h"
#include "core/vector_field.h"

#include "filters/gaussian_lowpass.h"


namespace openpiv::filter {
    using namespace openpiv::core;

    /**
     * @brief Perform a variance filter based on difference of two Gaussian kernels.
     * 
     * @param src The input image to filter.
     * @param out The output of the filtered image.
     * @param sigma1 The standard deviation of the first Gaussian kernel.
     * @param sigma2 The standard deviation of the second Gaussian kernel.
     * @param truncate The cutoff for the gaussian kernel (radius = sigma * truncate).
     * @param clip Clip negative pixel intensities to zero.
     * @return void.
     *
     * @note Only floating point data types are allowed.
     * @note Image borders are handled using a mirror scheme (dcb|abcd|cba).
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
    void gaussian_variance(
        const ImageT<ContainedT>& src,
        ImageT<ContainedT>& out,
        ValueT sigma1,
        ValueT sigma2,
        ValueT truncate,
        bool clip
    );


    /**
     * @brief Perform a variance filter based on difference of two Gaussian kernels.
     * 
     * @param src The input image to filter.
     * @param sigma1 The standard deviation of the first Gaussian kernel.
     * @param sigma2 The standard deviation of the second Gaussian kernel.
     * @param truncate The cutoff for the gaussian kernel (radius = sigma * truncate).
     * @param clip Clip negative pixel intensities to zero.
     * @return The output of the filtered image.
     *
     * @note Only floating point data types are allowed.
     * @note Image borders are handled using a mirror scheme (dcb|abcd|cba).
     */
     template <
        template<typename> class ImageT,
        typename ContainedT,
        typename ValueT = typename ContainedT::value_t,
        typename ResultT = core::image<ContainedT>,
        typename = std::enable_if_t<
            is_imagetype_v<ImageT<ContainedT>> &&
            is_real_mono_pixeltype_v<ContainedT> &&
            std::is_floating_point_v<ValueT>
        >
    >
    ResultT gaussian_variance(
        const ImageT<ContainedT>& src,
        ValueT sigma1,
        ValueT sigma2,
        ValueT truncate,
        bool clip
    );

} // end of namespace


#include "filters/detail/gaussian_variance.impl.h"