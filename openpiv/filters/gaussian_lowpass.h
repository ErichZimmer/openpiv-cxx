#pragma once

#include <type_traits>

#include "core/image.h"
#include "core/image_type_traits.h"
#include "core/pixel_types.h"
#include "core/vector_field.h"


namespace openpiv::filter {
    using namespace openpiv::core;

    /**
     * @brief Perform a lowpass filter using a Gaussian kernel.
     * 
     * @param src The input image to filter.
     * @param out The output of the filtered image.
     * @param sigma The standard deviation of the Gaussian kernel.
     * @param truncate The cutoff for the gaussian kernel (radius = sigma * truncate).
     * @return void.
     *
     * @note Only floating point data types are allowed.
     * @note Image borders are handled using a mirror scheme (dcb|abcd|cba).
     * @note May throw runtime error if kernel size is found to be unrealistic.
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
    void gaussian_lowpass(
        const ImageT<ContainedT>& src,
        ImageT<ContainedT>& out,
        ValueT sigma,
        ValueT truncate
    );   

    /**
     * @brief Perform a lowpass filter using a Gaussian kernel.
     * 
     * @param src The input image to filter.
     * @param sigma The standard deviation of the Gaussian kernel.
     * @param truncate The cutoff for the gaussian kernel (radius = sigma * truncate).
     * @return The output of the filtered image.
     *
     * @note Only floating point data types are allowed.
     * @note Image borders are handled using a mirror scheme (dcb|abcd|cba).
     * @note May throw runtime error if kernel size is found to be unrealistic.
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
    ResultT gaussian_lowpass(
        const ImageT<ContainedT>& src,
        ValueT sigma,
        ValueT truncate
    );

    /**
     * @brief Perform a lowpass filter using a Gaussian kernel.
     * 
     * @param src The input image to filter.
     * @param out The output of the filtered image.
     * @param kernel_half_size The radius of the Gaussian kernel.
     * @return void.
     *
     * @note Only floating point data types are allowed.
     * @note Image borders are handled using a mirror scheme (dcb|abcd|cba).
     * @note May throw runtime error if kernel size is found to be unrealistic.
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
    void gaussian_lowpass(
        const ImageT<ContainedT>& src,
        ImageT<ContainedT>& out,
        uint32_t kernel_half_size
    );   

    /**
     * @brief Perform a lowpass filter using a Gaussian kernel.
     * 
     * @param src The input image to filter.
     * @param kernel_half_size The radius of the Gaussian kernel.
     * @return The output of the filtered image.
     *
     * @note Only floating point data types are allowed.
     * @note Image borders are handled using a mirror scheme (dcb|abcd|cba).
     * @note May throw runtime error if kernel size is found to be unrealistic.
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
    ResultT gaussian_lowpass(
        const ImageT<ContainedT>& src,
        uint32_t kernel_half_size
    );  

} // end of namespace


#include "filters/detail/gaussian_lowpass.impl.h"