#pragma once

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
     * @brief Perform a highpass filter using a Gaussian kernel.
     * 
     * @param src The input image to filter.
     * @param out The output of the filtered image.
     * @param sigma The standard deviation of the Gaussian kernel.
     * @param truncate The cutoff for the gaussian kernel (radius = sigma * truncate).
     * @return void.
     *
     * @note Only floating point data types are allowed.
     * @note Image borders are handled using a mirror scheme (dcb|abcd|cba).
     */
     template <
        template<typename> class ImageT,
        typename ContainedT,
        typename ValueT,
        typename
    >
    void gaussian_highpass(
        const ImageT<ContainedT>& src,
        ImageT<ContainedT>& out,
        ValueT sigma,
        ValueT truncate,
        bool clip
    ) {
        
        core::image<ContainedT> temp(src.size());

        gaussian_lowpass(
            src,
            temp,
            sigma,
            truncate
        );

        // Make sure out is zeroed out
        // core::fill(out, ContainedT(0));

        // Subtract the low passed image and clip negative values
        for (size_t i=0; i < out.pixel_count(); i++)
        {
            out[i] = src[i] - temp[i];
        }

        // Remove negative pixels
        if (clip)
        {
            ContainedT zero = ContainedT(0);

            for (size_t i=0; i < out.pixel_count(); i++)
            {
                out[i] = (out[i] > zero) * out[i];
            }
        }
    }

    template <
        template <typename> class ImageT,
        typename ContainedT,
        typename ValueT,
        typename ResultT,
        typename
    >
    ResultT gaussian_highpass(
        const ImageT<ContainedT>& src,
        ValueT sigma,
        ValueT truncate,
        bool clip
    ) {
        ResultT out{ src.size() };

        gaussian_highpass<core::image, ContainedT>(
            src,
            out,
            sigma,
            truncate,
            clip
        );

        return out;
    }

} // end of namespace