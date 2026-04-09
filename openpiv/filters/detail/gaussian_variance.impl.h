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

     template <
        template<typename> class ImageT,
        typename ContainedT,
        typename ValueT,
        typename 
    >
    void gaussian_variance(
        const ImageT<ContainedT>& src,
        ImageT<ContainedT>& out,
        ValueT sigma1,
        ValueT sigma2,
        ValueT truncate,
        bool clip
    ) {

        core::image<ContainedT> temp1{ src.size() };
        core::image<ContainedT> temp2{ src.size() };

        gaussian_lowpass(
            src,
            temp1,
            sigma2,
            truncate
        );

        temp2 = temp1 * temp1;

        // Make sure out is zeroed out
        core::fill(out, ContainedT(0));

        // Then temp2 array to out
        gaussian_lowpass(
            temp2,
            out,
            sigma1,
            truncate
        );

        // and back to temp2 and zero out
        out.swap(temp2);
        core::fill(out, ContainedT(0));

        // Square root temp2
        std::transform(
            temp2.begin(), temp2.end(), 
            temp2.begin(), 
            [](ValueT x) { return std::sqrt(x); }
        );

        // Divide temp1 by temp2 and clip at 0
        ContainedT zero = ContainedT(0);
        
        for (size_t i=0; i < temp1.pixel_count(); i++)
        {
            // check for zero denominator
            out[i] = (temp2[i] != zero) ? (temp1[i] / temp2[i]) : 0;
        }

        // Remove negative pixels
        if (clip)
        {
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
    ResultT gaussian_variance(
        const ImageT<ContainedT>& src,
        ValueT sigma1,
        ValueT sigma2,
        ValueT truncate,
        bool clip
    ) {
        ResultT out{ src.size() };

        gaussian_variance<core::image, ContainedT>(
            src,
            out,
            sigma1,
            sigma2,
            truncate,
            clip
        );

        return out;
    }

} // end of namespace