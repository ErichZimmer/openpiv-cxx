#pragma once

#include <type_traits>

#include "core/image.h"
#include "core/image_type_traits.h"
#include "core/pixel_types.h"
#include "core/vector_field.h"

#include "filters/convolve_sep_2d.h"
#include "filters/filter_kernels.h"

namespace openpiv::filter {
    using namespace openpiv::core;

    template <
        template<typename> class ImageT,
        typename ContainedT,
        typename ValueT,
        typename
    >
    void gaussian_lowpass(
        const ImageT<ContainedT>& src,
        ImageT<ContainedT>& out,
        ValueT sigma,
        ValueT truncate
    ) {
        auto kernel_x = generate_gaussian_kernel1d<ContainedT>(
            sigma,
            truncate
        );

        auto kernel_y = kernel_x;

        convolve_2d_sep<core::image, ContainedT>(
            src,
            out,
            kernel_x,
            kernel_y
        );
    }   

     template <
        template<typename> class ImageT,
        typename ContainedT,
        typename ValueT,
        typename ResultT,
        typename
    >
    ResultT gaussian_lowpass(
        const ImageT<ContainedT>& src,
        ValueT sigma,
        ValueT truncate
    ) {
        ResultT out{ src.size() };

        gaussian_lowpass<core::image, ContainedT>(
            src,
            out,
            sigma,
            truncate
        );

        return out;
    };  

} // end of namespace