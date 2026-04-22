#pragma once

#include <cstdint>
#include <type_traits>
#include <tuple>

#include "core/image.h"
#include "core/image_type_traits.h"
#include "core/image_utils.h"
#include "core/pixel_types.h"
#include "core/point.h"
#include "core/vector_field.h"
#include "core/grid.h"
#include "core/dll_export.h"

#include "interp/map_polynomial.h"
#include "interp/map_sinc.h"
#include "interp/map_lanczos.h"

#include "piv/piv_common.h"


namespace openpiv::piv
{
    using namespace openpiv::core;

    DLL_EXPORT core::grid_coords generate_fine_grid(
        const core::grid_coords& coarse_grid,
        const core::size fine_size
    );

    DLL_EXPORT std::tuple<core::grid_coords, core::grid_data> create_deformation_field(
        const core::grid_coords& coarse_grid,
        const core::grid_data& coarse_data,
        const core::size fine_size
    );

    DLL_EXPORT core::grid_coords create_deformation_forward(
        const core::grid_coords& coarse_grid,
        const core::grid_data& coarse_data,
        const core::size fine_size
    );

    DLL_EXPORT std::tuple<core::grid_coords, core::grid_coords> create_deformation_symmetric(
        const core::grid_coords& coarse_grid,
        const core::grid_data& coarse_data,
        const core::size fine_size
    );


    template < 
        template <typename> class ImageT,
        typename ContainedT,
        typename ValueT = typename ContainedT::value_t,
        typename OutT = image<g<ValueT>>,
        typename = typename std::enable_if_t<
            is_imagetype_v<ImageT<ContainedT>> &&
            is_real_mono_pixeltype_v<ContainedT> &&
            std::is_floating_point<ValueT>::value
        >
    >
    OutT sparse_to_dense(
        const core::image<ContainedT>& coarse_data,
        const core::grid_coords& fine_grid
    ) {
        // Allocate memory for interpolated values
        core::image<ContainedT> fine_data(fine_grid.size());

        // Interpolate coarse data onto fine grid
        interp::lagrange_interp2d<core::image, ContainedT>(
            coarse_data,
            fine_grid,
            fine_data,
            2, // 4x4 interpolation kernel
            1 // Only use a single thread
        );

        return fine_data;
    }

            
    template<
        template <typename> class ImageT,
        typename ContainedT,
        typename,
        typename
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
    ) {
        auto frame_a_deform = frame_a;
        auto frame_b_deform = frame_b;

        if (order == deform_order::FORWARD)
        {
            auto deform_forward = create_deformation_forward(
                coarse_grid,
                coarse_data,
                frame_a.size()
            );

            if (method == deform_method::LAGRANGE)
            {
                interp::lagrange_interp2d<core::image, ContainedT>(
                    frame_b,
                    deform_forward,
                    frame_b_deform,
                    k,
                    threads
                );
            }

            else if (method == deform_method::LANCZOS)
            {
                interp::lanczos_interp2d<core::image, ContainedT>(
                    frame_b,
                    deform_forward,
                    frame_b_deform,
                    k,
                    threads
                );
            }

            else
            {
                interp::sinc_interp2d<core::image, ContainedT>(
                    frame_b,
                    deform_forward,
                    frame_b_deform,
                    k,
                    threads
                );
            }
        }
        
        else
        {
            auto [deform_backward, deform_forward] = create_deformation_symmetric(
                coarse_grid,
                coarse_data,
                frame_a.size()
            );

            if (method == deform_method::LAGRANGE)
            {
                interp::lagrange_interp2d<core::image, ContainedT>(
                    frame_a,
                    deform_backward,
                    frame_a_deform,
                    k,
                    threads
                );

                interp::lagrange_interp2d<core::image, ContainedT>(
                    frame_b,
                    deform_forward,
                    frame_b_deform,
                    k,
                    threads
                );
            }

            else if (method == deform_method::LANCZOS)
            {
                interp::lanczos_interp2d<core::image, ContainedT>(
                    frame_a,
                    deform_backward,
                    frame_a_deform,
                    k,
                    threads
                );

                interp::lanczos_interp2d<core::image, ContainedT>(
                    frame_b,
                    deform_forward,
                    frame_b_deform,
                    k,
                    threads
                );
            }

            else
            {
                interp::sinc_interp2d<core::image, ContainedT>(
                    frame_a,
                    deform_backward,
                    frame_a_deform,
                    k,
                    threads
                );

                interp::sinc_interp2d<core::image, ContainedT>(
                    frame_b,
                    deform_forward,
                    frame_b_deform,
                    k,
                    threads
                );
            }
        }

        return {frame_a_deform, frame_b_deform};
    }

} // end of namespace