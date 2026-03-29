#pragma once

#include <cstdint>
#include <type_traits>
#include <tuple>

#include "core/image.h"
#include "core/image_type_traits.h"
#include "core/image_utils.h"
#include "core/pixel_types.h"
#include "core/vector_field.h"
#include "core/grid.h"
#include "algos/stats.h"
#include "interp/map_polynomial.h"

namespace openpiv::windef
{
    using namespace openpiv::core;
    using namespace openpiv::algos;

    core::grid_coords generate_fine_grid(
        const core::grid_coords& coarse_grid,
        const core::size fine_size
    ) {
        // Get the grid origin and spacing
        core::size origin = {
            coarse_grid[{0,0}][0],
            coarse_grid[{0,0}][1]
        };
        
        core::size spacing = {
            coarse_grid[{1,0}][0] - coarse_grid[{0,0}][0],
            coarse_grid[{0,1}][1] - coarse_grid[{0,0}][1],
        };

        // Create a fine grid
        auto fine_grid = core::grid_coords(fine_size);

        core::apply(
            fine_grid,
            [w=fine_grid.width(), h=fine_grid.height()]( auto i, auto) -> core::point2<double>
            {
                auto x = ( i % w );
                auto y = ( i / w );

                // shift to origin
                x = x - origin[0];
                y = y - origin[1];

                // scale to spacing
                x = x / spacing[0];
                y = y / spacing[1];
                
                return {x, y};
            }
        );

        return fine_grid;
    }


    template < 
        template <typename> class ImageT,
        typename ContainedT,
        typename ValueT = typename ContainedT::value_t,
        typename OutT = image<g<ValueT>>,
        typename = typename std::enable_if_t<
            is_imagetype_v<ImageT<ContainedT>> &&
            is_real_mono_pixeltype_v<ContainedT> &&
            std::is_same_v<ValueT,double>
        >
    >
    OutT sparse_to_dense(
        const core::image<ContainedT>& coarse_data,
        const core::grid_coords& fine_grid
    ) {
      

        // Allocate memory for interpolated values
        core::image<ContainedT> fine_data(fine_grid.size());

        // Interpolate coarse data onto fine grid
        interp::lagrange_interp2d(
            coarse_data,
            fine_grid,
            fine_data
            k=2 // 4x4 interpolation kernel
        );

        return fine_data;
    }


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
    std::tuple<core::grid_coords, core::grid_data> create_deformation_field(
        const core::grid_coords& coarse_grid,
        const core::grid_data& coarse_data,
        const core::size fine_size
    ) {
        auto fine_grid = generate_fine_grid(
            coarse_grid,
            fine_size
        );

        auto fine_data = core::grid_data(fine_size);

        fine_data.u = sparse_to_dense(
            coarse_data.u,
            fine_grid
        );

        fine_data.v = sparse_to_dense(
            coarse_data.v,
            fine_grid
        );

        return {fine_grid, fine_data};
    }


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
    core::grid_coords create_deformation_forward(
        const core::grid_coords& coarse_grid,
        const core::grid_data& coarse_data,
        const core::size fine_size
    ) {
        auto [fine_grid, fine_data] = create_deformation_field(
            coarse_grid,
            coarse_data,
            fine_size
        )

        core::apply(
            fine_grid,
            [w=fine_grid.width(), h=fine_grid.height(), &fine_data]( auto i, auto) -> core::point2<double>
            {
                auto x = ( i % w );
                auto y = ( i / w );

                x = x + fine_data[{x,y}].u;
                y = y + fine_data[{x,y}].v;
                
                return {x, y};
            }
        );

        return fine_grid;
    }


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
    std::tuple<core::grid_coords, core::grid_coords> create_deformation_symmetric(
        const core::grid_coords& coarse_grid,
        const core::grid_data& coarse_data,
        const core::size fine_size
    ) {
        auto [fine_grid, fine_data] = create_deformation_field(
            coarse_grid,
            coarse_data,
            fine_size
        )

        auto fine_grid_forward = fine_grid;
        auto fine_grid_reverse = fine_grid;

        core::apply(
            fine_grid_reverse,
            [w=fine_grid.width(), h=fine_grid.height(), &fine_data]( auto i, auto) -> core::point2<double>
            {
                auto x = ( i % w );
                auto y = ( i / w );

                x = x - (fine_data[{x,y}].u / ContainedT(2));
                y = y - (fine_data[{x,y}].v / ContainedT(2));
                
                return {x, y};
            }
        );

        core::apply(
            fine_grid_forward,
            [w=fine_grid.width(), h=fine_grid.height(), &fine_data]( auto i, auto) -> core::point2<double>
            {
                auto x = ( i % w );
                auto y = ( i / w );

                x = x + (fine_data[{x,y}].u / ContainedT(2));
                y = y + (fine_data[{x,y}].v / ContainedT(2));
                
                return {x, y};
            }
        );

        return {fine_grid_reverse, fine_grid_forward};
    }



} // end of namespace