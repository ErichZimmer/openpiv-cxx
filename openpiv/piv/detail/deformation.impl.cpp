#include "deformation.impl.h"

#include <cstdint>

#include "core/point.h"
#include "core/vector_field.h"
#include "core/grid.h"


namespace openpiv::piv
{
    using namespace openpiv::core;

    core::grid_coords generate_fine_grid(
        const core::grid_coords& coarse_grid,
        const core::size fine_size
    ) {
        // Get the grid origin and spacing so we can convert to pixel units
        core::point2<uint32_t> origin = {
            coarse_grid[{0,0}][0],
            coarse_grid[{0,0}][1]
        };
        
        core::point2<uint32_t> spacing = {
            coarse_grid[{1,0}][0] - coarse_grid[{0,0}][0],
            coarse_grid[{0,1}][1] - coarse_grid[{0,0}][1],
        };

        // Create a fine grid
        auto fine_grid = core::grid_coords(fine_size);

        // Iterate over field grid to make mappings in pixel unites
        for (uint32_t y_ind=0; y_ind<fine_size.height(); y_ind++)
        {
            for (uint32_t x_ind=0; x_ind<fine_size.width(); x_ind++)
            {
                double x = static_cast<double>(x_ind);
                double y = static_cast<double>(y_ind);

                // shift to origin
                x = x - origin[0];
                y = y - origin[1];

                // scale spacing to pixels
                x = x / spacing[0];
                y = y / spacing[1];

                fine_grid[{x_ind,y_ind}] = {x,y};
            }
        }

        // Requires c++20
//        core::apply(
//            fine_grid,
//            [w=fine_grid.width(), h=fine_grid.height(), origin, spacing]( auto i, auto) -> core::point2<double>
//            {
//                auto x = ( i % w );
//                auto y = ( i / w );
//
//                // shift to origin
//                x = x - origin[0];
//                y = y - origin[1];
//
//                // scale spacing to pixels
//                x = x / spacing[0];
//                y = y / spacing[1];
//                
//                return {x, y};
//            }
//        );

        return fine_grid;
    }


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

        fine_data.u = sparse_to_dense<core::image, core::g_f64>(
            coarse_data.u,
            fine_grid
        );

        fine_data.v = sparse_to_dense<core::image, core::g_f64>(
            coarse_data.v,
            fine_grid
        );

        return {fine_grid, fine_data};
    }


    core::grid_coords create_deformation_forward(
        const core::grid_coords& coarse_grid,
        const core::grid_data& coarse_data,
        const core::size fine_size
    ) {
        auto [fine_grid, fine_data] = create_deformation_field(
            coarse_grid,
            coarse_data,
            fine_size
        );

        // Overwrite fine grid values with update grid coords
        // We basically are matching this in NumPy syntax
        // dx, dy = np.meshgrid(np.arange(img_shape[1]), np.arange(img_shape[0]))
        // dx = dx + du # du is dense interpolation of u
        // dy = dy + dv # du is dense interpolation of v
        for (uint32_t y_ind=0; y_ind<fine_size.height(); y_ind++)
        {
            for (uint32_t x_ind=0; x_ind<fine_size.width(); x_ind++)
            {
                double x = static_cast<double>(x_ind);
                double y = static_cast<double>(y_ind);

                auto u = fine_data.u[{x_ind,y_ind}];
                auto v = fine_data.v[{x_ind,y_ind}];

                fine_grid[{x_ind,y_ind}] = {
                    x + u,
                    y + v
                };
            }
        }

        return fine_grid;
    }


    std::tuple<core::grid_coords, core::grid_coords> create_deformation_symmetric(
        const core::grid_coords& coarse_grid,
        const core::grid_data& coarse_data,
        const core::size fine_size
    ) {
        auto [fine_grid, fine_data] = create_deformation_field(
            coarse_grid,
            coarse_data,
            fine_size
        );

        auto fine_grid_forward = fine_grid;
        auto fine_grid_reverse = fine_grid;

        // Overwrite fine grid values with update grid coords
        // We basically are matching this in NumPy syntax
        // dx, dy = np.meshgrid(np.arange(img_shape[1]), np.arange(img_shape[0]))
        // dx = dx - du/2 # du is dense interpolation of u
        // dy = dy - dv/2 # du is dense interpolation of v
        for (uint32_t y_ind=0; y_ind<fine_size.height(); y_ind++)
        {
            for (uint32_t x_ind=0; x_ind<fine_size.width(); x_ind++)
            {
                double x = static_cast<double>(x_ind);
                double y = static_cast<double>(y_ind);

                auto u = fine_data.u[{x_ind,y_ind}];
                auto v = fine_data.v[{x_ind,y_ind}];

                fine_grid_reverse[{x_ind,y_ind}] = {
                    x - (u / 2),
                    y - (v / 2)
                };
                
                fine_grid_forward[{x_ind,y_ind}] = {
                    x + (u / 2),
                    y + (v / 2)
                };
            }
        }

        return {fine_grid_reverse, fine_grid_forward};
    }

} // end of namespace