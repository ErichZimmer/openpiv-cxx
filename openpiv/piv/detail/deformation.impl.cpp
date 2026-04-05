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

} // end of namespace