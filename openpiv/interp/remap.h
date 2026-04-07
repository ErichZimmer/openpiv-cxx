#pragma once

#include <vector>
#include <cmath>
#include <type_traits>

#include "core/vector_field.h"
#include "core/image.h"
#include "core/image_type_traits.h"
#include "core/pixel_types.h"

#include "interp/map_polynomial.h"


namespace openpiv::interp
{
    using namespace openpiv::core;

    template <
        template<typename> class ImageT,
        typename ContainedT,
        typename ValueT = typename ContainedT::value_t,
        typename = std::enable_if_t<
            is_imagetype_v<ImageT<ContainedT>> &&
            is_real_mono_pixeltype_v<ContainedT> &&
            std::is_floating_point<ValueT>::value
        >
    >
    void remap2d(
        const core::grid_coords& grid_coords,
        const core::image<ContainedT>& grid_data,
        const core::grid_coords& mappings,
        core::image<ContainedT>& out
    ) {
        // EPS accounts for floating point error
        // const double EPS = 1e-4;

        // Get the grid origin and spacing so we can convert to pixel units
        const core::point2<uint32_t> origin = {
            grid_coords[{0,0}][0],
            grid_coords[{0,0}][1]
        };
        
        const core::point2<uint32_t> spacing = {
            grid_coords[{1,0}][0] - grid_coords[{0,0}][0],
            grid_coords[{0,1}][1] - grid_coords[{0,0}][1],
        };

        // TODO: Make interpolation function take 1-D grid vector of grid for x and y
        // Make sure grid is equidistant for ALL values (incase user accidentaly changes something)
        /*
        for (uint32_t y = 0; y < grid_coords.height() - 1; y++)
        {
            for (uint32_t x = 0; x < grid_coords.width() - 1; x++)
            {
                const auto coord_x1 = grid_)coords[{x,y}];
                const auto coord_x2 = grid_)coords[{x+1,y}];
                const auto coord_y1 = grid_)coords[{x,y}];
                const auto coord_y2 = grid_)coords[{x,y+1}];

                core::point2<uint32_t> spacing_g = {
                     coords_x2[0] - coord_x1[0],
                     coords_y2[1] - coord_y1[1]
                };

                spacing_g = spacing_g - spacing;

                if ((std::abs(spacing_g[0]) > EPS) || (std::abs(spacing_g[1]) > EPS))
                {
                    // Do something/raise error


                    continue;
                }
            }
        }
        */

        // Map interpolation values in pixels (e.g., grid spacing --> pixel spacing)
        out.resize(mappings.size());

        auto mappings_px = mappings;

        for (uint32_t y = 0; y < mappings.height(); ++y)
        {
            for (uint32_t x = 0; x < mappings.width(); ++x)
            {
                double x_val = mappings[{x,y}][0];
                double y_val = mappings[{x,y}][1];

                // shift to origin
                x_val = x_val - origin[0];
                y_val = y_val - origin[1];

                // scale spacing to pixels
                x_val = x_val / spacing[0];
                y_val = y_val / spacing[1];

                mappings_px[{x,y}] = {
                    x_val,
                    y_val
                };
            }
        }

        // Now interpolate!
        interp::lagrange_interp2d<core::image, ContainedT>(
            grid_data,
            mappings_px,
            out,
            1, // 2x2 interpolation kernel
            1 // Only use a single thread
        );
    }

}// end of namespace