#pragma once

#include <vector>
#include <cmath>
#include <type_traits>

#include "core/vector_field.h"
#include "core/image.h"
#include "core/image_type_traits.h"
#include "core/pixel_types.h"
#include "interp_common.impl.h"


namespace openpiv::interp
{
    using namespace openpiv::core;

    inline std::vector<double> poly_weights(double r, int k);

    template <
        template<typename> class ImageT,
        typename ContainedT,
        typename ValueT = typename ContainedT::value_t,
        typename = std::enable_if_t<
            is_imagetype_v<ImageT<ContainedT>> &&
            is_real_mono_pixeltype_v<ContainedT> &&
            std::is_same_v<ValueT, double>
        >
    >
    void lagrange_interp2d(
        const core::image<ContainedT>& src,
        const core::grid_coords& mappings,
        core::image<ContainedT>& out,
        int kernel_half_size
    )
    {
        const int width  = src.width();
        const int height = src.height();
        
        const int kernel_full_size = 2 * kernel_half_size;

        out.resize(mappings.size());

        for (uint32_t y = 0; y < mappings.height(); ++y)
        {
            for (uint32_t x = 0; x < mappings.width(); ++x)
            {
                const auto& point_xy = mappings[{x,y}];

                double grid_coord_x = point_xy[0];
                double grid_coord_y = point_xy[1];

                double cell_ix = std::floor(grid_coord_x);
                double cell_iy = std::floor(grid_coord_y);

                double centered_offset_x = grid_coord_x - (cell_ix + 0.5);
                double centered_offset_y = grid_coord_y - (cell_iy + 0.5);

                std::vector<double> wx = poly_weights(centered_offset_x, kernel_half_size);
                std::vector<double> wy = poly_weights(centered_offset_y, kernel_half_size);

                double value = 0.0;

                for (int j = 0; j < kernel_full_size; ++j)
                {
                    int jj = mirror_index(cell_iy - kernel_half_size + 1 + j, height);

                    const ContainedT* row = src.line(jj);

                    for (int i = 0; i < kernel_full_size; ++i)
                    {
                        int ii = mirror_index(cell_ix - kernel_half_size + 1 + i, width);

                        value += row[ii] * wx[i] * wy[j];
                    }
                }

                out[{x,y}] = static_cast<ContainedT>(value);
            }
        }
    }
}