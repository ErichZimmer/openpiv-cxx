#include "multipass.impl.h"

#include <atomic>
#include <thread>
#include <exception>
#include <cmath>
#include <vector>
#include <array>
#include <tuple>

#include "algos/pocket_fft.h"
#include "algos/stats.h"

#include "core/enumerate.h"
#include "core/exception_builder.h"
#include "core/grid.h"
#include "core/image.h"
#include "core/pixel_types.h"
#include "core/image_utils.h"
#include "core/image_expression.h"
#include "core/stream_utils.h"
#include "core/vector.h"
#include "core/vector_field.h"

#include "threadpool.hpp"

#include "interp/remap.h"

#include "piv/deformation.h"
#include "piv/evaluation.h"

namespace openpiv::piv
{

    using namespace openpiv::core;
    using namespace openpiv::interp;

    // Basic cross-correlation
    std::tuple<core::grid_coords, core::grid_data> process_images_standard_multi(
        const ImageT& image_a,
        const ImageT& image_b,
        const core::grid_coords& old_coords,
        const core::grid_data&   old_data,
        std::array<uint32_t, 2> window_size,
        std::array<uint32_t, 2> overlap_size,
        int method,
        int order,
        int k,
        bool step,
        bool zero_pad,
        bool centered,
        bool limit_search,
        int32_t threads
    ) {
        // Don't bother with checks as process image standard takes care of most of them
        if (!step)
        {
            overlap_size[0] = window_size[0] - overlap_size[0];
            overlap_size[1] = window_size[1] - overlap_size[1];
        }

        // create a grid for processing
        auto ia_size = core::size{window_size[0], window_size[1]};
        auto grid = core::generate_cartesian_grid(
            image_b.size(), 
            ia_size, 
            overlap_size,
            centered
        );

        auto field_shape = core::generate_grid_shape(
            image_b.size(), 
            ia_size, 
            overlap_size
        );

        // Populate field_coords with grid
        auto field_grid = core::grid_coords(field_shape);
        auto field_data = core::grid_data(field_shape);

        for (uint32_t i=0; i <field_grid.pixel_count(); i++)
            field_grid[i] = grid[i].midpoint();


        // Interpolate old data onto new grid
        interp::remap2d<core::image, ContainerT>(
            old_coords,
            old_data.u,
            field_grid,
            field_data.u
        );

        interp::remap2d<core::image, ContainerT>(
            old_coords,
            old_data.v,
            field_grid,
            field_data.v
        );

        // Now deform the images
        auto [frame_a, frame_b] = piv::deform_images<core::image, ContainerT>(
            image_a,
            image_b,
            field_grid,
            field_data,
            method,
            order,
            k
        );

        return piv::process_images_standard(
            frame_a,
            frame_b,
            window_size,
            overlap_size,
            step,
            zero_pad,
            centered,
            limit_search,
            threads
        );
    }

} // end of namespace