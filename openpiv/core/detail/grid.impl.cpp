#include "core/grid.h"

#include <vector>

// openpiv
#include "core/util.h"
#include "core/size.h"
#include "core/rect.h"


namespace openpiv::core {

    core::size
    generate_grid_shape( const core::size& image_size,
                          const core::size& interrogation_size,
                          double percentage_offset )
    {
        if ( percentage_offset < 0.0 || percentage_offset > 1.0 )
            core::exception_builder<std::runtime_error>() << "offsets must be between 0.0 and 1.0";

        return generate_grid_shape(
            image_size,
            interrogation_size,
            { (uint32_t)(interrogation_size.width()  * percentage_offset),
              (uint32_t)(interrogation_size.height() * percentage_offset) } );
    }

    core::size
    generate_grid_shape( const core::size& image_size,
                          const core::size& interrogation_size,
                          std::array<uint32_t, 2> offsets )
    {
        if ( image_size.area() == 0 )
            core::exception_builder<std::runtime_error>() << "image size must be non-zero";

        if ( interrogation_size.area() == 0 )
            core::exception_builder<std::runtime_error>() << "interrogation size must be non-zero";

        if ( offsets[0] == 0 || offsets[1] == 0 )
            core::exception_builder<std::runtime_error>() << "offsets must be non-zero";

        if ( interrogation_size.width() > image_size.width() ||
             interrogation_size.height() > image_size.height() )
            core::exception_builder<std::runtime_error>() << "interrogation size is bigger than image";

        auto [x_offset, y_offset] = offsets;

        auto x_count = 1 + (image_size.width()  - interrogation_size.width())  / x_offset;
        auto y_count = 1 + (image_size.height() - interrogation_size.height()) / y_offset;

        core::size field_shape(x_count, y_count);

        return field_shape;
    }

    std::vector<core::rect>
    generate_cartesian_grid( const core::size& image_size,
                             const core::size& interrogation_size,
                             double percentage_offset,
                             bool centered )
    {
        if ( percentage_offset < 0.0 || percentage_offset > 1.0 )
            core::exception_builder<std::runtime_error>() << "offsets must be between 0.0 and 1.0";

        return generate_cartesian_grid(
            image_size,
            interrogation_size,
            { (uint32_t)(interrogation_size.width()  * percentage_offset),
              (uint32_t)(interrogation_size.height() * percentage_offset) },
            centered );
    }

    std::vector<core::rect>
    generate_cartesian_grid( const core::size& image_size,
                             const core::size& interrogation_size,
                             std::array<uint32_t, 2> offsets,
                             bool centered )
    {
        if ( image_size.area() == 0 )
            core::exception_builder<std::runtime_error>() << "image size must be non-zero";

        if ( interrogation_size.area() == 0 )
            core::exception_builder<std::runtime_error>() << "interrogation size must be non-zero";

        if ( offsets[0] == 0 || offsets[1] == 0 )
            core::exception_builder<std::runtime_error>() << "offsets must be non-zero";

        if ( interrogation_size.width() > image_size.width() ||
             interrogation_size.height() > image_size.height() )
            core::exception_builder<std::runtime_error>() << "interrogation size is bigger than image";

        auto [x_offset, y_offset] = offsets;
        auto field_shape = generate_grid_shape(
            image_size,
            interrogation_size,
            offsets
        );

        auto x_count = field_shape.width();
        auto y_count = field_shape.height();

        auto x_start_offset = centered ? (image_size.width()  - interrogation_size.width()  - (x_offset * (x_count-1)))/2 : 0;
        auto y_start_offset = centered ? (image_size.height() - interrogation_size.height() - (y_offset * (y_count-1)))/2 : 0;

        std::vector<core::rect> result;
        for ( size_t y=0; y < y_count; ++y )
        {
            for ( size_t x=0; x < x_count; ++x )
            {
                rect r{ {x_start_offset + x * x_offset, y_start_offset + y * y_offset}, interrogation_size };
                result.emplace_back( std::move(r) );
            }
        }

        return result;
    }

}
