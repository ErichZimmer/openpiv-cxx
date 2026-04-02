#pragma once

#include <cstdint>

// local
#include "core/image.h"
#include "core/size.h"
#include "core/point.h"


namespace openpiv::core {

    template < typename T>
    struct vector_point_data
    {
        core::image<core::g<T>> u;
        core::image<core::g<T>> v;
        core::image<core::g<T>> s2n;
        core::image<core::g<T>> peak;

        vector_point_data(uint32_t w, uint32_t h)
        : u(w, h), v(w, h), s2n(w, h), peak{w, h}{}

        vector_point_data(core::size s)
        : u(s), v(s), s2n(s), peak(s) {}
    };

    using grid_data = vector_point_data<double>;
    using grid_coords  = core::image<core::point2<double>>;

}
