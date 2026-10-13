#pragma once

#include <cstdint>

// local
#include "core/image.h"
#include "core/size.h"
#include "core/point.h"


namespace openpiv::core {

    enum class FLAG : uint8_t
    {
        VALID = 0,
        INVALID,
        SUBSTITUTED,
        INTERPOLATED,
        MASKED
    };

    template < typename ContainedT>
    struct vector_point_data
    {
        core::image<ContainedT> u;
        core::image<ContainedT> v;
        core::image<ContainedT> s2n;
        core::image<ContainedT> p2p;
        core::image<ContainedT> peak;
        core::image<ContainedT> u2;
        core::image<ContainedT> v2;
        core::image<ContainedT> s2n2;
        core::image<ContainedT> peak2;
        core::image<core::g_8> flag;

        vector_point_data(uint32_t w, uint32_t h)
        : u(w, h), v(w, h), s2n(w, h), p2p(w, h), peak{w, h}, u2(w, h), v2(w, h), s2n2(w, h),  peak2{w, h}, flag{w, h}{}

        vector_point_data(core::size s)
        : u(s), v(s), s2n(s), p2p(s), peak(s), u2(s), v2(s), s2n2(s), peak2(s), flag(s) {}
    };

    template < typename ContainedT>
    struct deformation_point_data
    {
        core::image<ContainedT> u;
        core::image<ContainedT> v;

        deformation_point_data(uint32_t w, uint32_t h)
        : u(w, h), v(w, h){}

        deformation_point_data(core::size s)
        : u(s), v(s) {}
    };

    using grid_data_t = core::g_f64;
    using grid_coords_t = double;
    using grid_data = vector_point_data<grid_data_t>;
    using grid_coords = core::image<core::point2<grid_coords_t>>;
    using deform_data = deformation_point_data<grid_data_t>;

}
