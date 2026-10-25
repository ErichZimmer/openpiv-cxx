#pragma once

#include <cstdint>
#include <limits>

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
        core::image<ContainedT> u1;
        core::image<ContainedT> v1;
        core::image<ContainedT> s2n1;
        core::image<ContainedT> peak1;
        core::image<ContainedT> u2;
        core::image<ContainedT> v2;
        core::image<ContainedT> s2n2;
        core::image<ContainedT> peak2;
        core::image<ContainedT> u3;
        core::image<ContainedT> v3;
        core::image<ContainedT> s2n3;
        core::image<ContainedT> peak3;
        core::image<core::g_8> peak_ind;
        core::image<core::g_8> flag;

        vector_point_data(uint32_t w, uint32_t h)
        : vector_point_data(core::size{w, h}) {}

        vector_point_data(core::size s)
        : u(s), v(s), s2n(s), p2p(s), peak(s)
        , u1(   s, std::numeric_limits<typename ContainedT::value_t>::quiet_NaN())
        , v1(   s, std::numeric_limits<typename ContainedT::value_t>::quiet_NaN())
        , s2n1( s, std::numeric_limits<typename ContainedT::value_t>::quiet_NaN())
        , peak1(s, std::numeric_limits<typename ContainedT::value_t>::quiet_NaN())
        , u2(   s, std::numeric_limits<typename ContainedT::value_t>::quiet_NaN())
        , v2(   s, std::numeric_limits<typename ContainedT::value_t>::quiet_NaN())
        , s2n2( s, std::numeric_limits<typename ContainedT::value_t>::quiet_NaN())
        , peak2(s, std::numeric_limits<typename ContainedT::value_t>::quiet_NaN())
        , u3(   s, std::numeric_limits<typename ContainedT::value_t>::quiet_NaN())
        , v3(   s, std::numeric_limits<typename ContainedT::value_t>::quiet_NaN())
        , s2n3( s, std::numeric_limits<typename ContainedT::value_t>::quiet_NaN())
        , peak3(s, std::numeric_limits<typename ContainedT::value_t>::quiet_NaN())
        , peak_ind(s), flag(s) {}
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
