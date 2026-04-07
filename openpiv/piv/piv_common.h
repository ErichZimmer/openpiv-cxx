#pragma once


namespace openpiv::piv {

    enum class deform_method {
        SINC,
        LAGRANGE,
        LANCZOS
    };

    enum class deform_order {
        FORWARD,
        SYMMETRIC
    };

} // end of namespace