#include "map_lanczos.impl.h"

#include <cstdint>
#include <cmath>

#include "map_sinc.impl.h"

namespace openpiv::interp
{
    double lanczos_weighting(double x, int32_t a)
    {
        double ax = std::abs(x);

        if (ax >= static_cast<double>(a))
            return 0.0;
        
        // L(x) = sinc(x) * sinc(x / a)
        return sinc_weighting(x) * sinc_weighting(x / static_cast<double>(a));
    }
}