#include "map_sinc.impl.h"

#include <cstdint>
#include <cmath>

namespace openpiv::interp
{
    double sinc_weighting(double x)
    {
        if (x == 0.0)
            return 1.0;

        constexpr double pi = 3.14159265358979323846;
        double t = pi * x;
        
        return std::sin(t) / t;
    }
}