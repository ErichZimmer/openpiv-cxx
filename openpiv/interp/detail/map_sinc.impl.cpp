#include "map_sinc.impl.h"

#include <cstdlib>


namespace openpiv::interp
{
    // TODO: optimize using lookup tables (LUT) descritized to 256 with linear interp
    double sinc_weighting(double x)
    {
        if (x == 0.0)
            return 1.0;

        constexpr double pi = 3.14159265358979323846;
        double t = pi * x;
        
        return std::sin(t) / t;
    }
}