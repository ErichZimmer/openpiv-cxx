#include "map_polynomial.impl.h"

#include <vector>
#include <cstdint>
#include <stdexcept>

namespace openpiv::interp
{
    namespace
    {
        void stencil_2(double r, std::vector<double>& w)
        {
            w[0] = 0.5 - r;
            w[1] = 0.5 + r;
        }

        void stencil_4(double r, std::vector<double>& w)
        {
            const double q = r * r;
            const double a = 0.25 - q;
            const double b = 2.25 - q;

            const double s0 = b / 2.0;
            const double s1 = -a / 6.0;

            w[0] = (1.5 - r) * s1;
            w[1] = (0.5 - r) * s0;
            w[2] = (0.5 + r) * s0;
            w[3] = (1.5 + r) * s1;
        }

        void stencil_6(double r, std::vector<double>& w)
        {
            const double q = r * r;
            const double a = 0.25 - q;
            const double b = 2.25 - q;
            const double c = 6.25 - q;

            const double s0 = b * c / 12.0;
            const double s1 = -a * c / 24.0;
            const double s2 = a * b / 120.0;

            w[0] = (2.5 - r) * s2;
            w[1] = (1.5 - r) * s1;
            w[2] = (0.5 - r) * s0;
            w[3] = (0.5 + r) * s0;
            w[4] = (1.5 + r) * s1;
            w[5] = (2.5 + r) * s2;
        }

        void stencil_8(double r, std::vector<double>& w)
        {
            const double q = r * r;
            const double a = 0.25 - q;
            const double b = 2.25 - q;
            const double c = 6.25 - q;
            const double d = 12.25 - q;
            const double ab = a * b;
            const double cd = c * d;

            const double s0 = b * cd / 144.0;
            const double s1 = -a * cd / 240.0;
            const double s2 = ab * d / 720.0;
            const double s3 = -ab * c / 5040.0;

            w[0] = (3.5 - r) * s3;
            w[1] = (2.5 - r) * s2;
            w[2] = (1.5 - r) * s1;
            w[3] = (0.5 - r) * s0;
            w[4] = (0.5 + r) * s0;
            w[5] = (1.5 + r) * s1;
            w[6] = (2.5 + r) * s2;
            w[7] = (3.5 + r) * s3;
        }

        void stencil_10(double r, std::vector<double>& w)
        {
            const double q = r * r;
            const double a = 0.25 - q;
            const double b = 2.25 - q;
            const double c = 6.25 - q;
            const double d = 12.25 - q;
            const double e = 20.25 - q;
            const double ab = a * b;
            const double de = d * e;
            const double cde = c * de;
            const double abc = ab * c;

            const double s0 = b * cde / 2880.0;
            const double s1 = -a * cde / 4320.0;
            const double s2 = ab * de / 10080.0;
            const double s3 = -abc * e / 40320.0;
            const double s4 = abc * d / 362880.0;

            w[0] = (4.5 - r) * s4;
            w[1] = (3.5 - r) * s3;
            w[2] = (2.5 - r) * s2;
            w[3] = (1.5 - r) * s1;
            w[4] = (0.5 - r) * s0;
            w[5] = (0.5 + r) * s0;
            w[6] = (1.5 + r) * s1;
            w[7] = (2.5 + r) * s2;
            w[8] = (3.5 + r) * s3;
            w[9] = (4.5 + r) * s4;
        }
    }

    void poly_weights(
        double r, 
        int32_t k, 
        std::vector<double>& w
    ) {
        if (k < 1 || k > 5)
            core::exception_builder<std::runtime_error>()
                << "kernel_half_size must be between 1 and 5";

        if (w.size() < static_cast<size_t>(2 * k))
            core::exception_builder<std::runtime_error>()
                << "weight buffer must contain at least " << 2 * k << " elements";

        switch (k)
        {
            case 1: stencil_2(r, w); break;
            case 2: stencil_4(r, w); break;
            case 3: stencil_6(r, w); break;
            case 4: stencil_8(r, w); break;
            case 5: stencil_10(r, w); break;
        }
    }

} // end of namespace
