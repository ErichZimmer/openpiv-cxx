#include "map_polynomial.impl.h"

#include <vector>

namespace openpiv::interp
{
    // TODO: optimize using lookup tables (LUT) descritized to 256 with linear interp
    std::vector<double> poly_weights(double r, int k)
    {
        const int n = 2 * k;

        std::vector<double> w(n);

        for (int i = 0; i < n; ++i)
        {
            double xi = i - (k - 0.5);

            double wi = 1.0;

            for (int j = 0; j < n; ++j)
            {
                if (j == i) continue;

                double xj = j - (k - 0.5);

                wi *= (r - xj) / (xi - xj);
            }

            w[i] = wi;
        }

        return w;
    }

} // end of namespace