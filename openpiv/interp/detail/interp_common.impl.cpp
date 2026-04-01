#include "interp_common.impl.h"

#include <cstdint>
#include <cstdlib>


namespace openpiv::interp
{
    int32_t mirror_index(int32_t i, int32_t n)
    {
        const int32_t period = 2*n - 2;   // reflection period
        int32_t x = i % period;           // wrap into period
        x += (x < 0) * period;             // fix negative mod (branchless)

        return n - 1 - std::abs(x - (n - 1));
    }
}