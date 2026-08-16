#include "algos/duccfft_backend.h"
#include <ducc_fft.h>

namespace openpiv::algos
{
    const char* duccfft_simd_backend() noexcept
    {
        return ducc_fft::backend_name();
    }

} // openpiv::algos