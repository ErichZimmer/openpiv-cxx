#pragma once

#include <cstdint>
#include <cmath>
#include <numeric>
#include <vector>


namespace openpiv::filter
{

    template <
        typename ContainedT,
        typename ValueT,
        typename
    >
    std::vector<ValueT> generate_gaussian_kernel1d(
         ValueT sigma,
         ValueT truncate
    ) {
        // Calculate kernel size
        size_t kernel_half_size = static_cast<size_t>(sigma * truncate);
        size_t kernel_full_size = 2*kernel_half_size + 1;

        std::vector<ValueT> weights(kernel_full_size);

        // Calculate 1D weights
        ValueT sigma2 = sigma * sigma;

        for (size_t i=0; i < kernel_full_size; i++)
        {
            ValueT x = static_cast<ValueT>(i) - kernel_half_size;

            weights[i] = 2 * std::exp(-0.5 * (x*x) / sigma2);
        }

        // Normalize the weights
        ValueT sum = std::accumulate(weights.begin(), weights.end(), ValueT(0));

        for (auto& val : weights)
        {
            val /= sum;
        }
        
        return weights;
    }

} // end of namespace