#pragma once

#include <cstdint>
#include <cmath>
#include <numeric>
#include <vector>

#include "core/exception_builder.h"


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
        const uint32_t MAX_RADIUS = 128u; 

        // Calculate kernel size
        size_t kernel_half_size = static_cast<size_t>(sigma * truncate);
        size_t kernel_full_size = 2*kernel_half_size + 1;

        // Make sure kernel_half_size is sane
        if (kernel_half_size > MAX_RADIUS)
            core::exception_builder<std::runtime_error>() << "Large kernel size found indicating possible overflow";

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


    template <
        typename ContainedT,
        typename ValueT,
        typename
    >
    std::vector<ValueT> generate_gaussian_kernel1d(
         uint32_t kernel_half_size
    ) {
        const uint32_t MAX_RADIUS = 128u; 

        // Make sure kernel_half_size is sane
        if (kernel_half_size > MAX_RADIUS)
            core::exception_builder<std::runtime_error>() << "Large kernel size found indicating possible overflow";

        // Calculate kernel size
        uint32_t kernel_full_size = 2*kernel_half_size + 1;

        std::vector<ValueT> weights(kernel_full_size);

        // Calculate 1D weights
        for (uint32_t i=0; i < kernel_full_size; i++)
        {
            ValueT x = static_cast<ValueT>(i) - static_cast<ValueT>(kernel_half_size);

            weights[i] = std::exp( -((x*x) / static_cast<ValueT>(kernel_half_size)) );
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