#pragma once

// std
#include <cstdint>
#include <vector>
#include <algorithm>

// openpiv
#include "core/image.h"
#include "core/image_type_traits.h"
#include "core/pixel_types.h"
#include "core/exception_builder.h"

#include "filter_common.impl.h"


namespace openpiv::filter{

    using namespace openpiv::core;


    template <
        template <typename> class ImageT,
        typename ContainedT,
        typename ValueT,
        typename
    >
    void convolve_2d_sep(
        const ImageT<ContainedT>& src,
        ImageT<ContainedT>& out,
        const std::vector<ValueT>& kernel_x_orig,
        const std::vector<ValueT>& kernel_y_orig
    ) {
        // NOTE: This is a 3 year old implementation from my private repo. This needs to be redone!!!

        // Validate kernel size
        if (kernel_x_orig.size() != kernel_y_orig.size())
            core::exception_builder<std::runtime_error>() << "kernels must be of same size";

        // Validate image size
        if (src.size() != out.size())
            core::exception_builder<std::runtime_error>() << "images must be the same size";

        // Image dims
        uint32_t size_x = src.width();
        uint32_t size_y = src.height();
        uint32_t numel_elem = size_x * size_y;
        uint32_t kSize = kernel_x_orig.size();

        // Convolution buffers
        std::vector<ValueT> temp_img(numel_elem, ValueT(0));
        std::vector<ValueT> temp_row(size_x, ValueT(0));

        // Get center of kernel and kernel offsets
        uint32_t kCenter = kSize >> 1;
        uint32_t endIndex = size_x - kCenter;

        // offsets for borders
        // current offset is mirror (e.g., bcd|abcd|cba)
        std::vector<uint32_t> offsets(kCenter * kSize, 0);

        for (int32_t i = 0; i < static_cast<int32_t>(kCenter); ++i)
        {
            for (int32_t j = 0; j < static_cast<int32_t>(kSize); ++j)
            {
                int32_t index = j - static_cast<int32_t>(kCenter) + i;
                
                // Set the high bounds to max image size since we only care about the left border boundary
                offsets[i * kSize + j] = mirror_index(index, size_x);
            }
        }

        // reverse kernel due to convolution theorem
        auto kernel_x = kernel_x_orig;
        auto kernel_y = kernel_y_orig;

        std::reverse(kernel_x.begin(), kernel_x.end());
        std::reverse(kernel_y.begin(), kernel_y.end());
        
        // Convolution is split into left/top, middle, and right/bottom cases
        // We start with horizontal 1D convolutions, then perform vertical 1D convolutions.
        // This allows a convolution to be performed in approximately O(2N) instead of O(N^2).

        uint32_t off = 0;

        // left border case
        for (uint32_t j = 0; j < size_y; ++j)
        {
            for (uint32_t i = 0; i < kCenter; ++i)
            {
                uint32_t idx = j * size_x + i;
                uint32_t kInd = i * kSize; // get correct offsets

                for (uint32_t k = 0; k < kSize; ++k)
                    temp_img[idx] += src[idx + offsets[kInd + k] - i] * kernel_x[k];
            }
        }

        // center case
        for (uint32_t j = 0; j < size_y; ++j)
        {
            for (uint32_t i = kCenter; i < endIndex; ++i)
            {
                uint32_t idx = j * size_x + i;

                for (uint32_t k = 0; k < kSize; k++)
                    temp_img[idx] += src[idx - kCenter + k] * kernel_x[k];
            }
        }

        // right border case
        for (uint32_t j = 0; j < size_y; ++j)
        {
            off = kCenter - 1;
            for (uint32_t i = endIndex; i < size_x; ++i)
            {
                uint32_t idx = j * size_x + i;
                uint32_t kInd = off * kSize;

                for (uint32_t k = 0, m = kSize - 1; k < kSize; ++k, --m)
                    temp_img[idx] += src[idx - offsets[kInd + m] + off] * kernel_x[k];
                
                off--;
            }
        }

        // Now perform convolution in vertical direction.
        endIndex = size_y - kCenter;

        // top border case
        for (uint32_t j = 0; j < kCenter; ++j)
        {
            uint32_t idx = 0;
            uint32_t kInd = j * kSize;

            for (uint32_t k = 0; k < kSize; ++k)
            {
                uint32_t row = offsets[kInd + k];
                
                for (uint32_t i = 0; i < size_x; ++i)
                {
                    idx = row * size_x + i;
                    temp_row[i] += temp_img[idx] * kernel_y[k];
                }
            }

            // Copy temp_row to output image
            for (uint32_t i = 0; i < size_x; ++i)
            {
                idx = j * size_x + i;
                out[idx] = temp_row[i];
                temp_row[i] = 0.f;
            }
        }

        // center case
        for (uint32_t j = kCenter; j < endIndex; ++j)
        {
            uint32_t idx = 0;
            for (uint32_t k = 0; k < kSize; ++k)
            {
                uint32_t row = j - kCenter + k;
                for (uint32_t i = 0; i < size_x; ++i)
                {
                    idx = row * size_x + i;
                    temp_row[i] += temp_img[idx] * kernel_y[k];
                }
            }

            for (uint32_t i = 0; i < size_x; ++i)
            {
                idx = j * size_x + i;
                out[idx] = temp_row[i];
                temp_row[i] = 0.f;
            }
        }

        // bottom border case
        off = kCenter - 1;
        for (uint32_t j = endIndex; j < size_y; ++j)
        {
            uint32_t idx = 0;
            uint32_t kInd = off * kSize;

            for (uint32_t k = 0, m = kSize - 1; k < kSize; ++k, --m)
            {
                uint32_t row = size_y - offsets[kInd + m] - 1;
                for (uint32_t i = 0; i < size_x; ++i)
                {
                    idx = row * size_x + i;
                    temp_row[i] += temp_img[idx] * kernel_y[k];
                }
            }
            off--;

            // Copy temp_row to output image
            for (uint32_t i = 0; i < size_x; ++i)
            {
                idx = j * size_x + i;
                out[idx] = temp_row[i];
                temp_row[i] = 0.f;
            }
        }
    }

} // end of namespace