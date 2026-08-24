#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>

#include "core/image.h"
#include "core/image_type_traits.h"
#include "core/image_utils.h"
#include "core/pixel_types.h"
#include "core/size.h"

namespace openpiv::piv::detail
{

    template <
        template <typename> class ImageT,
        typename ContainedT,
        typename ValueT = typename ContainedT::value_t,
        typename = typename std::enable_if_t<
            core::is_imagetype_v<ImageT<ContainedT>> &&
            core::is_real_mono_pixeltype_v<ContainedT> &&
            std::is_floating_point<ValueT>::value
        >
    >
    void squared_integral_array(
        const ImageT<ContainedT>& input,
        std::vector<ValueT>& int_arr
    ) {
        const uint32_t width = input.width();
        const uint32_t height = input.height();
        const uint32_t stride = width + 1;

        int_arr.resize((height + 1) * stride);
        std::fill(int_arr.begin(), int_arr.end(), ValueT{0});

        for (uint32_t y = 0; y < height; ++y)
        {
            ValueT row_sum = ValueT{0};

            for (uint32_t x = 0; x < width; ++x)
            {
                const ValueT value = static_cast<ValueT>(input[{x, y}].v);

                row_sum += value * value;

                int_arr[(y + 1) * stride + (x + 1)] = int_arr[y * stride + (x + 1)] + row_sum;
            }
        }
    }

    template <typename ValueT>
    ValueT integral_area(
        const std::vector<ValueT>& int_arr,
        uint32_t stride,
        uint32_t x0,
        uint32_t y0,
        uint32_t x1,
        uint32_t y1
    ) {
        return std::max(
            ValueT{0},
              int_arr[y1 * stride + x1]
            - int_arr[y0 * stride + x1]
            - int_arr[y1 * stride + x0]
            + int_arr[y0 * stride + x0]
        );
    }

    template <
        template <typename> class ImageT,
        typename ContainedT,
        typename ValueT = typename ContainedT::value_t,
        typename = typename std::enable_if_t<
            core::is_imagetype_v<ImageT<ContainedT>> &&
            core::is_real_mono_pixeltype_v<ContainedT> &&
            std::is_floating_point<ValueT>::value
        >
    >
    ValueT squared_energy(
        const ImageT<ContainedT>& input
    ) {
        ValueT energy = ValueT{0};

        for (size_t i = 0; i < input.pixel_count(); ++i)
        {
            const ValueT value = static_cast<ValueT>(input[i].v);
            energy += value * value;
        }

        return energy;
    }

    template <
        template <typename> class ImageT,
        typename ContainedT,
        typename ValueT = typename ContainedT::value_t,
        typename = typename std::enable_if_t<
            core::is_imagetype_v<ImageT<ContainedT>> &&
            core::is_real_mono_pixeltype_v<ContainedT> &&
            std::is_floating_point<ValueT>::value
        >
    >
    void normalize_nsqecc(
        ImageT<ContainedT>& correlation,
        const ImageT<ContainedT>& interrogation,
        const std::vector<ValueT>& search_int_arr,
        const core::size& support_size
    ) {
        const uint32_t width = correlation.width();
        const uint32_t height = correlation.height();
        const uint32_t stride = width + 1;

        const ValueT interrogation_energy = squared_energy(interrogation);

        const ValueT search_energy = integral_area(
            search_int_arr,
            stride,
            0,
            0,
            width,
            height
        );

        if (
            interrogation_energy <= std::numeric_limits<ValueT>::epsilon() ||
            search_energy <= std::numeric_limits<ValueT>::epsilon()
        ) {
            core::fill(
                correlation,
                ContainedT(ValueT{0})
            );

            return;
        }

        const ValueT fft_scale = static_cast<ValueT>(correlation.pixel_count());

        const int32_t center_x = static_cast<int32_t>(width / 2);
        const int32_t center_y = static_cast<int32_t>(height / 2);

        const int32_t support_width = static_cast<int32_t>(support_size.width());
        const int32_t support_height = static_cast<int32_t>(support_size.height());

        const int32_t support_x0 = (static_cast<int32_t>(width) - support_width) / 2;
        const int32_t support_y0 = (static_cast<int32_t>(height) - support_height) / 2;

        const uint32_t valid_x0 = static_cast<uint32_t>(support_x0);
        const uint32_t valid_y0 = static_cast<uint32_t>(support_y0);
        const uint32_t valid_x1 = valid_x0 + support_size.width();
        const uint32_t valid_y1 = valid_y0 + support_size.height();

        for (uint32_t y = valid_y0; y < valid_y1; ++y)
        {
            for (uint32_t x = valid_x0; x < valid_x1; ++x)
            {
                const int32_t dx = static_cast<int32_t>(x) - center_x;
                const int32_t dy = static_cast<int32_t>(y) - center_y;

                const int32_t x0 = support_x0 + dx;
                const int32_t y0 = support_y0 + dy;
                const int32_t x1 = x0 + support_width;
                const int32_t y1 = y0 + support_height;

                const ValueT local_search_energy = integral_area(
                    search_int_arr,
                    stride,
                    static_cast<uint32_t>(x0),
                    static_cast<uint32_t>(y0),
                    static_cast<uint32_t>(x1),
                    static_cast<uint32_t>(y1)
                );

                if ( local_search_energy <= std::numeric_limits<ValueT>::epsilon() )
                {
                    correlation[{x, y}] = ContainedT(ValueT{0});

                    continue;
                }

                const ValueT cross_correlation =
                    static_cast<ValueT>(
                        correlation[{x, y}].v
                    ) / fft_scale;

                const ValueT squared_difference = std::max(
                    ValueT{0},
                    interrogation_energy
                        + local_search_energy
                        - ValueT{2} * cross_correlation
                );

                const ValueT denominator = std::sqrt(interrogation_energy * local_search_energy);
                const ValueT normalized_error = squared_difference / denominator;

                correlation[{x, y}] = ContainedT(ValueT{1} / (ValueT{1} + normalized_error));
            }
        }
    }

} // namespace openpiv::piv::detail