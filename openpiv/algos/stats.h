
#pragma once

// std
#include <tuple>
#include <cmath>

// local
#include "core/image.h"
#include "core/image_type_traits.h"

namespace openpiv::algos {

    using namespace core;

    template < template<typename> class ImageT,
               typename ContainedT,
               typename = typename std::enable_if_t< is_imagetype_v<ImageT<ContainedT>> >
               >
    std::tuple< ContainedT, ContainedT >
    find_image_range( const ImageT<ContainedT>& im )
    {
        ContainedT min, max;
        auto p = std::cbegin( im );
        auto e = std::cend( im );
        min = max = *p++;
        while ( p != e )
        {
            min = *p < min ? *p : min;
            max = *p > max ? *p : max;
            ++p;
        }

        return std::make_tuple( min, max );
    }

    template < template<typename> class ImageT,
               typename ContainedT,
               typename = typename std::enable_if_t<
                       is_imagetype_v<ImageT<ContainedT>> &&
                       is_real_mono_pixeltype_v<ContainedT>
                       >
              >
    std::tuple< ContainedT, ContainedT >
    find_mean_std( const ImageT<ContainedT>& im )
    {
        if (im.pixel_count() == 0)
            return std::make_tuple(ContainedT(0), ContainedT(0));

        double mean, std_temp, val;
        auto p = std::cbegin( im );
        auto e = std::cend( im );
        mean = std_temp = val = 0;

        while ( p != e )
        {
            val = static_cast<double>(*p);
            mean = mean + val;
            std_temp = std_temp + (val*val);
            ++p;
        }

        auto num_pixels = im.pixel_count();
        mean = mean / static_cast<double>(num_pixels);

        double var = (std_temp / num_pixels) - (mean * mean);

        // Guard against tiny negatives due to FP error; also convert to double for sqrt
        if (var < 0.0 && var > -std::numeric_limits<double>::epsilon())
            var = 0.0;

        double stdev = std::sqrt(std::max(0.0, var));

        return std::make_tuple( ContainedT(mean), ContainedT(stdev) );
    }

}
