
#pragma once

// std
#include <tuple>
#include <cmath>
#include <algorithm>

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
    ContainedT
    find_mean( const ImageT<ContainedT>& im )
    {
        if (im.pixel_count() == 0)
            return ContainedT(0);

        double mean;
        auto p = std::cbegin( im );
        auto e = std::cend( im );
        mean = 0;

        while ( p != e )
        {
            mean = mean + static_cast<double>(*p);
            ++p;
        }

        auto num_pixels = im.pixel_count();

        mean = mean / static_cast<double>(num_pixels);

        return ContainedT(mean);
    }


    template < 
        template<typename> class ImageT,
        typename ContainedT,
        typename ValueT = typename ContainedT::value_t,
        typename = typename std::enable_if_t<
            is_imagetype_v<ImageT<ContainedT>> &&
            is_real_mono_pixeltype_v<ContainedT> &&
            std::is_floating_point_v<ValueT>
        >
    >
    std::tuple< ContainedT, ContainedT >
    find_nanmean_nanstd( const ImageT<ContainedT>& im )
    {
        if (im.pixel_count() == 0)
            return std::make_tuple(ContainedT(0), ContainedT(0));

        ValueT mean, std_temp, val;
        auto p = std::cbegin( im );
        auto e = std::cend( im );
        mean = std_temp = val = 0;

        while ( p != e )
        {
            val = static_cast<ValueT>(*p);

            if (std::isfinite(val))
            {
                mean = mean + val;
                std_temp = std_temp + (val*val);
            }

            ++p;
        }

        auto num_pixels = im.pixel_count();
        mean = mean / static_cast<ValueT>(num_pixels);

        ValueT var = (std_temp / num_pixels) - (mean * mean);

        // Guard against tiny negatives due to FP error
        if (var < 0.0 && var > -std::numeric_limits<ValueT>::epsilon())
            var = 0.0;

        ValueT stdev = std::sqrt(std::max(0.0, var));

        return std::make_tuple( ContainedT(mean), ContainedT(stdev) );
    }

    template < 
        template<typename> class ImageT,
        typename ContainedT,
        typename ValueT = typename ContainedT::value_t,
        typename = typename std::enable_if_t<
            is_imagetype_v<ImageT<ContainedT>> &&
            is_real_mono_pixeltype_v<ContainedT> &&
            std::is_floating_point_v<ValueT>
        >
    >
    bool 
    has_nan( const ImageT<ContainedT>& im )
    {
        auto p = std::cbegin( im );
        auto e = std::cend( im );
        ValueT val{ 0 };

        while ( p != e )
        {
            val = static_cast<ValueT>(*p);
            if ( !std::isfinite(val) )
                return true;

            ++p;
        }

        return false;
    }


    template < 
        template<typename> class ImageT,
        typename ContainedT,
        typename ValueT = typename ContainedT::value_t,
        typename ResultT = core::image<g_8>,
        typename = typename std::enable_if_t<
            is_imagetype_v<ImageT<ContainedT>> &&
            is_real_mono_pixeltype_v<ContainedT> &&
            std::is_floating_point_v<ValueT>
        >
    >
    ResultT is_nan( const ImageT<ContainedT>& im )
    {
        ResultT nan_flag{ im.size() };
        
        for (uint32_t i=0; i<im.pixel_count(); i++)
        {
            if ( !std::isfinite(im[i]) )
                nan_flag[i] = 1;
        }

        return nan_flag;
    }

    // Median of 1D vector
    // Note, this is pass by copy since we directly manipulate arr
    template<typename ValueT>
    ValueT median(std::vector<ValueT> arr)
    {
        if (arr.empty())
            return ValueT(0);
        if (arr.size() == 1)
            return arr[0];

        // Quick sort to find middle value
        std::size_t n = arr.size()/2;
        std::nth_element(arr.begin(), arr.begin()+n, arr.end());
        ValueT mid = arr[n];

        // For odd, we can simply select the middle
        if (arr.size()%2==1)
            return mid;

        // Get average of left side maximum and previous middle for even-sized vectors
        ValueT left_max = *std::max_element(arr.begin(), arr.begin()+n);

        return ValueT(0.5) * (mid + left_max);
    }

    // Wrapper for core::image
    template <
        template<typename> class ImageT,
        typename ContainedT,
        typename ValueT = typename ContainedT::value_t,
        typename = typename std::enable_if_t< is_imagetype_v<ImageT<ContainedT>> >
    >
    ValueT median(const ImageT<ContainedT>& im)
    {
        const auto n = im.pixel_count();

        if (n == 0)
            return ContainedT(0);
        if (n == 1)
            return im[0];

        std::vector<ValueT> vals;
        vals.reserve(n);

        for (std::size_t i = 0; i < n; ++i)
            vals.push_back(static_cast<ValueT>(im[i]));

        return ContainedT( median(std::move(vals)) ); // Call the std::vector implementation
    }

}
