#pragma once

#include <cstdint>
#include <cmath>
#include <vector>

#include "core/exception_builder.h"

namespace openpiv::interp
{
   using namespace openpiv::core;

   int32_t mirror_index(int32_t i, int32_t n);

   struct lut_v
   {
      std::vector<std::vector<double>> rows;
      double r_min;
      double r_max;
   };

   template <typename Func>
   lut_v make_lut(Func func, size_t lut_size, double r_min, double r_max)
   {
      if (lut_size == 0) 
         core::exception_builder<std::runtime_error>() << "lut_size must be non-zero";

      // Probe to determine row length
      auto sample_vec = func(r_min);

      std::vector<double> sample;
      sample.reserve(sample_vec.size());

      for (auto v : sample_vec)
         sample.push_back(static_cast<double>(v));

      size_t n = sample.size();
      if (n == 0) 
         core::exception_builder<std::runtime_error>() << "lut values must be non-zero in length";

      lut_v lut;
      lut.rows.resize(lut_size, std::vector<double>(n));
      lut.r_min = r_min;
      lut.r_max = r_max;

      for (size_t i = 0; i < lut_size; ++i)
      {
         double t = (lut_size > 1) ? (static_cast<double>(i) / static_cast<double>(lut_size - 1)) : 0.0;
         double r = r_min + t * (r_max - r_min);

         auto row_vec = func(r);

         if (row_vec.size() != n)
            core::exception_builder<std::runtime_error>() << "func must return fixed-length vectors";

         for (size_t k = 0; k < n; ++k)
            lut.rows[i][k] = static_cast<double>(row_vec[k]);
      }

      return lut;
   }

   inline std::vector<double> get_weights(double r, const lut_v& lut)
   {
      if (lut.rows.empty())
         return {};

      size_t n_rows = lut.rows.size();
      size_t row_len = lut.rows[0].size();

      double rmin = lut.r_min;
      double rmax = lut.r_max;

      // Clamp out of bounds
      if (r <= rmin) return lut.rows.front();
      if (r >= rmax) return lut.rows.back();

      // Get position in lookup table
      double pos = (r - rmin) * static_cast<double>(n_rows - 1) / (rmax - rmin);

      size_t i0 = static_cast<size_t>(std::floor(pos));
      size_t i1 = std::min(i0 + 1, n_rows - 1);

      double t = pos - static_cast<double>(i0);

      const auto& row0 = lut.rows[i0];
      const auto& row1 = lut.rows[i1];

      // linear interpolation to lower descritization error
      std::vector<double> out(row_len);
      for (size_t k = 0; k < row_len; ++k)
         out[k] = (1.0 - t) * row0[k] + t * row1[k];

      return out;
   }
}