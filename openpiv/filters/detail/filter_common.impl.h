#pragma once

#include <cstdint>
#include <cmath>
#include <vector>

#include "core/exception_builder.h"

namespace openpiv::filter
{
   using namespace openpiv::core;

   int32_t mirror_index(int32_t i, int32_t n);

}