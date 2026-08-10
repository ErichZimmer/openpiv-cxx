#include "ducc_fft.h"

// These target-independent DUCC declarations and all standard-library headers
// are parsed once, outside the per-target namespaces.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <limits>
#include <memory>
#include <new>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <typeindex>
#include <typeinfo>
#include <utility>
#include <vector>

// Keep only target-independent DUCC infrastructure in ::ducc0.  Every
// template which depends on the SIMD value type (Cmplx, roots, plans and FFT
// kernels) is included by ducc_fft_target-inl.h inside HWY_NAMESPACE.
#include "ducc0/infra/aligned_array.h"
#include "ducc0/infra/error_handling.h"
#include "ducc0/infra/mav.h"
#include "ducc0/infra/misc_utils.h"
#include "ducc0/infra/threading.h"
#include "ducc0/infra/useful_macros.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "ducc_fft_target-inl.h"
#include "hwy/foreach_target.h"

// foreach_target handles every non-static target; this is the static target.
#include "ducc_fft_target-inl.h"

#include "hwy/highway.h"

#if HWY_ONCE
namespace ducc_fft_detail {
	HWY_EXPORT(BackendName);
	HWY_EXPORT(C2CF32);
	HWY_EXPORT(C2CF64);
	HWY_EXPORT(R2CF32);
	HWY_EXPORT(R2CF64);
	HWY_EXPORT(C2RF32);
	HWY_EXPORT(C2RF64);

}  // namespace ducc_fft_detail

namespace ducc_fft {
	namespace {
		void validate_call(
			const void* input,
			const void* output,
			const shape_t& shape, 
			const shape_t& axes
		) {
			if (input == nullptr || output == nullptr)
				throw std::invalid_argument("ducc_fft: input and output must be non-null");
			if (shape.empty())
				throw std::invalid_argument("ducc_fft: shape must not be empty");
			if (axes.empty())
				throw std::invalid_argument("ducc_fft: axes must not be empty");

			std::size_t elements = 1;
			for (const std::size_t extent : shape) {
				if (extent == 0)
					throw std::invalid_argument("ducc_fft: shape extents must be positive");
				if (elements > std::numeric_limits<std::size_t>::max() / extent)
					throw std::overflow_error("ducc_fft: shape product overflows size_t");
				elements *= extent;
			}

			for (const std::size_t axis : axes)
				if (axis >= shape.size())
					throw std::invalid_argument("ducc_fft: transform axis is out of range");
		}

	}  // namespace

	const char* backend_name() noexcept {
		return HWY_DYNAMIC_DISPATCH(ducc_fft_detail::BackendName)();
	}

	void c2c(
		const std::complex<float>* input, 
		std::complex<float>* output,
		const shape_t& shape, 
		const shape_t& axes,
		bool forward,
		float scale
	) {
		validate_call(
			input, 
			output, 
			shape, 
			axes
		);
		HWY_DYNAMIC_DISPATCH(ducc_fft_detail::C2CF32)(
			input,
			output, 
			shape, 
			axes,
			forward, 
			scale
		);
	}

	void c2c(
		const std::complex<double>* input, 
		std::complex<double>* output,
		const shape_t& shape, 
		const shape_t& axes, 
		bool forward,
		double scale
	) {
		validate_call(
			input,
			output, 
			shape, 
			axes
		);
		HWY_DYNAMIC_DISPATCH(ducc_fft_detail::C2CF64)(
			input, 
			output,
			shape, 
			axes,
			forward, 
			scale
		);
	}

	void r2c(
		const float* input, 
		std::complex<float>* output,
		const shape_t& shape, 
		const shape_t& axes, 
		bool forward,
		float scale
	) {
		validate_call(
			input, 
			output, 
			shape, 
			axes
		);
		HWY_DYNAMIC_DISPATCH(ducc_fft_detail::R2CF32)(
			input, 
			output, 
			shape, 
			axes,
			forward, 
			scale
		);
	}

	void r2c(
		const double* input, 
		std::complex<double>* output,
		const shape_t& shape, 
		const shape_t& axes, 
		bool forward,
		double scale
	) {
		validate_call(
			input, 
			output, 
			shape, 
			axes
		);
		HWY_DYNAMIC_DISPATCH(ducc_fft_detail::R2CF64)(
			input, 
			output, 
			shape, 
			axes,
			forward, 
			scale
		);
	}

	void c2r(
		const std::complex<float>* input, 
		float* output,
		const shape_t& real_shape, 
		const shape_t& axes, 
		bool forward,
		float scale
	) {
		validate_call(
			input, 
			output, 
			real_shape, 
			axes
		);
		HWY_DYNAMIC_DISPATCH(ducc_fft_detail::C2RF32)(
			input, 
			output, 
			real_shape,
			axes, 
			forward, 
			scale
		);
	}

	void c2r(
		const std::complex<double>* input, 
		double* output,
		const shape_t& real_shape, 
		const shape_t& axes, 
		bool forward,
		double scale
	) {
		validate_call(
			input, 
			output, 
			real_shape, 
			axes
		);
		HWY_DYNAMIC_DISPATCH(ducc_fft_detail::C2RF64)(
			input, 
			output, 
			real_shape,
			axes, 
			forward, 
			scale
		);
	}

}  // namespace ducc_fft
#endif  // HWY_ONCE
