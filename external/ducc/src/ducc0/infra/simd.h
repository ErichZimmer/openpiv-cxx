/* SPDX-License-Identifier: BSD-3-Clause OR GPL-2.0-or-later */

// Highway re-includes target-specific headers while toggling
// HWY_TARGET_TOGGLE. Match that convention so this header can participate in
// foreach_target compilation as well as an ordinary single-target build.
#if defined(DUCC0_SIMD_H) == defined(HWY_TARGET_TOGGLE)
#ifdef DUCC0_SIMD_H
#undef DUCC0_SIMD_H
#else
#define DUCC0_SIMD_H
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <type_traits>

#include "hwy/highway.h"

namespace ducc0 {
namespace detail_simd {

namespace hn = hwy::HWY_NAMESPACE;

struct element_aligned_tag {};

template <typename T>
constexpr inline bool vectorizable = false;

#if !defined(DUCC0_NO_SIMD) && HWY_TARGET != HWY_SCALAR
template <>
constexpr inline bool vectorizable<float> = true;
#if HWY_HAVE_FLOAT64
template <>
constexpr inline bool vectorizable<double> = true;
#endif
#endif

template <typename T>
constexpr std::size_t native_lane_count()
{
	if constexpr (!vectorizable<T>) {
		return 1;
	} else {
		// DUCC's FFT kernels intentionally top out at AVX2-sized batches. The
		// fixed cap is also suitable for NEON and SVE targets.
		constexpr std::size_t cap = std::is_same_v<T, double> ? 4 : 8;
		constexpr std::size_t available = HWY_LANES(T);

		return available < cap ? available : cap;
	}
}

template <typename T, std::size_t len>
constexpr inline bool simd_exists =
	vectorizable<T> && len > 1 && ((len & (len - 1)) == 0) &&
	len <= native_lane_count<T>();

template <typename T, std::size_t len>
class vtp;

template <typename T, std::size_t len>
class vmask_ 
{
	private:
	// Highway documents that LoadMaskBits may read at least eight bytes even
	// when only one byte contains useful bits (notably on SVE).
	static constexpr std::size_t mask_bytes =
		((len + 7) / 8) < 8 ? 8 : ((len + 7) / 8);
	std::array<std::uint8_t, mask_bytes> bits_{};

	constexpr bool get(std::size_t i) const 
	{
		return (bits_[i / 8] & static_cast<std::uint8_t>(1u << (i & 7))) != 0;
	}

	constexpr void set(std::size_t i, bool value) 
	{
		const auto bit = static_cast<std::uint8_t>(1u << (i & 7));
		if (value)
			bits_[i / 8] |= bit;
		else
			bits_[i / 8] &= static_cast<std::uint8_t>(~bit);
	}

	template <typename, std::size_t>
	friend class vtp;

	public:
	vmask_() = default;
	vmask_(const vmask_&) = default;
	vmask_& operator=(const vmask_&) = default;

	bool none() const 
	{
		for (std::size_t i = 0; i < len; ++i)
		{
			if (get(i))
				return false;
		}

		return true;
	}

	bool any() const { 
		return !none(); 
	}

	bool all() const 
	{
		for (std::size_t i = 0; i < len; ++i)
		{
			if (!get(i))
				return false;
		}

		return true;
	}

	vmask_ operator&(const vmask_& other) const 
	{
		vmask_ result;
		for (std::size_t i = 0; i < mask_bytes; ++i)
			result.bits_[i] = static_cast<std::uint8_t>(bits_[i] & other.bits_[i]);

		return result;
	}

	vmask_& operator&=(const vmask_& other) 
	{
		for (std::size_t i = 0; i < mask_bytes; ++i)
			bits_[i] &= other.bits_[i];

		return *this;
	}

	vmask_ operator|(const vmask_& other) const 
	{
		vmask_ result;
		for (std::size_t i = 0; i < mask_bytes; ++i)
			result.bits_[i] = static_cast<std::uint8_t>(bits_[i] | other.bits_[i]);

		return result;
	}

	vmask_& operator|=(const vmask_& other) 
	{
		for (std::size_t i = 0; i < mask_bytes; ++i)
			bits_[i] |= other.bits_[i];

		return *this;
	}
};

template <typename T, std::size_t len>
class alignas(sizeof(T) * len) vtp 
{
	static_assert(len != 0, "a SIMD value must contain at least one lane");
	static_assert((len & (len - 1)) == 0,
		"DUCC SIMD lane counts must be powers of two");
	static_assert((sizeof(T) * len & (sizeof(T) * len - 1)) == 0,
		"DUCC SIMD storage size must be a power of two");

	public:
	using value_type = T;
	using Tm = vmask_<T, len>;

	static constexpr std::size_t size()
	{
		return len;
	}

	private:
	std::array<T, len> lanes_;

	template <class HighwayOp, class ScalarOp>
	vtp binary(const vtp& other, HighwayOp hop, ScalarOp sop) const 
	{
		vtp result;
		if constexpr (simd_exists<T, len>) {
			const hn::CappedTag<T, len> d;
			const std::size_t step = hn::Lanes(d);
			for (std::size_t i = 0; i < len; i += step) {
				const auto a = hn::LoadU(d, lanes_.data() + i);
				const auto b = hn::LoadU(d, other.lanes_.data() + i);
				hn::StoreU(hop(a, b), d, result.lanes_.data() + i);
			}
		} else {
			for (std::size_t i = 0; i < len; ++i)
				result.lanes_[i] = sop(lanes_[i], other.lanes_[i]);
		}
		return result;
	}

	template <class HighwayOp, class ScalarOp>
	vtp unary(HighwayOp hop, ScalarOp sop) const 
	{
		vtp result;
		if constexpr (simd_exists<T, len>) {
			const hn::CappedTag<T, len> d;
			const std::size_t step = hn::Lanes(d);
			for (std::size_t i = 0; i < len; i += step) {
				const auto value = hn::LoadU(d, lanes_.data() + i);
				hn::StoreU(hop(value), d, result.lanes_.data() + i);
			}
		} else {
			for (std::size_t i = 0; i < len; ++i)
				result.lanes_[i] = sop(lanes_[i]);
		}
		return result;
	}

	template <class HighwayOp, class ScalarOp>
	Tm compare(const vtp& other, HighwayOp hop, ScalarOp sop) const 
	{
		Tm result;
		if constexpr (simd_exists<T, len>) {
			const hn::CappedTag<T, len> d;
			const std::size_t step = hn::Lanes(d);
			for (std::size_t i = 0; i < len; i += step) {
				std::array<std::uint8_t, 8> chunk_bits{};
				const auto a = hn::LoadU(d, lanes_.data() + i);
				const auto b = hn::LoadU(d, other.lanes_.data() + i);
				(void)hn::StoreMaskBits(d, hop(a, b), chunk_bits.data());
				for (std::size_t lane = 0; lane < step; ++lane) {
					const bool value =
						(chunk_bits[lane / 8] &
						static_cast<std::uint8_t>(1u << (lane & 7))) != 0;
					result.set(i + lane, value);
				}
			}
		} else {
			for (std::size_t i = 0; i < len; ++i)
				result.set(i, sop(lanes_[i], other.lanes_[i]));
		}
		return result;
	}

	static vtp select(const Tm& mask, const vtp& yes, const vtp& no) 
	{
		vtp result;
		if constexpr (simd_exists<T, len>) 
		{
			const hn::CappedTag<T, len> d;
			const std::size_t step = hn::Lanes(d);
			for (std::size_t i = 0; i < len; i += step) 
			{
				std::array<std::uint8_t, 8> chunk_bits{};
				for (std::size_t lane = 0; lane < step; ++lane) 
				{
					if (mask.get(i + lane))
						chunk_bits[lane / 8] |=
							static_cast<std::uint8_t>(1u << (lane & 7));
				}
				const auto m = hn::LoadMaskBits(d, chunk_bits.data());
				const auto y = hn::LoadU(d, yes.lanes_.data() + i);
				const auto n = hn::LoadU(d, no.lanes_.data() + i);
				hn::StoreU(hn::IfThenElse(m, y, n), d, result.lanes_.data() + i);
			}
		} 
		else 
		{
			for (std::size_t i = 0; i < len; ++i)
				result.lanes_[i] = mask.get(i) ? yes.lanes_[i] : no.lanes_[i];
		}
		return result;
	}

	public:
	vtp() = default;
	vtp(const vtp&) = default;
	vtp& operator=(const vtp&) = default;

	vtp(T value) 
	{
		if constexpr (simd_exists<T, len>) 
		{
			const hn::CappedTag<T, len> d;
			const std::size_t step = hn::Lanes(d);
			const auto values = hn::Set(d, value);
			for (std::size_t i = 0; i < len; i += step)
				hn::StoreU(values, d, lanes_.data() + i);
		} 
		else 
		{
			lanes_.fill(value);
		}
	}

	vtp& operator=(T value) 
	{
		if constexpr (simd_exists<T, len>) 
		{
			const hn::CappedTag<T, len> d;
			const std::size_t step = hn::Lanes(d);
			const auto values = hn::Set(d, value);
			for (std::size_t i = 0; i < len; i += step)
				hn::StoreU(values, d, lanes_.data() + i);
		} 
		else 
		{
			lanes_.fill(value);
		}

		return *this;
	}

	vtp(const T* ptr, element_aligned_tag) 
	{
		std::memcpy(lanes_.data(), ptr, sizeof(T) * len);
	}

	void copy_to(T* ptr, element_aligned_tag) const 
	{
		std::memcpy(ptr, lanes_.data(), sizeof(T) * len);
	}

	vtp operator-() const 
	{
		return unary([](auto v) { return hn::Neg(v); },
			[](T v) { return static_cast<T>(-v); });
	}

	vtp operator+(const vtp& other) const 
	{
		return binary(other, [](auto a, auto b) { return hn::Add(a, b); },
			[](T a, T b) { return static_cast<T>(a + b); });
	}

	vtp operator-(const vtp& other) const 
	{
		return binary(other, [](auto a, auto b) { return hn::Sub(a, b); },
			[](T a, T b) { return static_cast<T>(a - b); });
	}

	vtp operator*(const vtp& other) const 
	{
		return binary(other, [](auto a, auto b) { return hn::Mul(a, b); },
			[](T a, T b) { return static_cast<T>(a * b); });
	}

	vtp operator/(const vtp& other) const 
	{
		return binary(other, [](auto a, auto b) { return hn::Div(a, b); },
			[](T a, T b) { return static_cast<T>(a / b); });
	}

	vtp& operator+=(const vtp& other) { return *this = *this + other; }
	vtp& operator-=(const vtp& other) { return *this = *this - other; }
	vtp& operator*=(const vtp& other) { return *this = *this * other; }
	vtp& operator/=(const vtp& other) { return *this = *this / other; }

	vtp abs() const 
	{
		return unary([](auto v) { return hn::Abs(v); },
			[](T v) {
				using std::abs;
				return static_cast<T>(abs(v));
			});
	}

	vtp sqrt() const 
	{
		return unary([](auto v) { return hn::Sqrt(v); },
			[](T v) {
				using std::sqrt;
				return static_cast<T>(sqrt(v));
			});
	}

	vtp max(const vtp& other) const 
	{
		return binary(other, [](auto a, auto b) { return hn::Max(a, b); },
			[](T a, T b) { return std::max(a, b); });
	}

	vtp min(const vtp& other) const 
	{
		return binary(other, [](auto a, auto b) { return hn::Min(a, b); },
			[](T a, T b) { return std::min(a, b); });
	}

	Tm operator>(const vtp& other) const 
	{
		return compare(other, [](auto a, auto b) { return hn::Gt(a, b); },
			[](T a, T b) { return a > b; });
	}

	Tm operator>=(const vtp& other) const 
	{
		return compare(other, [](auto a, auto b) { return hn::Ge(a, b); },
			[](T a, T b) { return a >= b; });
	}

	Tm operator<(const vtp& other) const 
	{
		return compare(other, [](auto a, auto b) { return hn::Lt(a, b); },
			[](T a, T b) { return a < b; });
	}

	Tm operator<=(const vtp& other) const 
	{
		return compare(other, [](auto a, auto b) { return hn::Le(a, b); },
			[](T a, T b) { return a <= b; });
	}

	Tm operator==(const vtp& other) const 
	{
		return compare(other, [](auto a, auto b) { return hn::Eq(a, b); },
			[](T a, T b) { return a == b; });
	}

	Tm operator!=(const vtp& other) const 
	{
		return compare(other, [](auto a, auto b) { return hn::Ne(a, b); },
			[](T a, T b) { return a != b; });
	}

	static vtp blend(const Tm& mask, const vtp& a, const vtp& b) {
		return select(mask, a, b);
	}

	class reference 
	{
		private:
			vtp& value_;
			std::size_t index_;

		public:
			reference(vtp& value, std::size_t index)
				: value_(value), index_(index) {}

			reference& operator=(T value) {
				value_.lanes_[index_] = value;
				return *this;
			}

			reference& operator*=(T value) {
				value_.lanes_[index_] *= value;
				return *this;
			}

			operator T() const { return value_.lanes_[index_]; }
	};

	void Set(std::size_t i, T value) { lanes_[i] = value; }
	reference operator[](std::size_t i) { return reference(*this, i); }
	T operator[](std::size_t i) const { return lanes_[i]; }

	class where_expr 
	{
		private:
			vtp& value_;
			Tm mask_;

		public:
			where_expr(Tm mask, vtp& value) : value_(value), mask_(mask) {}

			where_expr& operator=(const vtp& other) 
			{
				value_ = select(mask_, other, value_);
				return *this;
			}

			where_expr& operator*=(const vtp& other) 
			{
				value_ = select(mask_, value_ * other, value_);
				return *this;
			}

			where_expr& operator+=(const vtp& other) 
			{
				value_ = select(mask_, value_ + other, value_);
				return *this;
			}

			where_expr& operator-=(const vtp& other) 
			{
				value_ = select(mask_, value_ - other, value_);
				return *this;
			}
	};
};

template <typename T, std::size_t len>
inline vtp<T, len> abs(vtp<T, len> value) 
{
	return value.abs();
}

template <typename T, std::size_t len>
inline typename vtp<T, len>::where_expr where(
	typename vtp<T, len>::Tm mask, 
	vtp<T, len>& value
) {
	return typename vtp<T, len>::where_expr(mask, value);
}

template <typename T0, typename T, std::size_t len>
inline vtp<T, len> operator*(
	T0 a, 
	vtp<T, len> b
) {
	return b * a;
}

template <typename T, std::size_t len>
inline vtp<T, len> operator+(
	T a, 
	vtp<T, len> b
) {
	return b + a;
}

template <typename T, std::size_t len>
inline vtp<T, len> operator-(
	T a, 
	vtp<T, len> b
) {
	return vtp<T, len>(a) - b;
}

template <typename T, std::size_t len>
inline vtp<T, len> max(
	vtp<T, len> a, 
	vtp<T, len> b
) {
	return a.max(b);
}

template <typename T, std::size_t len>
inline vtp<T, len> min(
	vtp<T, len> a, 
	vtp<T, len> b
) {
	return a.min(b);
}

template <typename T, std::size_t len>
inline vtp<T, len> sqrt(vtp<T, len> value) 
{
	return value.sqrt();
}

template <typename T, std::size_t len>
inline bool none_of(const vmask_<T, len>& mask) 
{
	return mask.none();
}

template <typename T, std::size_t len>
inline bool any_of(const vmask_<T, len>& mask) 
{
	return mask.any();
}

template <typename T, std::size_t len>
inline bool all_of(const vmask_<T, len>& mask) 
{
	return mask.all();
}

template <typename T, std::size_t len>
inline vtp<T, len> blend(
	const vmask_<T, len>& mask, 
	const vtp<T, len>& a,
	const vtp<T, len>& b
) {
	return vtp<T, len>::blend(mask, a, b);
}

template <typename Op, typename T, std::size_t len>
inline T reduce(
	const vtp<T, len>& value, 
	Op op
) {
	T result = value[0];
	for (std::size_t i = 1; i < len; ++i)
		result = op(result, value[i]);

	return result;
}

template <typename Func, typename T, std::size_t len>
inline vtp<T, len> apply(
	vtp<T, len> input,
	Func func
) {
	vtp<T, len> result;
	for (std::size_t i = 0; i < len; ++i)
		result[i] = func(input[i]);

	return result;
}

template <typename T>
using native_simd = vtp<T, native_lane_count<T>()>;

template <typename T, int len>
struct simd_select 
{
	using type = vtp<T, static_cast<std::size_t>(len)>;
};

template <typename T, std::size_t len>
inline vtp<T, len> sin(vtp<T, len> input) 
{
	return apply(input, [](T value) 
	{
		using std::sin;
		return static_cast<T>(sin(value));
	});
}

template <typename T, std::size_t len>
inline vtp<T, len> cos(vtp<T, len> input)
{
	return apply(input, [](T value) 
	{
		using std::cos;
		return static_cast<T>(cos(value));
	});
}

template <typename Tsimd>
inline Tsimd loadu(const typename Tsimd::value_type* ptr) 
{
	return Tsimd(ptr, element_aligned_tag());
}

template <typename Tsimd>
inline void storeu(
	Tsimd value, 
	typename Tsimd::value_type* ptr
) {
	value.copy_to(ptr, element_aligned_tag());
}

}  // namespace detail_simd

using detail_simd::all_of;
using detail_simd::any_of;
using detail_simd::blend;
using detail_simd::element_aligned_tag;
using detail_simd::loadu;
using detail_simd::native_simd;
using detail_simd::none_of;
using detail_simd::simd_exists;
using detail_simd::simd_select;
using detail_simd::storeu;
using detail_simd::vectorizable;
using detail_simd::where;

template <typename Tsimd>
inline void unaligned_add(
	typename Tsimd::value_type* ptr, 
	Tsimd value
) {
	storeu(loadu<Tsimd>(ptr) + value, ptr);
}

}  // namespace ducc0

#endif  // DUCC0_SIMD_H toggle