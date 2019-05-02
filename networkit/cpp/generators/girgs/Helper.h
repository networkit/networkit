#pragma once

namespace girgs {

// build multiplication stack
namespace detail {
    template<unsigned D, typename T>
    struct pow_helper;

    template<typename T>
    struct pow_helper<0, T> {
        static constexpr T pow(const T) { return {1.0}; }
    };

    template<typename T>
    struct pow_helper<1, T> {
        static constexpr T pow(const T x) { return x; }
    };

    template<typename T>
    struct pow_helper<2, T> {
        static constexpr T pow(const T x) { return x * x; }
    };

    template<unsigned D, typename T>
    struct pow_helper {
        static constexpr T pow(const T x) {
            return pow_helper<D / 2, T>::pow(x) * pow_helper<(D + 1) / 2, T>::pow(x);
        }
    };
}

/// computes x to the power of D where D is a compile-time constant
template <unsigned D, typename T>
constexpr T pow_to_the(const T x) {
    return detail::pow_helper<D, T>::pow(x);
}

}
