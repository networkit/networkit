#pragma once

#include <cstdint>

namespace girgs {

template<unsigned D, typename T = uint32_t, unsigned I = 0>
struct BitPattern {
    constexpr static T setEveryDthBit(int i) {
        return (i < 8*sizeof(T)) ? (1llu << i) | setEveryDthBit(i + D) : 0;
    }

    constexpr static unsigned kBits = 8 * sizeof(T);
    constexpr static T kEveryDthBit = setEveryDthBit(0);
    constexpr static T kDBits = (kBits == D) ? T{-1} : ((T{1} << D) - 1);
};
}

// Load implementations
#include <girgs/BitManipulationGeneric.inl>

#ifdef USE_BMI2
    #include <girgs/BitManipulationBMI2.inl>
#endif

namespace girgs {
#ifdef USE_BMI2
    template <unsigned D>
    using BitManipulation = BitManipulationDetails::BMI2::Implementation<D>;
#else
    template <unsigned D>
    using BitManipulation = BitManipulationDetails::Generic::Implementation<D>;
#endif
}
