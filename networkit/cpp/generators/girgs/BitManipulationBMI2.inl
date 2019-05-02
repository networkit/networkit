#pragma once

#include <immintrin.h>

namespace girgs {
namespace BitManipulationDetails {
namespace BMI2 {
template <unsigned D>
struct Implementation {
    static constexpr unsigned kDimensions = D;

    static uint32_t deposit(const std::array<uint32_t, D>& coords) noexcept {
        uint32_t result = 0;

        constexpr auto mask = BitPattern<D, uint32_t>::kEveryDthBit;
        for(unsigned i=0; i < D; ++i) {
            result |= _pdep_u32(coords[i], mask << i);
        }

        return result;
    }

    static std::array<uint32_t, kDimensions> extract(uint32_t cell) noexcept {
        std::array<uint32_t, D> result;

        for(int i = 0; i < D; ++i)
            result[i] = _pext_u32(cell, BitPattern<D, uint32_t>::kEveryDthBit << i);

        return result;
    }

    static std::string name() {
        return "BMI2";
    }
};

} // namespace BMI2
} // namespace BitManipulationDetails
} // namespace girgs
