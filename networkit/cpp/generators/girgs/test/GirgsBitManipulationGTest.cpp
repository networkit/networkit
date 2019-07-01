#include <algorithm>
#include <array>
#include <random>
#include <sstream>

#include <gtest/gtest.h>
#include "../BitManipulation.hpp"

namespace NetworKit {

template <typename T>
class GirgsBitManipulationGTest : public ::testing::Test {
public:
    using Implementation = T;
};

using Implementations = ::testing::Types<
#ifdef USE_BMI2
    girgs::BitManipulationDetails::BMI2::Implementation<1>,
    girgs::BitManipulationDetails::BMI2::Implementation<2>,
    girgs::BitManipulationDetails::BMI2::Implementation<3>,
    girgs::BitManipulationDetails::BMI2::Implementation<4>,
    girgs::BitManipulationDetails::BMI2::Implementation<5>,
#endif
    girgs::BitManipulationDetails::Generic::Implementation<1>,
    girgs::BitManipulationDetails::Generic::Implementation<2>,
    girgs::BitManipulationDetails::Generic::Implementation<3>,
    girgs::BitManipulationDetails::Generic::Implementation<4>,
    girgs::BitManipulationDetails::Generic::Implementation<5>
>;

TYPED_TEST_CASE(GirgsBitManipulationGTest, Implementations);

#define COMMON_DEFS \
    using Impl = typename TestFixture::Implementation; \
    constexpr auto D = Impl::kDimensions;              \
    using Coordinate = std::array<uint32_t, D>;        \
    Coordinate coord;                                  \
    std::fill_n(coord.begin(), D, 0);                  \
    std::default_random_engine prng(1 + D);

template <unsigned D>
static uint32_t ReferenceDeposite(std::array<uint32_t, D> coord) {
    uint32_t res = 0;

    for(int i = 0; i < 32; ++i) {
        const auto d = i % D;
        res |= ((coord[d] >> (i / D)) & 1) << i;
    }

    return res;
}

template <unsigned D>
static std::array<uint32_t, D> ReferenceExtract(uint32_t x) {
    std::array<uint32_t, D> res;
    std::fill_n(res.begin(), D, 0);

    for(int i = 0; i < 32; ++i) {
        const auto d = i % D;
        const auto k = i / D;
        res[d] |= ((x >> i) & 1) << k;
    }

    return res;
}


TYPED_TEST(GirgsBitManipulationGTest, DepositeSingle) {
    COMMON_DEFS

    for(unsigned d = 0; d < D; d++) {
        const auto bits = 32 / D + (d < 32u % D);
        const auto allowedBits = girgs::BitPattern<D>::kEveryDthBit << d;
        std::uniform_int_distribution<unsigned> distr(0, static_cast<uint32_t>((1llu << bits) - 1));

        for(int i=0; i < 1000; i++) {
            coord[d] = distr(prng);
            const auto tested = Impl::deposit(coord);
            const auto ref = ReferenceDeposite<D>(coord);

            ASSERT_EQ(tested, tested & allowedBits);
            ASSERT_EQ(tested, ref);
        }

        coord[d] = 0;
    }
}

TYPED_TEST(GirgsBitManipulationGTest, DepositeAll) {
    COMMON_DEFS

    std::uniform_int_distribution<unsigned> dim_distr(0, D - 1);
    for(int i=0; i < 10000; i++) {
        const auto d = dim_distr(prng);
        const auto bits = 32 / D + (d < 32 % D);
        std::uniform_int_distribution<unsigned> distr(0, static_cast<uint32_t>((1llu << bits) - 1));

        coord[d] = distr(prng);
        const auto tested = Impl::deposit(coord);
        const auto ref = ReferenceDeposite<D>(coord);

        ASSERT_EQ(tested, ref);
    }
}

TYPED_TEST(GirgsBitManipulationGTest, ExtractSingle) {
    COMMON_DEFS

    for(unsigned d = 0; d < D; d++) {
        const auto allowedBits = girgs::BitPattern<D>::kEveryDthBit << d;
        std::uniform_int_distribution<unsigned> distr;

        for(unsigned i=0; i < 1000u; i++) {
            const auto cell = distr(prng) & allowedBits;

            const auto tested = Impl::extract(cell);
            const auto ref = ReferenceExtract<D>(cell);

            for(unsigned j = 0; j < D; ++j) {
                if (j == d) continue;
                ASSERT_EQ(tested[j], 0);
            }

            ASSERT_EQ(tested, ref);
        }
    }
}

TYPED_TEST(GirgsBitManipulationGTest, ExtractAll) {
    COMMON_DEFS

    std::uniform_int_distribution<unsigned> distr;
    for(unsigned d = 0; d < D; d++) {
        for(unsigned i=0; i < 1000; i++) {
            const auto cell = distr(prng);

            const auto tested = Impl::extract(cell);
            const auto ref = ReferenceExtract<D>(cell);

            ASSERT_EQ(tested, ref);
        }
    }
}

TYPED_TEST(GirgsBitManipulationGTest, XRef) {
    COMMON_DEFS

    std::uniform_int_distribution<unsigned> distr;
    for(unsigned i=0; i < 10000u; i++) {
        const auto cell = distr(prng);

        const auto extr = Impl::extract(cell);
        const auto depo = Impl::deposit(extr);

        ASSERT_EQ(cell, depo);
    }
}

} // namespace NetworKit

