#pragma once
#include <array>
#include <cassert>
#include <cstdint>

#ifdef USE_BMI2
#include <immintrin.h>
#endif    


namespace girgs {
namespace BitManipulationDetails {
namespace Generic {

template <unsigned Bits, unsigned Space, unsigned Block>
struct ExtractionHelper {
private:
    static_assert( (Block + Space) * (Bits / Block) - Space + (Bits % Block ? (Space + Bits % Block) : 0) <= 32,
                   "Pattern does not fit into 32 bits");

    constexpr static uint32_t nBits(int n) {
        return n <= 0 ? 0 : ((1u << n) - 1);
    }

    constexpr static uint32_t compile_pattern(unsigned i = 0, unsigned b = Bits) {
        return i > 31 ? 0 : (
            (b ? nBits(b < Block ? b : Block) << i : 0) // own contribution
            | compile_pattern(i + Block + Space, b > Block ? b - Block : 0) // recursion
        );
    }

public:
    static constexpr uint32_t kPattern32 = compile_pattern();
    static constexpr uint64_t kPattern64 = (static_cast<uint64_t>(kPattern32) << 32) | kPattern32;

    template<typename T>
    static constexpr auto shift(const T x) -> typename std::enable_if<sizeof(T) == 4, uint32_t>::type {
        return recurse<uint32_t>(x & kPattern32);
    }

    template<typename T>
    static constexpr auto shift(const T x) -> typename std::enable_if<sizeof(T) == 8, uint64_t>::type {
        return recurse<uint64_t>(x & kPattern64);
    }

private:
    static constexpr bool kDoneAfter = (Block >= Bits);

    template<typename T>
    static constexpr auto recurse(T x) -> typename std::enable_if<kDoneAfter, T>::type {
        return x;
    }

    template<typename T>
    static constexpr auto recurse(T x) -> typename std::enable_if<!kDoneAfter, T>::type {
        return ExtractionHelper<Bits, 2*Space, 2*Block>::template shift<T>( x | (x >> Space) );
    }
};

static_assert(ExtractionHelper<16, 1, 1>::kPattern32 == 0x55555555, "Unittest failed");
static_assert(ExtractionHelper<15, 1, 1>::kPattern32 == 0x15555555, "Unittest failed");
static_assert(ExtractionHelper<11, 2, 1>::kPattern32 == 0x49249249, "Unittest failed");
static_assert(ExtractionHelper<8,  3, 1>::kPattern32 == 0x11111111, "Unittest failed");
static_assert(ExtractionHelper<16, 2, 2>::kPattern32 == 0x33333333, "Unittest failed");
static_assert(ExtractionHelper<11, 4, 2>::kPattern32 == 0x430c30c3, "Unittest failed");

static_assert(ExtractionHelper<16, 1, 1>::shift(0x55555555u) == 0xffff, "Unittest failed");
static_assert(ExtractionHelper<16, 1, 1>::shift(0xaaaaaaaau) == 0x0000, "Unittest failed");
static_assert(ExtractionHelper<16, 1, 1>::shift(0x0000ffffu) == 0x00ff, "Unittest failed");
static_assert(ExtractionHelper<16, 1, 1>::shift(0x5555555500000000u) == 0xffff00000000, "Unittest failed");
static_assert(ExtractionHelper<16, 1, 1>::shift(0xaaaaaaaa00000000u) == 0x000000000000, "Unittest failed");
static_assert(ExtractionHelper<16, 1, 1>::shift(0x0000ffff00000000u) == 0x00ff00000000, "Unittest failed");

// Generic Implementation of Extractting
template <size_t D>
struct Extract {
    static std::array<uint32_t, D> extract(uint32_t cell) {
        std::array<uint32_t, D> result;

        using Extractor = ExtractionHelper<(32 + D - 1) / D, D-1, 1>;

        if (D % 2) {
            // if D is odd, shifts are slightly more complicated, hence
            // we cannot process two coordinates in parallel
            for(int i = 0; i < D; ++i)
                result[i] = Extractor::shift(cell >> i);
        } else {
            for(int i = 0; i < D - 1; i += 2) {
                const auto twoWords = (cell >> i) | (static_cast<uint64_t>(cell >> (i + 1)) << 32);
                const auto shifted = Extractor::shift(twoWords);
                result[i  ] = shifted & 0xffffffff;
                result[i+1] = shifted >> 32;
            }
        }

        return result;
    }
};

template <>
struct Extract<1> {
    static std::array<uint32_t, 1> extract(uint32_t cell) {
        std::array<uint32_t, 1> result;
        result[0] = cell;
        return result;
    }
};

// Generic Implementation of Interleaving (special case for certain D below)
template <size_t D>
struct Deposit {
    static uint32_t deposit(const std::array<uint32_t, D>& coords) {
        if (D == 1)
            return coords.front();

        unsigned int result = 0u;
        unsigned int bit = 0;

        for(auto l = 0u; l*D < 32 + D; l++) {
            for(auto d = 0u; d != D; d++) {
                result |= ((coords[d] >> l) & 1) << bit++;
            }
        }

        return result;
    }
};

template <>
struct Deposit<2> {
    static uint32_t deposit(const std::array<uint32_t, 2>& coords) {
#ifndef NDEBUG
        for(auto x : coords) assert(x <= 0xffff);
#endif

        uint64_t z = coords[0] | (static_cast<uint64_t>(coords[1]) << 32);

        z = (z | (z << 8)) & 0x00FF00FF00FF00FF;
        z = (z | (z << 4)) & 0x0F0F0F0F0F0F0F0F;
        z = (z | (z << 2)) & 0x3333333333333333;
        z = (z | (z << 1)) & 0x5555555555555555;

        z = z | (z >> 31);

        return static_cast<uint32_t>(z);
    }
};

template <>
struct Deposit<3> {
    static uint32_t deposit(const std::array<uint32_t, 3>& coords) {
        assert(coords[0] <= 0x7ff);
        assert(coords[1] <= 0x7ff);
        assert(coords[2] <= 0x3ff);

        uint32_t x = coords[0];
        uint32_t y = coords[1];
        uint32_t z = coords[2];

        x = (x | x << 16) & 0x070000ff;
        y = (y | y << 16) & 0x070000ff;
        z = (z | z << 16) & 0x070000ff;
        x = (x | x <<  8) & 0x0700f00f;
        y = (y | y <<  8) & 0x0700f00f;
        z = (z | z <<  8) & 0x0700f00f;
        x = (x | x <<  4) & 0xc30c30c3;
        y = (y | y <<  4) & 0xc30c30c3;
        z = (z | z <<  4) & 0xc30c30c3;
        x = (x | x <<  2) & 0x49249249;
        y = (y | y <<  2) & 0x49249249;
        z = (z | z <<  2) & 0x49249249;

        return x | (y << 1) | (z << 2);
    }
};

template <>
struct Deposit<4> {
    static uint32_t deposit(const std::array<uint32_t, 4>& coords) {
#ifndef NDEBUG
        for(auto x : coords) assert(x <= 0xff);
#endif

        uint64_t x = coords[0] | (static_cast<uint64_t>(coords[1]) << 32);
        uint64_t y = coords[2] | (static_cast<uint64_t>(coords[3]) << 32);

        x = (x | (x << 12)) & 0x000F000F000F000F;
        y = (y | (y << 12)) & 0x000F000F000F000F;
        x = (x | (x <<  6)) & 0x0303030303030303;
        y = (y | (y <<  6)) & 0x0303030303030303;
        x = (x | (x <<  3)) & 0x1111111111111111;
        y = (y | (y <<  3)) & 0x1111111111111111;

        x |= x >> 31;
        y |= y >> 31;

        return static_cast<uint32_t>(x | (y << 2));
    }
};

template <unsigned D>
struct Implementation {
    static constexpr unsigned kDimensions = D;

    static uint32_t deposit(const std::array<uint32_t, D>& coords) {
        return Deposit<D>::deposit(coords);
    }

    static std::array<uint32_t, kDimensions> extract(uint32_t cell) {
        return Extract<D>::extract(cell);
    }

    static std::string name() {
        return "Generic";
    }
};

} // namespace Generic
} // namespace BitManipulationDetails
} // namespace girgs
