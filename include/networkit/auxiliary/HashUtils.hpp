#ifndef NETWORKIT_AUXILIARY_HASH_UTILS_HPP_
#define NETWORKIT_AUXILIARY_HASH_UTILS_HPP_

#include <utility>

namespace Aux {

template <typename T>
void hashCombine(std::size_t &seed, T const &v) {
    // The Boost Software License, Version 1.0 applies to this function.
    // See https://www.boost.org/LICENSE_1_0.txt
    // https://www.boost.org/doc/libs/1_75_0/doc/html/hash/reference.html#boost.hash_combine
    seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

// capacity: must be power of 2
// T shall be unsigned.
template <std::unsigned_integral T>
T fastMod(T const idx, size_t const capacity) {
    return idx & (capacity - 1);
}

// Thanks to David Stafford for his further research on Austin Appleby's
// MurmurHash3 specifically for input values with low entropy such as it is the
// case for our dynamic hashtable that shall store vertex ids - and are more
// close to counting numbers than random ones.
//
// The one we are using here is Mix13.
//
// Please go to
// http://zimbry.blogspot.com/2011/09/better-bit-mixing-improving-on.html for
// further reading.
//
// Also thanks to Thomas Mueller for his answer on good integer hashing:
// https://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
//
// The final function definition comes from Sebastiano Vigna:
// https://xorshift.di.unimi.it/splitmix64.c
inline uint64_t hash64(uint64_t x) {
    constexpr uint64_t op_a = 0xbf58476d1ce4e5b9;
    constexpr uint64_t op_b = 0x94d049bb133111eb;

    x = (x ^ (x >> 30)) * op_a;
    x = (x ^ (x >> 27)) * op_b;
    x = x ^ (x >> 31);
    return x;
}

struct PairHash {
    template <typename A, typename B>
    std::size_t operator()(const std::pair<A, B> &pair) const {
        std::size_t seed = 0;
        hashCombine(seed, pair.first);
        hashCombine(seed, pair.second);
        return seed;
    }
};

} // namespace Aux

#endif // NETWORKIT_AUXILIARY_HASH_UTILS_HPP_
