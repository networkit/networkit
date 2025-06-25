#ifndef NETWORKIT_AUXILIARY_ATOMIC_UTILS_HPP_
#define NETWORKIT_AUXILIARY_ATOMIC_UTILS_HPP_

#include <atomic>

namespace Aux {
inline uint32_t unsetBitAtomically(std::atomic_uint32_t &bitmask, size_t position) {
    uint32_t current_bitmask = bitmask.load();
    uint32_t desired_bitmask = current_bitmask & ~(1u << position);

    while (!bitmask.compare_exchange_weak(current_bitmask, desired_bitmask)) {
        current_bitmask = bitmask.load();
        desired_bitmask = current_bitmask & ~(1u << position);
    }

    return desired_bitmask;
}

inline uint32_t setBitAtomically(std::atomic_uint32_t &bitmask, size_t position) {
    uint32_t current_bitmask = bitmask.load();
    uint32_t desired_bitmask = current_bitmask | (1u << position);

    while (!bitmask.compare_exchange_weak(current_bitmask, desired_bitmask)) {
        current_bitmask = bitmask.load();
        desired_bitmask = current_bitmask | (1u << position);
    }

    return desired_bitmask;
}

template <typename T>
void swapAtomicsNonAtomically(std::atomic<T> &a, std::atomic<T> &b) noexcept {
    T const a_temp = a.load();
    a.store(b);
    b.store(a_temp);
}

template <typename T>
T incrementAtomically(std::atomic<T> &value, T increment = T(1)) {
    T current_value = value.load();
    T incremented_value = current_value + increment;

    while (!value.compare_exchange_weak(current_value, incremented_value)) {
        incremented_value = current_value + increment;
    }

    return incremented_value;
}
} // namespace Aux

#endif // NETWORKIT_AUXILIARY_ATOMIC_UTILS_HPP_
