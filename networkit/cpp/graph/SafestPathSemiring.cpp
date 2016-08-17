/*
 * SafestPathSemiring.cpp
 */

#include "SafestPathSemiring.h"

#include <algorithm>

SafestPathSemiring::SafestPathSemiring() {
    this->value = zero;
}

SafestPathSemiring::SafestPathSemiring(double v) {
    if (v >= 0 && v <= 1) {
        this->value = v;
    } else {
        this->value = 0;
    }
}

SafestPathSemiring& SafestPathSemiring::operator=(SafestPathSemiring x) {
    if (this == &x) {
        return *this;
    }
    if (x.value >= 0 && x.value <= 1) {
        value = x.value;
    }
    return *this;
}

SafestPathSemiring SafestPathSemiring::operator+(SafestPathSemiring& x) {
    return std::max(this->value, x.value);
}

SafestPathSemiring SafestPathSemiring::operator+(SafestPathSemiring x) {
    return std::max(this->value, x.value);
}

SafestPathSemiring SafestPathSemiring::operator*(SafestPathSemiring& x) {
    return this->value * x.value;
}

bool SafestPathSemiring::lessThan(const SafestPathSemiring& x) const {
    return this->value < x.value;
}

bool SafestPathSemiring::equals(const SafestPathSemiring& x) const {
    return this->value == x.value;
}

std::ostream& operator<<(std::ostream& os, const SafestPathSemiring& s) {
    os << s.value;
    return os;
}

