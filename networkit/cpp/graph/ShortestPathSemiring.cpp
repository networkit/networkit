/*
 * ShortestPathSemiring.cpp
 */

#include "ShortestPathSemiring.h"

ShortestPathSemiring::ShortestPathSemiring() {
   this->value = zero;
}

ShortestPathSemiring::ShortestPathSemiring(double v) {
    this->value = v;
}

ShortestPathSemiring& ShortestPathSemiring::operator=(ShortestPathSemiring x) {
    if (this == &x) {
        return *this;
    }
    value = x.value;
    return *this;
}

ShortestPathSemiring ShortestPathSemiring::operator+(ShortestPathSemiring& x) {
    return std::min(this->value, x.value);
}

ShortestPathSemiring ShortestPathSemiring::operator+(ShortestPathSemiring x) {
    return std::min(this->value, x.value);
}

ShortestPathSemiring ShortestPathSemiring::operator*(ShortestPathSemiring& x) {
    return this->value + x.value;
}

bool ShortestPathSemiring::lessThan(const ShortestPathSemiring& x) const {
    return this->value > x.value;
}

bool ShortestPathSemiring::equals(const ShortestPathSemiring& x) const {
    return this->value == x.value;
}

std::ostream& operator<<(std::ostream& os, const ShortestPathSemiring& s) {
    os << s.value;
    return os;
}

