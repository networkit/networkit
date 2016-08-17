/*
 * TimeVaryingSemiring.cpp
*/

#include "TimeVaryingSemiring.h"
#include <algorithm>

std::vector<std::vector<double>> TimeVaryingSemiring::zero = std::vector<std::vector<double>>();
//std::vector<std::vector<double>> TimeVaryingSemiring::zero(1, std::vector<double>(4));
std::vector<std::vector<double>> TimeVaryingSemiring::one {{ std::numeric_limits<double>::max(),
                                                            0,
                                                            std::numeric_limits<double>::max(),
                                                            0 }};

TimeVaryingSemiring::TimeVaryingSemiring() {
    zero.clear();
    this->value = zero;
}

TimeVaryingSemiring::TimeVaryingSemiring(std::vector<std::vector<double>> v) {
    for (u_int i = 0; i < v.size(); i++) {
        if (v[i].size() != 4) {
            this->value = std::vector<std::vector<double>>();
            return;
        }
        for (u_int j = 0; j < v[i].size(); j++) {
            if (v[i][j] < 0) {
                this->value = std::vector<std::vector<double>>();
                return;
            }
        }

    }
    this->value = v;
}

TimeVaryingSemiring& TimeVaryingSemiring::operator=(TimeVaryingSemiring x) {
    if (this == &x) {
        return *this;
    }
    for (u_int i = 0; i < x.value.size(); i++) {
        if (x.value[i].size() != 4) {
            this->value = std::vector<std::vector<double>>();
            return *this;
        }
        for (u_int j = 0; j < x.value[i].size(); j++) {
            if (x.value[i][j] < 0) {
                this->value = std::vector<std::vector<double>>();
                return *this;
            }
        }

    }
    this->value = x.value;
    return *this;
}

TimeVaryingSemiring TimeVaryingSemiring::operator+(TimeVaryingSemiring x) {
    std::vector<std::vector<double>> v = std::vector<std::vector<double>>();
    for (u_int i = 0; i < this->value.size(); i++) {
        v.push_back(this->value[i]);
    }
    for (u_int i = 0; i < x.value.size(); i++) {
        v.push_back(x.value[i]);
    }
    return v;
}

TimeVaryingSemiring TimeVaryingSemiring::operator+(TimeVaryingSemiring& x) {
    return (x + this->value);
}

TimeVaryingSemiring TimeVaryingSemiring::operator*(TimeVaryingSemiring& x) {

    if (this->equals(zero) || x.equals(zero)) {
        return zero;
    }

    if (this->equals(one)) {
        return x.value;
    }
    if (x.equals(one)) {
        return this->value;
    }
    std::vector<std::vector<double>> v = std::vector<std::vector<double>>();
    for (u_int i = 0; i < this->value.size(); i++) {
        for (u_int j = 0; j < x.value.size(); j++) {
            if (!(this->value[i][1] > this->value[i][2] ||
                    x.value[j][1] > x.value[j][2] ||
                    x.value[j][3] < this->value[i][0])) {
                v.push_back({std::max(this->value[i][0], x.value[j][0]),
                        this->value[i][1],
                        x.value[j][2],
                        std::min(this->value[i][3], x.value[j][3])});
            }
        }
    }
    return v;
}

bool TimeVaryingSemiring::equals(const TimeVaryingSemiring& x) const {
    return this->value == x.value;
}

std::ostream& operator<<(std::ostream& os, const TimeVaryingSemiring& x) {
    os << "{";
    for (u_int i = 0; i < x.value.size(); i++) {
        if (i == 0) {
            os << "[ ";
        } else {
            os << ", [ ";
        }
        for (u_int j = 0; j < x.value[i].size(); j++) {
            if (j == 0) {
                os << x.value[i][j];
            } else {
                os << ", " << x.value[i][j];
            }
        }
        os << " ]";
    }
    os << "}";
    return os;
}

