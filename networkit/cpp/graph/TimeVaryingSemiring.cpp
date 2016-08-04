/*
 * TimeVaryingSemiring.cpp
 */

#include "TimeVaryingSemiring.h"
#include <algorithm>

TimeVaryingSemiring::TimeVaryingSemiring() {
    this->interval = std::vector<std::vector<double>>();
}

TimeVaryingSemiring::TimeVaryingSemiring(std::vector<std::vector<double>> v) {
    for (u_int i = 0; i < v.size(); i++) {
        if (v[i].size() != 4) {
            this->interval = std::vector<std::vector<double>>();
            return;
        }
        for (u_int j = 0; j < v[i].size(); j++) {
            if (v[i][j] < 0) {
                this->interval = std::vector<std::vector<double>>();
                return;
            }
        }

    }
    this->interval = v;
}

// TimeVaryingSemiring::TimeVaryingSemiring(vector<double> v) {
    

TimeVaryingSemiring& TimeVaryingSemiring::operator=(TimeVaryingSemiring& x) {
    if (this == &x) {
        return *this;
    }
    for (u_int i = 0; i < x.interval.size(); i++) {
        if (x.interval[i].size() != 4) {
            this->interval = std::vector<std::vector<double>>();
            return *this;
        }
        for (u_int j = 0; j < x.interval[i].size(); j++) {
            if (x.interval[i][j] < 0) {
                this->interval = std::vector<std::vector<double>>();
                return *this;
            }
        }

    }
    this->interval = x.interval;
    return *this;
}

TimeVaryingSemiring& TimeVaryingSemiring::operator=(TimeVaryingSemiring x) {
    if (this == &x) {
        return *this;
    }
    for (u_int i = 0; i < x.interval.size(); i++) {
        if (x.interval[i].size() != 4) {
            this->interval = std::vector<std::vector<double>>();
            return *this;
        }
        for (u_int j = 0; j < x.interval[i].size(); j++) {
            if (x.interval[i][j] < 0) {
                this->interval = std::vector<std::vector<double>>();
                return *this;
            }
        }

    }
    this->interval = x.interval;
    return *this;
}

TimeVaryingSemiring TimeVaryingSemiring::operator+(TimeVaryingSemiring& x) {
    std::vector<std::vector<double>> v = std::vector<std::vector<double>>();
    for (u_int i = 0; i < this->interval.size(); i++) {
        v.push_back(this->interval[i]);
    }
    for (u_int i = 0; i < x.interval.size(); i++) {
        v.push_back(x.interval[i]);
    }
    return v;
}

TimeVaryingSemiring TimeVaryingSemiring::operator*(TimeVaryingSemiring& x) {
    if (this->interval[0][0] == std::numeric_limits<double>::max() &&
            this->interval[0][1] == 0.0 &&
            this->interval[0][2] == std::numeric_limits<double>::max() &&
            this->interval[0][3] == 0.0) {
        return x.interval;
    }
    if (x.interval[0][0] == std::numeric_limits<double>::max() &&
            x.interval[0][1] == 0.0 &&
            x.interval[0][2] == std::numeric_limits<double>::max() &&
            x.interval[0][3] == 0.0) {
        return this->interval;
    }
    std::vector<std::vector<double>> v = std::vector<std::vector<double>>();
    for (u_int i = 0; i < this->interval.size(); i++) {
        for (u_int j = 0; j < x.interval.size(); j++) {
            if (!(this->interval[i][1] > this->interval[i][2] ||
                    x.interval[j][1] > x.interval[j][2] ||
                    x.interval[j][3] < this->interval[i][0])) {
                v.push_back({std::max(this->interval[i][0], x.interval[j][0]),
                        this->interval[i][1],
                        x.interval[j][2],
                        std::min(this->interval[i][3], x.interval[j][3])});
            }
        }
    }
    return v;
}

std::ostream& operator<<(std::ostream& os, const TimeVaryingSemiring& s) {
    os << "{";
    for (u_int i = 0; i < s.interval.size(); i++) {
        if (i == 0) {
            os << "[ ";
        } else {
            os << ", [ ";
        }
        for (u_int j = 0; j < s.interval[i].size(); j++) {
            if (j == 0) {
                os << s.interval[i][j];
            } else {
                os << ", " << s.interval[i][j];
            }
        }
        os << " ]";
    }
    os << "}";
    return os;
}

