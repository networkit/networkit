/*
 * ShortestPathSemiring.h
 */

#ifndef SHORTESTPATHSEMIRING_H_
#define SHORTESTPATHSEMIRING_H_

#include "SelectiveSemiring.h"

#include <iostream>
#include <algorithm>

class ShortestPathSemiring : public SelectiveSemiring<ShortestPathSemiring> {

private:

    double value;
    static constexpr double one = 0;
    static constexpr double zero = std::numeric_limits<double>::max();

public:

    ShortestPathSemiring();

    ShortestPathSemiring(double v);

    ShortestPathSemiring& operator=(ShortestPathSemiring x);

    ShortestPathSemiring operator+(ShortestPathSemiring& x);

    ShortestPathSemiring operator+(ShortestPathSemiring x);

    ShortestPathSemiring operator*(ShortestPathSemiring& x);

    bool lessThan(const ShortestPathSemiring& x) const;

    bool equals(const ShortestPathSemiring& x) const;

    friend std::ostream& operator<<(std::ostream& os, const ShortestPathSemiring& s);

    static ShortestPathSemiring getZero() {
        return zero;
    }

    static ShortestPathSemiring getOne() {
        return one;
    }
};

#endif /* SHORTESTPATHSEMIRING_H_ */
