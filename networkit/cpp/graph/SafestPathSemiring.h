/*
 * SafestPathSemiring.h
 */

#ifndef SAFESTPATHSEMIRING_H_
#define SAFESTPATHSEMIRING_H_

#include "SelectiveSemiring.h"

#include <iostream>

class SafestPathSemiring : public SelectiveSemiring<SafestPathSemiring> {

private:

    double value;
    static constexpr double zero = 0;
    static constexpr double one = 1;

public:

    SafestPathSemiring();

    SafestPathSemiring(double v);

    SafestPathSemiring& operator=(SafestPathSemiring x);

    SafestPathSemiring operator+(SafestPathSemiring& x);

    SafestPathSemiring operator+(SafestPathSemiring x);

    SafestPathSemiring operator*(SafestPathSemiring& x);

    bool lessThan(const SafestPathSemiring& x) const;

    bool equals(const SafestPathSemiring& x) const;

    friend std::ostream& operator<<(std::ostream& os, const SafestPathSemiring& s);

    static SafestPathSemiring getZero() {
        return zero;
    }

    static SafestPathSemiring getOne() {
        return one;
    }
};

#endif /* SAFESTPATHSEMIRING_H_ */
