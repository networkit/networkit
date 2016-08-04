/*
 * SafestPathSemiring.h
 */

#ifndef SAFESTPATHSEMIRING_H_
#define SAFESTPATHSEMIRING_H_

#include <iostream>

class SafestPathSemiring {

private:

    double value;
    static constexpr double zero = 1;
    static constexpr double infinity = 0;

public:

    SafestPathSemiring();

    SafestPathSemiring(double v);

    SafestPathSemiring& operator=(SafestPathSemiring x);

    SafestPathSemiring operator+(SafestPathSemiring& x);

    SafestPathSemiring operator+(SafestPathSemiring x);

    SafestPathSemiring operator*(SafestPathSemiring& x);

    friend inline bool operator<(const SafestPathSemiring& x, const SafestPathSemiring& y) { return x.value < y.value; }
    friend inline bool operator> (const SafestPathSemiring& x, const SafestPathSemiring& y) { return y < x; }
    friend inline bool operator<= (const SafestPathSemiring& x, const SafestPathSemiring& y) { return !(x < y); }
    friend inline bool operator>= (const SafestPathSemiring& x, const SafestPathSemiring& y) { return !(x > y); }
    friend inline bool operator== (const SafestPathSemiring& x, const SafestPathSemiring& y) { return x.value == y.value; }
    friend inline bool operator!= (const SafestPathSemiring& x, const SafestPathSemiring& y) { return !(x == y); }

    friend std::ostream& operator<<(std::ostream& os, const SafestPathSemiring& s);

    static SafestPathSemiring getZero();

    static SafestPathSemiring getInfinity();
};

#endif /* SAFESTPATHSEMIRING_H_ */
