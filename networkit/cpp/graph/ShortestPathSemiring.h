/*
 * ShortestPathSemiring.h
 */

#ifndef SHORTESTPATHSEMIRING_H_
#define SHORTESTPATHSEMIRING_H_

#include <iostream>
#include <algorithm>

class ShortestPathSemiring {

private:

    double value;
    static constexpr double zero = 0;
    static constexpr double infinity = std::numeric_limits<double>::max();

public:

    ShortestPathSemiring();

    ShortestPathSemiring(double v);

    ShortestPathSemiring& operator=(ShortestPathSemiring x);

    ShortestPathSemiring operator+(ShortestPathSemiring& x);

    ShortestPathSemiring operator+(ShortestPathSemiring x);

    ShortestPathSemiring operator*(ShortestPathSemiring& x);

    friend inline bool operator<(const ShortestPathSemiring& x, const ShortestPathSemiring& y) { return x.value < y.value; }
    friend inline bool operator> (const ShortestPathSemiring& x, const ShortestPathSemiring& y) { return y < x; }
    friend inline bool operator<= (const ShortestPathSemiring& x, const ShortestPathSemiring& y) { return !(x < y); }
    friend inline bool operator>= (const ShortestPathSemiring& x, const ShortestPathSemiring& y) { return !(x > y); }
    friend inline bool operator== (const ShortestPathSemiring& x, const ShortestPathSemiring& y) { return x.value == y.value; }
    friend inline bool operator!= (const ShortestPathSemiring& x, const ShortestPathSemiring& y) { return !(x == y); }

    friend std::ostream& operator<<(std::ostream& os, const ShortestPathSemiring& s);

    static ShortestPathSemiring getZero();

    static ShortestPathSemiring getInfinity();
};

#endif /* SHORTESTPATHSEMIRING_H_ */
