/*
 * SelectiveSemiring.h
 */

#ifndef SELECTIVESEMIRING_H_
#define SELECTIVESEMIRING_H_

#include "Semiring.h"

#include <iostream>

template <class T>
class SelectiveSemiring : public Semiring<T> {

public:

    SelectiveSemiring() {
    }

    SelectiveSemiring(T v) {
    }

    virtual bool lessThan(const T& x) const;

    friend inline bool operator<(const T& x, const T& y) { return x.lessThan(y); }

    friend inline bool operator> (const T& x, const T& y) { return y < x; }

    friend inline bool operator<= (const T& x, const T& y) { return !(x < y); }

    friend inline bool operator>= (const T& x, const T& y) { return !(x > y); }

};

#endif /* SELECTIVESEMIRNG_H_ */
