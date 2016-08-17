/*
 * Semiring.h
 */

#ifndef SEMIRING_H_
#define SEMIRING_H_

#include <iostream>

template <class T>
class Semiring {

public:

    Semiring() {
    }

    Semiring(T v) {
    }

    virtual T& operator=(T x);

    virtual T operator+(T& x);

    virtual T operator+(T x);

    virtual T operator*(T& x);

    virtual bool equals(const T& x) const;

    friend inline bool operator== (const T& x, const T& y) { return x.equals(y); }

    friend inline bool operator!= (const T& x, const T& y) { return !(x == y); }

    friend std::ostream& operator<<(std::ostream& os, const T& x);

    static T getZero();

    static T getOne();

    //virtual T getValue();


};

#endif /* SEMIRNG_H_ */
