/*
 * TimeVaryingSemiring.h
 */

#ifndef TIMEVARYINGSEMIRING_H_
#define TIMEVARYINGSEMIRING_H_

#include "Semiring.h"

#include <iostream>
#include <vector>

class TimeVaryingSemiring : public Semiring<TimeVaryingSemiring> {

protected:

    std::vector<std::vector<double>> interval;
    std::vector<std::vector<double>> value;
    static std::vector<std::vector<double>> zero;
    static std::vector<std::vector<double>> one;

public:

    TimeVaryingSemiring();

    TimeVaryingSemiring(std::vector<std::vector<double>> i);

    TimeVaryingSemiring& operator=(TimeVaryingSemiring x);
    
    TimeVaryingSemiring operator+(TimeVaryingSemiring& x);

    TimeVaryingSemiring operator+(TimeVaryingSemiring x);

    TimeVaryingSemiring operator*(TimeVaryingSemiring& x);

    bool equals(const TimeVaryingSemiring& x) const;

    friend std::ostream& operator<<(std::ostream& os, const TimeVaryingSemiring& x);

    static TimeVaryingSemiring getZero() {
        return zero;
    }

    static TimeVaryingSemiring getOne() {
        return one;
    }
};

#endif /* TIMEVARYINGSEMIRING_H_ */
