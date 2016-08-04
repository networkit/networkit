/*
 * TimeVaryingSemiring.h
 */

#ifndef TIMEVARYINGSEMIRING_H_
#define TIMEVARYINGSEMIRING_H_

#include <iostream>
#include <vector>

class TimeVaryingSemiring {

private:

    std::vector<std::vector<double>> interval;

public:
    TimeVaryingSemiring();

    TimeVaryingSemiring(std::vector<std::vector<double>> i);

    TimeVaryingSemiring& operator=(TimeVaryingSemiring& x);

    TimeVaryingSemiring& operator=(TimeVaryingSemiring x);
    
    TimeVaryingSemiring operator+(TimeVaryingSemiring& x);

    TimeVaryingSemiring operator*(TimeVaryingSemiring& x);

    friend std::ostream& operator<<(std::ostream& os, const TimeVaryingSemiring& s);
};

#endif /* TIMEVARYINGSEMIRING_H_ */
