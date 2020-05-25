/*
 * StaticDegreeSequenceGenerator.cpp
 *
 *  Created on: 24.02.2014
 *      Author: Henning
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/generators/StaticDegreeSequenceGenerator.hpp>

namespace NetworKit {

StaticDegreeSequenceGenerator::StaticDegreeSequenceGenerator(const std::vector<count> &sequence):
        seq(sequence), realizable(UNKNOWN)
{

}

bool StaticDegreeSequenceGenerator::getRealizable() const {
    return realizable;
}


bool StaticDegreeSequenceGenerator::isRealizable() {
    DEBUG("check if sequence is realizable");
    count n = seq.size();

    /* First inequality. */
    count deg_sum = 0;
    for (count i = 0; i < n; ++i) {
        if (seq[i] >= n) {
            realizable = NO;
            DEBUG("not realizable: ", seq[i], ", n: ", n);
            return false;
        }
        deg_sum += seq[i];
    }

    if (deg_sum % 2 != 0) {
        DEBUG("not realizable, degree sum not even!");
        realizable = NO;
        return false;
    }

    /* Second inequality. */
    // this inequality needs a sorted sequence
    std::vector<count> sortedSeq = seq;
    Aux::Parallel::sort(sortedSeq.begin(), sortedSeq.end(), std::greater<count>());

    deg_sum = 0;
    for (count j = 0; j < n; ++j) {
        deg_sum += sortedSeq[j];

        /* sum of min(deg(i), j) for i from j + 1 to n - 1. */
        count min_deg_sum = 0;
        for (count i = j + 1; i < n; ++i) {
            min_deg_sum += std::min(sortedSeq[i], j + 1);
        }

        if (deg_sum > (j + 1) * j + min_deg_sum) {
            DEBUG("not realizable");
            realizable = NO;
            return false;
        }
    }

    realizable = YES;
    return true;
}


} /* namespace NetworKit */
