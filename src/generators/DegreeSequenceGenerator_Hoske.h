/*
 * DegreeSequenceGenerator_Hoske.h
 *
 *  Created on: 05-12-2013
 *      Author: dhoske
 */

#ifndef DEGREESEQUENCEGENERATOR_HOSKE_H_
#define DEGREESEQUENCEGENERATOR_HOSKE_H_

#include "StaticGraphGenerator.h"
#include <cmath>

namespace NetworKit {

/**
 * Generates a graph with a given degree sequence.
 */
class DegreeSequenceGenerator_Hoske: public NetworKit::StaticGraphGenerator {
protected:
    std::vector<count> deg_seq; /* Desired degree sequence. */
    count n;                    /* Number of nodes. */

public:
    /** nullary constructor needed for Python Shell - do not use this to construct instance */
    DegreeSequenceGenerator_Hoske() {
    };

    /** Constructor that takes a degree sequence. */
    DegreeSequenceGenerator_Hoske(const std::vector<count>& degree_sequence)
        : deg_seq(degree_sequence), n(count(degree_sequence.size())) {
    };

    /** Default destrutor. */
    virtual ~DegreeSequenceGenerator_Hoske() override = default;

    /** Generates a graph with the given degree sequence. */
    virtual Graph generate() override;
};

} /* namespace NetworKit */
#endif /* DEGREESEQUENCEGENERATOR_HOSKE_H_ */
