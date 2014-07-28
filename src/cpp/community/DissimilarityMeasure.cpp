/*
 * DissimilarityMeasure.cpp
 *
 *  Created on: 19.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "DissimilarityMeasure.h"

namespace NetworKit {

DissimilarityMeasure::~DissimilarityMeasure() {
}

double DissimilarityMeasure::getDissimilarity(Graph &G, Cover &first, Cover &second) {
	throw std::runtime_error("Not implemented");
}

} /* namespace NetworKit */
