/*
 * HavelHakimiGenerator.cpp
 *
 *  Created on: Dec 10, 2013
 *      Author: Henning
 */

#include "HavelHakimiGenerator.h"

namespace NetworKit {

HavelHakimiGenerator::HavelHakimiGenerator(const std::vector<count>& sequence): seq(sequence),
		realizable(false) {

}

HavelHakimiGenerator::~HavelHakimiGenerator() {

}

// TODO: Students, please rename class and implement this method
bool HavelHakimiGenerator::isRealizable() {
	realizable = false;
	count n = seq.size();


	// check first inequation


	// check second inequation


	realizable = true;
	return realizable;
}


// TODO: Students, please rename class and implement this method
Graph HavelHakimiGenerator::generate() {
	Graph G;

	if (! realizable) {
		WARN("Degree sequence not realizable or not checked for realizability yet! Will return empty graph!")
	}
	else {
		DEBUG("Degree sequence is realizable, continue with generation algorithm.")


		// TODO
	}

	return G;
}

} /* namespace NetworKit */
