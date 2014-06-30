/*
 * StaticDegreeSequenceGenerator.cpp
 *
 *  Created on: 24.02.2014
 *      Author: Henning
 */

#include "StaticDegreeSequenceGenerator.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

StaticDegreeSequenceGenerator::StaticDegreeSequenceGenerator(const std::vector<unsigned int>& sequence):
		seq(sequence), realizable(UNKNOWN)
{

}

StaticDegreeSequenceGenerator::~StaticDegreeSequenceGenerator() {

}


bool StaticDegreeSequenceGenerator::getRealizable() const {
	return realizable;
}


bool StaticDegreeSequenceGenerator::isRealizable() {
	DEBUG("check if sequence is realizable");
	count n = seq.size();

	/* First inequality. */
	unsigned int deg_sum = 0;
	for (count i = 0; i < n; ++i) {
		if (seq[i] < 0 || seq[i] >= n) {
			realizable = NO;
			DEBUG("not realizable: ", seq[i], ", n: ", n);
			return false;
		}
		deg_sum += seq[i];
	}

	if (deg_sum % 2 != 0) {
		DEBUG("not realizable");
		realizable = NO;
		return false;
	}

	/* Second inequality. */
	deg_sum = 0;
	for (unsigned int j = 0; j < n; ++j) {
		deg_sum += seq[j];

		/* sum of min(deg(i), j) for i from j + 1 to n - 1. */
		unsigned int min_deg_sum = 0;
		for (unsigned int i = j + 1; i < n; ++i) {
			min_deg_sum += std::min(seq[i], j + 1);
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
