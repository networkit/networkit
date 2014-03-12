/*
 * HavelHakimiGenerator.h
 *
 *  Created on: Dec 10, 2013
 *      Author: Henning
 */

#ifndef HAVELHAKIMIGENERATOR_H_
#define HAVELHAKIMIGENERATOR_H_

#include <vector>
#include <stack>

#include "../graph/Graph.h"
#include "../auxiliary/ShellList.h"
#include "StaticDegreeSequenceGenerator.h"

namespace NetworKit {

/**
 * 	Havel-Hakimi algorithm for generating a graph according to a given degree sequence.
 * 	The sequence, if it is realizable, is reconstructed exactly.
 *  The resulting graph usually has a high clustering coefficient.
 *  Construction runs in linear time O(m). However, the test if a sequence is realizable
 *  is quadratic in the sequence length.
 */
class HavelHakimiGenerator: public NetworKit::StaticDegreeSequenceGenerator  {
protected:

public:
	/**
	 * @param[in] sequence Degree sequence to realize. Must be non-increasing.
	 * @param[in] skipTest If true, the test if the sequence is realizable is skipped.
	 *            Default value is false. Set ONLY to true if you are certain that the
	 *            sequence is realizable.
	 */
	HavelHakimiGenerator(const std::vector<unsigned int>& sequence, bool skipTest = false);
	virtual ~HavelHakimiGenerator() = default;

	/**
	 * Generates degree sequence seq (if it is realizable).
	 * @return Empty graph if graph is not realizable, otherwise graph with degree sequence seq.
	 */
	Graph generate() override;
};


} /* namespace NetworKit */
#endif /* HAVELHAKIMIGENERATOR_H_ */
