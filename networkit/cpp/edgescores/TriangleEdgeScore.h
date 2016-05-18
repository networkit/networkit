/*
 * TriangleEdgeScore.h
 *
 *  Created on: 22.05.2014
 *      Author: Gerd Lindner, Michael Hamann
 */

#ifndef TRIANGLE_COUNTER_H_
#define TRIANGLE_COUNTER_H_

#include "EdgeScore.h"

namespace NetworKit {

/**
 * A parallel triangle counting implementation based on ideas in [0].
 *
 * This implementation avoids less work but needs therefore also less checks and
 * is (apart from a fast initialization) parallel without any locks. In experiments
 * this implementation seems to be fast both in non-parallel as well as in parallel
 * settings. With only one thread its performance is similar to the sequential
 * ChibaNishizekiTriangleEdgeScore in NetworKit.
 *
 * [0] Triangle Listing Algorithms: Back from the Diversion
 * Mark Ortmann and Ulrik Brandes                                                                          *
 * 2014 Proceedings of the Sixteenth Workshop on Algorithm Engineering and Experiments (ALENEX). 2014, 1-8
 */
class TriangleEdgeScore : public EdgeScore<count> {

public:

	TriangleEdgeScore(const Graph& G);
	virtual count score(edgeid eid) override;
	virtual count score(node u, node v) override;
	virtual void run() override;

};

} /* namespace NetworKit */

#endif /* TRIANGLE_COUNTER_H_ */
