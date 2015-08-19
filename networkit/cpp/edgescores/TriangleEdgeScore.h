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
 * An implementation of the triangle counting algorithm by Ortmann et al.
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
