/*
 * ChibaNishizekiTriangleEdgeScore.h
 *
 *  Created on: 22.05.2014
 *      Author: Gerd Lindner
 */

#ifndef CHIBANISHIZEKITRIANGLEEDGESCORE_H_
#define CHIBANISHIZEKITRIANGLEEDGESCORE_H_

#include "../graph/Graph.h"
#include "EdgeScore.h"

namespace NetworKit {

/**
 * An implementation of the triangle counting algorithm by Chiba/Nishizeki.
 *
 * @deprecated Use TriangleEdgeScore instead which is parallelized and has a similar performance even in the sequential case.
 */
class ChibaNishizekiTriangleEdgeScore : public EdgeScore<count> {

public:

	ChibaNishizekiTriangleEdgeScore(const Graph& G);
	virtual count score(edgeid eid) override;
	virtual count score(node u, node v) override;
	virtual void run() override;
};

} /* namespace NetworKit */

#endif /* CHIBANISHIZEKITRIANGLEEDGESCORE_H_ */
