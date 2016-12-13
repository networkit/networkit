/*
 * ChibaNishizekiQuadrangleEdgeScore.h
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann, Gerd Lindner
 */

#ifndef CHIBANISHIZEKIQUADRANGLEEDGESCORE_H
#define CHIBANISHIZEKIQUADRANGLEEDGESCORE_H

#include "EdgeScore.h"

namespace NetworKit {

class ChibaNishizekiQuadrangleEdgeScore : public EdgeScore<count> {

public:
	ChibaNishizekiQuadrangleEdgeScore(const Graph& G);
	virtual count score(edgeid eid) override;
	virtual count score(node u, node v) override;
	virtual void run() override;
};

} // namespace NetworKit

#endif // CHIBANISHIZEKIQUADRANGLEEDGESCORE_H
