/*
 * MultilevelLayouter.h
 *
 *  Created on: 27.01.2014
 *      Author: Henning
 */

#ifndef MULTILEVELLAYOUTER_H_
#define MULTILEVELLAYOUTER_H_

#include "Layouter.h"
#include "FruchtermanReingold.h"
#include "../coarsening/Contracter.h"
#include "../coarsening/MatchingContracter.h"
#include "../matching/ParallelMatcher.h"
#include "../matching/Matching.h"

namespace NetworKit {


class MultilevelLayouter: public NetworKit::Layouter {
protected:
	static const count N_THRSH;

public:
	MultilevelLayouter(Point<float> bottomLeft, Point<float> topRight);
	virtual ~MultilevelLayouter();

	virtual void draw(Graph& G);
	virtual void drawWrapper(Graph& G, count level);

	virtual void prolongCoordinates(Graph& Gcon, Graph& G);
};

} /* namespace NetworKit */
#endif /* MULTILEVELLAYOUTER_H_ */
