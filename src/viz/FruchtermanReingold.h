/*
 * ForceDirected.h
 *
 *  Created on: Apr 11, 2013
 *      Author: Henning
 */

#ifndef FORCEDIRECTED_H_
#define FORCEDIRECTED_H_

#include "Layouter.h"
#include "Point.h"
#include "PostscriptWriter.h"

#include <vector>
#include <cmath>


namespace NetworKit {

/**
 * Fruchterman-Reingold graph drawing algorithm. We mostly follow
 * the description in Stephen G. Kobourov: Spring Embedders and Force
 * Directed Graph Drawing Algorithms.
 */
class FruchtermanReingold: public NetworKit::Layouter {
private:
	const count MAX_ITER = 250;

public:

	/**
	 * Constructor. DO NOT use for creating objects, only needed for
	 * Python shell.
	 */
	FruchtermanReingold() {}

	/**
	 * Constructor.
	 * @param[in] bottomLeft Coordinate of point in bottom/left corner
	 * @param[in] topRight Coordinate of point in top/right corner
	 */
	FruchtermanReingold(Point<float> bottomLeft, Point<float> topRight);

	virtual ~FruchtermanReingold();

	/**
	 * Assigns coordinates to vertices in graph @a g
	 */
	virtual void draw(Graph& g);
};

} /* namespace NetworKit */
#endif /* FORCEDIRECTED_H_ */
