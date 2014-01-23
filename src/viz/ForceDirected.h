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
#include <vector>
#include <cmath>

namespace NetworKit {

/**
 * TODO: class documentation
 */
class FruchtermanReingold: public NetworKit::Layouter {
private:
	const count MAX_ITER = 200;

public:

	/**
	 * Constructor
	 */
	FruchtermanReingold(); // nullary constructor needed for Python shell

	/**
	 * @param bottomLeft Coordinate of point in bottom/left corner
	 * @param topRight Coordinate of point in top/right corner
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
