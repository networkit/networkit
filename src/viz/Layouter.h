/*
 * Layouter.h
 *
 *  Created on: Apr 11, 2013
 *      Author: Henning
 */

#ifndef LAYOUTER_H_
#define LAYOUTER_H_

#include "../graph/Graph.h"
#include "../auxiliary/RandomProbability.h"

namespace NetworKit {

class Layouter {
protected:
	Point<float> bottomLeft;
	Point<float> topRight;
	std::vector<Point<float> > layout;

public:
	Layouter(); // nullary constructor needed for Python shell
	Layouter(Point<float> bottomLeft, Point<float> topRight);
	virtual ~Layouter();

	virtual void draw(Graph& g) = 0;

	virtual void randomInitCoordinates(Graph& g);
};

} /* namespace NetworKit */
#endif /* LAYOUTER_H_ */
