/*
 * Layouter.h
 *
 *  Created on: Apr 11, 2013
 *      Author: Henning
 */

#ifndef LAYOUTER_H_
#define LAYOUTER_H_

#include "../graph/Graph.h"

namespace NetworKit {

class Layouter {
public:
	Layouter();
	virtual ~Layouter();

	virtual void draw(Graph& g) = 0;
};

} /* namespace NetworKit */
#endif /* LAYOUTER_H_ */
