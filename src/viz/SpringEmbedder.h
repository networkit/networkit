/*
 * SpringEmbedder.h
 *
 *  Created on: Apr 11, 2013
 *      Author: Henning
 */

#ifndef SPRINGEMBEDDER_H_
#define SPRINGEMBEDDER_H_

#include "Layouter.h"

namespace NetworKit {

const count MAX_ITER = 200; // TODO: maybe better as object variable

class SpringEmbedder: public NetworKit::Layouter {
protected:

public:

	SpringEmbedder(); // nullary constructor needed for Python shell

	SpringEmbedder(Point<float> bottomLeft, Point<float> topRight);

	virtual ~SpringEmbedder();

	virtual void draw(Graph& g) = 0;
};

} /* namespace NetworKit */
#endif /* SPRINGEMBEDDER_H_ */
