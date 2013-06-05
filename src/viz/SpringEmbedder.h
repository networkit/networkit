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

const count MAX_ITER = 1; // TODO: maybe better as object variable

class SpringEmbedder: public NetworKit::Layouter {
protected:
	float forceCorrection; //<! value for scaling the forces

public:
	SpringEmbedder(Point<float> bottomLeft, Point<float> topRight);
	virtual ~SpringEmbedder();

	virtual void draw(Graph& g) = 0;
};

} /* namespace NetworKit */
#endif /* SPRINGEMBEDDER_H_ */
