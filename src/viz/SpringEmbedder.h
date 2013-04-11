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

class SpringEmbedder: public NetworKit::Layouter {
public:
	SpringEmbedder();
	virtual ~SpringEmbedder();

	virtual void draw(Graph& g) = 0;
};

} /* namespace NetworKit */
#endif /* SPRINGEMBEDDER_H_ */
