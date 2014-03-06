/*
 * Betweenness.h
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#ifndef BETWEENNESS_H_
#define BETWEENNESS_H_

#include "Centrality.h"

namespace NetworKit {

class Betweenness: public NetworKit::Centrality {
public:
	Betweenness(const Graph& G, bool normalized=false);

	void run() override;

};

} /* namespace NetworKit */

#endif /* BETWEENNESS_H_ */
