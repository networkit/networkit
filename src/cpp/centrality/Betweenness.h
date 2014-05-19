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
	Betweenness(const IGraph& G, bool normalized=false);

	void run(bool runUnweightedInParallel);

	void run() override { run(false); }

};

} /* namespace NetworKit */

#endif /* BETWEENNESS_H_ */
