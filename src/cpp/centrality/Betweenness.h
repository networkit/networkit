/*
 * Betweenness.h
 *
 *  Created on: 19.02.2014
 *      Author: hm
 */

#ifndef BETWEENNESS_H_
#define BETWEENNESS_H_

#include "Centrality.h"

namespace NetworKit {

class Betweenness: public NetworKit::Centrality {
public:
	Betweenness(const Graph& G, bool normalized=false);

	void run(bool runUnweightedInParallel);

	void run() override { run(false); }

};

} /* namespace NetworKit */

#endif /* BETWEENNESS_H_ */
