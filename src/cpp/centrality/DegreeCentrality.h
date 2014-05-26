/*
 * DegreeCentrality.h
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#ifndef DEGREECENTRALITY_H_
#define DEGREECENTRALITY_H_

#include "Centrality.h"

namespace NetworKit {

/** 
 * Node centrality index which ranks nodes by their degree.
 * Optional normalization by maximum degree.
 */
class DegreeCentrality: public NetworKit::Centrality {
public:
	DegreeCentrality(const Graph& G, bool normalized=false);

	void run() override;
};

} /* namespace NetworKit */

#endif /* DEGREECENTRALITY_H_ */
