/*
 * TNodeDistance.h
 *
 *  Created on: 24.06.2013
 *      Author: cls
 */

#ifndef TNODEDISTANCE_H_
#define TNODEDISTANCE_H_

#include "../base/Parameters.h"
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup distmeasures
 */
class TNodeDistance {

public:

	TNodeDistance(const Graph& G);

	/** Default destructor */
	virtual ~TNodeDistance();

	void initialize(const Parameters& param);

	double distance(node u, node v);

protected:

	const Graph& G;
};

} /* namespace NetworKit */
#endif /* TNODEDISTANCE_H_ */
