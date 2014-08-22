/*
 * TNeighborhoodDistance.h
 *
 *  Created on: 24.06.2013
 *      Author: cls
 */

#ifndef TNEIGHBORHOODDISTANCE_H_
#define TNEIGHBORHOODDISTANCE_H_

#include "TNodeDistance.h"

namespace NetworKit {

/**
 * @ingroup distmeasures
 */
class TNeighborhoodDistance: public NetworKit::TNodeDistance {
public:

	TNeighborhoodDistance(const Graph& G);

	/** Default destructor */
	virtual ~TNeighborhoodDistance();

	void initialize(const Parameters& param);

	double distance(node u, node v);
};

} /* namespace NetworKit */
#endif /* TNEIGHBORHOODDISTANCE_H_ */
