/*
 * ApproximateClusteringCoefficient.h
 *
 *  Created on: 16.11.2013
 *      Author: gbrueckner
 */

#ifndef APPROXIMATECLUSTERINGCOEFFICIENT_H_
#define APPROXIMATECLUSTERINGCOEFFICIENT_H_

#include "../graph/Graph.h"

namespace NetworKit {

// TODO: is class necessary?
class ApproximateClusteringCoefficient {

public:

	ApproximateClusteringCoefficient();

	virtual ~ApproximateClusteringCoefficient();

	virtual double calculate(Graph& G, int k);
};

} /* namespace NetworKit */
#endif /* APPROXIMATECLUSTERINGCOEFFICIENT_H_ */
