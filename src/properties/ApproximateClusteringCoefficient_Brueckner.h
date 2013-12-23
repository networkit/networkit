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

class ApproximateClusteringCoefficient_Brueckner {

public:

	ApproximateClusteringCoefficient_Brueckner();

	virtual ~ApproximateClusteringCoefficient_Brueckner();

	virtual double calculate(Graph& G, int k);
};

} /* namespace NetworKit */
#endif /* APPROXIMATECLUSTERINGCOEFFICIENT_H_ */
