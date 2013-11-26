/*
 * ApproximateClusteringCoefficient_Ritter.h
 *
 */

#ifndef APPROXIMATE_CLUSTERING_COEFFICIENT_RITTER_H_
#define APPROXIMATE_CLUSTERING_COEFFICIENT_RITTER_H_

#include "../graph/Graph.h"

namespace NetworKit {

class ApproximateClusteringCoefficient_Ritter {
public:

	ApproximateClusteringCoefficient_Ritter();

	virtual ~ApproximateClusteringCoefficient_Ritter();

	virtual double calculate(const Graph& G, count k);
	
private:
	index indexOfBiggestValueLoEq(std::vector<count> values, count key);
};

} /* namespace NetworKit */
#endif /* APPROXIMATE_CLUSTERING_COEFFICIENT_RITTER_H_ */
