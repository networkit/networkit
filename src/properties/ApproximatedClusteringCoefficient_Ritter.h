/*
 * ApproximatedClusteringCoefficient_Ritter.h
 *
 */

#ifndef APPROXIMATED_CLUSTERING_COEFFICIENT_RITTER_H_
#define APPROXIMATED_CLUSTERING_COEFFICIENT_RITTER_H_

#include "../graph/Graph.h"

namespace NetworKit {

class ApproximatedClusteringCoefficient_Ritter {
public:

	ApproximatedClusteringCoefficient_Ritter();

	virtual ~ApproximatedClusteringCoefficient_Ritter();

	virtual double calculate(Graph& G, count k);
	
private:
	index indexOfBiggestValueLoEq(count* array, count size, count key);
	count randomNode(count* cumulated_propabilities, count n, count propability_sum);
};

} /* namespace NetworKit */
#endif /* APPROXIMATED_CLUSTERING_COEFFICIENT_RITTER_H_ */
