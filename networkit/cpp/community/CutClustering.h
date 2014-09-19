#ifndef CUTCLUSTERING_H
#define CUTCLUSTERING_H

#include "CommunityDetectionAlgorithm.h"

namespace NetworKit {

/**
 * Cut clustering algorithm as defined in
 * Flake, Gary William; Tarjan, Robert E.; Tsioutsiouliklis, Kostas. Graph Clustering and Minimum Cut Trees.
 * Internet Mathematics 1 (2003), no. 4, 385--408.
 */
class CutClustering : public CommunityDetectionAlgorithm {
public:
	/**
	 * Initialize cut clustering algorithm with parameter alpha.
	 *
	 * A value of 0 gives a clustering with one cluster with all nodes,
	 * A value that equals to the largest edge weight gives singleton clusters.
	 *
	 * @param alpha The parameter for the cut clustering
	 */
	CutClustering(edgeweight alpha);

	/**
	 * Apply algorithm to graph
	 * @return partition of the node set
	 */
	virtual Partition run(const Graph& G) override;

	/**
	 * @return string representation of algorithm and parameters.
	 */
	virtual std::string toString() const override;
private:
	edgeweight alpha;
};

} // namespace NetworKit

#endif // CUTCLUSTERING_H
