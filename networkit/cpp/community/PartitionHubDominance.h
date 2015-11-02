#ifndef PARTITIONHUBDOMINANCE_H
#define PARTITIONHUBDOMINANCE_H

#include "LocalPartitionEvaluation.h"

namespace NetworKit {

/**
 * A quality measure that measures the dominance of hubs in clusters. The hub dominance of a single
 * cluster is defined as the maximum cluster-internal degree of a node in that cluster divided by
 * the maximum cluster-internal degree, i.e. the number of nodes in the cluster minus one. The
 * value for all clusters is defined as the average of all clusters.
 * Strictly speaking this is not a quality measure as this is rather dependent on the type of the
 * considered graph, for more information see
 * Lancichinetti A, Kivelä M, Saramäki J, Fortunato S (2010)
 * Characterizing the Community Structure of Complex Networks
 * PLoS ONE 5(8): e11976. doi: 10.1371/journal.pone.0011976
 * http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0011976
 */
class PartitionHubDominance : public LocalPartitionEvaluation {
public:
	using LocalPartitionEvaluation::LocalPartitionEvaluation;

	/**
	 * Execute the algorithm.
	 */
	virtual void run() override;

	/**
	 * @return false - small values are not better, large values indicate better cluster cohesion.
	 */
	virtual bool isSmallBetter() const override { return false; }

	/**
	 * @return false - this implementation is not paralle.
	 */
	virtual bool isParallel() const override { return false; }
};

}

#endif // PARTITIONHUBDOMINANCE_H
