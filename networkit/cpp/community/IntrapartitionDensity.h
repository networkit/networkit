#ifndef INTRAPARTITIONDENSITY_H
#define INTRAPARTITIONDENSITY_H

#include "LocalPartitionEvaluation.h"

namespace NetworKit {

/**
 * The intra-cluster density of a partition is defined as the number of existing edges divided by the number of possible edges.
 * The global value is the sum of all existing intra-cluster edges divided by the sum of all possible intra-cluster edges.
 */
class IntrapartitionDensity : public LocalPartitionEvaluation {
public:
	using LocalPartitionEvaluation::LocalPartitionEvaluation;

	/**
	 * Execute the algorithm.
	 */
	virtual void run() override;

	/**
	 * Get the global intra-cluster density.
	 *
	 * @return The global intra-cluster density.
	 */
	double getGlobal() const { assureFinished(); return globalValue; };

	/**
	 * @return false - this algorithm is not parallel.
	 */
	bool isParallel() const override { return false; }

	/**
	 * This value should be high in a good clustering.
	 * @return false - high values are better than small values.
	 */
	bool isSmallBetter() const override { return false; }
protected:
	double globalValue;
};

}

#endif // INTRAPARTITIONDENSITY_H
