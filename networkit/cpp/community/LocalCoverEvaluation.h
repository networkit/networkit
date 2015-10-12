#ifndef LOCALCOVEREVALUATION_H
#define LOCALCOVEREVALUATION_H

#include "LocalCommunityEvaluation.h"
#include "../structures/Cover.h"
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * Virtual base class of all evaluation methods for a single Cover which is based on the evaluation of single clusters.
 * This is the base class for Covers.
 */
class LocalCoverEvaluation : public LocalCommunityEvaluation {
public:
	/**
	 * Initialize the cover evaluation method.
	 *
	 * @param G The graph on which the evaluation shall be performed
	 * @param C The cover that shall be evaluated.
	 */
	LocalCoverEvaluation(const Graph &G, const Cover &C);
protected:
	const Graph &G;
	const Cover &C;
};

}

#endif // LOCALCOVEREVALUATION_H
