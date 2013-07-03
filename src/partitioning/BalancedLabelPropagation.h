/*
 * BalancedLabelPropagation.h
 *
 *  Created on: Jun 19, 2013
 *      Author: Henning
 */

#ifndef BALANCEDLABELPROPAGATION_H_
#define BALANCEDLABELPROPAGATION_H_

#include "../clustering/ClusteringGenerator.h"
#include "../clustering/EdgeCut.h"
#include "../io/ClusteringWriter.h"
#include "../clustering/Clustering.h"
#include "../auxiliary/RandomProbability.h"
#include "../coarsening/MatchingContracter.h"
#include "../matching/ParallelMatcher.h"
#include "../coarsening/ClusteringProjector.h"
#include "BalancedPartitioner.h"
#include <vector>
#include <algorithm>

namespace NetworKit {

class BalancedLabelPropagation: public BalancedPartitioner {
protected:
	double exponent;

public:
	BalancedLabelPropagation(double exponent);
	virtual ~BalancedLabelPropagation();

	virtual Clustering run(Graph& graph, count numBlocks);
	virtual Clustering& rerun(Graph& graph, count numBlocks, Clustering& partition);
//	virtual Clustering multilevelRun(Graph& graph, count numParts);
//	virtual Clustering& multilevelRerun(Graph& graph, count numParts, Clustering& partition);

	virtual Clustering& postsmooth(Graph& graph, count numBlocks, Clustering& partition);

	void setExponent(double exponent) {
		this->exponent = exponent;
	}
};

} /* namespace NetworKit */
#endif /* BALANCEDLABELPROPAGATION_H_ */
