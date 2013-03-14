/*
 * LabelPropagation.h
 *
 *  Created on: 07.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef LABELPROPAGATION_H_
#define LABELPROPAGATION_H_

#include "Clusterer.h"
#include "../base/Clustering.h"

#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <random>
#include <omp.h>

#include "../../Globals.h"
#include "../../aux/Log.h"
#include "../../aux/ProgressMeter.h"
#include "../../aux/Timer.h"
#include "../../aux/RandomInteger.h"
#include "../../graph/NodeMap.h"
#include "../../base/IndexMap.h"
#include "../../io/GraphIO.h"


namespace EnsembleClustering {

/**
 *     As described in Ovelgoenne et al: An Ensemble Learning Strategy for Graph Clustering
        Raghavan et al. proposed a label propagation algorithm for graph clustering.
        This algorithm initializes every vertex of a graph with a unique label. Then, in iterative
        sweeps over the set of vertices the vertex labels are updated. A vertex gets the label
        that the maximum number of its neighbors have. The procedure is stopped when every vertex
        has the label that at least half of its neighbors have.
 *
 *
 */
class LabelPropagation: public EnsembleClustering::Clusterer {

protected:

	count updateThreshold = 0;

public:

	LabelPropagation(count theta = 0);

	virtual ~LabelPropagation();

	/**
	 * Run the label propagation clustering algorithm.
	 *
	 * @param[in]	G	input graph
	 * @return			clustering
	 */
	virtual Clustering run(Graph& G);

	/**
	 * @return string representation of algorithm and parameters.
	 */
	virtual std::string toString() const;


	/**
	 * The algorithm runs until a number of nodes less than
	 * the threshold is updated.
	 *
	 */
	virtual void setUpdateThreshold(count th);
};

} /* namespace EnsembleClustering */
#endif /* LABELPROPAGATION_H_ */
