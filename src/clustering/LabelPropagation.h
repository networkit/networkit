/*
 * LabelPropagation.h
 *
 *  Created on: 07.12.2012
 *      Author: cls
 */

#ifndef LABELPROPAGATION_H_
#define LABELPROPAGATION_H_

#include "Clusterer.h"
#include "Clustering.h"

#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <stdexcept>

#include "../graph/NodeMap.h"
#include "../base/IndexMap.h"
#include "../aux/log.h"

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

public:

	LabelPropagation();

	virtual ~LabelPropagation();

	virtual Clustering run(Graph& G);
};

} /* namespace EnsembleClustering */
#endif /* LABELPROPAGATION_H_ */
