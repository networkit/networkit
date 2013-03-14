/*
 * GraphContraction.h
 *
 *  Created on: 25.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef GRAPHCONTRACTION_H_
#define GRAPHCONTRACTION_H_

#include "../graph/NodeMap.h"

namespace EnsembleClustering {

class GraphContraction {

protected:

	Graph fine;	//!< the fine graph
	Graph coarse;	//!< the coarse graph
	NodeMap<node> fineToCoarse;	//!< maps each node in the fine graph to a supernode in the coarse graph

public:

	GraphContraction(Graph& fine, Graph& coarse, NodeMap<node>& fineToCoarse);

	virtual ~GraphContraction();

	virtual Graph& getFineGraph();

	virtual Graph& getCoarseGraph();

	virtual NodeMap<node>& getFineToCoarseMap();

};

} /* namespace EnsembleClustering */
#endif /* GRAPHCONTRACTION_H_ */
