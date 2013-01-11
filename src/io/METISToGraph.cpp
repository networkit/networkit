/*
 * METISToGraph.cpp
 *
 *  Created on: 27.11.2012
 *      Author: cls
 */

#include "METISToGraph.h"

#include <string>
#include <vector>

#include "../graph/Graph.h"
#include "GraphFromAdjacencies.h"
#include "METISParser.h"

namespace EnsembleClustering {

METISToGraph::METISToGraph() {
	// TODO Auto-generated constructor stub

}

METISToGraph::~METISToGraph() {
	// TODO Auto-generated destructor stub
}

Graph METISToGraph::read(std::string graphPath) {

	METISParser parser;



	parser.open(graphPath);
	std::pair<int64_t, int64_t> header = parser.getHeader();

	GraphFromAdjacencies builder;

	while (parser.hasNext()) {
		builder.addAdjacencies(parser.getNext());
	}

	parser.close();

	return builder.getGraph();
}

} /* namespace EnsembleClustering */
