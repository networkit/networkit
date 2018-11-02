/*
 * GraphIO.h
 *
 *  Created on: 09.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef GRAPHIO_H_
#define GRAPHIO_H_

#include <string>
#include <iostream>
#include <fstream>

#include "../graph/Graph.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

/**
 * @ingroup io
 */
class GraphIO {
public:
	/**
	 * Writes graph to text file in edge list format.
	 * Keep in mind that isolated nodes are ignored.
	 *
	 * @param[in]	G	graph
	 * @param[in]	path	file path
	 *
	 * Edge list format:
	 * 		for each edge {u, v}:
	 * 			 write line "u v"
	 */
	virtual void writeEdgeList(Graph& G, std::string path);


	/**
	 * Writes graph to text file in adjacency list format.
	 *
	 * @param[in]	G	graph
	 * @param[in]	path	file path
	 *
	 * Adjacency list format:
	 * 		for each node v:
	 * 			write "v"
	 * 			write " x" for each edge {v, x}
	 * 			end line
	 */
	virtual void writeAdjacencyList(Graph& G, std::string path);


};

} /* namespace NetworKit */
#endif /* GRAPHIO_H_ */
