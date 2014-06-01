/*
 * GraphData.h
 *
 *  Created on: 01.06.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef GRAPHDATA_H_
#define GRAPHDATA_H_

#include <vector>

#include "../Globals.h"

namespace NetworKit {

namespace graph_impl {

// data structures for special graph classes, the BasicGraph class will inherit this private as needed
struct UnweightedData {
	count getMemoryUsage() const { return 0; }
	void shrinkToFit() {}
	void addNode() {}
};

struct WeightedData {
	std::vector< std::vector<edgeweight> > edgeWeights;

	count getMemoryUsage() const;
	void shrinkToFit();
	void addNode();
};

struct UndirectedData {
	UndirectedData(count n = 0) :
		deg(n, 0),
		adja(n)
	{}
	std::vector<count> deg;
	std::vector< std::vector<node> > adja;

	count getMemoryUsage() const;
	void shrinkToFit();
	void addNode();
	index indexInEdgeArray(node u, node v) const;
};

struct DirectedData {
	DirectedData(count n = 0) :
		inDeg(n, 0),
		outDeg(n, 0),
		inEdges(n),
		outEdges(n)
	{}
	std::vector<count> inDeg;
	std::vector<count> outDeg;
	std::vector< std::vector<node> > inEdges;
	std::vector< std::vector<node> > outEdges;

	count getMemoryUsage() const;
	void shrinkToFit();
	void addNode();
	index indexInEdgeArray(node u, node v) const { return indexInOutEdgeArray(u, v); }
	index indexInInEdgeArray(node u, node v) const;
	index indexInOutEdgeArray(node u, node v) const;
};

} /* namespace graph_impl */

} /* namespace NetworKit */

#endif /* GRAPHDATA_H_ */
