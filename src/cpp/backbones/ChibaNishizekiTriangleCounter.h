/*
 * ChibaNishizekiTriangleCounter.h
 *
 *  Created on: 22.05.2014
 *      Author: Gerd Lindner
 */

#ifndef CHIBANISHIZEKI_H_
#define CHIBANISHIZEKI_H_

#include "../graph/Graph.h"
#include "TriangleCounter.h"
#include <unordered_map>
#include <utility> //for pair
#include <set>

namespace NetworKit {

typedef std::vector<std::pair<std::pair<node,node>, count>> edgeCountSet; //TEMPORARY return type, for debug/output purposes only.

/** 
 * An implementation of the triangle counting algorithm by Chiba/Nishizeki.
 */
class ChibaNishizekiTriangleCounter : public TriangleCounter {

public:

	edgeCountMap triangleCounts(const Graph& graph);
	edgeCountSet triangleCountsDebug(const Graph& graph); //TEMPORARY

private:
	void triangleFound(edgeCountMap& map, const node& u, const node& v, const node& w);
	void removeNode(Graph& graph, const node& u);
};

} /* namespace NetworKit */
#endif /* CHIBANISHIZEKI_H_ */
