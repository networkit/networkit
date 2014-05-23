/*
 * ChibaNishizekiTriangleCounter.h
 *
 *  Created on: 22.05.2014
 *      Author: Gerd Lindner
 */

#ifndef CHIBANISHIZEKI_H_
#define CHIBANISHIZEKI_H_

#include "../graph/Graph.h"
#include <unordered_map>
#include <utility> //for pair
#include <set>

//TODO: this should not be defined here (stolen from ClusteringCut.h)

namespace std
{
  template<typename S, typename T> struct hash<pair<S, T>>
  {
    inline size_t operator()(const pair<S, T> & v) const
    {
      return v.first + v.second;
    }
  };
}

namespace NetworKit {

typedef std::unordered_map<std::pair<node,node>, count> edgeCountMap; //map: edge --> triangle count
typedef std::vector<std::pair<std::pair<node,node>, count>> edgeCountSet; //TEMPORARY return type, for debug/output purposes only.

/** 
 * An implementation of the triangle counting algorithm by Chiba/Nishizeki.
 */
class ChibaNishizekiTriangleCounter {

public:
	ChibaNishizekiTriangleCounter();

	~ChibaNishizekiTriangleCounter();

	edgeCountMap triangleCounts(const Graph& graph);
	edgeCountSet triangleCountsDebug(const Graph& graph); //TEMPORARY

private:
	void triangleFound(edgeCountMap& map, const node& u, const node& v, const node& w);
	void removeNode(Graph& graph, const node& u);
};

} /* namespace NetworKit */
#endif /* CHIBANISHIZEKI_H_ */
