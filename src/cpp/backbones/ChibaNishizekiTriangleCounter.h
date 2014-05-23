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
template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std
{
  template<typename S, typename T> struct hash<pair<S, T>>
  {
    inline size_t operator()(const pair<S, T> & v) const
    {
      size_t seed = 0;
      ::hash_combine(seed, v.first);
      ::hash_combine(seed, v.second);
      return seed;
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

	edgeCountSet triangleCounts(const Graph& graph);

private:
	void triangleFound(const edgeCountMap& map, const node& u, const node& v, const node& w);
	void removeNode(Graph& graph, const node& u);
};

} /* namespace NetworKit */
#endif /* CHIBANISHIZEKI_H_ */
