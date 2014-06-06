/*
 * TriangleCounter.h
 *
 *  Created on: 22.05.2014
 *      Author: Gerd Lindner
 */

#ifndef TRIANGLECOUNTER_H_
#define TRIANGLECOUNTER_H_

#include "../graph/Graph.h"
#include <unordered_map>

namespace NetworKit {

/**
 * Represents an undirected edge.
 */
struct uEdge
{
	node lowNode;
	node highNode;

	uEdge(node a, node b) {
		if (a < b) {
			lowNode = a;
			highNode = b;
		} else {
			lowNode = b;
			highNode = a;
		}
	}

	bool operator==(const uEdge& rhs) const {
		return lowNode == rhs.lowNode && highNode == rhs.highNode;
	}

	/**
	 * We define an ordering on uEdge: Sort by lowNode, then by highNode.
	 */
	bool operator<(const uEdge& other) const {
		return lowNode < other.lowNode || (lowNode == other.lowNode && highNode < other.highNode);
	}
	bool operator>(const uEdge& other) const {
		return lowNode > other.lowNode || (lowNode == other.lowNode && highNode > other.highNode);
	}
};

typedef std::unordered_map<uEdge, count> edgeCountMap; //map: edge --> triangle count

/** 
 * Abstract base class for per-edge triangle counting algorithms.
 */
class TriangleCounter {

public:
	/**
	 * Calculates triangle counts for the edges of the given graph.
	 */
	virtual edgeCountMap triangleCounts(const Graph& g) = 0;
};

} /* namespace NetworKit */



namespace std
{

template <class T>
inline void hash_combine(std::size_t & s, const T & v)
{
  std::hash<T> h;
  s^= h(v) + 0x9e3779b9 + (s<< 6) + (s>> 2);
}

template<typename S, typename T> struct hash<pair<S, T>>
  {
    inline size_t operator()(const pair<S, T> & v) const
    {
      return v.first + v.second;
    }
  };

template <>
struct hash<NetworKit::uEdge>
{
    std::size_t operator()(const NetworKit::uEdge& edge) const
    {
    	std::size_t res = 0;
    	hash_combine(res, edge.lowNode);
    	hash_combine(res, edge.highNode);
        return res;
    }
};
} /* namespace std */

#endif /* TRIANGLECOUNTER_H_ */
