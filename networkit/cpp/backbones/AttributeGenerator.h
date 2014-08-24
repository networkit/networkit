/*
 * AttributeGenerator.h
 *
 *  Created on: 22.07.2014
 *      Author: Gerd Lindner
 */

#ifndef ATTRIBUTEGENERATOR_H_
#define ATTRIBUTEGENERATOR_H_

#include "../graph/Graph.h"
#include <unordered_map>
#include <vector>

namespace NetworKit {

//TODO: he following will later be realized using graph edge attributes.

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
	 * We define a lexicographic ordering on the edges.
	 */
	bool operator<(const uEdge& other) const {
		return std::tie(lowNode, highNode) < std::tie(other.lowNode, other.highNode);
	}

	bool operator>(const uEdge& other) const {
		return std::tie(lowNode, highNode) > std::tie(other.lowNode, other.highNode);
	}
};

} /* namespace NetworKit */



namespace std {

template <class T>
inline void hash_combine(std::size_t & s, const T & v) {
	std::hash<T> h;
	s^= h(v) + 0x9e3779b9 + (s<< 6) + (s>> 2);
}

template<typename S, typename T> struct hash<pair<S, T>> {
inline size_t operator()(const pair<S, T> & v) const {
	return v.first + v.second;
}
};

template <>
struct hash<NetworKit::uEdge> {
	std::size_t operator()(const NetworKit::uEdge& edge) const {
    	std::size_t res = 0;
    	hash_combine(res, edge.lowNode);
    	hash_combine(res, edge.highNode);
        return res;
    }
};

} /* namespace std */

namespace NetworKit {
/**
 * Abstract base class for graph attribute generator. It takes a graph (weighted or unweighted)
 * and calculates a graph attribute from the input graph.
 */
template<typename TInput, typename TOutput>
class AttributeGenerator {

public:

	/**
	 * Calculates an edge attribute for the edges of the given graph.
	 * (Possibly under consideration of the given attribute).
	 */
	virtual std::vector<TOutput> getAttribute(const Graph& g, const std::vector<TInput>& attribute) = 0;

	virtual ~AttributeGenerator() = default;

	std::vector<TOutput>* _getAttribute(Graph& g, std::vector<TInput>& attribute) {
		return new std::vector<TOutput>{std::move(getAttribute(g, attribute))};
	};

};

}


#endif /* ATTRIBUTEGENERATOR_H_ */
