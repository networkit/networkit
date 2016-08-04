/*
 * APP.h
 *
 *  Created on: 15.04.2014
 *      Author: cls
 */

#ifndef APP_H_
#define APP_H_

#include <set>
#include <stack>

#include "Graph.h"
#include "../base/Algorithm.h"


namespace NetworKit {

/**
 * @ingroup graph
 * Abstract base class for single-source algebraic path algorithms.
 */
template <class T>
class APP: public Algorithm {

public:

	/**
	 * Creates the APP class for @a G and source @a s.
	 *
	 * @param G The graph.
	 * @param s The source node.
	 */
    APP(const Graph& G, node s, bool storePaths=true, bool storeStack=false, node target=none, std::vector<T> edgeWeights={}) : Algorithm(), G(G), source(s), target(target), storePaths(storePaths), storeStack(storeStack), edgeWeights(edgeWeights) {
}
	virtual ~APP() = default;

	/** Computes the shortest paths from the source to all other nodes. */
	virtual void run() = 0;

	/**
	 * Returns a vector of weighted distances from the source node, i.e. the
 	 * length of the shortest path from the source node to any other node.
 	 *
 	 * @param moveOut If set to true, the container will be moved out of the class instead of copying it; default=true.
 	 * @return The weighted distances from the source node to any other node in the graph.
	 */
    std::vector<T> getDistances(bool moveOut=true) {
        return (moveOut)?std::move(distances):distances;
    }

	/**
	 * Returns the distance from the source node to @a t.
	 * @param  t Target node.
	 * @return The distance from source to target node @a t.
	 */
	T distance(node t) const {
        return distances[t];
    }

	/**
	 * Returns the number of shortest paths between the source node and @a t.
	 * @param  t Target node.
	 * @return The number of shortest paths between source and @a t.
	 */
	bigfloat numberOfPaths(node t) const;

	/**
	 * Returns the number of shortest paths between the source node and @a t
	 * as a double value. Workaround for Cython
	 * @param  t Target node.
	 * @return The number of shortest paths between source and @a t.
	 */
	double _numberOfPaths(node t) const;

	/**
	 * Returns the predecessor nodes of @a t on all shortest paths from source to @a t.
	 * @param t Target node.
	 * @return The predecessors of @a t on all shortest paths from source to @a t.
	 */
	std::vector<node> getPredecessors(node t) const;

	/**
	 * Returns a shortest path from source to @a t and an empty path if source and @a t are not connected.
	 *
	 * @param t Target node.
	 * @param forward If @c true (default) the path is directed from source to @a t, otherwise the path is reversed.
	 * @return A shortest path from source to @a t or an empty path.
	 */
    std::vector<node> getPath(node t, bool forward) const {
        if (! storePaths) {
            throw std::runtime_error("paths have not been stored");
        }
        std::vector<node> path;
        if (previous[t].empty()) { // t is not reachable from source
            WARN("there is no path from ", source, " to ", t);
            return path;
        }
        node v = t;
        while (v != source) {
            path.push_back(v);
            v = previous[v].front();
        }
        path.push_back(source);

        if (forward) {
            std::reverse(path.begin(), path.end());
        }
        return path;
    }

	/**
	 * Returns all shortest paths from source to @a t and an empty set if source and @a t are not connected.
	 *
	 * @param t Target node.
	 * @param forward If @c true (default) the path is directed from source to @a t, otherwise the path is reversed.
	 * @return All shortest paths from source node to target node @a t.
	 */
    std::set<std::vector<node>> getPaths(node t, bool forward) const {

        std::set<std::vector<node>> paths;
        if (previous[t].empty()) { // t is not reachable from source
            WARN("there is no path from ", source, " to ", t);
            return paths;
        }

        std::function<void (std::vector<node> suffix, node v) > trace = [&](std::vector<node> suffix, node v) {
            // base case
            suffix.push_back(v);
            if (v == source) {
                paths.insert(suffix);
                return;
            }
            for (node u : previous[v]) {
                trace(suffix, u);
            }
        };

        std::vector<node> emptyPath;
        trace(emptyPath, t);

        if (forward) {
            std::set<std::vector<node>> paths1;
            for (std::vector<node> path : paths) {
                std::reverse(std::begin(path), std::end(path));
                paths1.insert(path);
            }
            return paths1;
        }

	return paths;
}

	/* Returns the number of shortest paths to node t.*/
	bigfloat getNumberOfPaths(node t) const;

	/**
	* Returns a stack of nodes ordered in decreasing distance from the source
	*
	* @param moveOut If set to true, the container will be moved out of the class instead of copying it; default=true.
	* @return stack of nodes
	*/
    std::vector<node> getStack(bool moveOut=false) {
        if (!storeStack) {
            throw std::runtime_error("stack has not been stored");
        }
        return (moveOut)?std::move(stack):stack;
    }

protected:

	const Graph& G;
	const node source;
	node target;
	std::vector<T> distances;
	std::vector<std::vector<node> > previous; // predecessors on shortest path
	std::vector<bigfloat> npaths;
    std::vector<T> edgeWeights;

	std::vector<node> stack;

	bool storePaths;		//!< if true, paths are reconstructable and the number of paths is stored
	bool storeStack;		//!< if true, store a stack of nodes ordered in decreasing distance from the source
};

template <class T>
inline bigfloat APP<T>::numberOfPaths(node t) const {
	if (! storePaths) {
		throw std::runtime_error("number of paths have not been stored");
	}
	return npaths[t];
}

template <class T>
inline double APP<T>::_numberOfPaths(node t) const {
	if (! storePaths) {
		throw std::runtime_error("number of paths have not been stored");
	}
	bigfloat limit = std::numeric_limits<double>::max();
	if (npaths[t] > limit) {
		throw std::overflow_error("number of paths do not fit into a double");
	}
	double res;
	npaths[t].ToDouble(res);
	return res;
}

template <class T>
inline std::vector<node> APP<T>::getPredecessors(node t) const {
	if (! storePaths) {
		throw std::runtime_error("predecessors have not been stored");
	}
	return previous[t];
}

template <class T>
inline bigfloat APP<T>::getNumberOfPaths(node t) const {
	return npaths[t];
}

} /* namespace NetworKit */

#endif /* APP_H_ */
