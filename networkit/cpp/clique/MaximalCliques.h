#ifndef MAXIMAL_CLIQUES_H_
#define MAXIMAL_CLIQUES_H_

#include "../graph/Graph.h"
#include "../base/Algorithm.h"
#include <functional>

namespace NetworKit {

/**
 * Algorithm for listing all maximal cliques.
 *
 * The implementation is based on the "hybrid" algorithm described in
 *
 * Eppstein, D., & Strash, D. (2011).
 * Listing All Maximal Cliques in Large Sparse Real-World Graphs.
 * In P. M. Pardalos & S. Rebennack (Eds.),
 * Experimental Algorithms (pp. 364–375). Springer Berlin Heidelberg.
 * Retrieved from http://link.springer.com/chapter/10.1007/978-3-642-20662-7_31
 *
 * The running time of this algorithm should be in @f$O(d^2 \cdot n \cdot 3^{d/3})@f$
 * where @f$d@f$ is the degeneracy of the graph, i.e., the maximum core number.
 * The running time in practive depends on the structure of the graph. In
 * particular for complex networks it is usually quite fast, even graphs with
 * millions of edges can usually be processed in less than a minute.
 */
class MaximalCliques : public Algorithm {

public:
	/**
	 * Construct the maximal cliques algorithm with the given graph.
	 *
	 * If the @a maximumOnly argument is set, the algorithm will only store
	 * the clique of maximum size. Further, this enables some optimizations
	 * to skip smaller cliques more efficiently leading to a reduced
	 * running time.
	 *
	 * @param G The graph to list the cliques for.
	 * @param maximumOnly If only a maximum clique shall be found.
	 */
	MaximalCliques(const Graph& G, bool maximumOnly = false);

	/**
	 * Construct the maximal cliques algorithm with the given graph and a callback.
	 *
	 * The callback is called once for each found clique with a reference to the clique.
	 * Note that the reference is to an internal object, the callback should not assume that
	 * this reference is still valid after it returned.
	 *
	 * @param G The graph to list cliques for
	 * @param callback The callback to call for each clique.
	 */
	MaximalCliques(const Graph& G, std::function<void(const std::vector<node>&)> callback);

	/**
	 * Execute the maximal clique listing algorithm.
	 */
	void run() override;

	/**
	 * Return all found cliques unless a callback was given.
	 *
	 * This method will throw if a callback was given and thus the cliques were not stored.
	 * If only the maximum clique was stored, it will return exactly one clique unless the graph
	 * is empty.
	 *
	 * @return a vector of cliques, each being represented as a vector of nodes.
	 */
	const std::vector<std::vector<node>>& getCliques() const;

protected:
	const Graph& G;

	std::vector<std::vector<node>> result;

	std::function<void(const std::vector<node>&)> callback;
	bool maximumOnly;
};

}

#endif /* MAXIMAL_CLIQUES_H_ */
