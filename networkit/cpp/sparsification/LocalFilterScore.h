/*
 * LocalFilterScore.h
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann, Gerd Lindner
 */

#ifndef LOCALLOGSCORE_H
#define LOCALLOGSCORE_H

#include "../edgescores/EdgeScore.h"
#include "../auxiliary/Parallel.h"

namespace NetworKit {

/**
 * Local filtering edge scoring. Edges with high score are more important.
 *
 * Edges are ranked locally, the top d^e (logarithmic, default) or 1+e*(d-1) edges (non-logarithmic) are kept.
 * For equal attribute values, neighbors of low degree are preferred.
 * If bothRequired is set (default: false), both neighbors need to indicate that they want to keep the edge.
 */
template<typename InType>
class LocalFilterScore : public EdgeScore<double> {

public:
	/**
	 * Initialize the local edge filtering score.
	 *
	 * @param G The graph for which the score shall be.
	 * @param attribute The input attribute according to which the edges shall be fitlered locally.
	 * @param logarithmic If the score shall be logarithmic in the rank (then d^e edges are kept). Linear otherwise.
	 * @param bothRequired if both neighbors need to indicate that they want to keep an edge (default: one suffices).
	 */
	LocalFilterScore(const Graph& G, const std::vector< InType > &attribute, bool logarithmic = true, bool bothRequired = false) :
		EdgeScore<double>(G), attribute(attribute), bothRequired(bothRequired), logarithmic(logarithmic) {}

	/**
	 * Execute the algorithm.
	 */
	virtual void run() {
		if (!G.hasEdgeIds()) {
			throw std::runtime_error("edges have not been indexed - call indexEdges first");
		}

		/*
		* For each edge, we calculate the minimum required sparsification exponent e
		* such that the edge is contained in the sparse graph.
		*/

		std::vector<double> sparsificationExp(G.upperEdgeIdBound(), (bothRequired ? 1.0 : .0));
		count n = G.numberOfNodes();

		G.balancedParallelForNodes([&](node i) {
			count d = G.degree(i);

			/*
			 * The top d^e edges (sorted by similarity in descending order)
			 * are to be kept in the sparse graph.
			 */

			std::vector<std::tuple<double, count, index> > neighbors;
			neighbors.reserve(d);
			G.forNeighborsOf(i, [&](node _i, node j, edgeid eid) {
				neighbors.emplace_back(attribute[eid], n - d, eid); // if in doubt, prefer links to low-degree nodes
			});
			Aux::Parallel::sort(neighbors.begin(), neighbors.end(), std::greater<std::tuple<double, count, index> >());

			count rank = 1;

			#pragma omp critical // each value is set twice, the value can be wrong if the wrong thread wins
			for (auto it : neighbors) {
				edgeid eid = std::get<2>(it);

				double e = 1.0;

				if (d > 1) {
					if (logarithmic) {
						e = 1.0 - log(rank) / log(d);
					} else {
						e = 1.0 - (rank-1) * 1.0 / (d - 1); // Keep top 1 + e * (d-1) edges
					}
				}

				if ((e < sparsificationExp[eid]) == bothRequired) {
					sparsificationExp[eid] = e; // do not always write in order to avoid cache synchronization
				}

				rank++;
			}

		});

		scoreData = std::move(sparsificationExp);
		hasRun = true;
	}

	virtual double score(node u, node v) {
		throw std::runtime_error("Not implemented: Use scores() instead.");
	}

	virtual double score(edgeid eid) {
		throw std::runtime_error("Not implemented: Use scores() instead.");
	}

private:
	const std::vector<InType>& attribute;
	bool bothRequired;
	bool logarithmic;

};

} // namespace NetworKit

#endif // LOCALLOGSCORE_H
