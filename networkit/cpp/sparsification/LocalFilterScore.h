/*
 * LocalFilterScore.h
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann, Gerd Lindner
 */

#ifndef LOCALLOGSCORE_H
#define LOCALLOGSCORE_H

#include "../edgescores/EdgeScore.h"

namespace NetworKit {

template<typename InType>
class LocalFilterScore : public EdgeScore<double> {

public:
	LocalFilterScore(const Graph& G, const std::vector< InType > &attribute, bool logarithmic = true, bool bothRequired = false) :
		EdgeScore<double>(G), attribute(attribute), bothRequired(bothRequired), logarithmic(logarithmic) {}

	virtual void run() {
		if (!G.hasEdgeIds()) {
			throw std::runtime_error("edges have not been indexed - call indexEdges first");
		}

		/*
		* For each edge, we calculate the minimum required sparsification exponent e
		* such that the edge is contained in the sparse graph.
		*/

		std::vector<double> sparsificationExp(G.upperEdgeIdBound(), (bothRequired ? .0 : 1.0));
		count n = G.numberOfNodes();

		G.parallelForNodes([&](node i) {
			count d = G.degree(i);

			/*
			 * The top d^e edges (sorted by similarity in descending order)
			 * are to be kept in the sparse graph.
			 */

			std::vector<std::tuple<double, count, index> > neighbors;
			G.forNeighborsOf(i, [&](node _i, node j, edgeid eid) {
				neighbors.emplace_back(attribute[eid], n - G.degree(j), eid);
			});
			std::sort(neighbors.begin(), neighbors.end(), std::greater<std::tuple<double, count, index> >());

			count rank = 1;

			#pragma omp critical // each value is set twice, the value can be wrong if the wrong thread wins
			for (auto it : neighbors) {
				edgeid eid = std::get<2>(it);

				double e = 0.0;

				if (d > 1) {
					if (logarithmic) {
						e = log(rank) / log(d);
					} else {
						e = (rank-1) * 1.0 / (d - 1); // Keep top 1 + e * (d-1) edges
					}
				}

				if ((e > sparsificationExp[eid]) == bothRequired) {
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
