/*
 * LocalFilterAttributizer.h
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann, Gerd Lindner
 */

#ifndef LOCALLOGATTRIBUTIZER_H
#define LOCALLOGATTRIBUTIZER_H

#include "../edgescores/EdgeAttribute.h"

namespace NetworKit {

template<typename InType>
class LocalFilterAttributizer : public EdgeAttribute<double> {

public:
	LocalFilterAttributizer(const Graph &graph, const std::vector< InType > &attribute, bool logarithmic = true, bool bothRequired = false) :
		graph(graph), attribute(attribute), bothRequired(bothRequired), logarithmic(logarithmic) {}

	virtual std::vector< double > getAttribute() {
		if (!graph.hasEdgeIds()) {
			throw std::runtime_error("edges have not been indexed - call indexEdges first");
		}

		/*
		* For each edge, we calculate the minimum required sparsification exponent e
		* such that the edge is contained in the backbone.
		*/

		std::vector<double> sparsificationExp(graph.upperEdgeIdBound(), (bothRequired ? .0 : 1.0));

		count n = graph.numberOfNodes();

		graph.parallelForNodes([&](node i) {
			count d = graph.degree(i);

			/*
			 * The top d^e edges (sorted by similarity in descending order)
			 * are to be kept in the backbone.
			 */

			std::vector<std::tuple<double, count, index> > neighbors;
			graph.forNeighborsOf(i, [&](node _i, node j, edgeid eid) {
				neighbors.emplace_back(attribute[eid], n - graph.degree(j), eid);
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

		return sparsificationExp;

	}

private:
	const Graph& graph;
	const std::vector<InType>& attribute;
	bool bothRequired;
	bool logarithmic;

};

} // namespace NetworKit

#endif // LOCALLOGATTRIBUTIZER_H
